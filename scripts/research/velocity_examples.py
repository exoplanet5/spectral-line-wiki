#!/usr/bin/env python
"""
velocity_examples.py -- velocity reference frames for precise spectroscopy.

Computes (never invents) the numbers used in the "wavelength calibration /
velocity frames" section of the Spectral Line Wiki:

  1. Observer terms: Earth rotation, Earth orbit, Earth-Moon barycentre wobble,
     Sun about the solar-system barycentre (SSB), light-time, relativistic
     terms of the barycentric correction (Wright & Eastman 2014, eq. 28),
     timing sensitivity and exposure-midpoint curvature.
  2. A check of astropy's SkyCoord.radial_velocity_correction against an
     independent implementation of WE14 eq. 28 driven directly by the JPL
     DE440s kernel through jplephem, and builtin-vs-DE440s ephemeris.
  3. Solar-motion / LSR / Galactic-rotation / Local-Group / CMB apex vectors
     from the literature (values typed from the cited papers, see SOURCES),
     converted with astropy to Galactic (U,V,W) and ICRS (RA, Dec) apices.
  4. Worked example: target ICRS RA=10h, Dec=+20d observed from Mauna Kea
     (Keck site) at 2026-03-20T10:00:00 UTC; projected corrections for every
     frame, Monte-Carlo propagated apex uncertainties, astropy frame cross-checks,
     Mould et al. (2000) multi-attractor model for a hypothetical galaxy,
     exact multiplicative CMB-frame redshift (Davis+2011; Carr+2022).
  5. Cosmology-scale secular terms: SSB acceleration (Gaia EDR3) projected on
     the target and Sandage-Loeb redshift drift for astropy's Planck18.

Outputs: a Markdown note and a JSON file with every computed value.

Usage
  ~/.venvs/astro313/bin/python velocity_examples.py \
      --ephem ~/JPLDE/de440s.bsp \
      --out-md   data/notes/velocity_frames.md \
      --out-json data/notes/velocity_frames.json

All literature inputs are in the LITERATURE dict below with a citation key;
everything else is computed here and labelled "derived".
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import math
import warnings

import numpy as np
import astropy
import astropy.units as u
import astropy.constants as const
from astropy.time import Time
from astropy.coordinates import (
    SkyCoord, EarthLocation, ICRS, Galactic, FK4, LSRK, LSRD, LSR, AltAz,
    Galactocentric, CartesianRepresentation, CartesianDifferential, GCRS, TETE,
    UnitSphericalRepresentation, solar_system_ephemeris,
    galactocentric_frame_defaults, get_body_barycentric_posvel,
)
from jplephem.spk import SPK

C_MS = const.c.to_value(u.m / u.s)          # 299792458 m/s (exact)
C_KMS = C_MS / 1e3
AU_M = const.au.to_value(u.m)                 # IAU 2012 exact
K_PM = 4.740470463533348                      # km/s per (mas/yr * kpc) = 1 au/yr in km/s (derived: au/(Julian yr))

# ---------------------------------------------------------------------------
# Literature inputs (typed from the papers; see the SOURCES section of the md)
# l, b in degrees (Galactic), V in km/s; uvw = (U toward GC, V toward l=90, W toward NGP)
# ---------------------------------------------------------------------------
LITERATURE = {
    # ---- Earth constants
    "omega_earth": dict(value=7.292115e-5, unit="rad/s", src="IERS Conventions (2010) TN36 Table 1.2 (GRS80) / WGS84: nominal mean angular velocity"),
    "a_wgs84": dict(value=6378137.0, unit="m", src="IERS TN36 Table 1.2 / WGS84 equatorial radius"),
    "earth_orbit_factsheet": dict(mean=29.78, max=30.29, min=29.29, unit="km/s", src="NASA NSSDCA Earth Fact Sheet"),
    # ---- Solar motion w.r.t. LSR (peculiar / standard)
    "LSRK_Gordon1975": dict(V=20.0, ra_b1900=270.0, dec_b1900=30.0, src="Gordon 1975 (Meth. Exp. Phys. 12C, sec 6.1.5); astropy LSRK; SLALIB sla_RVLSRK"),
    "LSRD_Delhaye1965": dict(uvw=(9.0, 12.0, 7.0), src="Delhaye 1965 (Galactic Structure, Stars & Stellar Systems vol 5); astropy LSRD; SLALIB sla_RVLSRD"),
    "LSR_Schoenrich2010": dict(uvw=(11.1, 12.24, 7.25), stat=((0.75, 0.69), (0.47, 0.47), (0.36, 0.37)), sys=(1.0, 2.0, 0.5),
                               src="Schoenrich, Binney & Dehnen 2010 MNRAS 403, 1829 (abstract); astropy LSR / GalacticLSR default"),
    "LSR_DehnenBinney1998": dict(uvw=(10.00, 5.25, 7.17), err=(0.36, 0.62, 0.38), src="Dehnen & Binney 1998 MNRAS 298, 387 (published abstract; arXiv v1 abstract gives V0=5.23)"),
    "LSR_BHG2016": dict(uvw=(10.0, 11.0, 7.0), err=(1.0, 2.0, 0.5), src="Bland-Hawthorn & Gerhard 2016 ARA&A 54, 529 (adopted)"),
    "LSR_Ding2019_GaiaDR2": dict(uvw=(8.63, 4.76, 7.26), err=(0.64, 0.49, 0.36), src="Ding, Zhu & Liu 2019 RAA 19, 068 (Gaia DR2)"),
    # ---- Galactic rotation / Sun about Galactic centre
    "IAU1985": dict(R0=8.5, Theta0=220.0, src="Kerr & Lynden-Bell 1986 MNRAS 221, 1023 (IAU Commission 33, 1985 Nov 21, Delhi)"),
    "NED_GSR_RC3": dict(V=232.3, l=87.8, b=1.7, src="NED velocity-correction help (RC3, de Vaucouleurs et al. 1991)"),
    "astropy_v4.0": dict(galcen_v_sun=(12.9, 245.6, 7.78), R0=8.122, z_sun=20.8, galcen_radec=(266.4051, -28.936175),
                         src="astropy Galactocentric 'v4.0' defaults: Drimmel & Poggio 2018 RNAAS 2, 210; GRAVITY 2018 A&A 615, L15; Reid & Brunthaler 2004; Bennett & Bovy 2019"),
    "astropy_pre-v4.0": dict(galcen_v_sun=(11.1, 232.24, 7.25), R0=8.3, z_sun=27.0, src="astropy 'pre-v4.0' defaults (Schoenrich+2010; Bovy 2015; Gillessen+2009; Chen+2001)"),
    "Reid2019": dict(R0=8.15, R0_err=0.15, Theta0=236.0, Theta0_err=7.0, Theta0_plus_Vsun=247.0, err=4.0, Omega=30.32, Omega_err=0.27, Zsun_pc=5.5, Zsun_err=5.8,
                     src="Reid et al. 2019 ApJ 885, 131 (abstract)"),
    "BHG2016_galactic": dict(R0=8.2, R0_err=0.1, Theta0=238.0, Theta0_err=15.0, Vg_sun=248.0, Vg_sun_err=3.0, z0_pc=25.0, z0_err=5.0,
                             src="Bland-Hawthorn & Gerhard 2016 ARA&A 54, 529"),
    "GRAVITY2018": dict(R0_pc=8122.0, err_pc=31.0, src="GRAVITY Collab. 2018 A&A 615, L15 (table; astropy v4.0 R0)"),
    "GRAVITY2019": dict(R0_pc=8178.0, stat=13.0, sys=22.0, src="GRAVITY Collab. 2019 A&A 625, L10 (abstract)"),
    "GRAVITY2021": dict(R0_pc=8275.0, stat=9.0, sys=33.0, src="GRAVITY Collab. 2021 A&A 647, A59 (eq. 38)"),
    "GRAVITY2022": dict(R0_pc=8277.0, stat=9.0, sys=33.0, src="GRAVITY Collab. 2022 A&A 657, L12 (sec 4; sys from 2021)"),
    "SgrA_pm_RB2020": dict(mu_l=-6.411, mu_l_err=0.008, mu_b=-0.219, mu_b_err=0.007, unit="mas/yr", src="Reid & Brunthaler 2020 ApJ 892, 39 (abstract)"),
    "SgrA_pm_RB2004": dict(mu=6.379, mu_err=0.024, src="Reid & Brunthaler 2004 ApJ 616, 872 (abstract)"),
    # ---- Sun relative to Local Group barycentre
    "LG_RC2_300sinl": dict(V=300.0, l=90.0, b=0.0, src="RC2 (de Vaucouleurs, de Vaucouleurs & Corwin 1976) 'IAU' 300 sin l cos b; quoted by Mould+2000, Courteau & vdB 1999"),
    "LG_YTS1977": dict(xyz=(-79.0, 296.0, -36.0), src="Yahil, Tammann & Sandage 1977 ApJ 217, 903; formula as given in Mould+2000 eq. A1"),
    "LG_KM1996": dict(V=316.0, V_err=5.0, l=93.0, l_err=2.0, b=-4.0, b_err=2.0, src="Karachentsev & Makarov 1996 AJ 111, 794 (abstract); used by NED"),
    "LG_CvdB1999": dict(V=306.0, V_err=18.0, l=99.0, l_err=5.0, b=-4.0, b_err=4.0, src="Courteau & van den Bergh 1999 AJ 118, 337 (abstract)"),
    "LG_Diaz2014": dict(V=299.0, V_err=15.0, l=98.4, l_err=3.6, b=-5.9, b_err=3.0, src="Diaz et al. 2014 MNRAS 443, 1688 (abstract); Planck 2018 I Table 3"),
    # ---- CMB dipole
    "CMB_Planck2018": dict(V=369.82, V_err=0.11, l=264.021, l_err=0.011, b=48.253, b_err=0.005, beta=1.23357e-3, beta_err=0.00036e-3,
                           T_uK=3362.08, T_err=0.99, radec=(167.942, -6.944), radec_err=0.007,
                           src="Planck Collab. 2020 (2018 results I) A&A 641, A1, Table 2/3 & sec 3"),
    "CMB_Fixsen1996_NED": dict(V=371.0, l=264.14, l_err95=0.30, b=48.26, b_err95=0.30, T_mK=3.372, T_err95=0.007, T0=2.728, T0_err95=0.004,
                               src="Fixsen et al. 1996 ApJ 473, 576 (abstract, 95% CL); V=371.0 as used by NED"),
    "CMB_Kogut1993_dipole": dict(T_mK=3.365, T_err=0.027, l=264.4, l_err=0.3, b=48.4, b_err=0.5, src="Kogut et al. 1993 ApJ 419, 1 (abstract)"),
    "LG_CMB_Kogut1993": dict(V=627.0, V_err=22.0, l=276.0, l_err=3.0, b=30.0, b_err=3.0, src="Kogut et al. 1993 ApJ 419, 1 (abstract)"),
    "LG_CMB_Planck2018": dict(V=620.0, V_err=15.0, l=271.9, l_err=2.0, b=29.6, b_err=1.4, src="Planck 2018 I Table 3 (from Diaz+2014 Sun-LG)"),
    "GC_CMB_Planck2018": dict(V=565.0, V_err=5.0, l=265.76, l_err=0.20, b=28.38, b_err=0.28, src="Planck 2018 I Table 3 (Schoenrich+2010 + McMillan 2011)"),
    "Sun_LSR_Planck2018": dict(V=17.9, V_err=2.0, l=48.0, l_err=7.0, b=23.0, b_err=4.0, src="Planck 2018 I Table 3 (from Schoenrich+2010)"),
    "LSR_GC_Planck2018": dict(V=239.0, V_err=5.0, l=90.0, b=0.0, src="Planck 2018 I Table 3 (McMillan 2011)"),
    "M31_vdM2012": dict(V_rad=-109.3, V_rad_err=4.4, V_tan=17.0, V_tan_1sigma_upper=34.3, M_LG=3.17e12, M_LG_err=0.57e12, src="van der Marel et al. 2012 ApJ 753, 8 (abstract)"),
    "LocalSheet_Tully2008": dict(V=631.0, components=(259.0, 185.0, 455.0), src="Tully et al. 2008 ApJ 676, 184 (abstract; Local-Sheet frame)"),
    # ---- Large-scale flow model (Mould et al. 2000, Table A1; RA/Dec B1950; NED corrections)
    "Mould2000": dict(gamma=2.0, attractors={
        "Virgo": dict(ra_b1950="12h28m19s", dec_b1950="+12d40m00s", VH=1035.0, VLG=957.0, Vfid=200.0, radius=10.0, vrange=(600.0, 2300.0)),
        "GreatAttractor": dict(ra_b1950="13h20m00s", dec_b1950="-44d00m00s", VH=4600.0, VLG=4380.0, Vfid=400.0, radius=10.0, vrange=(2600.0, 6600.0)),
        "Shapley": dict(ra_b1950="13h30m00s", dec_b1950="-31d00m00s", VH=13800.0, VLG=13600.0, Vfid=85.0, radius=12.0, vrange=(10000.0, 16000.0)),
    }, src="Mould et al. 2000 ApJ 529, 786, eqs. 1-2, A1-A2, Table A1; NED objresult help (sign/declination corrections; Virgo vrange 600 per NED, -600 in arXiv preprint)"),
    # ---- secular / cosmological
    "SSB_accel_GaiaEDR3": dict(a=2.32e-10, a_err=0.16e-10, ra=269.1, ra_err=5.4, dec=-31.6, dec_err=4.1, src="Gaia Collab., Klioner et al. 2021 A&A 649, A9 (abstract)"),
    "Sun_grav_redshift_LD2003": dict(photosphere=636.486, err=0.024, at_1au=633.5, convective_blueshift=-400.0,
                                     src="Lindegren & Dravins 2003 A&A 401, 1185, sec 3.1-3.2 and footnote 3"),
}

# ---------------------------------------------------------------------------
# Sources (bibcodes validated against DOIs through the ADS link_gateway; see check_bibcodes.py)
# ---------------------------------------------------------------------------
SOURCES = [
 ("Wright, J. T. & Eastman, J. D. 2014, PASP 126, 838-852, 'Barycentric Corrections at 1 cm/s for precise Doppler velocities'", "2014PASP..126..838W", "10.1086/678541", "https://arxiv.org/abs/1409.4774", "zB formula (eq. 2, 10, 11, 27, 28), precision budget"),
 ("Lindegren, L. & Dravins, D. 2003, A&A 401, 1185-1201, 'The fundamental definition of radial velocity' (includes IAU 2000 Resolutions C1, C2)", "2003A&A...401.1185L", "10.1051/0004-6361:20030181", "https://arxiv.org/abs/astro-ph/0302522", "barycentric RV measure, solar gravitational redshift"),
 ("Kanodia, S. & Wright, J. T. 2018, RNAAS 2, 4 (barycorrpy)", "2018RNAAS...2....4K", "10.3847/2515-5172/aaa4b7", "https://arxiv.org/abs/1801.01634", "barycorrpy precision statement"),
 ("Tronsgaard, R., Buchhave, L. A., Wright, J. T., Eastman, J. D., Blackman, R. T. 2019, MNRAS 489, 2395-2402", "2019MNRAS.489.2395T", "10.1093/mnras/stz2181", "https://arxiv.org/abs/1908.00991", "photon-weighted midpoint second-order error"),
 ("Blackman, R. T., Szymkowiak, A. E., Fischer, D. A., Jurgenson, C. A. 2017, ApJ 837, 18", "2017ApJ...837...18B", "10.3847/1538-4357/aa5ead", "https://arxiv.org/abs/1703.00012", "chromatic barycentric correction"),
 ("Blackman, R. T., Ong, J. M. J., Fischer, D. A. 2019, AJ 158, 40", "2019AJ....158...40B", "10.3847/1538-3881/ab24c3", "https://arxiv.org/abs/1906.01653", "chromatic barycentric correction measured with EXPRES"),
 ("Eastman, J., Siverd, R., Gaudi, B. S. 2010, PASP 122, 935-946", "2010PASP..122..935E", "10.1086/655938", "", "BJD time reduction (cited via WE14)"),
 ("Petit, G. & Luzum, B. (eds.) 2010, IERS Conventions (2010), IERS Technical Note 36, Table 1.2", "", "", "https://iers-conventions.obspm.fr/content/tn36.pdf", "omega = 7.292115e-5 rad/s, a_E = 6378137 m"),
 ("NASA NSSDCA Earth Fact Sheet (D. R. Williams)", "", "", "https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html", "Earth orbital speed mean/max/min"),
 ("Gordon 1975, Methods of Experimental Physics vol. 12C (ed. Meeks), sec. 6.1.5(.2), p. 281 - as cited in astropy LSRK and SLALIB sla_RVLSRK", "", "", "https://docs.astropy.org/en/stable/api/astropy.coordinates.LSRK.html", "LSRK definition 20 km/s to 18h,+30 (1900)"),
 ("Delhaye, J. 1965, in Galactic Structure (Stars and Stellar Systems V), eds. Blaauw & Schmidt, p. 61", "1965gast.book...61D", "", "https://irtfweb.ifa.hawaii.edu/~tcs3/tcs3/Misc/slalib_html/node173.html", "LSRD (9,12,7) km/s"),
 ("Wallace, P. T. 2000, SLALIB Starlink User Note 67 (sla_RVLSRK, sla_RVLSRD)", "", "", "https://irtfweb.ifa.hawaii.edu/~tcs3/tcs3/Misc/slalib_html/node174.html", "LSR definitions, sign convention"),
 ("Schoenrich, R., Binney, J., Dehnen, W. 2010, MNRAS 403, 1829-1833", "2010MNRAS.403.1829S", "10.1111/j.1365-2966.2010.16253.x", "https://arxiv.org/abs/0912.3693", "(U,V,W)sun = (11.1, 12.24, 7.25)"),
 ("Dehnen, W. & Binney, J. J. 1998, MNRAS 298, 387-394", "1998MNRAS.298..387D", "10.1046/j.1365-8711.1998.01600.x", "https://arxiv.org/abs/astro-ph/9710077", "(U,V,W) = (10.00, 5.25, 7.17)"),
 ("Ding, P.-J., Zhu, Z., Liu, J.-C. 2019, RAA 19, 068 (Gaia DR2 LSR)", "2019RAA....19...68D", "10.1088/1674-4527/19/5/68", "", "(U,V,W) = (8.63, 4.76, 7.26)"),
 ("Reid, M. J. et al. 2009, ApJ 700, 137-148", "2009ApJ...700..137R", "10.1088/0004-637X/700/1/137", "https://arxiv.org/abs/0902.3913", "standard solar motion (10.3, 15.3, 7.7)"),
 ("Kerr, F. J. & Lynden-Bell, D. 1986, MNRAS 221, 1023-1038", "1986MNRAS.221.1023K", "10.1093/mnras/221.4.1023", "", "IAU 1985 R0 = 8.5 kpc, Theta0 = 220 km/s"),
 ("GRAVITY Collaboration 2018, A&A 615, L15", "2018A&A...615L..15G", "10.1051/0004-6361/201833718", "https://arxiv.org/abs/1807.09409", "R0 = 8122 +- 31 pc (astropy v4.0)"),
 ("GRAVITY Collaboration 2019, A&A 625, L10", "2019A&A...625L..10G", "10.1051/0004-6361/201935656", "https://arxiv.org/abs/1904.05721", "R0 = 8178 +- 13 +- 22 pc"),
 ("GRAVITY Collaboration 2021, A&A 647, A59", "2021A&A...647A..59G", "10.1051/0004-6361/202040208", "https://arxiv.org/abs/2101.12098", "R0 = 8275 +- 9 +- 33 pc"),
 ("GRAVITY Collaboration 2022, A&A 657, L12", "2022A&A...657L..12G", "10.1051/0004-6361/202142465", "https://arxiv.org/abs/2112.07478", "R0 = 8277 +- 9 pc"),
 ("Reid, M. J. et al. 2019, ApJ 885, 131", "2019ApJ...885..131R", "10.3847/1538-4357/ab4a11", "https://arxiv.org/abs/1910.03357", "R0 = 8.15, Theta0 = 236 +- 7, Theta0+Vsun = 247 +- 4"),
 ("Reid, M. J. & Brunthaler, A. 2004, ApJ 616, 872-884", "2004ApJ...616..872R", "10.1086/424960", "https://arxiv.org/abs/astro-ph/0408107", "Sgr A* proper motion; astropy galcen_coord"),
 ("Reid, M. J. & Brunthaler, A. 2020, ApJ 892, 39", "2020ApJ...892...39R", "10.3847/1538-4357/ab76cd", "https://arxiv.org/abs/2001.04386", "Sgr A* pm -6.411 +- 0.008, -0.219 +- 0.007 mas/yr"),
 ("Bland-Hawthorn, J. & Gerhard, O. 2016, ARA&A 54, 529-596", "2016ARA&A..54..529B", "10.1146/annurev-astro-081915-023441", "https://arxiv.org/abs/1602.07702", "R0 = 8.2 +- 0.1, Theta0 = 238 +- 15, Vg = 248 +- 3, z0 = 25 +- 5"),
 ("Drimmel, R. & Poggio, E. 2018, RNAAS 2, 210", "2018RNAAS...2..210D", "10.3847/2515-5172/aaef8b", "", "astropy v4.0 galcen_v_sun"),
 ("Bennett, M. & Bovy, J. 2019, MNRAS 482, 1417-1425", "2019MNRAS.482.1417B", "10.1093/mnras/sty2813", "https://arxiv.org/abs/1809.03507", "z_sun = 20.8 +- 0.3 pc"),
 ("Gaia Collaboration, Drimmel, R. et al. 2023, A&A 674, A37", "2023A&A...674A..37G", "10.1051/0004-6361/202243797", "https://arxiv.org/abs/2206.06207", "Gaia DR3 disc kinematics (context; no solar-motion value used)"),
 ("Yahil, A., Tammann, G. A., Sandage, A. 1977, ApJ 217, 903", "1977ApJ...217..903Y", "10.1086/155636", "", "LG correction -79 cos l cos b + 296 sin l cos b - 36 sin b"),
 ("Karachentsev, I. D. & Makarov, D. A. 1996, AJ 111, 794", "1996AJ....111..794K", "10.1086/117825", "", "Sun-LG 316 +- 5 km/s to (93, -4)"),
 ("Courteau, S. & van den Bergh, S. 1999, AJ 118, 337-345", "1999AJ....118..337C", "10.1086/300942", "https://arxiv.org/abs/astro-ph/9903298", "Sun-LG 306 +- 18 to (99, -4); RC2 300 sin l cos b; Delhaye LSR"),
 ("Diaz, J. D., Koposov, S. E., Irwin, M., Belokurov, V., Evans, W. 2014, MNRAS 443, 1688-1703", "2014MNRAS.443.1688D", "10.1093/mnras/stu1210", "https://arxiv.org/abs/1405.3662", "Sun-LG 299 +- 15 to (98.4, -5.9)"),
 ("van der Marel, R. P. et al. 2012, ApJ 753, 8 (M31 velocity vector II)", "2012ApJ...753....8V", "10.1088/0004-637X/753/1/8", "https://arxiv.org/abs/1205.6864", "M31 radial -109.3 +- 4.4, tangential 17.0 (<34.3) km/s"),
 ("Planck Collaboration 2020, Planck 2018 results I, A&A 641, A1", "2020A&A...641A...1P", "10.1051/0004-6361/201833880", "https://arxiv.org/abs/1807.06205", "solar dipole 369.82 +- 0.11 km/s to (264.021, 48.253); Table 3"),
 ("Planck Collaboration 2020, Planck 2018 results VI, A&A 641, A6", "2020A&A...641A...6P", "10.1051/0004-6361/201833910", "", "Planck18 cosmology (astropy)"),
 ("Fixsen, D. J. et al. 1996, ApJ 473, 576-587", "1996ApJ...473..576F", "10.1086/178173", "https://arxiv.org/abs/astro-ph/9605054", "COBE/FIRAS dipole 3.372 mK to (264.14, 48.26); NED 3K frame"),
 ("Kogut, A. et al. 1993, ApJ 419, 1", "1993ApJ...419....1K", "10.1086/173453", "https://arxiv.org/abs/astro-ph/9312056", "COBE DMR dipole; LG-CMB 627 +- 22 to (276, 30)"),
 ("Tully, R. B. et al. 2008, ApJ 676, 184-205", "2008ApJ...676..184T", "10.1086/527428", "https://arxiv.org/abs/0705.4139", "Local Void / Virgo / large-scale decomposition 259, 185, 455 -> 631 km/s"),
 ("Mould, J. R. et al. 2000, ApJ 529, 786-794", "2000ApJ...529..786M", "10.1086/308304", "https://arxiv.org/abs/astro-ph/9909260", "multi-attractor flow model eqs. 1-2, A1-A2, Table A1"),
 ("NASA/IPAC Extragalactic Database: velocity-correction help and object-result help", "", "", "https://ned.ipac.caltech.edu/help/velc_help.html ; https://ned.ipac.caltech.edu/help/objresult_help.html", "NED apex vectors and Mould-model corrections"),
 ("Tully, R. B. et al. 2023, ApJ 944, 94 (Cosmicflows-4)", "2023ApJ...944...94T", "10.3847/1538-4357/ac94d8", "https://arxiv.org/abs/2209.11238", "peculiar-velocity catalogue"),
 ("Davis, T. M. et al. 2011, ApJ 741, 67", "2011ApJ...741...67D", "10.1088/0004-637X/741/1/67", "https://arxiv.org/abs/1012.2912", "(1+z) = (1+zbar)(1+z_pec,sun), z_sun = -v.n/c"),
 ("Carr, A. et al. 2022, PASA 39, e046", "2022PASA...39...46C", "10.1017/pasa.2022.41", "https://arxiv.org/abs/2112.01471", "exact multiplicative heliocentric->CMB conversion"),
 ("Gaia Collaboration, Klioner, S. A. et al. 2021, A&A 649, A9", "2021A&A...649A...9G", "10.1051/0004-6361/202039734", "https://arxiv.org/abs/2012.02036", "SSB acceleration 2.32e-10 m/s^2"),
 ("Sandage, A. 1962, ApJ 136, 319", "1962ApJ...136..319S", "10.1086/147385", "", "redshift drift"),
 ("Loeb, A. 1998, ApJ 499, L111-L114", "1998ApJ...499L.111L", "10.1086/311375", "", "redshift drift"),
 ("astropy 7.2.0 source: coordinates/sky_coordinate.py (radial_velocity_correction), earth.py (gravitational_redshift), builtin_frames/lsr.py, builtin_frames/galactocentric.py", "", "", "https://github.com/astropy/astropy", "implementation details verified locally"),
 ("astropy issue #18550 'A possible LSRK documentation error' (docstring wording vs code)", "", "", "https://github.com/astropy/astropy/issues/18550", "LSR docstring ambiguity; code follows Gordon 1975"),
]

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def gal_xyz_from_lbV(l, b, V):
    l, b = np.radians(l), np.radians(b)
    return V * np.array([np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)])


def lbV_from_gal_xyz(xyz):
    x, y, z = xyz
    V = math.sqrt(x * x + y * y + z * z)
    return V, math.degrees(math.atan2(y, x)) % 360.0, math.degrees(math.asin(z / V))


def gal_dir_to_icrs(l, b):
    c = SkyCoord(l=l * u.deg, b=b * u.deg, frame="galactic").icrs
    return c.ra.deg, c.dec.deg


def proj_lb(l, b, V, l0, b0):
    """Projected correction V*cos(theta) toward target (l0, b0)  [same units as V]."""
    l, b, l0, b0 = map(np.radians, (l, b, l0, b0))
    return V * (np.cos(b) * np.cos(b0) * np.cos(l - l0) + np.sin(b) * np.sin(b0))


def fmt_ra(ra_deg):
    return SkyCoord(ra=ra_deg * u.deg, dec=0 * u.deg).ra.to_string(unit=u.hourangle, sep="hms", precision=2, pad=True)


def fmt_dec(dec_deg):
    return SkyCoord(ra=0 * u.deg, dec=dec_deg * u.deg).dec.to_string(sep="dms", precision=1, alwayssign=True, pad=True)


class DE:
    """Minimal DE440s access via jplephem (positions km, velocities km/day, TDB)."""

    def __init__(self, path):
        self.k = SPK.open(path)

    def _pv(self, a, b, t):
        p, v = self.k[a, b].compute_and_differentiate(t.tdb.jd1, t.tdb.jd2)
        return p.T * 1e3, v.T * 1e3 / 86400.0   # m, m/s ; shape (N,3) or (3,)

    def body(self, name, t):
        if name == "emb":
            return self._pv(0, 3, t)
        if name == "earth":
            p1, v1 = self._pv(0, 3, t); p2, v2 = self._pv(3, 399, t); return p1 + p2, v1 + v2
        if name == "moon":
            p1, v1 = self._pv(0, 3, t); p2, v2 = self._pv(3, 301, t); return p1 + p2, v1 + v2
        if name == "sun":
            return self._pv(0, 10, t)
        if name == "jupiter":
            return self._pv(0, 5, t)
        raise KeyError(name)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ephem", required=True, help="path to JPL DE kernel (de440s.bsp)")
    ap.add_argument("--out-md", required=True)
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--nmc", type=int, default=20000)
    args = ap.parse_args()

    R = {"meta": dict(generated_utc=_dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
                      astropy=astropy.__version__, numpy=np.__version__, ephemeris=args.ephem,
                      script="scripts/research/velocity_examples.py"),
         "literature": LITERATURE}
    de = DE(args.ephem)
    rng = np.random.default_rng(20260320)
    ephem_url = "file://" + args.ephem

    # ======================================================================
    # 0. velocity <-> wavelength scale (why 1e-3 Angstrom)
    # ======================================================================
    scale = {}
    for lam in (1216.0, 3000.0, 5000.0, 6563.0, 10000.0, 21000.0):
        scale[str(lam)] = dict(mA_per_ms=lam / C_MS * 1e3, ms_per_mA=1e-3 / lam * C_MS)
    R["wavelength_velocity_scale"] = scale

    # ======================================================================
    # 1. Earth rotation
    # ======================================================================
    om = LITERATURE["omega_earth"]["value"]; a = LITERATURE["a_wgs84"]["value"]
    keck = EarthLocation.from_geodetic(lon=-155.47833333 * u.deg, lat=19.82833333 * u.deg, height=4160.0 * u.m)  # astropy site registry 'Keck Observatory'
    xyz = np.array([keck.x.to_value(u.m), keck.y.to_value(u.m), keck.z.to_value(u.m)])
    rho_axis = math.hypot(xyz[0], xyz[1])
    R["earth_rotation"] = dict(v_equator_ms=om * a, keck_axis_distance_m=rho_axis, v_keck_ms=om * rho_axis,
                               accel_equator_ms2=om * om * a, dv_per_s_equator=om * om * a,
                               note="v = omega * a (WGS84); Keck: distance from rotation axis from WGS84 geodetic->ITRS")

    # ======================================================================
    # 2. Ephemeris statistics (DE440s via jplephem; independent of astropy ephemeris code)
    # ======================================================================
    t_year = Time("2026-01-01T00:00:00", scale="tdb") + np.arange(0, 365.25 * 24, 1.0) * u.hour
    pE, vE = de.body("earth", t_year); pS, vS = de.body("sun", t_year); pB, vB = de.body("emb", t_year)
    sE = np.linalg.norm(vE, axis=1) / 1e3
    sEh = np.linalg.norm(vE - vS, axis=1) / 1e3
    sBh = np.linalg.norm(vB - vS, axis=1) / 1e3
    wob = np.linalg.norm(vE - vB, axis=1)
    accE = np.linalg.norm(np.gradient(vE, 3600.0, axis=0), axis=1)
    R["earth_orbit_2026"] = dict(
        earth_bary_speed_kms=dict(min=sE.min(), max=sE.max(), mean=sE.mean()),
        earth_helio_speed_kms=dict(min=sEh.min(), max=sEh.max(), mean=sEh.mean()),
        emb_helio_speed_kms=dict(min=sBh.min(), max=sBh.max(), mean=sBh.mean()),
        earth_about_emb_ms=dict(min=wob.min(), max=wob.max(), mean=wob.mean()),
        earth_accel_ms2=dict(min=accE.min(), max=accE.max()),
        light_time_1au_s=AU_M / C_MS,
        dv_over_1au_lighttime_ms=accE.mean() * AU_M / C_MS)

    t_long = Time("1990-01-01T00:00:00", scale="tdb") + np.arange(0, 50 * 365.25, 1.0) * u.day
    pS2, vS2 = de.body("sun", t_long); pJ2, vJ2 = de.body("jupiter", t_long)
    sS = np.linalg.norm(vS2, axis=1); dS = np.linalg.norm(pS2, axis=1)
    vJh = np.linalg.norm(vJ2 - vS2, axis=1)
    mJ = (const.GM_jup / const.GM_sun).decompose().value
    R["sun_about_ssb_1990_2040"] = dict(
        speed_ms=dict(min=sS.min(), max=sS.max(), mean=sS.mean()),
        distance_Rsun=dict(min=dS.min() / const.R_sun.to_value(u.m), max=dS.max() / const.R_sun.to_value(u.m)),
        max_light_time_s=dS.max() / C_MS,
        jupiter_only_reflex_ms=dict(mean=(vJh * mJ / (1 + mJ)).mean(), min=(vJh * mJ / (1 + mJ)).min(), max=(vJh * mJ / (1 + mJ)).max()),
        note="Jupiter-only reflex = |v_Jup,helio| * GM_J/(GM_sun+GM_J) with astropy GM_jup (IAU 2015 nominal)")

    # ======================================================================
    # 3. Worked-example geometry
    # ======================================================================
    tobs = Time("2026-03-20T10:00:00", scale="utc")
    tgt = SkyCoord(ra=150.0 * u.deg, dec=20.0 * u.deg, frame="icrs")
    n_icrs = unit(tgt.cartesian.xyz.value)
    gal = tgt.galactic
    l0, b0 = gal.l.deg, gal.b.deg
    altaz = tgt.transform_to(AltAz(obstime=tobs, location=keck))
    lst = tobs.sidereal_time("apparent", longitude=keck.lon)
    app = tgt.transform_to(TETE(obstime=tobs, location=keck))   # apparent (true equator & equinox, topocentric)
    ha = (lst - app.ra).wrap_at(180 * u.deg)
    ecl = tgt.barycentrictrueecliptic
    R["example_geometry"] = dict(time_utc=tobs.isot, time_tdb=tobs.tdb.isot, jd_utc=tobs.jd, ut1_utc_s=float(tobs.delta_ut1_utc),
                                 site="Keck Observatory, Mauna Kea (astropy site registry): lon -155.47833 deg, lat +19.82833 deg, h 4160 m",
                                 ra_deg=150.0, dec_deg=20.0, l_deg=l0, b_deg=b0, ecl_lon_deg=ecl.lon.deg, ecl_lat_deg=ecl.lat.deg,
                                 alt_deg=altaz.alt.deg, az_deg=altaz.az.deg, airmass=float(altaz.secz),
                                 lst_apparent_h=lst.hour, hour_angle_h=ha.hour, ra_app_deg=app.ra.deg, dec_app_deg=app.dec.deg)

    # ======================================================================
    # 4. Barycentric correction: astropy vs independent WE14 eq.28
    # ======================================================================
    def astropy_rvc(t, kind="barycentric", eph=ephem_url):
        with solar_system_ephemeris.set(eph):
            return tgt.radial_velocity_correction(kind=kind, obstime=t, location=keck).to_value(u.m / u.s)

    def we14(t, n=n_icrs, include=("sun", "jupiter", "moon", "earth"), return_parts=False):
        """Independent zB (WE14 eq. 28 without stellar terms, Shapiro or zLT): DE440s via jplephem."""
        obs_p, obs_v = keck.get_gcrs_posvel(t)
        op = obs_p.xyz.to_value(u.m).T; ov = obs_v.xyz.to_value(u.m / u.s).T
        pE_, vE_ = de.body("earth", t)
        beta = (vE_ + ov) / C_MS
        gam = 1.0 / np.sqrt(1.0 - np.sum(beta * beta, axis=-1))
        pos_obs = pE_ + op
        GM = {"sun": const.GM_sun.value, "jupiter": const.GM_jup.value,
              "moon": (const.G * 7.34767309e22 * u.kg).to_value(u.m ** 3 / u.s ** 2), "earth": const.GM_earth.value}
        zgr = 0.0; parts = {}
        for bname in include:
            if bname == "earth":
                d = np.linalg.norm(op, axis=-1)
            else:
                pb, _ = de.body(bname, t); d = np.linalg.norm(pos_obs - pb, axis=-1)
            parts[bname] = -GM[bname] / (C_MS * d)          # m/s (astropy sign convention)
            zgr = zgr + parts[bname] / C_MS
        bn = np.sum(beta * n, axis=-1)
        zB = gam * (1.0 + bn) / (1.0 + zgr) - 1.0
        if return_parts:
            return zB * C_MS, dict(beta_dot_n_ms=bn * C_MS, gamma_minus1_ms=(gam - 1) * C_MS,
                                   grav_ms={k: float(np.atleast_1d(v)[0]) for k, v in parts.items()},
                                   v_obs_gcrs=ov, v_earth=vE_, pos_obs=pos_obs, op=op)
        return zB * C_MS

    vb_ast = float(astropy_rvc(tobs))
    vb_ast_builtin = float(astropy_rvc(tobs, eph="builtin"))
    vh_ast = float(astropy_rvc(tobs, kind="heliocentric"))
    vb_ind, parts = we14(tobs, return_parts=True)
    vb_ind = float(vb_ind)

    # decomposition (Galilean projections)
    pS_t, vS_t = de.body("sun", tobs); pB_t, vB_t = de.body("emb", tobs); pE_t, vE_t = de.body("earth", tobs)
    v_rot = float(np.dot(parts["v_obs_gcrs"], n_icrs))
    v_emb = float(np.dot(vB_t, n_icrs))
    v_wob = float(np.dot(vE_t - vB_t, n_icrs))
    v_sun = float(np.dot(vS_t, n_icrs))
    galilean = float(parts["beta_dot_n_ms"])
    # own rotation estimate: omega_nominal * (z_CIP x r_obs) in GCRS
    zc = SkyCoord(CartesianRepresentation(0, 0, 6.4e6, unit=u.m), frame="itrs", obstime=tobs).transform_to(GCRS(obstime=tobs))
    zhat = unit(zc.cartesian.xyz.to_value(u.m))
    v_rot_own = float(np.dot(om * np.cross(zhat, parts["op"]), n_icrs))
    v_rot_formula = float(-om * rho_axis * math.cos(app.dec.radian) * math.sin(ha.radian))
    # light-time / Roemer
    roemer_s = float(np.dot(parts["pos_obs"], n_icrs) / C_MS)
    dv_roemer = float(astropy_rvc(tobs + roemer_s * u.s) - vb_ast)
    # time sensitivity, exposure curvature (uniform flux, 1800 s)
    dv_dt = float((astropy_rvc(tobs + 30 * u.s) - astropy_rvc(tobs - 30 * u.s)) / 60.0)
    texp = tobs + np.linspace(-900, 900, 181) * u.s
    vexp = astropy_rvc(texp)
    curv_1800 = float(np.mean(vexp) - vb_ast)
    texp2 = tobs + np.linspace(-1800, 1800, 361) * u.s
    curv_3600 = float(np.mean(astropy_rvc(texp2)) - vb_ast)
    # Shapiro (WE14 eq. 27), x = observer - Sun, rho_hat = n
    x = parts["pos_obs"] - pS_t; xn = np.linalg.norm(x); xh = x / xn
    beta_E = (vE_t + parts["v_obs_gcrs"]) / C_MS
    zsd = -2 * const.GM_sun.value / C_MS ** 2 * np.dot(beta_E, n_icrs - np.dot(xh, n_icrs) * xh) / (xn * (1 + np.dot(xh, n_icrs)))
    elong = math.degrees(math.acos(np.dot(-xh, n_icrs)))
    R["example_barycentric"] = dict(
        astropy_bary_de440s_ms=vb_ast, astropy_bary_builtin_ms=vb_ast_builtin, astropy_helio_ms=vh_ast,
        independent_we14_de440s_ms=vb_ind, diff_astropy_minus_independent_ms=vb_ast - vb_ind,
        diff_builtin_minus_de440s_ms=vb_ast_builtin - vb_ast,
        bary_minus_helio_ms=vb_ast - vh_ast,
        decomposition_ms=dict(rotation_observer_gcrs=v_rot, rotation_own_omega_cross_r=v_rot_own,
                              rotation_simple_formula=v_rot_formula,
                              emb_orbital=v_emb, earth_about_emb=v_wob, galilean_total=galilean,
                              sun_about_ssb=v_sun,
                              relativistic_remainder=vb_ind - galilean,
                              gamma_minus1=float(parts["gamma_minus1_ms"]), grav_redshift_terms=parts["grav_ms"]),
        roemer_delay_s=roemer_s, change_if_evaluated_at_t_plus_roemer_ms=dv_roemer,
        dvcorr_dt_ms_per_s=dv_dt,
        exposure_mean_minus_midpoint_ms=dict(texp_1800s=curv_1800, texp_3600s=curv_3600),
        shapiro_we14_eq27_ms=float(zsd * C_MS), sun_elongation_deg=elong,
        cross_term_examples_ms={str(v): vb_ast * v * 1e3 / C_MS for v in (10.0, 50.0, 300.0)},
    )

    # ---- year-long comparison grid (every 6 h, 2026, at Keck, fixed target)
    tg = Time("2026-01-01T00:00:00", scale="utc") + np.arange(0, 365 * 4) * 6 * u.hour
    a_de = astropy_rvc(tg); a_bi = astropy_rvc(tg, eph="builtin"); ind = we14(tg)
    pEg, vEg = de.body("earth", tg)
    obs_p, obs_v = keck.get_gcrs_posvel(tg)
    gal_only = np.sum((vEg + obs_v.xyz.to_value(u.m / u.s).T) * n_icrs, axis=1)
    R["year_grid_2026"] = dict(
        n=len(tg), bary_corr_range_ms=dict(min=float(a_de.min()), max=float(a_de.max())),
        astropy_minus_independent_ms=dict(max_abs=float(np.abs(a_de - ind).max()), rms=float(np.sqrt(np.mean((a_de - ind) ** 2)))),
        builtin_minus_de440s_ms=dict(max_abs=float(np.abs(a_bi - a_de).max()), rms=float(np.sqrt(np.mean((a_bi - a_de) ** 2)))),
        full_minus_galilean_ms=dict(min=float((a_de - gal_only).min()), max=float((a_de - gal_only).max())),
    )

    # ======================================================================
    # 5. Frame apex table
    # ======================================================================
    frames = {}

    def add_uvw(key, uvw, sig=None, label=None, kind="uvw"):
        V, l, b = lbV_from_gal_xyz(uvw)
        ra, dec = gal_dir_to_icrs(l, b)
        corr = float(np.dot(uvw, unit(gal_xyz_from_lbV(l0, b0, 1.0))))
        d = dict(label=label or key, uvw=list(map(float, uvw)), V=V, l=l, b=b, ra=ra, dec=dec, ra_hms=fmt_ra(ra), dec_dms=fmt_dec(dec),
                 corr_kms=corr, src=LITERATURE.get(key, {}).get("src", ""))
        if sig is not None:
            s = np.array(sig, float)
            smp = np.asarray(uvw)[None, :] + rng.standard_normal((args.nmc, 3)) * s[None, :]
            d["corr_sigma_kms"] = float(np.std(smp @ unit(gal_xyz_from_lbV(l0, b0, 1.0))))
            Vs = np.linalg.norm(smp, axis=1); d["V_sigma"] = float(np.std(Vs))
        frames[key] = d

    def add_lbv(key, V, l, b, sV=None, sl=None, sb=None, label=None):
        uvw = gal_xyz_from_lbV(l, b, V)
        ra, dec = gal_dir_to_icrs(l, b)
        d = dict(label=label or key, uvw=list(map(float, uvw)), V=V, l=l, b=b, ra=ra, dec=dec, ra_hms=fmt_ra(ra), dec_dms=fmt_dec(dec),
                 corr_kms=float(proj_lb(l, b, V, l0, b0)), src=LITERATURE.get(key, {}).get("src", ""))
        if sV is not None:
            Vs = V + rng.standard_normal(args.nmc) * sV
            ls = l + rng.standard_normal(args.nmc) * (sl or 0.0)
            bs = b + rng.standard_normal(args.nmc) * (sb or 0.0)
            d["corr_sigma_kms"] = float(np.std(proj_lb(ls, bs, Vs, l0, b0)))
        frames[key] = d

    # LSRK from its B1900 definition
    lk = LITERATURE["LSRK_Gordon1975"]
    apex_k = SkyCoord(ra=lk["ra_b1900"] * u.deg, dec=lk["dec_b1900"] * u.deg, frame=FK4(equinox="B1900"))
    apex_k_icrs = apex_k.icrs; apex_k_gal = apex_k.galactic
    uvw_k = gal_xyz_from_lbV(apex_k_gal.l.deg, apex_k_gal.b.deg, lk["V"])
    add_uvw("LSRK_Gordon1975", uvw_k, label="LSRK (standard solar motion, radio convention)")
    v_icrs_k = lk["V"] * unit(apex_k_icrs.cartesian.xyz.value)
    frames["LSRK_Gordon1975"]["icrs_xyz_kms"] = v_icrs_k.tolist()
    frames["LSRK_Gordon1975"]["astropy_hardcoded_V_OFFSET_LSRK"] = [0.28999706839034606, -17.317264789717928, 10.00141199546947]
    add_uvw("LSRD_Delhaye1965", LITERATURE["LSRD_Delhaye1965"]["uvw"], label="LSRD (Delhaye 1965 peculiar solar motion)")
    sb = LITERATURE["LSR_Schoenrich2010"]
    sig_s = [math.hypot(0.5 * (a_ + b_), sy) for (a_, b_), sy in zip(sb["stat"], sb["sys"])]
    add_uvw("LSR_Schoenrich2010", sb["uvw"], sig=sig_s, label="LSR (Schoenrich, Binney & Dehnen 2010; stat+sys)")
    add_uvw("LSR_DehnenBinney1998", LITERATURE["LSR_DehnenBinney1998"]["uvw"], sig=LITERATURE["LSR_DehnenBinney1998"]["err"], label="LSR (Dehnen & Binney 1998)")
    add_uvw("LSR_BHG2016", LITERATURE["LSR_BHG2016"]["uvw"], sig=LITERATURE["LSR_BHG2016"]["err"], label="LSR (Bland-Hawthorn & Gerhard 2016 adopted)")
    add_uvw("LSR_Ding2019_GaiaDR2", LITERATURE["LSR_Ding2019_GaiaDR2"]["uvw"], sig=LITERATURE["LSR_Ding2019_GaiaDR2"]["err"], label="LSR (Ding+2019, Gaia DR2)")
    # Galactic standard of rest variants
    add_uvw("GSR_IAU1985_Delhaye", (9.0, 12.0 + 220.0, 7.0), label="GSR: IAU 1985 Theta0=220 + Delhaye (9,12,7)  [derived]")
    frames["GSR_IAU1985_Delhaye"]["src"] = "derived: Kerr & Lynden-Bell 1986 + Delhaye 1965 (as in Courteau & vdB 1999 / RC3)"
    add_lbv("NED_GSR_RC3", 232.3, 87.8, 1.7, label="GSR as used by NED (RC3)")
    add_uvw("astropy_v4.0", LITERATURE["astropy_v4.0"]["galcen_v_sun"], label="GSR: astropy Galactocentric v4.0 default (approx. Galactic axes)")
    add_uvw("astropy_pre-v4.0", LITERATURE["astropy_pre-v4.0"]["galcen_v_sun"], label="GSR: astropy pre-v4.0 (approx. Galactic axes)")
    r19 = LITERATURE["Reid2019"]
    add_uvw("GSR_Reid2019_composite", (11.1, r19["Theta0_plus_Vsun"], 7.25), sig=(1.2, 4.0, 0.6),
            label="GSR: Reid+2019 Theta0+Vsun=247+-4 with Schoenrich U,W  [derived composite]")
    frames["GSR_Reid2019_composite"]["src"] = "derived: Reid+2019 ApJ 885, 131 (V) + Schoenrich+2010 (U, W)"
    # From Sgr A* proper motion x GRAVITY R0 (derived)
    pm = LITERATURE["SgrA_pm_RB2020"]
    for key, R0pc, sR0 in (("GRAVITY2019", 8178.0, math.hypot(13, 22)), ("GRAVITY2021", 8275.0, math.hypot(9, 33))):
        Vy = K_PM * abs(pm["mu_l"]) * R0pc / 1e3
        Vz = K_PM * abs(pm["mu_b"]) * R0pc / 1e3
        sVy = Vy * math.hypot(pm["mu_l_err"] / pm["mu_l"], sR0 / R0pc)
        sVz = Vz * math.hypot(pm["mu_b_err"] / pm["mu_b"], sR0 / R0pc)
        k = f"GSR_SgrApm_{key}"
        add_uvw(k, (11.1, Vy, Vz), sig=(1.2, sVy, sVz), label=f"GSR: Sgr A* pm (Reid & Brunthaler 2020) x R0 ({key}) + U=11.1  [derived]")
        frames[k]["src"] = f"derived: V_phi=4.74047*|mu_l|*R0, V_z=4.74047*|mu_b|*R0 with mu from Reid & Brunthaler 2020, R0 from {key}; U from Schoenrich+2010"
        frames[k]["Vy_kms"] = Vy; frames[k]["Vz_kms"] = Vz; frames[k]["sVy"] = sVy
    # Local Group
    add_lbv("LG_RC2_300sinl", 300.0, 90.0, 0.0, label="LG: RC2/'IAU' 300 sin l cos b")
    add_uvw("LG_YTS1977", LITERATURE["LG_YTS1977"]["xyz"], label="LG: Yahil, Tammann & Sandage 1977")
    for k in ("LG_KM1996", "LG_CvdB1999", "LG_Diaz2014"):
        d = LITERATURE[k]
        add_lbv(k, d["V"], d["l"], d["b"], d["V_err"], d["l_err"], d["b_err"], label={"LG_KM1996": "LG: Karachentsev & Makarov 1996 (NED)",
                                                                                       "LG_CvdB1999": "LG: Courteau & van den Bergh 1999",
                                                                                       "LG_Diaz2014": "LG: Diaz+2014 (Planck 2018 I Table 3)"}[k])
    # CMB
    p = LITERATURE["CMB_Planck2018"]
    add_lbv("CMB_Planck2018", p["V"], p["l"], p["b"], p["V_err"], p["l_err"], p["b_err"], label="CMB: Planck 2018 (HFI) solar dipole")
    f = LITERATURE["CMB_Fixsen1996_NED"]
    add_lbv("CMB_Fixsen1996_NED", f["V"], f["l"], f["b"], 1.0, f["l_err95"] / 2, f["b_err95"] / 2, label="CMB: Fixsen+1996 (NED's 3K frame)")
    frames["CMB_Fixsen1996_NED"]["V_from_dT_over_T0_derived"] = C_KMS * f["T_mK"] * 1e-3 / f["T0"]
    kg = LITERATURE["CMB_Kogut1993_dipole"]
    frames["CMB_Fixsen1996_NED"]["note"] = "sigma_V=1 km/s and 95%->1sigma halving are assumptions for the MC only"
    for k in ("LG_CMB_Kogut1993", "LG_CMB_Planck2018", "GC_CMB_Planck2018", "Sun_LSR_Planck2018"):
        d = LITERATURE[k]
        add_lbv(k, d["V"], d["l"], d["b"], d["V_err"], d.get("l_err"), d.get("b_err"), label=k)
    # cross-check: Sun->CMB = Sun->LG + LG->CMB (vector sum, Planck Table 3 inputs)
    v_sun_lg = gal_xyz_from_lbV(98.4, -5.9, 299.0); v_lg_cmb = gal_xyz_from_lbV(271.9, 29.6, 620.0)
    R["vector_closure_check"] = dict(sun_lg_plus_lg_cmb=lbV_from_gal_xyz(v_sun_lg + v_lg_cmb),
                                     planck_sun_cmb=(369.82, 264.021, 48.253),
                                     note="Planck 2018 I Table 3 says LG-CMB was derived by non-relativistic subtraction; this re-adds them")
    # Planck RA/Dec check
    pr = gal_dir_to_icrs(p["l"], p["b"])
    R["planck_apex_radec_check"] = dict(computed=pr, published=p["radec"])
    R["frames"] = frames
    lsr_keys = ["LSRK_Gordon1975", "LSRD_Delhaye1965", "LSR_Schoenrich2010", "LSR_DehnenBinney1998", "LSR_BHG2016", "LSR_Ding2019_GaiaDR2"]
    R["lsr_definition_differences_kms"] = {f"{a_}-{b_}": float(np.linalg.norm(np.subtract(frames[a_]["uvw"], frames[b_]["uvw"])))
                                           for i, a_ in enumerate(lsr_keys) for b_ in lsr_keys[i + 1:]}
    R["kogut1993_v_from_dT_with_FIRAS_T0_kms"] = C_KMS * 3.365e-3 / 2.728
    R["sources"] = [dict(citation=c_, bibcode=b_, doi=d_, url=u_, used_for=f_) for c_, b_, d_, u_, f_ in SOURCES]

    # ---- astropy frame-transform cross-checks at the target (object at rest w.r.t. SSB)
    obj = SkyCoord(ra=150 * u.deg, dec=20 * u.deg, distance=1e6 * u.kpc, pm_ra_cosdec=0 * u.mas / u.yr, pm_dec=0 * u.mas / u.yr,
                   radial_velocity=0 * u.km / u.s, frame="icrs")
    xc = {}
    xc["LSRK"] = float(obj.transform_to(LSRK()).radial_velocity.to_value(u.km / u.s))
    xc["LSRD"] = float(obj.transform_to(LSRD()).radial_velocity.to_value(u.km / u.s))
    xc["LSR_default"] = float(obj.transform_to(LSR()).radial_velocity.to_value(u.km / u.s))
    for name in ("v4.0", "pre-v4.0"):
        with galactocentric_frame_defaults.set(name):
            gc = Galactocentric()
            o = obj.transform_to(gc)
            sun = SkyCoord(ra=150 * u.deg, dec=20 * u.deg, distance=1e-12 * u.kpc, frame="icrs").transform_to(gc)
            los = unit((o.cartesian.xyz - sun.cartesian.xyz).to_value(u.kpc))
            vgc = o.velocity.d_xyz.to_value(u.km / u.s)
            xc[f"GSR_astropy_{name}"] = float(np.dot(vgc, los))
            xc[f"galcen_v_sun_{name}"] = vgc.tolist()
    R["astropy_frame_crosscheck_kms"] = xc

    # ======================================================================
    # 6. Hypothetical galaxy: CMB frame exact vs additive, Mould+2000 flow model
    # ======================================================================
    cz_hel = 3000.0  # km/s, hypothetical (optical convention, cz)
    z_hel = cz_hel / C_KMS
    vproj = frames["CMB_Planck2018"]["corr_kms"]
    # observer moving with beta relative to CMB; source at rest in CMB frame:
    # 1+z_sun = 1/(gamma*(1+beta.n)) exact SR; first order z_sun = -v.n/c (Davis+2011 eq.5)
    beta = p["V"] / C_KMS; cos_t = vproj / p["V"]
    zsun_1st = -vproj / C_KMS
    zsun_exact = 1.0 / ((1.0 / math.sqrt(1 - beta ** 2)) * (1 + beta * cos_t)) - 1.0
    z_cmb_add = z_hel - zsun_1st
    z_cmb_mult = (1 + z_hel) / (1 + zsun_1st) - 1
    z_cmb_mult_exact = (1 + z_hel) / (1 + zsun_exact) - 1
    R["example_galaxy_cmb"] = dict(cz_hel_kms=cz_hel, z_hel=z_hel, sun_cmb_projection_kms=vproj,
                                   z_sun_first_order=zsun_1st, z_sun_exact_SR=zsun_exact,
                                   z_cmb_additive=z_cmb_add, z_cmb_multiplicative=z_cmb_mult, z_cmb_mult_exactSR=z_cmb_mult_exact,
                                   cz_cmb_additive_kms=z_cmb_add * C_KMS, cz_cmb_multiplicative_kms=z_cmb_mult * C_KMS,
                                   difference_kms=(z_cmb_mult - z_cmb_add) * C_KMS,
                                   at_z1_difference=dict(z_hel=1.0, dz=((2.0) / (1 + zsun_1st) - 1) - (1.0 - zsun_1st)))

    # Mould et al. 2000 model (NED implementation: V_LG via Karachentsev & Makarov; + infall terms)
    md = LITERATURE["Mould2000"]
    v_lg_km = cz_hel + frames["LG_KM1996"]["corr_kms"]
    v_lg_yts = cz_hel + frames["LG_YTS1977"]["corr_kms"]
    mould = {}
    for which, Vo in (("NED_KM1996", v_lg_km), ("Mould_YTS1977", v_lg_yts)):
        tot = Vo; terms = {}
        for name, at in md["attractors"].items():
            ac = SkyCoord(at["ra_b1950"], at["dec_b1950"], frame=FK4(equinox="B1950")).icrs
            cth = float(np.dot(unit(ac.cartesian.xyz.value), n_icrs))
            Va = at["VLG"]; Vf = at["Vfid"]; g = md["gamma"]
            roa = math.sqrt(Vo ** 2 + Va ** 2 - 2 * Vo * Va * cth)
            vin = Vf * cth + Vf * ((Vo - Va * cth) / roa) * (roa / Va) ** (1 - g)
            sep = math.degrees(math.acos(cth))
            in_core = (sep < at["radius"]) and (at["vrange"][0] <= cz_hel <= at["vrange"][1])
            terms[name] = dict(theta_deg=sep, r_oa_kms=roa, V_infall_kms=vin, in_core=in_core,
                               icrs=(ac.ra.deg, ac.dec.deg))
            tot += vin
        mould[which] = dict(V_LG_kms=Vo, terms=terms, V_cosmic_kms=tot)
    R["example_galaxy_mould2000"] = mould

    # ======================================================================
    # 7. Secular / cosmological drifts
    # ======================================================================
    sa = LITERATURE["SSB_accel_GaiaEDR3"]
    a_dir = SkyCoord(ra=sa["ra"] * u.deg, dec=sa["dec"] * u.deg).cartesian.xyz.value
    yr = 365.25 * 86400.0
    R["ssb_acceleration"] = dict(a_ms2=sa["a"], drift_mm_s_per_yr_max=sa["a"] * yr * 1e3,
                                 drift_mm_s_per_yr_target=float(-sa["a"] * yr * 1e3 * np.dot(a_dir, n_icrs)),
                                 note="d(v_bary)/dt = -a.n for a distant source at rest in the extragalactic frame (derived)")
    from astropy.cosmology import Planck18
    H0 = Planck18.H0.to_value(1 / u.s)
    sl = {}
    for zz in (0.5, 1.0, 2.0, 3.0, 4.0, 5.0):
        zdot = (1 + zz) * H0 - Planck18.H(zz).to_value(1 / u.s)   # 1/s
        vdot = C_MS * zdot / (1 + zz) * yr * 100.0               # cm/s per yr
        sl[str(zz)] = dict(zdot_per_yr=zdot * yr, vdot_cm_s_per_yr=vdot, dv_cm_s_in_20yr=vdot * 20)
    R["sandage_loeb_planck18"] = sl
    # perspective (secular) acceleration, hypothetical star mu=1"/yr at 10 pc
    mu_as = 1.0; dpc = 10.0
    vt = K_PM * mu_as * 1e3 * dpc / 1e3  # km/s
    R["perspective_accel_hypothetical"] = dict(mu_arcsec_yr=mu_as, d_pc=dpc, v_t_kms=vt,
                                               dvr_dt_m_s_per_yr=(vt * 1e3) ** 2 / (dpc * 3.0856775814913673e16) * yr)
    # solar gravitational redshift recomputed from IAU 2015 nominal values
    R["solar_grav_redshift_derived"] = dict(
        photosphere_ms=(const.GM_sun / (const.R_sun * const.c)).to_value(u.m / u.s),
        at_1au_ms=(const.GM_sun / (const.R_sun * const.c) - const.GM_sun / (const.au * const.c)).to_value(u.m / u.s),
        note="astropy IAU 2015 nominal GM_sun and R_sun; LD2003 published 636.486+-0.024 (R=6.95508e8 m)")

    # ======================================================================
    # write outputs
    # ======================================================================
    def default(o):
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
        return str(o)

    with open(args.out_json, "w") as fh:
        json.dump(R, fh, indent=1, default=default)
    write_md(R, args.out_md)
    print("wrote", args.out_json, args.out_md)


# ---------------------------------------------------------------------------
# Markdown writer
# ---------------------------------------------------------------------------

def write_md(R, path):
    E = R["example_barycentric"]; G = R["example_geometry"]; F = R["frames"]; Y = R["year_grid_2026"]
    O = R["earth_orbit_2026"]; S = R["sun_about_ssb_1990_2040"]; ER = R["earth_rotation"]; XC = R["astropy_frame_crosscheck_kms"]
    L = []
    w = L.append
    w("# Velocity reference frames and corrections for precise spectroscopy")
    w("")
    w(f"Generated {R['meta']['generated_utc']} by `scripts/research/velocity_examples.py` (astropy {R['meta']['astropy']}, ephemeris {R['meta']['ephemeris']}).")
    w("All literature inputs are typed from the cited papers (see Sources); every other number is **derived** by the script and is labelled so. Machine-readable values: `data/notes/velocity_frames.json`.")
    w("")
    w("## 0. Velocity <-> wavelength scale (why 1e-3 Angstrom)")
    w("")
    w("Doppler: dlambda/lambda = v/c (first order). Derived from c = 299 792 458 m/s (exact):")
    w("")
    w("| lambda (A) | 1 m/s = ... mA | 1 mA = ... m/s |")
    w("|---|---|---|")
    for lam, d in R["wavelength_velocity_scale"].items():
        w(f"| {float(lam):.0f} | {d['mA_per_ms']:.5f} | {d['ms_per_mA']:.2f} |")
    w("")
    k5 = R["wavelength_velocity_scale"]["5000.0"]["mA_per_ms"]
    rel = E["decomposition_ms"]["relativistic_remainder"]
    dl = R["lsr_definition_differences_kms"]
    w(f"So 1e-3 A (1 mA) corresponds to {R['wavelength_velocity_scale']['21000.0']['ms_per_mA']:.0f}-{R['wavelength_velocity_scale']['1216.0']['ms_per_mA']:.0f} m/s between 2.1 um and Ly-alpha. At 5000 A (derived): "
      f"Earth rotation ({ER['v_equator_ms']:.1f} m/s, equator) = {ER['v_equator_ms']*k5:.2f} mA; Earth orbit ({O['earth_bary_speed_kms']['mean']:.2f} km/s) = {O['earth_bary_speed_kms']['mean']*1e3*k5/1e3:.3f} A; "
      f"Sun about the SSB (<= {S['speed_ms']['max']:.1f} m/s) = {S['speed_ms']['max']*k5:.2f} mA; Earth about the EMB (~{O['earth_about_emb_ms']['mean']:.1f} m/s) = {O['earth_about_emb_ms']['mean']*k5:.2f} mA; "
      f"relativistic terms of the barycentric correction ({rel:.2f} m/s here) = {rel*k5:.3f} mA; 1 m/s = {k5*1e-3:.2e} A; 1 cm/s = {k5*1e-5:.2e} A. "
      f"Different LSR conventions differ by {min(dl.values()):.1f}-{max(dl.values()):.1f} km/s in vector terms (= {min(dl.values())*1e3*k5:.0f}-{max(dl.values())*1e3*k5:.0f} mA). "
      "Hence at 1 mA precision Earth rotation and every frame-definition choice are resolved (so the topocentric/barycentric/LSR bookkeeping is mandatory), "
      "while the Sun's SSB reflex, the Earth-EMB wobble and the relativistic terms stay below one unit of the last digit; m/s-class (exoplanet) work needs ~1e-5 A wavelengths and the full relativistic correction.")
    w("")
    w("## 1. Observer (topocentric -> barycentric) terms")
    w("")
    w(f"* **Earth rotation.** omega = {LITERATURE['omega_earth']['value']:.6e} rad/s (IERS Conventions 2010, TN36 Table 1.2, GRS80/WGS84 nominal mean angular velocity), "
      f"a = 6 378 137 m -> v_eq = omega a = **{ER['v_equator_ms']:.2f} m/s** (derived). At Keck/Mauna Kea (distance from the rotation axis {ER['keck_axis_distance_m']/1e3:.1f} km) v = {ER['v_keck_ms']:.2f} m/s.")
    w("  Line-of-sight projection: the observer's velocity component *toward* a target at hour angle H, declination dec is v_toward = -omega rho cos(phi') cos(dec) sin(H) "
      "(rho cos phi' = distance from the axis; H > 0 west). The correction to *add* to the observed RV is +v_toward, i.e. the observed RV contains "
      "+omega rho cos phi' cos dec sin H (receding for setting targets).")
    w(f"  Rate of change: up to omega^2 a = {ER['accel_equator_ms2']*100:.2f} cm/s per second (derived; WE14 quote 3.4 cm/s/s), so the exposure time must be known to <~0.25 s for 1 cm/s (WE14 sec 6.2.2).")
    w(f"* **Earth orbit** (DE440s, 2026, derived): Earth-centre barycentric speed {O['earth_bary_speed_kms']['min']:.3f}-{O['earth_bary_speed_kms']['max']:.3f} km/s (mean {O['earth_bary_speed_kms']['mean']:.3f}); "
      f"EMB heliocentric {O['emb_helio_speed_kms']['min']:.3f}-{O['emb_helio_speed_kms']['max']:.3f} km/s. NASA Earth Fact Sheet: mean 29.78, max 30.29, min 29.29 km/s.")
    w(f"* **Earth about the Earth-Moon barycentre** (DE440s, 2026, derived): {O['earth_about_emb_ms']['min']:.2f}-{O['earth_about_emb_ms']['max']:.2f} m/s (mean {O['earth_about_emb_ms']['mean']:.2f}).")
    w(f"* **Sun about the SSB** (DE440s 1990-2040, derived): {S['speed_ms']['min']:.2f}-{S['speed_ms']['max']:.2f} m/s; Jupiter-only reflex {S['jupiter_only_reflex_ms']['mean']:.2f} m/s (range {S['jupiter_only_reflex_ms']['min']:.2f}-{S['jupiter_only_reflex_ms']['max']:.2f}). "
      f"Sun-SSB distance {S['distance_Rsun']['min']:.2f}-{S['distance_Rsun']['max']:.2f} R_sun (<= {S['max_light_time_s']:.2f} s light time). "
      "This is the heliocentric-minus-barycentric velocity difference (plus relativistic terms, which astropy applies only to the barycentric kind).")
    w(f"* **Light time.** 1 au / c = {O['light_time_1au_s']:.3f} s (derived). The Earth's acceleration ({O['earth_accel_ms2']['min']*1e3:.3f}-{O['earth_accel_ms2']['max']*1e3:.3f} mm/s^2) changes its velocity by "
      f"~{O['dv_over_1au_lighttime_ms']:.2f} m/s in that time, so the barycentric *velocity* must be evaluated at the observation epoch (UTC->TT/TDB), not at BJD; "
      "BJD (Roemer-corrected) is the time tag for the RV time series (WE14 sec 6.2.2; Eastman et al. 2010).")
    w("* **Relativistic barycentric correction** (Wright & Eastman 2014, PASP 126, 838, eq. 28):")
    w("")
    w("  `z_B = gamma_obs * (1 + beta_obs . rho_hat) / (1 + z_GR) * (1 + beta_s . r0_hat)/(1 + beta_s . rho_hat) - 1 - z_LT - z_SD`, "
      "with `z_GR = -sum_i G M_i / (c^2 d_i)` (WE14 eq. 2; Sun and Earth matter at 1 cm/s; ~3 m/s), and the correction applied **multiplicatively**: "
      "`(1 + z_true) = (1 + z_meas)(1 + z_B)` (WE14 eq. 10), i.e. `v_true ~= v_meas + c z_B + v_meas z_B` (cross term ~ c beta^2 ~ 3 m/s, needed at 1 m/s).")
    w("  WE14 budget: Galilean beta.n adequate at 10 m/s; multiplicative cross term at 1 m/s; precession matters at 1 m/s; parallax-proper-motion term at 10 cm/s for nearby stars; "
      "GR (~3 m/s, varying ~3 cm/s with orbital eccentricity) and time dilation (variation ~3 cm/s) at 1 cm/s; nutation 3 cm/s; UT1-UTC 2 cm/s; star coordinates to 69 mas; "
      "Shapiro delay <= 1.5 mm/s for tau Ceti; their code agrees with TEMPO2 to 0.037 cm/s RMS (0.24 cm/s peak-to-peak).")
    w("* **IAU 2000 Resolution C1** (text in Lindegren & Dravins 2003, A&A 401, 1185): report the *barycentric radial-velocity measure* c z_B, with z_B = (lambda - lambda0)/lambda0 "
      "for a hypothetical observer at zero gravitational potential at rest at the SSB, after correcting for solar-system gravitational effects and the observer's motion; "
      "epoch = barycentric time of light arrival. It is *not* a physical velocity. (LD03's z_B is the source's barycentric redshift; WE14's z_B is the correction factor - different quantities with the same symbol.)")
    w(f"* **Solar gravitational redshift** (relevant to solar Fraunhofer lines): GM/(R c) = 636.486 +- 0.024 m/s at the photosphere, 633.5 m/s for light received at 1 au (LD03 sec 3.1); "
      f"photospheric lines are blueshifted ~400 m/s by convection (LD03 sec 3.2). Recomputed with IAU 2015 nominal GM_sun, R_sun: {R['solar_grav_redshift_derived']['photosphere_ms']:.2f} m/s (derived).")
    w("* **Exposure midpoint.** Use the photon-weighted (exposure-meter) midpoint; geometric midpoints can differ by >250 ms (WE14). Even with the photon-weighted midpoint, "
      "the non-linearity of v_bary(t) gives a typical 10 cm/s, sometimes >1 m/s, second-order error (Tronsgaard et al. 2019, MNRAS 489, 2395); chromatic (airmass-dependent) "
      "flux weighting adds >= several cm/s (Blackman et al. 2017, ApJ 837, 18) and up to 1 m/s in EXPRES data (Blackman, Ong & Fischer 2019, AJ 158, 40).")
    w("")
    w("### astropy `SkyCoord.radial_velocity_correction` (verified in astropy 7.2 source)")
    w("")
    w("* Returns the velocity to **add** to an observed RV (`kind='barycentric'` default, or `'heliocentric'`). Barycentric: "
      "`zb = gamma_obs (1 + beta_obs . n)/(1 + gr/c)`, `vcorr = (zb - 1) c`, i.e. WE14 eq. 28 'neglecting Shapiro delay and effects of the star's own motion' (source comment); "
      "`beta_obs` = (SSB->geocentre velocity from the selected ephemeris + observer GCRS velocity)/c; `gr` = `EarthLocation.gravitational_redshift` = sum of -GM/(c d) for Sun, Jupiter, Moon and Earth (Earth at the site's geocentric radius).")
    w("* If the SkyCoord carries full space motion (distance, pm, rv) it multiplies by (1 + beta_s . r0)/(1 + beta_s . n) (WE14 parallax/proper-motion term) and, when obstime is set, applies space motion; distances enable a parallax-corrected direction.")
    w("* Heliocentric kind = plain Galilean projection (no gamma, no GR) - astropy says not to use it below m/s.")
    w("* The docstring states the algorithm is 'sufficient to perform corrections at the mm/s level' and that the result is in the optical convention; apply as `rv + vcorr + rv*vcorr/c`.")
    w("* Default ephemeris is 'builtin' (ERFA epv00); set `solar_system_ephemeris` to a JPL kernel for best accuracy. Omitted vs WE14: Shapiro delay z_SD, light-travel term z_LT, no exposure averaging.")
    w(f"* Year-2026 test (Keck, target below, every 6 h, N={Y['n']}): astropy(DE440s) minus an independent WE14 eq.-28 implementation driven by jplephem/DE440s: max |diff| = {Y['astropy_minus_independent_ms']['max_abs']*1e3:.4f} mm/s "
      f"(RMS {Y['astropy_minus_independent_ms']['rms']*1e3:.4f} mm/s); astropy 'builtin' minus DE440s: max {Y['builtin_minus_de440s_ms']['max_abs']*100:.3f} cm/s (RMS {Y['builtin_minus_de440s_ms']['rms']*100:.3f} cm/s); "
      f"full minus Galilean: {Y['full_minus_galilean_ms']['min']:.3f} to {Y['full_minus_galilean_ms']['max']:.3f} m/s. (Both implementations share astropy's observer GCRS velocity/Earth orientation; this checks formula + ephemeris handling, not Earth orientation.)")
    w("* astropy LSRK/LSRD/LSR/Galactocentric velocity changes are Galilean vector offsets (AffineTransform), not relativistic compositions.")
    w("")
    w("## 2-5. Frame apices (Sun's velocity relative to each frame)")
    w("")
    w(f"Correction to the barycentric (heliocentric) RV toward target (l, b): **v_frame = v_bary + V [cos b cos b_a cos(l - l_a) + sin b sin b_a]** (first order; Galilean). "
      f"Apex RA/Dec (ICRS ~ J2000) and (U,V,W) are **derived** with astropy from the published (V, l, b) or (U, V, W). Last two columns: projected correction for the worked-example target "
      f"(l = {G['l_deg']:.4f}, b = {G['b_deg']:.4f}) and its 1-sigma from Monte-Carlo of the published errors.")
    w("")
    w("| Frame | V (km/s) | l (deg) | b (deg) | RA (ICRS) | Dec (ICRS) | (U, V, W) km/s | corr @ target (km/s) | sigma | Source |")
    w("|---|---|---|---|---|---|---|---|---|---|")
    for k, d in F.items():
        uvw = ", ".join(f"{x:.2f}" for x in d["uvw"])
        sig = f"{d['corr_sigma_kms']:.2f}" if "corr_sigma_kms" in d else "-"
        w(f"| {d['label']} | {d['V']:.2f} | {d['l']:.3f} | {d['b']:.3f} | {d['ra_hms']} ({d['ra']:.3f}) | {d['dec_dms']} ({d['dec']:.3f}) | {uvw} | {d['corr_kms']:+.3f} | {sig} | {d['src']} |")
    w("")
    k = F["LSRK_Gordon1975"]
    w(f"* LSRK: 20 km/s toward RA 18h, Dec +30 (B1900) -> ICRS apex {k['ra_hms']} {k['dec_dms']} (derived via astropy FK4(B1900)->ICRS); ICRS vector {', '.join(f'{x:.5f}' for x in k['icrs_xyz_kms'])} km/s, "
      f"identical to astropy's hard-coded V_OFFSET_LSRK; Galactic (U,V,W) = ({', '.join(f'{x:.2f}' for x in k['uvw'])}) km/s (Reid et al. 2009, ApJ 700, 137, Table 1 note: 10.3, 15.3, 7.7 km/s).")
    w(f"* NED 'GSR' (232.3 km/s to l=87.8, b=+1.7, RC3) equals IAU-1985 Theta0 = 220 km/s plus Delhaye (9,12,7): derived apex {F['GSR_IAU1985_Delhaye']['V']:.2f} km/s to ({F['GSR_IAU1985_Delhaye']['l']:.2f}, {F['GSR_IAU1985_Delhaye']['b']:.2f}).")
    w(f"* YTS77 formula V_LG = V_H - 79 cos l cos b + 296 sin l cos b - 36 sin b is a {F['LG_YTS1977']['V']:.1f} km/s apex toward (l, b) = ({F['LG_YTS1977']['l']:.1f}, {F['LG_YTS1977']['b']:.1f}) (derived).")
    for gk in ("GSR_SgrApm_GRAVITY2019", "GSR_SgrApm_GRAVITY2021"):
        d = F[gk]
        w(f"* {d['label']}: V_phi = {d['Vy_kms']:.2f} +- {d['sVy']:.2f} km/s, V_z = {d['Vz_kms']:.2f} km/s (derived; 4.74047 km/s per mas/yr per kpc).")
    w(f"* Vector closure: Sun->LG (Diaz+2014) + LG->CMB (Planck Table 3) = {R['vector_closure_check']['sun_lg_plus_lg_cmb'][0]:.1f} km/s toward ({R['vector_closure_check']['sun_lg_plus_lg_cmb'][1]:.2f}, {R['vector_closure_check']['sun_lg_plus_lg_cmb'][2]:.2f}) vs Planck Sun->CMB 369.82 toward (264.021, 48.253).")
    w(f"* Planck apex RA/Dec computed {R['planck_apex_radec_check']['computed'][0]:.3f}, {R['planck_apex_radec_check']['computed'][1]:.3f} vs published 167.942, -6.944 (J2000).")
    w(f"* Fixsen 1996: c dT/T0 = {F['CMB_Fixsen1996_NED']['V_from_dT_over_T0_derived']:.2f} km/s from 3.372 mK / 2.728 K (derived); NED adopts 371.0 km/s.")
    w("* astropy Galactocentric defaults (astropy 7.2: 'latest' = 'v4.0'): galcen RA, Dec = 266.4051, -28.936175 (Reid & Brunthaler 2004), R0 = 8.122 kpc (GRAVITY 2018), "
      "galcen_v_sun = (12.9, 245.6, 7.78) km/s (Drimmel & Poggio 2018), z_sun = 20.8 pc (Bennett & Bovy 2019), roll = 0. 'pre-v4.0': 8.3 kpc, (11.1, 232.24, 7.25), 27 pc.")
    w(f"* Galactic constants: IAU 1985 R0 = 8.5 kpc, Theta0 = 220 km/s (Kerr & Lynden-Bell 1986); GRAVITY R0 = 8178 +- 13 +- 22 pc (2019), 8275 +- 9 +- 33 pc (2021), 8277 +- 9 (+-33 sys) pc (2022); "
      "Reid+2019 R0 = 8.15 +- 0.15 kpc, Theta0 = 236 +- 7, Theta0 + V_sun = 247 +- 4 km/s, Omega = 30.32 +- 0.27 km/s/kpc; Bland-Hawthorn & Gerhard 2016: R0 = 8.2 +- 0.1 kpc, Theta0 = 238 +- 15, V_g,sun = 248 +- 3 km/s, z0 = 25 +- 5 pc; "
      "Sgr A* apparent motion -6.411 +- 0.008 mas/yr along the plane, -0.219 +- 0.007 mas/yr toward the NGP (Reid & Brunthaler 2020).")
    w("")
    w("astropy frame-transform cross-check at the target (object at rest w.r.t. SSB, 1 Gpc): "
      + ", ".join(f"{k} {v:+.4f}" for k, v in XC.items() if not k.startswith("galcen")) + " km/s (compare the projection column).")
    w("")
    w("## 6. Worked example: RA = 10h, Dec = +20 deg, Keck/Mauna Kea, 2026-03-20T10:00:00 UTC")
    w("")
    w(f"Geometry (derived): l = {G['l_deg']:.4f}, b = {G['b_deg']:.4f}; ecliptic lon/lat {G['ecl_lon_deg']:.3f}/{G['ecl_lat_deg']:.3f}; altitude {G['alt_deg']:.2f} deg, azimuth {G['az_deg']:.2f}, airmass {G['airmass']:.3f}; "
      f"apparent LST {G['lst_apparent_h']:.4f} h, apparent (TETE) RA/Dec {G['ra_app_deg']:.4f}/{G['dec_app_deg']:.4f}, hour angle {G['hour_angle_h']:+.4f} h; TDB {G['time_tdb']}; UT1-UTC {G['ut1_utc_s']:.4f} s (IERS-A).")
    w("")
    D = E["decomposition_ms"]
    w("| term | m/s |")
    w("|---|---|")
    w(f"| Barycentric correction, astropy + DE440s (to add) | **{E['astropy_bary_de440s_ms']:+.4f}** |")
    w(f"| same, independent WE14 eq.28 (jplephem DE440s) | {E['independent_we14_de440s_ms']:+.4f} (diff {E['diff_astropy_minus_independent_ms']*1e3:+.4f} mm/s) |")
    w(f"| same, astropy builtin ephemeris (ERFA epv00) | {E['astropy_bary_builtin_ms']:+.4f} (diff {E['diff_builtin_minus_de440s_ms']*100:+.3f} cm/s) |")
    w(f"| Heliocentric correction, astropy (Galilean) | {E['astropy_helio_ms']:+.4f} |")
    w(f"| barycentric - heliocentric | {E['bary_minus_helio_ms']:+.4f} |")
    w(f"| -- Earth rotation (observer GCRS velocity . n) | {D['rotation_observer_gcrs']:+.4f} |")
    w(f"| -- rotation, omega_nominal (z_CIP x r) . n | {D['rotation_own_omega_cross_r']:+.4f} |")
    w(f"| -- rotation, -omega rho cos(dec) sin(H) (apparent dec, H) | {D['rotation_simple_formula']:+.4f} |")
    w(f"| -- EMB barycentric orbital motion . n | {D['emb_orbital']:+.4f} |")
    w(f"| -- Earth about EMB . n | {D['earth_about_emb']:+.4f} |")
    w(f"| -- Galilean total beta.n c | {D['galilean_total']:+.4f} |")
    w(f"| -- relativistic remainder (gamma, GR, cross) | {D['relativistic_remainder']:+.4f} |")
    w(f"| ---- (gamma - 1) c | {D['gamma_minus1']:+.4f} |")
    w(f"| ---- GR terms -GM/(c d): " + ", ".join(f"{k} {v:+.5f}" for k, v in D['grav_redshift_terms'].items()) + " | |")
    w(f"| Sun about SSB . n (helio vs bary) | {D['sun_about_ssb']:+.4f} |")
    w(f"| Shapiro term (WE14 eq. 27; Sun elongation {E['sun_elongation_deg']:.1f} deg) | {E['shapiro_we14_eq27_ms']*1e3:+.4f} mm/s |")
    w(f"| Roemer delay r_obs.n/c | {E['roemer_delay_s']:+.3f} s |")
    w(f"| correction change if evaluated at t + Roemer delay (wrong epoch; dominated by Earth rotation) | {E['change_if_evaluated_at_t_plus_roemer_ms']:+.4f} |")
    w(f"| d(vcorr)/dt | {E['dvcorr_dt_ms_per_s']*100:+.4f} cm/s per s |")
    w(f"| exposure mean - midpoint, 1800 s uniform flux | {E['exposure_mean_minus_midpoint_ms']['texp_1800s']*100:+.4f} cm/s |")
    w(f"| exposure mean - midpoint, 3600 s uniform flux | {E['exposure_mean_minus_midpoint_ms']['texp_3600s']*100:+.4f} cm/s |")
    w(f"| cross term v_meas z_B for v_meas = 10 / 50 / 300 km/s | " + " / ".join(f"{v:+.3f}" for v in E['cross_term_examples_ms'].values()) + " |")
    w("")
    w(f"Over 2026 the barycentric correction for this target from Keck spans {Y['bary_corr_range_ms']['min']/1e3:+.3f} to {Y['bary_corr_range_ms']['max']/1e3:+.3f} km/s (derived).")
    w("")
    w("Frame corrections for the same target (add to the barycentric RV; first-order Galilean; km/s):")
    w("")
    w("| frame | correction (km/s) | 1-sigma |")
    w("|---|---|---|")
    for k, d in F.items():
        sig = f"{d['corr_sigma_kms']:.2f}" if "corr_sigma_kms" in d else "(conventional)"
        w(f"| {d['label']} | {d['corr_kms']:+.3f} | {sig} |")
    w("")
    vb = E['astropy_bary_de440s_ms'] / 1e3
    w(f"Total for an observed topocentric RV v_obs (km/s): v_LSRK = v_obs + {vb:+.4f} + {F['LSRK_Gordon1975']['corr_kms']:+.3f} (+ v_obs*vcorr/c); "
      f"v_GSR(astropy v4.0) = v_obs + {vb:+.4f} + {F['astropy_v4.0']['corr_kms']:+.3f}; v_LG(KM96) = v_obs + {vb:+.4f} + {F['LG_KM1996']['corr_kms']:+.3f}; "
      f"v_CMB(Planck) = v_obs + {vb:+.4f} + {F['CMB_Planck2018']['corr_kms']:+.3f} (first order).")
    w("")
    g = R["example_galaxy_cmb"]
    w(f"**Hypothetical galaxy, cz_hel = {g['cz_hel_kms']:.0f} km/s** in this direction (illustration): Sun-CMB projection {g['sun_cmb_projection_kms']:+.3f} km/s. "
      f"Additive cz_CMB = {g['cz_cmb_additive_kms']:.3f} km/s; exact multiplicative (1+z_hel) = (1+z_CMB)(1+z_sun) (Davis+2011 eq. 6; Carr+2022 eq. 2-3) gives {g['cz_cmb_multiplicative_kms']:.3f} km/s "
      f"(difference {g['difference_kms']:+.3f} km/s); with exact SR z_sun = 1/[gamma(1+beta cos theta)] - 1: cz_CMB = {g['z_cmb_mult_exactSR']*C_KMS:.3f} km/s. At z_hel = 1 the additive error is {g['at_z1_difference']['dz']:.2e} in z.")
    w("")
    for which, m in R["example_galaxy_mould2000"].items():
        w(f"Mould+2000 multi-attractor model ({which} LG step): V_LG = {m['V_LG_kms']:.1f}; " +
          "; ".join(f"{n}: theta {t['theta_deg']:.1f} deg, V_infall {t['V_infall_kms']:+.1f}" for n, t in m['terms'].items()) +
          f" -> V_cosmic = {m['V_cosmic_kms']:.1f} km/s (derived; gamma = 2; NED sign convention 'all positive').")
    w("")
    w("## 7. Secular and cosmology-scale terms")
    w("")
    sa = R["ssb_acceleration"]
    w(f"* SSB acceleration (Gaia EDR3; Klioner+2021): (2.32 +- 0.16)e-10 m/s^2 toward RA 269.1, Dec -31.6 -> up to {sa['drift_mm_s_per_yr_max']:.2f} mm/s per year in apparent RVs of distant sources; "
      f"d(v_bary)/dt = -a.n on the target {sa['drift_mm_s_per_yr_target']:+.2f} mm/s/yr (derived). The 'barycentric' frame is not inertial at this level.")
    w("* Sandage-Loeb redshift drift: dz/dt0 = (1+z) H0 - H(z) (Sandage 1962; Loeb 1998), dv = c dz/(1+z); astropy Planck18 (Planck 2018 VI, TT,TE,EE+lowE+lensing+BAO), derived:")
    w("")
    w("| z | dv/dt (cm/s/yr) | over 20 yr (cm/s) |")
    w("|---|---|---|")
    for zz, d in R["sandage_loeb_planck18"].items():
        w(f"| {zz} | {d['vdot_cm_s_per_yr']:+.3f} | {d['dv_cm_s_in_20yr']:+.2f} |")
    w("")
    pa = R["perspective_accel_hypothetical"]
    w(f"* Perspective (secular) acceleration dv_r/dt = v_t^2/d: hypothetical star with mu = 1 arcsec/yr at 10 pc (v_t = {pa['v_t_kms']:.2f} km/s) -> {pa['dvr_dt_m_s_per_yr']:.3f} m/s/yr (derived).")
    w("* Peculiar-velocity flows: Virgo infall, Great Attractor, Shapley (Mould+2000 model used by NED), Local Void / Local Sheet decomposition 259 + 185 + 455 km/s -> 631 km/s (Tully+2008); "
      "Cosmicflows-4: 55 877 galaxy distances in 38 065 groups (Tully+2023).")
    w("")
    w("## Additional literature values")
    w("")
    w(f"* M31 relative to the Milky Way (van der Marel et al. 2012, ApJ 753, 8): V_rad = -109.3 +- 4.4 km/s, V_tan = 17.0 km/s (< 34.3 km/s at 1 sigma); M_LG = (3.17 +- 0.57)e12 Msun.")
    w(f"* COBE DMR (Kogut et al. 1993): dipole 3.365 +- 0.027 mK toward (264.4 +- 0.3, 48.4 +- 0.5); with FIRAS T0 = 2.728 K, c dT/T0 = {R['kogut1993_v_from_dT_with_FIRAS_T0_kms']:.1f} km/s (derived). LG-CMB 627 +- 22 km/s toward (276 +- 3, 30 +- 3).")
    w("* Planck 2018 I Table 3: Sun-LSR 17.9 +- 2.0 to (48 +- 7, 23 +- 4) [Schoenrich+2010]; LSR-GC 239 +- 5 to (90, 0) [McMillan 2011]; GC-CMB 565 +- 5 to (265.76, 28.38); Sun-LG 299 +- 15 to (98.4, -5.9) [Diaz+2014]; LG-CMB 620 +- 15 to (271.9, 29.6) (non-relativistic vector subtraction). Dipole beta = (1.23357 +- 0.00036)e-3; amplitude 3362.08 +- 0.99 uK; apex RA, Dec = 167.942, -6.944 (+-0.007) J2000.")
    w("* NED predefined apices (velc_help): heliocentric->GSR 232.3 km/s to (87.8, +1.7) [RC3]; ->LG 316 to (93, -4) [Karachentsev & Makarov 1996]; ->3K 371.0 to (264.14, +48.26) [Fixsen+1996]. NED error model: 4% (GSR), 6% (LG), 7% (3K), 7% (flow-model) of the correction added in quadrature.")
    w("* Mould+2000 Table A1 (RA/Dec B1950): Virgo 12h28m19s +12d40m, V_H 1035, V_LG 957, V_fid 200, radius 10 deg, 600-2300 km/s; Great Attractor 13h20m -44d00m, 4600, 4380, 400, 10 deg, 2600-6600; Shapley 13h30m -31d00m, 13800, 13600, 85, 12 deg, 10000-16000; gamma = 2. "
      "Eq. 1: V_infall = V_fid cos(theta) + V_fid [(V_o - V_a cos theta)/r_oa] (r_oa/V_a)^(1-gamma), r_oa = sqrt(V_o^2 + V_a^2 - 2 V_o V_a cos theta). NED: Dec signs of GA and Shapley negative and all signs in eq. A2 positive.")
    w("")
    w("## Sign conventions")
    w("")
    w("* astropy radial_velocity_correction: positive = add to observed RV. Frame apex V toward (l_a, b_a) is the Sun's velocity relative to that frame; v_frame = v_helio + V cos(theta).")
    w("* SLALIB sla_RVLSRK / sla_RVLSRD / sla_RVEROT return the component that is 'positive when the Sun [observer] is receding from the given point' - check sign before combining with astropy numbers.")
    w("* Galactic (U, V, W): U toward the Galactic centre (l=0), V toward l=90 (rotation), W toward the NGP. YTS77/NED formulas use the same axes.")
    w("")
    w("## Sources")
    w("")
    for src in R["sources"]:
        extra = "; ".join(x for x in (("bibcode " + src["bibcode"]) if src["bibcode"] else "", ("doi:" + src["doi"]) if src["doi"] else "", src["url"]) if x)
        w(f"* {src['citation']}. {extra}. Used for: {src['used_for']}.")
    w("")
    w("Bibcodes were validated (DOI match) through the public ADS link_gateway with `scripts/research/check_bibcodes.py`; journal metadata via Crossref (`scripts/research/crossref_check.py`).")
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")


if __name__ == "__main__":
    warnings.simplefilter("default")
    main()
