# Velocity reference frames and corrections for precise spectroscopy

Generated 2026-11-23T08:50:13+00:00 by `scripts/research/velocity_examples.py` (astropy 7.2.0, ephemeris ~/JPLDE/de440s.bsp).
All literature inputs are typed from the cited papers (see Sources); every other number is **derived** by the script and is labelled so. Machine-readable values: `data/notes/velocity_frames.json`.

## 0. Velocity <-> wavelength scale (why 1e-3 Angstrom)

Doppler: dlambda/lambda = v/c (first order). Derived from c = 299 792 458 m/s (exact):

| lambda (A) | 1 m/s = ... mA | 1 mA = ... m/s |
|---|---|---|
| 1216 | 0.00406 | 246.54 |
| 3000 | 0.01001 | 99.93 |
| 5000 | 0.01668 | 59.96 |
| 6563 | 0.02189 | 45.68 |
| 10000 | 0.03336 | 29.98 |
| 21000 | 0.07005 | 14.28 |

So 1e-3 A (1 mA) corresponds to 14-247 m/s between 2.1 um and Ly-alpha. At 5000 A (derived): Earth rotation (465.1 m/s, equator) = 7.76 mA; Earth orbit (29.78 km/s) = 0.497 A; Sun about the SSB (<= 16.1 m/s) = 0.27 mA; Earth about the EMB (~12.4 m/s) = 0.21 mA; relativistic terms of the barycentric correction (4.71 m/s here) = 0.079 mA; 1 m/s = 1.67e-05 A; 1 cm/s = 1.67e-07 A. Different LSR conventions differ by 1.4-10.7 km/s in vector terms (= 24-178 mA). Hence at 1 mA precision Earth rotation and every frame-definition choice are resolved (so the topocentric/barycentric/LSR bookkeeping is mandatory), while the Sun's SSB reflex, the Earth-EMB wobble and the relativistic terms stay below one unit of the last digit; m/s-class (exoplanet) work needs ~1e-5 A wavelengths and the full relativistic correction.

## 1. Observer (topocentric -> barycentric) terms

* **Earth rotation.** omega = 7.292115e-05 rad/s (IERS Conventions 2010, TN36 Table 1.2, GRS80/WGS84 nominal mean angular velocity), a = 6 378 137 m -> v_eq = omega a = **465.10 m/s** (derived). At Keck/Mauna Kea (distance from the rotation axis 6006.2 km) v = 437.98 m/s.
  Line-of-sight projection: the observer's velocity component *toward* a target at hour angle H, declination dec is v_toward = -omega rho cos(phi') cos(dec) sin(H) (rho cos phi' = distance from the axis; H > 0 west). The correction to *add* to the observed RV is +v_toward, i.e. the observed RV contains +omega rho cos phi' cos dec sin H (receding for setting targets).
  Rate of change: up to omega^2 a = 3.39 cm/s per second (derived; WE14 quote 3.4 cm/s/s), so the exposure time must be known to <~0.25 s for 1 cm/s (WE14 sec 6.2.2).
* **Earth orbit** (DE440s, 2026, derived): Earth-centre barycentric speed 29.294-30.278 km/s (mean 29.783); EMB heliocentric 29.293-30.286 km/s. NASA Earth Fact Sheet: mean 29.78, max 30.29, min 29.29 km/s.
* **Earth about the Earth-Moon barycentre** (DE440s, 2026, derived): 11.73-13.41 m/s (mean 12.43).
* **Sun about the SSB** (DE440s 1990-2040, derived): 8.52-16.08 m/s; Jupiter-only reflex 12.44 m/s (range 11.86-13.08). Sun-SSB distance 0.06-1.98 R_sun (<= 4.60 s light time). This is the heliocentric-minus-barycentric velocity difference (plus relativistic terms, which astropy applies only to the barycentric kind).
* **Light time.** 1 au / c = 499.005 s (derived). The Earth's acceleration (5.709-6.156 mm/s^2) changes its velocity by ~2.96 m/s in that time, so the barycentric *velocity* must be evaluated at the observation epoch (UTC->TT/TDB), not at BJD; BJD (Roemer-corrected) is the time tag for the RV time series (WE14 sec 6.2.2; Eastman et al. 2010).
* **Relativistic barycentric correction** (Wright & Eastman 2014, PASP 126, 838, eq. 28):

  `z_B = gamma_obs * (1 + beta_obs . rho_hat) / (1 + z_GR) * (1 + beta_s . r0_hat)/(1 + beta_s . rho_hat) - 1 - z_LT - z_SD`, with `z_GR = -sum_i G M_i / (c^2 d_i)` (WE14 eq. 2; Sun and Earth matter at 1 cm/s; ~3 m/s), and the correction applied **multiplicatively**: `(1 + z_true) = (1 + z_meas)(1 + z_B)` (WE14 eq. 10), i.e. `v_true ~= v_meas + c z_B + v_meas z_B` (cross term ~ c beta^2 ~ 3 m/s, needed at 1 m/s).
  WE14 budget: Galilean beta.n adequate at 10 m/s; multiplicative cross term at 1 m/s; precession matters at 1 m/s; parallax-proper-motion term at 10 cm/s for nearby stars; GR (~3 m/s, varying ~3 cm/s with orbital eccentricity) and time dilation (variation ~3 cm/s) at 1 cm/s; nutation 3 cm/s; UT1-UTC 2 cm/s; star coordinates to 69 mas; Shapiro delay <= 1.5 mm/s for tau Ceti; their code agrees with TEMPO2 to 0.037 cm/s RMS (0.24 cm/s peak-to-peak).
* **IAU 2000 Resolution C1** (text in Lindegren & Dravins 2003, A&A 401, 1185): report the *barycentric radial-velocity measure* c z_B, with z_B = (lambda - lambda0)/lambda0 for a hypothetical observer at zero gravitational potential at rest at the SSB, after correcting for solar-system gravitational effects and the observer's motion; epoch = barycentric time of light arrival. It is *not* a physical velocity. (LD03's z_B is the source's barycentric redshift; WE14's z_B is the correction factor - different quantities with the same symbol.)
* **Solar gravitational redshift** (relevant to solar Fraunhofer lines): GM/(R c) = 636.486 +- 0.024 m/s at the photosphere, 633.5 m/s for light received at 1 au (LD03 sec 3.1); photospheric lines are blueshifted ~400 m/s by convection (LD03 sec 3.2). Recomputed with IAU 2015 nominal GM_sun, R_sun: 636.31 m/s (derived).
* **Exposure midpoint.** Use the photon-weighted (exposure-meter) midpoint; geometric midpoints can differ by >250 ms (WE14). Even with the photon-weighted midpoint, the non-linearity of v_bary(t) gives a typical 10 cm/s, sometimes >1 m/s, second-order error (Tronsgaard et al. 2019, MNRAS 489, 2395); chromatic (airmass-dependent) flux weighting adds >= several cm/s (Blackman et al. 2017, ApJ 837, 18) and up to 1 m/s in EXPRES data (Blackman, Ong & Fischer 2019, AJ 158, 40).

### astropy `SkyCoord.radial_velocity_correction` (verified in astropy 7.2 source)

* Returns the velocity to **add** to an observed RV (`kind='barycentric'` default, or `'heliocentric'`). Barycentric: `zb = gamma_obs (1 + beta_obs . n)/(1 + gr/c)`, `vcorr = (zb - 1) c`, i.e. WE14 eq. 28 'neglecting Shapiro delay and effects of the star's own motion' (source comment); `beta_obs` = (SSB->geocentre velocity from the selected ephemeris + observer GCRS velocity)/c; `gr` = `EarthLocation.gravitational_redshift` = sum of -GM/(c d) for Sun, Jupiter, Moon and Earth (Earth at the site's geocentric radius).
* If the SkyCoord carries full space motion (distance, pm, rv) it multiplies by (1 + beta_s . r0)/(1 + beta_s . n) (WE14 parallax/proper-motion term) and, when obstime is set, applies space motion; distances enable a parallax-corrected direction.
* Heliocentric kind = plain Galilean projection (no gamma, no GR) - astropy says not to use it below m/s.
* The docstring states the algorithm is 'sufficient to perform corrections at the mm/s level' and that the result is in the optical convention; apply as `rv + vcorr + rv*vcorr/c`.
* Default ephemeris is 'builtin' (ERFA epv00); set `solar_system_ephemeris` to a JPL kernel for best accuracy. Omitted vs WE14: Shapiro delay z_SD, light-travel term z_LT, no exposure averaging.
* Year-2026 test (Keck, target below, every 6 h, N=1460): astropy(DE440s) minus an independent WE14 eq.-28 implementation driven by jplephem/DE440s: max |diff| = 0.1228 mm/s (RMS 0.0815 mm/s); astropy 'builtin' minus DE440s: max 0.243 cm/s (RMS 0.095 cm/s); full minus Galilean: 4.512 to 4.788 m/s. (Both implementations share astropy's observer GCRS velocity/Earth orientation; this checks formula + ephemeris handling, not Earth orientation.)
* astropy LSRK/LSRD/LSR/Galactocentric velocity changes are Galilean vector offsets (AffineTransform), not relativistic compositions.

## 2-5. Frame apices (Sun's velocity relative to each frame)

Correction to the barycentric (heliocentric) RV toward target (l, b): **v_frame = v_bary + V [cos b cos b_a cos(l - l_a) + sin b sin b_a]** (first order; Galilean). Apex RA/Dec (ICRS ~ J2000) and (U,V,W) are **derived** with astropy from the published (V, l, b) or (U, V, W). Last two columns: projected correction for the worked-example target (l = 213.8010, b = 50.2637) and its 1-sigma from Monte-Carlo of the published errors.

| Frame | V (km/s) | l (deg) | b (deg) | RA (ICRS) | Dec (ICRS) | (U, V, W) km/s | corr @ target (km/s) | sigma | Source |
|---|---|---|---|---|---|---|---|---|---|
| LSRK (standard solar motion, radio convention) | 20.00 | 56.157 | 22.765 | 18h03m50.27s (270.959) | +30d00m16.8s (30.005) | 10.27, 15.32, 7.74 | -4.952 | - | Gordon 1975 (Meth. Exp. Phys. 12C, sec 6.1.5); astropy LSRK; SLALIB sla_RVLSRK |
| LSRD (Delhaye 1965 peculiar solar motion) | 16.55 | 53.130 | 25.017 | 17h49m58.67s (267.494) | +28d07m03.9s (28.118) | 9.00, 12.00, 7.00 | -3.665 | - | Delhaye 1965 (Galactic Structure, Stars & Stellar Systems vol 5); astropy LSRD; SLALIB sla_RVLSRD |
| LSR (Schoenrich, Binney & Dehnen 2010; stat+sys) | 18.04 | 47.796 | 23.690 | 17h48m13.99s (267.058) | +23d05m47.6s (23.097) | 11.10, 12.24, 7.25 | -4.674 | 1.10 | Schoenrich, Binney & Dehnen 2010 MNRAS 403, 1829 (abstract); astropy LSR / GalacticLSR default |
| LSR (Dehnen & Binney 1998) | 13.38 | 27.699 | 32.409 | 16h46m03.19s (251.513) | +10d07m44.2s (10.129) | 10.00, 5.25, 7.17 | -1.665 | 0.41 | Dehnen & Binney 1998 MNRAS 298, 387 (published abstract; arXiv v1 abstract gives V0=5.23) |
| LSR (Bland-Hawthorn & Gerhard 2016 adopted) | 16.43 | 47.726 | 25.214 | 17h41m54.92s (265.479) | +23d34m09.5s (23.569) | 10.00, 11.00, 7.00 | -3.841 | 0.97 | Bland-Hawthorn & Gerhard 2016 ARA&A 54, 529 (adopted) |
| LSR (Ding+2019, Gaia DR2) | 12.24 | 28.880 | 36.377 | 16h33m05.20s (248.272) | +12d41m50.4s (12.697) | 8.63, 4.76, 7.26 | -0.694 | 0.47 | Ding, Zhu & Liu 2019 RAA 19, 068 (Gaia DR2) |
| GSR: IAU 1985 Theta0=220 + Delhaye (9,12,7)  [derived] | 232.28 | 87.778 | 1.727 | 20h55m26.77s (313.862) | +47d49m23.5s (47.823) | 9.00, 232.00, 7.00 | -81.903 | - | derived: Kerr & Lynden-Bell 1986 + Delhaye 1965 (as in Courteau & vdB 1999 / RC3) |
| GSR as used by NED (RC3) | 232.30 | 87.800 | 1.700 | 20h55m39.10s (313.913) | +47d49m20.0s (47.822) | 8.91, 232.03, 6.89 | -81.950 | - | NED velocity-correction help (RC3, de Vaucouleurs et al. 1991) |
| GSR: astropy Galactocentric v4.0 default (approx. Galactic axes) | 246.06 | 86.993 | 1.812 | 20h52m04.43s (313.018) | +47d16m34.3s (47.276) | 12.90, 245.60, 7.78 | -88.211 | - | astropy Galactocentric 'v4.0' defaults: Drimmel & Poggio 2018 RNAAS 2, 210; GRAVITY 2018 A&A 615, L15; Reid & Brunthaler 2004; Bennett & Bovy 2019 |
| GSR: astropy pre-v4.0 (approx. Galactic axes) | 232.62 | 87.264 | 1.786 | 20h53m12.70s (313.303) | +47d28m02.4s (47.467) | 11.10, 232.24, 7.25 | -82.911 | - | astropy 'pre-v4.0' defaults (Schoenrich+2010; Bovy 2015; Gillessen+2009; Chen+2001) |
| GSR: Reid+2019 Theta0+Vsun=247+-4 with Schoenrich U,W  [derived composite] | 247.36 | 87.427 | 1.680 | 20h54m18.91s (313.579) | +47d31m26.7s (47.524) | 11.10, 247.00, 7.25 | -88.160 | 1.63 | derived: Reid+2019 ApJ 885, 131 (V) + Schoenrich+2010 (U, W) |
| GSR: Sgr A* pm (Reid & Brunthaler 2020) x R0 (GRAVITY2019) + U=11.1  [derived] | 248.93 | 87.443 | 1.955 | 20h53m07.50s (313.281) | +47d42m46.1s (47.713) | 11.10, 248.54, 8.49 | -87.754 | 0.73 | derived: V_phi=4.74047*|mu_l|*R0, V_z=4.74047*|mu_b|*R0 with mu from Reid & Brunthaler 2020, R0 from GRAVITY2019; U from Schoenrich+2010 |
| GSR: Sgr A* pm (Reid & Brunthaler 2020) x R0 (GRAVITY2021) + U=11.1  [derived] | 251.88 | 87.473 | 1.955 | 20h53m14.32s (313.310) | +47d44m08.9s (47.736) | 11.10, 251.49, 8.59 | -88.725 | 0.77 | derived: V_phi=4.74047*|mu_l|*R0, V_z=4.74047*|mu_b|*R0 with mu from Reid & Brunthaler 2020, R0 from GRAVITY2021; U from Schoenrich+2010 |
| LG: RC2/'IAU' 300 sin l cos b | 300.00 | 90.000 | 0.000 | 21h12m01.05s (318.004) | +48d19m46.7s (48.330) | 0.00, 300.00, 0.00 | -106.687 | - | RC2 (de Vaucouleurs, de Vaucouleurs & Corwin 1976) 'IAU' 300 sin l cos b; quoted by Mould+2000, Courteau & vdB 1999 |
| LG: Yahil, Tammann & Sandage 1977 | 308.47 | 104.943 | -6.702 | 22h52m03.66s (343.015) | +51d57m03.7s (51.951) | -79.00, 296.00, -36.00 | -90.983 | - | Yahil, Tammann & Sandage 1977 ApJ 217, 903; formula as given in Mould+2000 eq. A1 |
| LG: Karachentsev & Makarov 1996 (NED) | 316.00 | 93.000 | -4.000 | 21h41m28.91s (325.370) | +47d31m36.6s (47.527) | -16.50, 314.80, -22.04 | -120.137 | 10.45 | Karachentsev & Makarov 1996 AJ 111, 794 (abstract); used by NED |
| LG: Courteau & van den Bergh 1999 | 306.00 | 99.000 | -4.000 | 22h10m24.04s (332.600) | +51d13m54.2s (51.232) | -47.75, 301.50, -21.35 | -98.268 | 22.75 | Courteau & van den Bergh 1999 AJ 118, 337 (abstract) |
| LG: Diaz+2014 (Planck 2018 I Table 3) | 299.00 | 98.400 | -5.900 | 22h14m07.87s (333.533) | +49d19m52.0s (49.331) | -43.45, 294.23, -30.73 | -105.189 | 16.55 | Diaz et al. 2014 MNRAS 443, 1688 (abstract); Planck 2018 I Table 3 |
| CMB: Planck 2018 (HFI) solar dipole | 369.82 | 264.021 | 48.253 | 11h11m46.06s (167.942) | -06d56m39.3s (-6.944) | -25.65, -244.90, 275.92 | +312.899 | 0.10 | Planck Collab. 2020 (2018 results I) A&A 641, A1, Table 2/3 & sec 3 |
| CMB: Fixsen+1996 (NED's 3K frame) | 371.00 | 264.140 | 48.260 | 11h12m02.85s (168.012) | -06d58m58.9s (-6.983) | -25.22, -245.70, 276.83 | +313.655 | 0.93 | Fixsen et al. 1996 ApJ 473, 576 (abstract, 95% CL); V=371.0 as used by NED |
| LG_CMB_Kogut1993 | 627.00 | 276.000 | 30.000 | 11h06m37.60s (166.657) | -27d20m00.5s (-27.333) | 56.76, -540.02, 313.50 | +402.975 | 27.40 | Kogut et al. 1993 ApJ 419, 1 (abstract) |
| LG_CMB_Planck2018 | 620.00 | 271.900 | 29.600 | 10h51m48.94s (162.954) | -25d58m02.4s (-25.967) | 17.87, -538.79, 306.24 | +417.612 | 16.30 | Planck 2018 I Table 3 (from Diaz+2014 Sun-LG) |
| GC_CMB_Planck2018 | 565.00 | 265.760 | 28.380 | 10h29m03.84s (157.266) | -24d01m44.0s (-24.029) | -36.75, -495.73, 268.55 | +402.335 | 3.95 | Planck 2018 I Table 3 (Schoenrich+2010 + McMillan 2011) |
| Sun_LSR_Planck2018 | 17.90 | 48.000 | 23.000 | 17h51m19.63s (267.832) | +23d01m32.0s (23.026) | 11.03, 12.24, 6.99 | -4.833 | 1.35 | Planck 2018 I Table 3 (from Schoenrich+2010) |

* LSRK: 20 km/s toward RA 18h, Dec +30 (B1900) -> ICRS apex 18h03m50.27s +30d00m16.8s (derived via astropy FK4(B1900)->ICRS); ICRS vector 0.29000, -17.31726, 10.00141 km/s, identical to astropy's hard-coded V_OFFSET_LSRK; Galactic (U,V,W) = (10.27, 15.32, 7.74) km/s (Reid et al. 2009, ApJ 700, 137, Table 1 note: 10.3, 15.3, 7.7 km/s).
* NED 'GSR' (232.3 km/s to l=87.8, b=+1.7, RC3) equals IAU-1985 Theta0 = 220 km/s plus Delhaye (9,12,7): derived apex 232.28 km/s to (87.78, 1.73).
* YTS77 formula V_LG = V_H - 79 cos l cos b + 296 sin l cos b - 36 sin b is a 308.5 km/s apex toward (l, b) = (104.9, -6.7) (derived).
* GSR: Sgr A* pm (Reid & Brunthaler 2020) x R0 (GRAVITY2019) + U=11.1  [derived]: V_phi = 248.54 +- 0.84 km/s, V_z = 8.49 km/s (derived; 4.74047 km/s per mas/yr per kpc).
* GSR: Sgr A* pm (Reid & Brunthaler 2020) x R0 (GRAVITY2021) + U=11.1  [derived]: V_phi = 251.49 +- 1.09 km/s, V_z = 8.59 km/s (derived; 4.74047 km/s per mas/yr per kpc).
* Vector closure: Sun->LG (Diaz+2014) + LG->CMB (Planck Table 3) = 369.3 km/s toward (264.03, 48.25) vs Planck Sun->CMB 369.82 toward (264.021, 48.253).
* Planck apex RA/Dec computed 167.942, -6.944 vs published 167.942, -6.944 (J2000).
* Fixsen 1996: c dT/T0 = 370.56 km/s from 3.372 mK / 2.728 K (derived); NED adopts 371.0 km/s.
* astropy Galactocentric defaults (astropy 7.2: 'latest' = 'v4.0'): galcen RA, Dec = 266.4051, -28.936175 (Reid & Brunthaler 2004), R0 = 8.122 kpc (GRAVITY 2018), galcen_v_sun = (12.9, 245.6, 7.78) km/s (Drimmel & Poggio 2018), z_sun = 20.8 pc (Bennett & Bovy 2019), roll = 0. 'pre-v4.0': 8.3 kpc, (11.1, 232.24, 7.25), 27 pc.
* Galactic constants: IAU 1985 R0 = 8.5 kpc, Theta0 = 220 km/s (Kerr & Lynden-Bell 1986); GRAVITY R0 = 8178 +- 13 +- 22 pc (2019), 8275 +- 9 +- 33 pc (2021), 8277 +- 9 (+-33 sys) pc (2022); Reid+2019 R0 = 8.15 +- 0.15 kpc, Theta0 = 236 +- 7, Theta0 + V_sun = 247 +- 4 km/s, Omega = 30.32 +- 0.27 km/s/kpc; Bland-Hawthorn & Gerhard 2016: R0 = 8.2 +- 0.1 kpc, Theta0 = 238 +- 15, V_g,sun = 248 +- 3 km/s, z0 = 25 +- 5 pc; Sgr A* apparent motion -6.411 +- 0.008 mas/yr along the plane, -0.219 +- 0.007 mas/yr toward the NGP (Reid & Brunthaler 2020).

astropy frame-transform cross-check at the target (object at rest w.r.t. SSB, 1 Gpc): LSRK -4.9518, LSRD -3.6654, LSR_default -4.6740, GSR_astropy_v4.0 -88.1752, GSR_astropy_pre-v4.0 -82.8711 km/s (compare the projection column).

## 6. Worked example: RA = 10h, Dec = +20 deg, Keck/Mauna Kea, 2026-03-20T10:00:00 UTC

Geometry (derived): l = 213.8010, b = 50.2637; ecliptic lon/lat 145.124/7.291; altitude 69.22 deg, azimuth 273.92, airmass 1.070; apparent LST 11.4984 h, apparent (TETE) RA/Dec 150.3697/19.8728, hour angle +1.4737 h; TDB 2026-03-20T10:01:09.186; UT1-UTC 0.0571 s (IERS-A).

| term | m/s |
|---|---|
| Barycentric correction, astropy + DE440s (to add) | **-16495.0646** |
| same, independent WE14 eq.28 (jplephem DE440s) | -16495.0647 (diff +0.1194 mm/s) |
| same, astropy builtin ephemeris (ERFA epv00) | -16495.0645 (diff +0.011 cm/s) |
| Heliocentric correction, astropy (Galilean) | -16490.8396 |
| barycentric - heliocentric | -4.2249 |
| -- Earth rotation (observer GCRS velocity . n) | -155.0368 |
| -- rotation, omega_nominal (z_CIP x r) . n | -155.0368 |
| -- rotation, -omega rho cos(dec) sin(H) (apparent dec, H) | -155.0051 |
| -- EMB barycentric orbital motion . n | -16334.2598 |
| -- Earth about EMB . n | -10.4813 |
| -- Galilean total beta.n c | -16499.7779 |
| -- relativistic remainder (gamma, GR, cross) | +4.7132 |
| ---- (gamma - 1) c | +1.5331 |
| ---- GR terms -GM/(c d): sun -2.97131, jupiter -0.00058, moon -0.00004, earth -0.20840 | |
| Sun about SSB . n (helio vs bary) | -8.9383 |
| Shapiro term (WE14 eq. 27; Sun elongation 145.0 deg) | +0.1844 mm/s |
| Roemer delay r_obs.n/c | +406.636 s |
| correction change if evaluated at t + Roemer delay (wrong epoch; dominated by Earth rotation) | -13.2463 |
| d(vcorr)/dt | -3.2747 cm/s per s |
| exposure mean - midpoint, 1800 s uniform flux | +11.2585 cm/s |
| exposure mean - midpoint, 3600 s uniform flux | +44.7576 cm/s |
| cross term v_meas z_B for v_meas = 10 / 50 / 300 km/s | -0.550 / -2.751 / -16.506 |

Over 2026 the barycentric correction for this target from Keck spans -29.625 to +30.261 km/s (derived).

Frame corrections for the same target (add to the barycentric RV; first-order Galilean; km/s):

| frame | correction (km/s) | 1-sigma |
|---|---|---|
| LSRK (standard solar motion, radio convention) | -4.952 | (conventional) |
| LSRD (Delhaye 1965 peculiar solar motion) | -3.665 | (conventional) |
| LSR (Schoenrich, Binney & Dehnen 2010; stat+sys) | -4.674 | 1.10 |
| LSR (Dehnen & Binney 1998) | -1.665 | 0.41 |
| LSR (Bland-Hawthorn & Gerhard 2016 adopted) | -3.841 | 0.97 |
| LSR (Ding+2019, Gaia DR2) | -0.694 | 0.47 |
| GSR: IAU 1985 Theta0=220 + Delhaye (9,12,7)  [derived] | -81.903 | (conventional) |
| GSR as used by NED (RC3) | -81.950 | (conventional) |
| GSR: astropy Galactocentric v4.0 default (approx. Galactic axes) | -88.211 | (conventional) |
| GSR: astropy pre-v4.0 (approx. Galactic axes) | -82.911 | (conventional) |
| GSR: Reid+2019 Theta0+Vsun=247+-4 with Schoenrich U,W  [derived composite] | -88.160 | 1.63 |
| GSR: Sgr A* pm (Reid & Brunthaler 2020) x R0 (GRAVITY2019) + U=11.1  [derived] | -87.754 | 0.73 |
| GSR: Sgr A* pm (Reid & Brunthaler 2020) x R0 (GRAVITY2021) + U=11.1  [derived] | -88.725 | 0.77 |
| LG: RC2/'IAU' 300 sin l cos b | -106.687 | (conventional) |
| LG: Yahil, Tammann & Sandage 1977 | -90.983 | (conventional) |
| LG: Karachentsev & Makarov 1996 (NED) | -120.137 | 10.45 |
| LG: Courteau & van den Bergh 1999 | -98.268 | 22.75 |
| LG: Diaz+2014 (Planck 2018 I Table 3) | -105.189 | 16.55 |
| CMB: Planck 2018 (HFI) solar dipole | +312.899 | 0.10 |
| CMB: Fixsen+1996 (NED's 3K frame) | +313.655 | 0.93 |
| LG_CMB_Kogut1993 | +402.975 | 27.40 |
| LG_CMB_Planck2018 | +417.612 | 16.30 |
| GC_CMB_Planck2018 | +402.335 | 3.95 |
| Sun_LSR_Planck2018 | -4.833 | 1.35 |

Total for an observed topocentric RV v_obs (km/s): v_LSRK = v_obs + -16.4951 + -4.952 (+ v_obs*vcorr/c); v_GSR(astropy v4.0) = v_obs + -16.4951 + -88.211; v_LG(KM96) = v_obs + -16.4951 + -120.137; v_CMB(Planck) = v_obs + -16.4951 + +312.899 (first order).

**Hypothetical galaxy, cz_hel = 3000 km/s** in this direction (illustration): Sun-CMB projection +312.899 km/s. Additive cz_CMB = 3312.899 km/s; exact multiplicative (1+z_hel) = (1+z_CMB)(1+z_sun) (Davis+2011 eq. 6; Carr+2022 eq. 2-3) gives 3316.361 km/s (difference +3.461 km/s); with exact SR z_sun = 1/[gamma(1+beta cos theta)] - 1: cz_CMB = 3316.261 km/s. At z_hel = 1 the additive error is 1.05e-03 in z.

Mould+2000 multi-attractor model (NED_KM1996 LG step): V_LG = 2879.9; Virgo: theta 36.9 deg, V_infall +244.2; GreatAttractor: theta 79.2 deg, V_infall +233.5; Shapley: theta 72.3 deg, V_infall +17.3 -> V_cosmic = 3374.8 km/s (derived; gamma = 2; NED sign convention 'all positive').
Mould+2000 multi-attractor model (Mould_YTS1977 LG step): V_LG = 2909.0; Virgo: theta 36.9 deg, V_infall +243.2; GreatAttractor: theta 79.2 deg, V_infall +234.9; Shapley: theta 72.3 deg, V_infall +17.5 -> V_cosmic = 3404.6 km/s (derived; gamma = 2; NED sign convention 'all positive').

## 7. Secular and cosmology-scale terms

* SSB acceleration (Gaia EDR3; Klioner+2021): (2.32 +- 0.16)e-10 m/s^2 toward RA 269.1, Dec -31.6 -> up to 7.32 mm/s per year in apparent RVs of distant sources; d(v_bary)/dt = -a.n on the target +4.16 mm/s/yr (derived). The 'barycentric' frame is not inertial at this level.
* Sandage-Loeb redshift drift: dz/dt0 = (1+z) H0 - H(z) (Sandage 1962; Loeb 1998), dv = c dz/(1+z); astropy Planck18 (Planck 2018 VI, TT,TE,EE+lowE+lensing+BAO), derived:

| z | dv/dt (cm/s/yr) | over 20 yr (cm/s) |
|---|---|---|
| 0.5 | +0.251 | +5.01 |
| 1.0 | +0.225 | +4.50 |
| 2.0 | -0.011 | -0.22 |
| 3.0 | -0.280 | -5.61 |
| 4.0 | -0.537 | -10.74 |
| 5.0 | -0.776 | -15.53 |

* Perspective (secular) acceleration dv_r/dt = v_t^2/d: hypothetical star with mu = 1 arcsec/yr at 10 pc (v_t = 47.40 km/s) -> 0.230 m/s/yr (derived).
* Peculiar-velocity flows: Virgo infall, Great Attractor, Shapley (Mould+2000 model used by NED), Local Void / Local Sheet decomposition 259 + 185 + 455 km/s -> 631 km/s (Tully+2008); Cosmicflows-4: 55 877 galaxy distances in 38 065 groups (Tully+2023).

## Additional literature values

* M31 relative to the Milky Way (van der Marel et al. 2012, ApJ 753, 8): V_rad = -109.3 +- 4.4 km/s, V_tan = 17.0 km/s (< 34.3 km/s at 1 sigma); M_LG = (3.17 +- 0.57)e12 Msun.
* COBE DMR (Kogut et al. 1993): dipole 3.365 +- 0.027 mK toward (264.4 +- 0.3, 48.4 +- 0.5); with FIRAS T0 = 2.728 K, c dT/T0 = 369.8 km/s (derived). LG-CMB 627 +- 22 km/s toward (276 +- 3, 30 +- 3).
* Planck 2018 I Table 3: Sun-LSR 17.9 +- 2.0 to (48 +- 7, 23 +- 4) [Schoenrich+2010]; LSR-GC 239 +- 5 to (90, 0) [McMillan 2011]; GC-CMB 565 +- 5 to (265.76, 28.38); Sun-LG 299 +- 15 to (98.4, -5.9) [Diaz+2014]; LG-CMB 620 +- 15 to (271.9, 29.6) (non-relativistic vector subtraction). Dipole beta = (1.23357 +- 0.00036)e-3; amplitude 3362.08 +- 0.99 uK; apex RA, Dec = 167.942, -6.944 (+-0.007) J2000.
* NED predefined apices (velc_help): heliocentric->GSR 232.3 km/s to (87.8, +1.7) [RC3]; ->LG 316 to (93, -4) [Karachentsev & Makarov 1996]; ->3K 371.0 to (264.14, +48.26) [Fixsen+1996]. NED error model: 4% (GSR), 6% (LG), 7% (3K), 7% (flow-model) of the correction added in quadrature.
* Mould+2000 Table A1 (RA/Dec B1950): Virgo 12h28m19s +12d40m, V_H 1035, V_LG 957, V_fid 200, radius 10 deg, 600-2300 km/s; Great Attractor 13h20m -44d00m, 4600, 4380, 400, 10 deg, 2600-6600; Shapley 13h30m -31d00m, 13800, 13600, 85, 12 deg, 10000-16000; gamma = 2. Eq. 1: V_infall = V_fid cos(theta) + V_fid [(V_o - V_a cos theta)/r_oa] (r_oa/V_a)^(1-gamma), r_oa = sqrt(V_o^2 + V_a^2 - 2 V_o V_a cos theta). NED: Dec signs of GA and Shapley negative and all signs in eq. A2 positive.

## Sign conventions

* astropy radial_velocity_correction: positive = add to observed RV. Frame apex V toward (l_a, b_a) is the Sun's velocity relative to that frame; v_frame = v_helio + V cos(theta).
* SLALIB sla_RVLSRK / sla_RVLSRD / sla_RVEROT return the component that is 'positive when the Sun [observer] is receding from the given point' - check sign before combining with astropy numbers.
* Galactic (U, V, W): U toward the Galactic centre (l=0), V toward l=90 (rotation), W toward the NGP. YTS77/NED formulas use the same axes.

## Sources

* Wright, J. T. & Eastman, J. D. 2014, PASP 126, 838-852, 'Barycentric Corrections at 1 cm/s for precise Doppler velocities'. bibcode 2014PASP..126..838W; doi:10.1086/678541; https://arxiv.org/abs/1409.4774. Used for: zB formula (eq. 2, 10, 11, 27, 28), precision budget.
* Lindegren, L. & Dravins, D. 2003, A&A 401, 1185-1201, 'The fundamental definition of radial velocity' (includes IAU 2000 Resolutions C1, C2). bibcode 2003A&A...401.1185L; doi:10.1051/0004-6361:20030181; https://arxiv.org/abs/astro-ph/0302522. Used for: barycentric RV measure, solar gravitational redshift.
* Kanodia, S. & Wright, J. T. 2018, RNAAS 2, 4 (barycorrpy). bibcode 2018RNAAS...2....4K; doi:10.3847/2515-5172/aaa4b7; https://arxiv.org/abs/1801.01634. Used for: barycorrpy precision statement.
* Tronsgaard, R., Buchhave, L. A., Wright, J. T., Eastman, J. D., Blackman, R. T. 2019, MNRAS 489, 2395-2402. bibcode 2019MNRAS.489.2395T; doi:10.1093/mnras/stz2181; https://arxiv.org/abs/1908.00991. Used for: photon-weighted midpoint second-order error.
* Blackman, R. T., Szymkowiak, A. E., Fischer, D. A., Jurgenson, C. A. 2017, ApJ 837, 18. bibcode 2017ApJ...837...18B; doi:10.3847/1538-4357/aa5ead; https://arxiv.org/abs/1703.00012. Used for: chromatic barycentric correction.
* Blackman, R. T., Ong, J. M. J., Fischer, D. A. 2019, AJ 158, 40. bibcode 2019AJ....158...40B; doi:10.3847/1538-3881/ab24c3; https://arxiv.org/abs/1906.01653. Used for: chromatic barycentric correction measured with EXPRES.
* Eastman, J., Siverd, R., Gaudi, B. S. 2010, PASP 122, 935-946. bibcode 2010PASP..122..935E; doi:10.1086/655938. Used for: BJD time reduction (cited via WE14).
* Petit, G. & Luzum, B. (eds.) 2010, IERS Conventions (2010), IERS Technical Note 36, Table 1.2. https://iers-conventions.obspm.fr/content/tn36.pdf. Used for: omega = 7.292115e-5 rad/s, a_E = 6378137 m.
* NASA NSSDCA Earth Fact Sheet (D. R. Williams). https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html. Used for: Earth orbital speed mean/max/min.
* Gordon 1975, Methods of Experimental Physics vol. 12C (ed. Meeks), sec. 6.1.5(.2), p. 281 - as cited in astropy LSRK and SLALIB sla_RVLSRK. https://docs.astropy.org/en/stable/api/astropy.coordinates.LSRK.html. Used for: LSRK definition 20 km/s to 18h,+30 (1900).
* Delhaye, J. 1965, in Galactic Structure (Stars and Stellar Systems V), eds. Blaauw & Schmidt, p. 61. bibcode 1965gast.book...61D; https://irtfweb.ifa.hawaii.edu/~tcs3/tcs3/Misc/slalib_html/node173.html. Used for: LSRD (9,12,7) km/s.
* Wallace, P. T. 2000, SLALIB Starlink User Note 67 (sla_RVLSRK, sla_RVLSRD). https://irtfweb.ifa.hawaii.edu/~tcs3/tcs3/Misc/slalib_html/node174.html. Used for: LSR definitions, sign convention.
* Schoenrich, R., Binney, J., Dehnen, W. 2010, MNRAS 403, 1829-1833. bibcode 2010MNRAS.403.1829S; doi:10.1111/j.1365-2966.2010.16253.x; https://arxiv.org/abs/0912.3693. Used for: (U,V,W)sun = (11.1, 12.24, 7.25).
* Dehnen, W. & Binney, J. J. 1998, MNRAS 298, 387-394. bibcode 1998MNRAS.298..387D; doi:10.1046/j.1365-8711.1998.01600.x; https://arxiv.org/abs/astro-ph/9710077. Used for: (U,V,W) = (10.00, 5.25, 7.17).
* Ding, P.-J., Zhu, Z., Liu, J.-C. 2019, RAA 19, 068 (Gaia DR2 LSR). bibcode 2019RAA....19...68D; doi:10.1088/1674-4527/19/5/68. Used for: (U,V,W) = (8.63, 4.76, 7.26).
* Reid, M. J. et al. 2009, ApJ 700, 137-148. bibcode 2009ApJ...700..137R; doi:10.1088/0004-637X/700/1/137; https://arxiv.org/abs/0902.3913. Used for: standard solar motion (10.3, 15.3, 7.7).
* Kerr, F. J. & Lynden-Bell, D. 1986, MNRAS 221, 1023-1038. bibcode 1986MNRAS.221.1023K; doi:10.1093/mnras/221.4.1023. Used for: IAU 1985 R0 = 8.5 kpc, Theta0 = 220 km/s.
* GRAVITY Collaboration 2018, A&A 615, L15. bibcode 2018A&A...615L..15G; doi:10.1051/0004-6361/201833718; https://arxiv.org/abs/1807.09409. Used for: R0 = 8122 +- 31 pc (astropy v4.0).
* GRAVITY Collaboration 2019, A&A 625, L10. bibcode 2019A&A...625L..10G; doi:10.1051/0004-6361/201935656; https://arxiv.org/abs/1904.05721. Used for: R0 = 8178 +- 13 +- 22 pc.
* GRAVITY Collaboration 2021, A&A 647, A59. bibcode 2021A&A...647A..59G; doi:10.1051/0004-6361/202040208; https://arxiv.org/abs/2101.12098. Used for: R0 = 8275 +- 9 +- 33 pc.
* GRAVITY Collaboration 2022, A&A 657, L12. bibcode 2022A&A...657L..12G; doi:10.1051/0004-6361/202142465; https://arxiv.org/abs/2112.07478. Used for: R0 = 8277 +- 9 pc.
* Reid, M. J. et al. 2019, ApJ 885, 131. bibcode 2019ApJ...885..131R; doi:10.3847/1538-4357/ab4a11; https://arxiv.org/abs/1910.03357. Used for: R0 = 8.15, Theta0 = 236 +- 7, Theta0+Vsun = 247 +- 4.
* Reid, M. J. & Brunthaler, A. 2004, ApJ 616, 872-884. bibcode 2004ApJ...616..872R; doi:10.1086/424960; https://arxiv.org/abs/astro-ph/0408107. Used for: Sgr A* proper motion; astropy galcen_coord.
* Reid, M. J. & Brunthaler, A. 2020, ApJ 892, 39. bibcode 2020ApJ...892...39R; doi:10.3847/1538-4357/ab76cd; https://arxiv.org/abs/2001.04386. Used for: Sgr A* pm -6.411 +- 0.008, -0.219 +- 0.007 mas/yr.
* Bland-Hawthorn, J. & Gerhard, O. 2016, ARA&A 54, 529-596. bibcode 2016ARA&A..54..529B; doi:10.1146/annurev-astro-081915-023441; https://arxiv.org/abs/1602.07702. Used for: R0 = 8.2 +- 0.1, Theta0 = 238 +- 15, Vg = 248 +- 3, z0 = 25 +- 5.
* Drimmel, R. & Poggio, E. 2018, RNAAS 2, 210. bibcode 2018RNAAS...2..210D; doi:10.3847/2515-5172/aaef8b. Used for: astropy v4.0 galcen_v_sun.
* Bennett, M. & Bovy, J. 2019, MNRAS 482, 1417-1425. bibcode 2019MNRAS.482.1417B; doi:10.1093/mnras/sty2813; https://arxiv.org/abs/1809.03507. Used for: z_sun = 20.8 +- 0.3 pc.
* Gaia Collaboration, Drimmel, R. et al. 2023, A&A 674, A37. bibcode 2023A&A...674A..37G; doi:10.1051/0004-6361/202243797; https://arxiv.org/abs/2206.06207. Used for: Gaia DR3 disc kinematics (context; no solar-motion value used).
* Yahil, A., Tammann, G. A., Sandage, A. 1977, ApJ 217, 903. bibcode 1977ApJ...217..903Y; doi:10.1086/155636. Used for: LG correction -79 cos l cos b + 296 sin l cos b - 36 sin b.
* Karachentsev, I. D. & Makarov, D. A. 1996, AJ 111, 794. bibcode 1996AJ....111..794K; doi:10.1086/117825. Used for: Sun-LG 316 +- 5 km/s to (93, -4).
* Courteau, S. & van den Bergh, S. 1999, AJ 118, 337-345. bibcode 1999AJ....118..337C; doi:10.1086/300942; https://arxiv.org/abs/astro-ph/9903298. Used for: Sun-LG 306 +- 18 to (99, -4); RC2 300 sin l cos b; Delhaye LSR.
* Diaz, J. D., Koposov, S. E., Irwin, M., Belokurov, V., Evans, W. 2014, MNRAS 443, 1688-1703. bibcode 2014MNRAS.443.1688D; doi:10.1093/mnras/stu1210; https://arxiv.org/abs/1405.3662. Used for: Sun-LG 299 +- 15 to (98.4, -5.9).
* van der Marel, R. P. et al. 2012, ApJ 753, 8 (M31 velocity vector II). bibcode 2012ApJ...753....8V; doi:10.1088/0004-637X/753/1/8; https://arxiv.org/abs/1205.6864. Used for: M31 radial -109.3 +- 4.4, tangential 17.0 (<34.3) km/s.
* Planck Collaboration 2020, Planck 2018 results I, A&A 641, A1. bibcode 2020A&A...641A...1P; doi:10.1051/0004-6361/201833880; https://arxiv.org/abs/1807.06205. Used for: solar dipole 369.82 +- 0.11 km/s to (264.021, 48.253); Table 3.
* Planck Collaboration 2020, Planck 2018 results VI, A&A 641, A6. bibcode 2020A&A...641A...6P; doi:10.1051/0004-6361/201833910. Used for: Planck18 cosmology (astropy).
* Fixsen, D. J. et al. 1996, ApJ 473, 576-587. bibcode 1996ApJ...473..576F; doi:10.1086/178173; https://arxiv.org/abs/astro-ph/9605054. Used for: COBE/FIRAS dipole 3.372 mK to (264.14, 48.26); NED 3K frame.
* Kogut, A. et al. 1993, ApJ 419, 1. bibcode 1993ApJ...419....1K; doi:10.1086/173453; https://arxiv.org/abs/astro-ph/9312056. Used for: COBE DMR dipole; LG-CMB 627 +- 22 to (276, 30).
* Tully, R. B. et al. 2008, ApJ 676, 184-205. bibcode 2008ApJ...676..184T; doi:10.1086/527428; https://arxiv.org/abs/0705.4139. Used for: Local Void / Virgo / large-scale decomposition 259, 185, 455 -> 631 km/s.
* Mould, J. R. et al. 2000, ApJ 529, 786-794. bibcode 2000ApJ...529..786M; doi:10.1086/308304; https://arxiv.org/abs/astro-ph/9909260. Used for: multi-attractor flow model eqs. 1-2, A1-A2, Table A1.
* NASA/IPAC Extragalactic Database: velocity-correction help and object-result help. https://ned.ipac.caltech.edu/help/velc_help.html ; https://ned.ipac.caltech.edu/help/objresult_help.html. Used for: NED apex vectors and Mould-model corrections.
* Tully, R. B. et al. 2023, ApJ 944, 94 (Cosmicflows-4). bibcode 2023ApJ...944...94T; doi:10.3847/1538-4357/ac94d8; https://arxiv.org/abs/2209.11238. Used for: peculiar-velocity catalogue.
* Davis, T. M. et al. 2011, ApJ 741, 67. bibcode 2011ApJ...741...67D; doi:10.1088/0004-637X/741/1/67; https://arxiv.org/abs/1012.2912. Used for: (1+z) = (1+zbar)(1+z_pec,sun), z_sun = -v.n/c.
* Carr, A. et al. 2022, PASA 39, e046. bibcode 2022PASA...39...46C; doi:10.1017/pasa.2022.41; https://arxiv.org/abs/2112.01471. Used for: exact multiplicative heliocentric->CMB conversion.
* Gaia Collaboration, Klioner, S. A. et al. 2021, A&A 649, A9. bibcode 2021A&A...649A...9G; doi:10.1051/0004-6361/202039734; https://arxiv.org/abs/2012.02036. Used for: SSB acceleration 2.32e-10 m/s^2.
* Sandage, A. 1962, ApJ 136, 319. bibcode 1962ApJ...136..319S; doi:10.1086/147385. Used for: redshift drift.
* Loeb, A. 1998, ApJ 499, L111-L114. bibcode 1998ApJ...499L.111L; doi:10.1086/311375. Used for: redshift drift.
* astropy 7.2.0 source: coordinates/sky_coordinate.py (radial_velocity_correction), earth.py (gravitational_redshift), builtin_frames/lsr.py, builtin_frames/galactocentric.py. https://github.com/astropy/astropy. Used for: implementation details verified locally.
* astropy issue #18550 'A possible LSRK documentation error' (docstring wording vs code). https://github.com/astropy/astropy/issues/18550. Used for: LSR docstring ambiguity; code follows Gordon 1975.

Bibcodes were validated (DOI match) through the public ADS link_gateway with `scripts/research/check_bibcodes.py`; journal metadata via Crossref (`scripts/research/crossref_check.py`).
