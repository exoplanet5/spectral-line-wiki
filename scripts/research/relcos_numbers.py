#!/usr/bin/env python
"""Worked numbers for the relativity / cosmology part of the wavelength-
calibration notes (data/notes/relativity_cosmology.md).

Every number here is DERIVED from a documented formula and documented
constants (IAU / IERS / CODATA / Planck 2018 via astropy).  No catalogue
values are read.  Measured literature values quoted in the notes are NOT
produced here; they come from the cited papers.

Usage:
    ~/.venvs/astro313/bin/python relcos_numbers.py OUT.json
"""
import json
import sys

import numpy as np
from astropy import units as u
from astropy.cosmology import Planck18, FlatLambdaCDM
from scipy.optimize import brentq

# ---------------------------------------------------------------- constants
# taken from astropy.constants (IAU 2015 Res. B3 nominal values, IAU 2012 B2 au)
# and ERFA (IAU 2000 B1.9 L_G, IAU 2006 B3 L_B) so nothing is typed by hand
import astropy.constants as K
import erfa
C = K.c.value                          # 299 792 458 m/s exact
GM_SUN = K.GM_sun.value                # 1.3271244e20 m^3 s^-2 (nominal)
R_SUN = K.R_sun.value                  # 6.957e8 m (nominal)
AU = K.au.value                        # 149 597 870 700 m
GM_EARTH = K.GM_earth.value            # 3.986004e14 m^3 s^-2 (nominal)
R_EARTH_EQ = K.R_earth.value           # 6.3781e6 m (nominal equatorial)
# Earth rotation rate from the Earth Rotation Angle rate (IAU 2000 Res. B1.8):
# ERA = 2pi(0.7790572732640 + 1.00273781191135448 Tu), Tu in UT1 days
OMEGA_EARTH = 2 * np.pi * 1.00273781191135448 / 86400.0
L_G = erfa.ELG                         # 6.969290134e-10 (TT vs TCG)
L_B = erfa.ELB                         # 1.550519768e-8 (TDB vs TCB)
L_C = 1 - (1 - L_B) / (1 - L_G)        # TCG vs TCB mean rate, from (1-L_B)=(1-L_C)(1-L_G)
G0 = 9.80665                           # m/s^2 standard gravity (CGPM 1901), for g*h/c estimate
PC = K.pc.value
MPC = 1e6 * PC
YEAR = 365.25 * 86400.0                # Julian year, s
V_CMB_DIPOLE = 369.82e3                # m/s, Planck 2018 I solar dipole (as quoted by Carr et al. 2022)
ARCSEC = np.pi / 180 / 3600

out = {}


def put(key, value, unit, how):
    out[key] = {"value": value, "unit": unit, "derived_from": how}
    if isinstance(value, float):
        print(f"{key:55s} {value: .6g} {unit:10s} <- {how}")
    else:
        print(f"{key:55s} {value} {unit:10s} <- {how}")


print("\n=== 1. Equivalent velocity of a 1e-3 A wavelength step  dv = c*dlam/lam")
for lam in [1215.67, 2000.0, 2796.35, 4000.0, 5000.0, 6562.8, 10000.0, 16000.0, 22000.0]:
    put(f"dv_per_1e-3A_at_{lam:g}A", C * 1e-3 / lam, "m/s", "c*1e-3/lambda")
for dv in [1000.0, 100.0, 60.0, 20.0, 10.0, 1.0, 0.1, 0.01]:
    put(f"dlam_at_5000A_for_{dv:g}m/s", 5000.0 * dv / C, "A", "lambda*dv/c, lambda=5000 A")
    put(f"dlam_at_2000A_for_{dv:g}m/s", 2000.0 * dv / C, "A", "lambda*dv/c, lambda=2000 A")

print("\n=== 2. Gravitational redshifts  v_g = GM/(R c)  (weak field)")
vg_sun_surf = GM_SUN / (R_SUN * C)
vg_sun_exact = C * (1.0 / np.sqrt(1 - 2 * GM_SUN / (R_SUN * C**2)) - 1.0)
put("vg_sun_surface_weakfield", vg_sun_surf, "m/s", "GM_sun/(R_sun c), IAU 2015 nominal")
put("vg_sun_surface_schwarzschild", vg_sun_exact, "m/s", "c[(1-2GM/Rc^2)^-1/2 - 1]")
dv_sun_pot_1au = GM_SUN / (AU * C)
dv_earth_pot = GM_EARTH / (R_EARTH_EQ * C)
dv_geoid = L_G * C
put("blueshift_observer_in_solar_potential_at_1AU", dv_sun_pot_1au, "m/s", "GM_sun/(1 au c)")
put("blueshift_observer_on_Earth_GM_over_R", dv_earth_pot, "m/s", "GM_E/(R_E c) (no rotation)")
put("blueshift_observer_on_geoid_LG_c", dv_geoid, "m/s", "L_G c = W0/c (incl. rotation)")
net_sun_at_earth = vg_sun_surf - dv_sun_pot_1au - dv_geoid
put("solar_gravitational_redshift_seen_on_Earth_static", net_sun_at_earth, "m/s",
    "GM/(R_sun c) - GM/(au c) - W0/c  (potentials only)")
v_orb = np.sqrt(GM_SUN / AU)
put("transverse_doppler_Earth_orbit_v2_over_2c", v_orb**2 / (2 * C), "m/s",
    "v^2/(2c), v=sqrt(GM_sun/au)=%.1f m/s" % v_orb)
put("solar_redshift_seen_on_Earth_incl_TD_of_observer", net_sun_at_earth - v_orb**2 / (2 * C), "m/s",
    "previous minus observer transverse-Doppler blueshift (Sun at right angles to orbit velocity)")
for h in [2400.0, 2635.0, 4200.0, 5050.0]:
    put(f"observatory_height_{h:g}m_extra_blueshift_gh_over_c", G0 * h / C, "m/s", "g h / c")
put("clock_rate_LC_times_c", L_C * C, "m/s", "mean TCG-TCB rate (observer orbit: GM/a + v^2/2)")
put("clock_rate_LB_times_c", L_B * C, "m/s", "mean TT/TDB-TCB rate: total observer-frame offset")
put("velocity_rescale_TDB_vs_TCB_at_30kms", 30e3 * L_B, "m/s", "v*L_B (unit rescaling only)")

def vg(m_msun, r_rsun):
    return GM_SUN * m_msun / (R_SUN * r_rsun * C)

print("  -- illustrative stars (inputs are illustrative M, R; v_g = 636.3 m/s * M/R)")
for name, m, r in [("Sun", 1.0, 1.0), ("K_dwarf_0.8M_0.8R", 0.8, 0.8), ("M_dwarf_0.2M_0.22R", 0.2, 0.22),
                   ("subgiant_1.2M_2R", 1.2, 2.0), ("giant_1.5M_10R", 1.5, 10.0),
                   ("supergiant_15M_500R", 15.0, 500.0),
                   ("WD_0.6M_0.0125R", 0.6, 0.0125), ("WD_1.0M_0.008R", 1.0, 0.008)]:
    put(f"vg_{name}", vg(m, r), "m/s", "GM/(Rc)")
# neutron star, strong field
m_ns, r_ns = 1.4, 12e3 / R_SUN
z_ns = 1 / np.sqrt(1 - 2 * GM_SUN * m_ns / (12e3 * C**2)) - 1
put("z_grav_neutron_star_1.4Msun_12km", z_ns, "", "(1-2GM/Rc^2)^-1/2 - 1")
# cluster
M_cl = 1e15
put("GM_over_Rc_cluster_1e15Msun_1Mpc", GM_SUN * M_cl / (MPC * C), "m/s",
    "GM/(Rc), M=1e15 Msun, R=1 Mpc (order of magnitude only)")

print("\n=== 3. Transverse Doppler / time dilation  v^2/(2c)")
for name, v in [("Earth_rotation_equator", OMEGA_EARTH * R_EARTH_EQ),
                ("Earth_orbit_mean", v_orb), ("star_vt_50kms", 50e3), ("star_vt_100kms", 100e3),
                ("halo_star_vt_300kms", 300e3), ("Sun_wrt_CMB_369.82kms", V_CMB_DIPOLE)]:
    put(f"TD_{name}", v**2 / (2 * C), "m/s", f"v^2/(2c), v={v:.1f} m/s")
put("Earth_rotation_speed_equator", OMEGA_EARTH * R_EARTH_EQ, "m/s", "omega_E * R_E")

print("\n=== 4. Aberration and direction errors in the barycentric projection")
put("max_error_using_apparent_dir_with_1st_order_formula", v_orb**2 / C, "m/s",
    "c*beta^2*sin^2(theta) <= v^2/c (n'.v - n.v)")
put("rv_error_per_arcsec_position_error", v_orb * ARCSEC, "m/s/arcsec", "v_orb * 1 arcsec")
put("rv_error_per_arcsec_rotation_term", OMEGA_EARTH * R_EARTH_EQ * ARCSEC, "m/s/arcsec",
    "v_rot * 1 arcsec")

print("\n=== 5. Solar Shapiro delay derivative  dv = (2GM/c^2) cot(eps/2) d(eps)/dt")
two_rs = 2 * GM_SUN / C**2
Omega = 2 * np.pi / (365.25636 * 86400)
put("2GM_sun_over_c2", two_rs, "m", "2GM/c^2")
for eps in [1.0, 5.0, 10.0, 30.0, 60.0, 90.0, 150.0]:
    put(f"shapiro_rv_max_elong_{eps:g}deg", two_rs * Omega / np.tan(np.radians(eps) / 2), "m/s",
        "Delta_t=-(2GM/c^3) ln(1-cos eps); |d eps/dt|<=2pi/sidereal yr")

print("\n=== 6. Velocity conventions for the same redshift")
for z in [0.001, 0.01, 0.023, 0.05, 0.1, 0.5, 1.0, 3.0]:
    vopt = C * z
    vrad = C * z / (1 + z)
    vrel = C * ((1 + z)**2 - 1) / ((1 + z)**2 + 1)
    put(f"z={z:g}: v_opt,v_rad,v_rel (km/s)", [round(vopt / 1e3, 3), round(vrad / 1e3, 3), round(vrel / 1e3, 3)],
        "km/s", "cz ; cz/(1+z) ; c((1+z)^2-1)/((1+z)^2+1)")

print("\n=== 7. Additive vs multiplicative redshift composition")
for zc in [0.01, 0.023, 0.05, 0.1, 0.3, 0.5, 1.0, 2.0]:
    for vp in [300e3, V_CMB_DIPOLE]:
        zp = vp / C     # first-order; exact SR line of sight: sqrt((1+b)/(1-b))-1
        zp_sr = np.sqrt((1 + vp / C) / (1 - vp / C)) - 1
        zmul = (1 + zc) * (1 + zp_sr) - 1
        zadd = zc + zp
        put(f"z_cos={zc:g}, v_pec={vp/1e3:.2f}km/s: z_mult - z_add", zmul - zadd, "",
            "(1+zc)(1+zp)-1 - (zc+v/c); zp exact SR")
        put(f"   ... as velocity c*dz/(1+z)", C * (zmul - zadd) / (1 + zmul), "m/s", "c dz/(1+z)")
# distance-modulus bias from a redshift error dz at fixed z (Planck18 d_L)
print("  -- distance modulus change dmu for dz=1e-4 and dz=1e-5 (Planck18)")
for z in [0.01, 0.023, 0.05, 0.1, 0.3, 0.5, 1.0]:
    for dz in [1e-4, 1e-5]:
        dmu = (Planck18.distmod(z + dz) - Planck18.distmod(z)).value
        put(f"dmu(z={z:g}, dz={dz:g})", dmu, "mag", "Planck18.distmod(z+dz)-distmod(z)")
# peculiar velocity scatter as distance error
for z in [0.01, 0.0233, 0.05, 0.1, 0.15]:
    for sv in [250e3, 150e3]:
        # Davis et al. 2011 eqs. A1 + A4 (empty-universe approximation):
        # sigma_z = (1+z) sigma_v/c ; sigma_mu = (5/ln10) (1+z)/(z(1+z/2)) sigma_z
        sz = (1 + z) * sv / C
        smu = 5 / np.log(10) * (1 + z) / (z * (1 + z / 2)) * sz
        put(f"sigma_mu_pec(z={z:g}, sigma_v={sv/1e3:.0f}km/s)", smu, "mag",
            "(5/ln10) (1+z)/(z(1+z/2)) (1+z) sigma_v/c  [Davis+2011 A1,A4]")
        put(f"  frac distance error", smu * np.log(10) / 5, "", "sigma_d/d = sigma_mu ln10/5")

print("\n=== 8. CMB-frame conversion magnitude")
for z in [0.01, 0.1, 0.5, 1.0]:
    # Sun moving towards apex with beta; for source exactly at apex
    b = V_CMB_DIPOLE / C
    g = 1 / np.sqrt(1 - b**2)
    for cth in [1.0, 0.0, -1.0]:
        zh = (1 + z) / (g * (1 + b * cth)) - 1         # observed (heliocentric) for z_CMB=z
        z_add_back = zh + b * cth                       # naive additive inversion
        put(f"z_CMB={z:g}, cos(theta)={cth:+g}: z_hel", zh, "", "(1+zC)/(gamma(1+beta cos th)) - 1")
        put(f"   additive z_CMB error", z_add_back - z, "", "z_hel + beta cos th - z_CMB")

print("\n=== 9. Redshift drift (Sandage-Loeb) for Planck18 (astropy)")
H0 = Planck18.H0.to(1 / u.s).value
put("H0_Planck18", Planck18.H0.value, "km/s/Mpc", "astropy Planck18")
put("H0_in_1/s", H0, "1/s", "")
put("c*H0 (velocity drift scale)", C * H0 * YEAR * 100, "cm/s/yr", "c H0")
drift = {}
for z in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]:
    Ez = Planck18.efunc(z)
    dzdt = H0 * ((1 + z) - Ez)                       # 1/s
    dvdt = C * H0 * (1 - Ez / (1 + z))               # m/s per s
    drift[z] = dict(
        E=Ez,
        dzdt_per_yr=dzdt * YEAR,
        dvdt_cm_s_per_yr=dvdt * YEAR * 100,
        dv_10yr_cm_s=dvdt * 10 * YEAR * 100,
        dv_20yr_cm_s=dvdt * 20 * YEAR * 100,
        dz_20yr=dzdt * 20 * YEAR,
        dlam_20yr_A_at_Lya=1215.67 * (1 + z) * dvdt * 20 * YEAR / C,
    )
    put(f"drift z={z:g}", {k: float(f"{v:.5g}") for k, v in drift[z].items()}, "",
        "dz/dt0=(1+z)H0-H(z); dv/dt=cH0[1-E(z)/(1+z)]")
# Liske et al. 2008 eq. 16 photon-noise scaling (all lines down to Ly-gamma):
# sigma_v = 1.35 (S/N/2370)^-1 (N_QSO/30)^-1/2 ((1+z)/5)^q cm/s, q=-1.7 (z<4), -0.9 (z>4)
def liske_sigma_v(sn, nq, zq):
    q = -1.7 if zq < 4 else -0.9
    return 1.35 * (sn / 2370.0) ** -1 * (nq / 30.0) ** -0.5 * ((1 + zq) / 5.0) ** q
for sn, nq, zq in [(2370, 30, 4.0), (3000, 30, 3.0), (3000, 10, 3.0), (10000, 1, 3.0), (1000, 30, 2.5)]:
    put(f"Liske08_sigma_v(S/N={sn},N={nq},z={zq})", liske_sigma_v(sn, nq, zq), "cm/s", "Liske+2008 eq.16")
zzero = brentq(lambda z: (1 + z) - Planck18.efunc(z), 0.5, 5)
put("drift_zero_crossing_z_Planck18", zzero, "", "(1+z)=E(z)")
zz = np.linspace(0.01, 3, 3000)
dv = C * H0 * (1 - Planck18.efunc(zz) / (1 + zz)) * YEAR * 100
put("drift_max_positive_z", float(zz[np.argmax(dv)]), "", "argmax dv/dt")
put("drift_max_positive_cm_s_per_decade", float(dv.max() * 10), "cm/s/decade", "")
# comparison cosmologies
for name, cos in [("EdS_H0=67.66", FlatLambdaCDM(H0=67.66, Om0=1.0)),
                  ("flat_noRad_Om0.3", FlatLambdaCDM(H0=70, Om0=0.3))]:
    for z in [1.0, 2.0, 3.0, 4.0]:
        h0 = cos.H0.to(1 / u.s).value
        put(f"drift {name} z={z:g}", C * h0 * (1 - cos.efunc(z) / (1 + z)) * 10 * YEAR * 100,
            "cm/s/decade", "cH0[1-E/(1+z)]*10yr")

print("\n=== 10. Recession velocity vs Doppler interpretations (Davis & Lineweaver 2004)")
for z in [0.1, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 1089.8]:
    Dc = Planck18.comoving_distance(z).to(u.m).value
    vrec = H0 * Dc
    vsr = C * ((1 + z)**2 - 1) / ((1 + z)**2 + 1)
    put(f"z={z:g}: v_rec(t0)/c, v_SR/c, cz/c", [round(float(vrec / C), 4), round(float(vsr / C), 4), z], "",
        "H0 D_C(z)/c ; SR Doppler ; z")
zH = brentq(lambda z: Planck18.comoving_distance(z).to(u.m).value * H0 - C, 0.5, 3)
put("z_where_v_rec_today_equals_c_Planck18", zH, "", "H0 D_C(z) = c")
# low-z cosmographic: v_rec = cz[1 - (1+q0)z/2 + ...]
q0 = Planck18.Om0 / 2 * (1 + 0) - Planck18.Ode0   # (radiation neglected)
put("q0_Planck18_approx", q0, "", "Om/2 - OL (radiation neglected)")
for z in [0.01, 0.0233, 0.05, 0.1, 0.15]:
    Dc = Planck18.comoving_distance(z).to(u.m).value
    put(f"z={z:g}: H0*D_C/(cz) - 1", H0 * Dc / (C * z) - 1, "", "fractional difference cz vs H0 D_C")
    dl = Planck18.luminosity_distance(z).to(u.m).value
    put(f"z={z:g}: H0*D_L/(cz) - 1", H0 * dl / (C * z) - 1, "", "fractional difference cz vs H0 D_L")

print("\n=== 11. Rest-frame interval scaling: observed line separation stretches by (1+z)")
for z in [0.5, 1.0, 2.0, 3.0]:
    put(f"z={z:g}: 1e-3 A rest-frame -> observed", 1e-3 * (1 + z), "A", "(1+z)*dlam_rest")

json.dump(out, open(sys.argv[1], "w"), indent=1, default=float)
print("\nwrote", sys.argv[1])
