#!/usr/bin/env python
"""Ephemeris-derived magnitudes for the barycentric / heliocentric chain,
using the shared JPL DE440s kernel (~/JPLDE/de440s.bsp) through astropy.

Derived quantities (no catalogue values):
  * |v_sun| w.r.t. the solar-system barycentre (heliocentric vs barycentric
    frames differ by up to this, projected on the line of sight)
  * |v_earth| w.r.t. SSB (max/min over the year)
  * |v_earth - v_EMB| (lunar reflex of the Earth)
  * the observer-clock factor  Phi_sun/c^2 + Phi_E/c^2 + v^2/(2c^2)  expressed
    in m/s and its annual modulation (cf. Lindegren & Dravins 2003 eq. 42)

Usage:  ~/.venvs/astro313/bin/python relcos_ephemeris.py OUT.json
"""
import os
import json
import sys

import numpy as np
from astropy import units as u
from astropy.coordinates import get_body_barycentric_posvel, solar_system_ephemeris
from astropy.time import Time

C = 299_792_458.0
GM_SUN = 1.3271244e20
L_G = 6.969290134e-10   # Earth (geoid) potential term, IAU 2000 B1.9

solar_system_ephemeris.set(os.path.expanduser("~/JPLDE/de440s.bsp"))

t = Time("2000-01-01", scale="tdb") + np.arange(0, 40 * 365.25, 2.0) * u.day
ps, vs = get_body_barycentric_posvel("sun", t)
pe, ve = get_body_barycentric_posvel("earth", t)
pm, vm = get_body_barycentric_posvel("moon", t)
vsun = vs.norm().to(u.m / u.s).value
vearth = ve.norm().to(u.m / u.s).value
# Earth-Moon barycentre (mass ratio from DE440 GM values used by astropy is not exposed; use EMRAT)
from jplephem.spk import SPK
_cm = SPK.open(os.path.expanduser("~/JPLDE/de440s.bsp")).comments()
EMRAT = float(_cm.split("EMRAT", 1)[1].split()[0].replace("D", "E"))   # read from kernel header
vemb = (ve * EMRAT + vm) / (EMRAT + 1)
dv_lunar = (ve - vemb).norm().to(u.m / u.s).value
r_es = (pe - ps).norm().to(u.m).value
clock = GM_SUN / r_es / C + L_G * C + vearth**2 / (2 * C)   # m/s equivalent, Sun+Earth potential + Earth TD

out = {
    "sun_wrt_SSB_speed_m_s": {"min": float(vsun.min()), "max": float(vsun.max()), "mean": float(vsun.mean()),
                              "span": "2000-2040 TDB, DE440s"},
    "earth_wrt_SSB_speed_m_s": {"min": float(vearth.min()), "max": float(vearth.max()), "mean": float(vearth.mean())},
    "earth_wrt_EMB_speed_m_s": {"min": float(dv_lunar.min()), "max": float(dv_lunar.max()),
                                "note": f"EMRAT={EMRAT!r} read from DE440s comment header"},
    "observer_clock_term_m_s": {"min": float(clock.min()), "max": float(clock.max()), "mean": float(clock.mean()),
                                "peak_to_peak": float(clock.max() - clock.min()),
                                "formula": "c*[GM_sun/(r c^2) + L_G + v_E^2/(2c^2)]  (Earth-centre speed, geoid potential)"},
}
for k, v in out.items():
    print(k, v)
json.dump(out, open(sys.argv[1], "w"), indent=1)
