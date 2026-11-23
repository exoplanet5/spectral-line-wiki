#!/usr/bin/env python
"""Derive OH X2Pi (Meinel) vibrational band origins from the HITRAN line list.

Input : HITRAN .par file for 16OH (global iso id 48) downloaded from https://hitran.org/lbl/api
Method: collect unique lower-state levels (Omega, v, J, e/f) with energies E'' (cm-1, relative to
        the lowest level); average the e/f Lambda-doublets; for every v fit the eigenvalues of the
        2x2 Hund's case (a) 2Pi Hamiltonian with parameters T_v, A_v, B_v, D_v, gamma_v
        (F1 = Omega 3/2 = lower eigenvalue for inverted OH, F2 = Omega 1/2 = upper) for J <= Jmax.
        Band origin nu0(v'-v'') = T_v' - T_v''.  vacuum wavelength = 1e8/nu0; air via Peck & Reeder (1972).
Cross-check: Huber & Herzberg (1979) / Maillard et al. (1976) G(v) polynomial (NIST WebBook), valid v <= 5.

usage: python -I oh_band_origins.py <hitran.par> <out.csv>
"""
import sys, csv
import numpy as np
from scipy.optimize import least_squares

par, out = sys.argv[1:3]
levels = {}
for ln in open(par):
    if len(ln) < 127:
        continue
    try:
        Elow = float(ln[45:55])
    except ValueError:
        continue
    gl = ln[82:97]
    loc = ln[112:127]
    st = gl.split()
    if len(st) != 2 or not st[0].startswith('X'):
        continue
    omega, v = st[0], int(st[1])
    # lower local quanta: Br(2) J''(F5.1) sym(2) F''(F5.1)
    try:
        J = float(loc[3:8]); sym = loc[8:10].strip()
    except ValueError:
        continue
    key = (omega, v, J, sym)
    levels.setdefault(key, []).append(Elow)

# merge hyperfine duplicates and e/f
lv = {}
for (om, v, J, sym), Es in levels.items():
    lv.setdefault((om, v, J), []).append(np.mean(Es))
rows = [(om, v, J, np.mean(E), len(E)) for (om, v, J), E in lv.items()]

def model(p, J, F1):
    """2x2 case-(a) 2Pi Hamiltonian (T, A, B, D, gamma), x = (J+1/2)^2 - 1.
    F1 (Omega=3/2 for inverted OH, A<0) = lower eigenvalue; F2 = upper eigenvalue."""
    T, A, B, D, g = p
    x = (J + 0.5) ** 2 - 1.0
    h33 = T + A / 2 + B * (x - 1) - D * ((x - 1) ** 2 + x)
    h11 = T - A / 2 + B * (x + 1) - D * ((x + 1) ** 2 + x) - g
    h13 = -np.sqrt(np.maximum(x, 0)) * (B - 2 * D * (x + 1) - g / 2)
    mean = (h33 + h11) / 2; half = np.sqrt(((h33 - h11) / 2) ** 2 + h13 ** 2)
    return np.where(F1, mean - half, mean + half)

res = {}
JMAX = 12.5
for v in range(0, 14):
    sel = [(om, J, E) for om, vv, J, E, n in rows if vv == v and J <= JMAX]
    if len(sel) < 8:
        continue
    J = np.array([s[1] for s in sel]); E = np.array([s[2] for s in sel])
    F1 = np.array([s[0] == 'X3/2' for s in sel])
    p0 = [E.min() + 100, -139.0, 18.5 - 0.7 * v, 1.9e-3, -0.1]
    fit = least_squares(lambda p: model(p, J, F1) - E, p0)
    rms = np.sqrt(np.mean(fit.fun ** 2))
    res[v] = dict(T=fit.x[0], A=fit.x[1], B=fit.x[2], D=fit.x[3], g=fit.x[4], rms=rms, n=len(sel))
    print(f'v={v:2d} n={len(sel):3d} T={fit.x[0]:12.4f} A={fit.x[1]:9.3f} B={fit.x[2]:8.4f} D={fit.x[3]:.3e} g={fit.x[4]:.4f} rms={rms:.4f} cm-1')

# Huber-Herzberg / Maillard 1976 polynomial (WebBook notes 21-22): valid v<=5
def G(v):
    x = v + 0.5
    return 3737.761 * x - 84.881 * x ** 2 + 0.5409 * x ** 3 - 0.02134 * x ** 4 - 0.00113 * x ** 5

def n_air(lam_vac_A):
    s2 = (1e4 / lam_vac_A) ** 2   # (1/um)^2
    return 1 + 1e-8 * (8060.51 + 2480990 / (132.274 - s2) + 17455.7 / (39.32957 - s2))

def vac2air(lv):
    # iterate: lambda_air = lambda_vac / n(lambda_vac) (Peck & Reeder index evaluated at vacuum wavelength as in NIST ASD)
    return lv / n_air(lv)

with open(out, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['band', 'dv', 'nu0_cm-1', 'nu0_fit_unc_cm-1', 'wl_vac_A', 'wl_air_A', 'nu0_HH_poly_cm-1', 'HH_minus_fit_cm-1', 'rms_upper', 'rms_lower'])
    for dv in range(2, 10):
        for vl in range(0, 10):
            vu = vl + dv
            if vu > 9 or vu not in res or vl not in res:
                continue
            nu = res[vu]['T'] - res[vl]['T']
            unc = np.hypot(res[vu]['rms'], res[vl]['rms'])
            lvac = 1e8 / nu
            hh = G(vu) - G(vl)
            w.writerow([f'{vu}-{vl}', dv, f'{nu:.3f}', f'{unc:.3f}', f'{lvac:.3f}', f'{vac2air(lvac):.3f}' if lvac < 20000 else '',
                        f'{hh:.2f}', f'{hh - nu:.2f}', f"{res[vu]['rms']:.4f}", f"{res[vl]['rms']:.4f}"])
print('written', out)
