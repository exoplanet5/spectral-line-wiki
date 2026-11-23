#!/usr/bin/env python
"""Compact Chebyshev representation of the Earth (geocentre) and Sun barycentric state from
JPL DE440s, for the in-browser barycentric-correction calculator.

Fits, per fixed-length segment in TDB, Chebyshev series to
  Earth velocity (km/s)   -- high precision (target < 0.1 mm/s)
  Earth position (km)     -- for BJD_TDB / Roemer delay (target < 1 km)
  Sun position  (km)      -- for heliocentric correction & solar potential
  Sun velocity  (km/s)
Coefficients are stored as integers in fixed units to keep the JSON small.

Output: data/ephem_cheb.json  (+ printed max fit errors versus DE440s on a dense grid)
"""
import json, os
import numpy as np
from jplephem.spk import SPK

ROOT = os.path.join(os.path.dirname(__file__), '..')
K = SPK.open(os.path.expanduser('~/JPLDE/de440s.bsp'))

JD0 = 2447892.5      # 1990-01-01 00:00 TDB
JD1 = 2476102.5      # 2067-03-?? (covers through 2066)


def earth(jd):
    p1, v1 = K[0, 3].compute_and_differentiate(jd)
    p2, v2 = K[3, 399].compute_and_differentiate(jd)
    return p1 + p2, (v1 + v2) / 86400.0          # km, km/s


def sun(jd):
    p, v = K[0, 10].compute_and_differentiate(jd)
    return p, v / 86400.0


def fit(func, comp, seg, deg, scale):
    """Return integer coefficient array [nseg][3][deg+1] in units of `scale`."""
    nseg = int(np.ceil((JD1 - JD0) / seg))
    n = deg + 1
    # Chebyshev-Gauss nodes
    k = np.arange(n)
    x = np.cos(np.pi * (k + 0.5) / n)
    T = np.polynomial.chebyshev.chebvander(x, deg)    # n x n
    Tinv = np.linalg.inv(T)
    out = []
    for i in range(nseg):
        a = JD0 + i * seg
        jd = a + (x + 1) * seg / 2
        val = func(jd)[comp]                         # 3 x n
        c = (Tinv @ val.T).T                         # 3 x n
        out.append(np.round(c / scale).astype(np.int64))
    return np.array(out)


def evaluate(coef, seg, scale, jd):
    i = np.floor((jd - JD0) / seg).astype(int)
    x = 2 * (jd - (JD0 + i * seg)) / seg - 1
    c = coef[i] * scale                              # m x 3 x n
    res = np.empty((len(jd), 3))
    for j in range(3):
        res[:, j] = [np.polynomial.chebyshev.chebval(xx, cc) for xx, cc in zip(x, c[:, j, :])]
    return res


def main():
    cfg = {
        'earth_vel': dict(func=earth, comp=1, seg=16.0, deg=12, scale=1e-7),   # 1e-7 km/s = 0.1 mm/s
        'earth_pos': dict(func=earth, comp=0, seg=16.0, deg=10, scale=1e-1),   # 100 m
        'sun_pos':   dict(func=sun,   comp=0, seg=256.0, deg=10, scale=1.0),
        'sun_vel':   dict(func=sun,   comp=1, seg=64.0, deg=10, scale=1e-8),
    }
    rng = np.random.default_rng(1)
    jd_test = np.sort(rng.uniform(JD0, JD1 - 1, 4000))
    out = {'source': 'JPL DE440s (Park et al. 2021, AJ 161, 105) via jplephem; ICRF/ICRS axes, '
                     'Solar-System-barycentric, TDB',
           'jd0_tdb': JD0, 'jd1_tdb': JD1}
    for name, c in cfg.items():
        coef = fit(c['func'], c['comp'], c['seg'], c['deg'], c['scale'])
        ev = evaluate(coef, c['seg'], c['scale'], jd_test)
        truth = np.array([c['func'](j)[c['comp']] for j in jd_test])
        err = np.abs(ev - truth).max()
        print(f'{name:10s} nseg={coef.shape[0]} deg={c["deg"]} max|err|={err:.3e} '
              f'({"km/s" if "vel" in name else "km"})')
        out[name] = {'seg_days': c['seg'], 'deg': c['deg'], 'scale': c['scale'],
                     'coef': coef.reshape(coef.shape[0], -1).tolist(), 'max_err': float(err)}
    # test vectors for the JS port
    tv = []
    for jd in np.linspace(JD0 + 10, JD1 - 10, 25):
        pe, ve = earth(jd); ps, vs = sun(jd)
        tv.append({'jd': jd, 'earth_pos': pe.tolist(), 'earth_vel': ve.tolist(),
                   'sun_pos': ps.tolist(), 'sun_vel': vs.tolist()})
    out['test_vectors'] = tv
    fn = os.path.join(ROOT, 'data', 'ephem_cheb.json')
    json.dump(out, open(fn, 'w'), separators=(',', ':'))
    print('wrote', fn, os.path.getsize(fn) / 1e3, 'kB')


if __name__ == '__main__':
    main()
