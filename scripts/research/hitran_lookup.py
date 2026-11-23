#!/usr/bin/env python
"""Parse HITRAN .par (160-char) files and list lines near given AIR wavelengths.

Usage: python hitran_lookup.py <hitran_dir> air_wl_A[:halfwidth_A] ...
Prints vacuum wavenumber, its HITRAN uncertainty code, vacuum and air (Peck & Reeder 1972)
wavelengths, line intensity S (cm/molecule at 296 K), pressure shift delta_air and quanta.
"""
import glob, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from airvac import vac_to_air, air_to_vac

IERR_NU = {0: '>=1 or unreported', 1: '0.1-1', 2: '0.01-0.1', 3: '0.001-0.01', 4: '1e-4-1e-3',
           5: '1e-5-1e-4', 6: '<1e-5', 7: '<1e-6', 8: '<1e-7', 9: '<1e-8'}


def parse_par(fn):
    out = []
    for ln in open(fn):
        if len(ln) < 160:
            continue
        out.append(dict(
            mol=int(ln[0:2]), iso=ln[2], nu=float(ln[3:15]), nu_str=ln[3:15].strip(), S=float(ln[15:25]),
            A=float(ln[25:35]), g_air=float(ln[35:40]), g_self=float(ln[40:45]), Elow=float(ln[45:55]),
            n_air=float(ln[55:59]), d_air=float(ln[59:67]), gq_up=ln[67:82].strip(), gq_lo=ln[82:97].strip(),
            lq_up=ln[97:112].strip(), lq_lo=ln[112:127].strip(), ierr=ln[127:133], iref=ln[133:145].strip(),
            file=os.path.basename(fn)))
    return out


def load(d):
    rows = []
    for fn in sorted(glob.glob(os.path.join(d, '*.par'))):
        rows += parse_par(fn)
    return rows


if __name__ == '__main__':
    rows = load(sys.argv[1])
    for arg in sys.argv[2:]:
        wl, hw = (arg.split(':') + ['0.08'])[:2]
        wl, hw = float(wl), float(hw)
        vac = float(air_to_vac(wl))
        print(f'== target air {wl:.3f} A (vac {vac:.3f}, sigma {1e8/vac:.4f} cm-1) +/-{hw}')
        sel = [r for r in rows if abs(1e8 / r['nu'] - vac) < hw]
        for r in sorted(sel, key=lambda r: -r['S'])[:8]:
            lv = 1e8 / r['nu']
            print(f"   nu={r['nu_str']:>12s} ierr_nu={r['ierr'][0]}({IERR_NU[int(r['ierr'][0])]}) vac={lv:.4f} air={float(vac_to_air(lv)):.4f} "
                  f"S={r['S']:.3e} dair={r['d_air']:+.6f} {r['gq_up']}<-{r['gq_lo']} | {r['lq_lo']} | iso{r['iso']} {r['file']}")
