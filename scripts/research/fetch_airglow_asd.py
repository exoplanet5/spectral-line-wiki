#!/usr/bin/env python
"""Download NIST ASD line data around each airglow target line.

Two queries per wavelength window:
  av2 : show_av=2 (vacuum < 2000 A, standard air 2000-20000 A, vacuum > 20000 A) + show_wn=1
  av3 : show_av=3 (vacuum, all wavelengths)   [ASD rejects show_wn together with show_av=3]
so that both ASD's own air value and ASD's own vacuum value are on disk (no local conversion).

Usage: python fetch_airglow_asd.py <outdir>
"""
import os, sys, time, urllib.parse, urllib.request
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from airglow_targets import T

OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)
BASE = 'https://physics.nist.gov/cgi-bin/ASD/lines1.pl?'
COMMON = dict(unit=0, format=2, line_out=0, output=0, show_obs_wl=1, show_calc_wl=1, unc_out=1,
              intens_out='on', allowed_out=1, forbid_out=1, no_spaces=1, enrg_out='on', J_out='on',
              term_out='on', conf_out='on', A_out=0, f_out='on', bibrefs=1, order_out=0,
              remove_js='on', page_size=15)

# merge targets into windows per spectrum
win = {}
for t in T:
    lo, hi = t['nominal'] - t['tol'] - 1.0, t['nominal'] + t['tol'] + 1.0
    win.setdefault(t['spec'], []).append([lo, hi])
jobs = []
for spec, ws in win.items():
    ws.sort()
    merged = [ws[0]]
    for lo, hi in ws[1:]:
        if lo <= merged[-1][1] + 2:
            merged[-1][1] = max(hi, merged[-1][1])
        else:
            merged.append([lo, hi])
    for lo, hi in merged:
        jobs.append((spec, lo, hi))

def get(q, fn):
    if os.path.exists(fn) and os.path.getsize(fn) > 100:
        return 'cached'
    url = BASE + urllib.parse.urlencode(q)
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'spectral-line-wiki/1.0 (research script)'})
            txt = urllib.request.urlopen(req, timeout=120).read().decode('utf-8', 'replace')
            break
        except Exception as e:
            print('retry', fn, e, file=sys.stderr); time.sleep(8 * (attempt + 1))
    else:
        return 'FAILED'
    if txt.lstrip().startswith('<'):
        open(fn + '.err.html', 'w').write(txt)
        return 'NON-CSV'
    open(fn, 'w').write(txt)
    with open(os.path.join(OUT, 'queries.log'), 'a') as f:
        f.write(f'{os.path.basename(fn)}\t{url}\n')
    time.sleep(1.0)
    return f'{txt.count(chr(10))} rows'

for spec, lo, hi in jobs:
    tag = f"{spec.replace(' ', '_')}_{lo:.2f}_{hi:.2f}"
    q2 = dict(COMMON, spectra=spec, low_w=f'{lo:.3f}', upp_w=f'{hi:.3f}', show_av=2, show_wn=1)
    if lo > 1990:
        vlo, vhi = lo * 1.00026 - 1.0, hi * 1.00030 + 1.0
    else:
        vlo, vhi = lo, hi
    q3 = dict(COMMON, spectra=spec, low_w=f'{vlo:.3f}', upp_w=f'{vhi:.3f}', show_av=3)
    r2 = get(q2, os.path.join(OUT, tag + '_av2.csv'))
    r3 = get(q3, os.path.join(OUT, tag + '_av3.csv'))
    print(spec, f'{lo:.2f}-{hi:.2f}', 'av2:', r2, '| av3:', r3, flush=True)
print('DONE')
