#!/usr/bin/env python
"""Fetch NIST Atomic Spectra Database (ASD) line lists for every spectrum in the Handbook,
over the Handbook's wavelength span, so that each Handbook line can be matched to the current
ASD observed / Ritz wavelength and its stated uncertainty.

Main query: ASD option show_av=4 => vacuum below 1850 A, standard air from 1850 A upward (including
beyond 2 um).  The Handbook quotes vacuum below 2000 A and air above, so the two conventions differ in
1850-2000 A.  Converting ASD's rounded air strings back to vacuum there changes the last digit in about
one case in five, so a second query per spectrum with Handbook vacuum lines in 1850-2000 A asks ASD for
vacuum wavelengths everywhere (show_av=3) over 1848-2002 A; match_asd.py uses those vacuum strings
directly.  Wavenumbers (vacuum, cm^-1) are requested as an independent cross-check.

Output: data/raw/nist_asd/<El>_<ion>.csv          (raw ASD CSV, show_av=4, one file per spectrum)
        data/raw/nist_asd/<El>_<ion>_vac1850.csv  (raw ASD CSV, show_av=3, 1848-2002 A vacuum)
Existing non-empty files (or an earlier non-CSV answer saved as .err.html) are not fetched again;
queries are spaced by 1.5 s.
"""
import json, os, sys, time, urllib.parse, urllib.request

ROOT = os.path.join(os.path.dirname(__file__), '..')
OUTD = os.path.join(ROOT, 'data', 'raw', 'nist_asd')
os.makedirs(OUTD, exist_ok=True)

d = json.load(open(os.path.join(ROOT, 'data', 'handbook.json')))
span, vac1850 = {}, set()
for l in d['lines']:
    k = (l['el'], l['ion'])
    lo, hi = span.get(k, (1e9, 0))
    span[k] = (min(lo, l['wl']), max(hi, l['wl']))
    if l['medium'] == 'vacuum' and 1850 <= l['wl'] < 2000:
        vac1850.add(k)

BASE = 'https://physics.nist.gov/cgi-bin/ASD/lines1.pl?'
COMMON = dict(unit=0, format=2, line_out=0, output=0, show_obs_wl=1, show_calc_wl=1, unc_out=1,
              show_av=4, show_wn=1, intens_out='on', allowed_out=1, forbid_out=1, no_spaces=1,
              enrg_out='on', J_out='on', term_out='on', conf_out='on', A_out=0, bibrefs=1,
              order_out=0, remove_js='on', page_size=15)


def fetch(fn, q, label):
    if (os.path.exists(fn) and os.path.getsize(fn) > 200) or os.path.exists(fn + '.err.html'):
        return False            # already fetched (an .err.html means ASD answered without a line list, e.g. Np II)
    url = BASE + urllib.parse.urlencode(q)
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'spectral-line-wiki/1.0 (research script)'})
            txt = urllib.request.urlopen(req, timeout=180).read().decode('utf-8', 'replace')
            break
        except Exception as e:
            print('retry', label, e, file=sys.stderr); time.sleep(10 * (attempt + 1))
    else:
        print('FAILED', label); return True
    if not txt.startswith('obs_wl') and not txt.startswith('ritz') and 'wl' not in txt[:60]:
        open(fn + '.err.html', 'w').write(txt)
        print(f'{label}: NON-CSV response ({len(txt)} bytes)', flush=True)
    else:
        open(fn, 'w').write(txt)
        print(f'{label}: {txt.count(chr(10))} rows', flush=True)
    time.sleep(1.5)
    return True


todo = sorted(span, key=lambda k: (k[0], k[1]))
for i, (el, ion) in enumerate(todo):
    lo, hi = span[(el, ion)]
    q = dict(COMMON, spectra=f'{el} {ion}', low_w=f'{max(lo - 2, 1):.3f}', upp_w=f'{hi + 2:.3f}')
    fetch(os.path.join(OUTD, f'{el}_{ion}.csv'), q, f'[{i+1}/{len(todo)}] {el} {ion} {lo:.1f}-{hi:.1f} A')
for i, (el, ion) in enumerate(sorted(vac1850)):
    q = dict(COMMON, spectra=f'{el} {ion}', low_w='1848', upp_w='2002', show_av=3)
    fetch(os.path.join(OUTD, f'{el}_{ion}_vac1850.csv'), q, f'[vac {i+1}/{len(vac1850)}] {el} {ion} 1848-2002 A vacuum')
print('DONE')
