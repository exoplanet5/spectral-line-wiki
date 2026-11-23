#!/usr/bin/env python
"""Fetch NIST ASD (v5.12) observed/Ritz AIR wavelengths (show_av=2) for every species that occurs
in the Allende Prieto & Garcia Lopez (1998) solar line catalogue, over 3900-8000 A, so that each
solar line can be paired with a laboratory wavelength.

Usage: python fetch_nist_solar_species.py <apgl98_table1.tsv> <outdir>
"""
import os, sys, time, urllib.parse, urllib.request

tsv, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
ions = set()
for ln in open(tsv):
    p = ln.rstrip('\n').split('\t')
    if len(p) >= 6 and p[0].strip().isdigit():
        ions.add(p[5].strip())
ROMAN = {'1': 'I', '2': 'II', '3': 'III'}
specs = []
for code in sorted(i for i in ions if i):
    el, st = code[:-1].strip(), code[-1]
    specs.append(f'{el.capitalize()} {ROMAN[st]}')
print(len(specs), specs)
BASE = 'https://physics.nist.gov/cgi-bin/ASD/lines1.pl?'
COMMON = dict(unit=0, format=2, line_out=0, output=0, show_obs_wl=1, show_calc_wl=1, unc_out=1,
              show_wn=1, intens_out='on', allowed_out=1, forbid_out=1, no_spaces=1, enrg_out='on',
              J_out='on', term_out='on', conf_out='on', A_out=0, bibrefs=1, order_out=0,
              remove_js='on', page_size=15, show_av=2, low_w='3900', upp_w='8000')
for sp in specs:
    fn = os.path.join(OUT, sp.replace(' ', '_') + '_3900_8000_air.csv')
    if os.path.exists(fn) and os.path.getsize(fn) > 200:
        continue
    url = BASE + urllib.parse.urlencode(dict(COMMON, spectra=sp))
    txt = None
    for k in range(4):
        try:
            txt = urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'spectral-line-wiki/1.0 (research script)'}), timeout=180).read().decode('utf-8', 'replace')
            break
        except Exception as e:
            print('retry', sp, e); time.sleep(8 * (k + 1))
    if txt is None:
        print('FAILED', sp); continue
    if txt.lstrip().startswith(('obs_wl', 'ritz_wl')):
        open(fn, 'w').write(f'# query: {url}\n' + txt)
        print(sp, txt.count('\n'), 'rows', flush=True)
    else:
        open(fn + '.err.html', 'w').write(txt)
        print(sp, 'NON-CSV', flush=True)
    time.sleep(1.5)
print('DONE')
