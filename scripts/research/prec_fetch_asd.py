#!/usr/bin/env python
"""Fetch small NIST ASD windows used in precision_rationale.md (isotope / fine-structure examples).
Usage: python prec_fetch_asd.py OUTDIR
Output CSVs are untrusted downloaded data; parse them with python -I.
show_av=2: vacuum < 2000 A, standard air 2000-20000 A, vacuum > 20000 A.
"""
import os, sys, time, urllib.parse, urllib.request

BASE = 'https://physics.nist.gov/cgi-bin/ASD/lines1.pl?'
COMMON = dict(unit=0, format=2, line_out=0, output=0, show_obs_wl=1, show_calc_wl=1, unc_out=1,
              show_av=2, show_wn=1, intens_out='on', allowed_out=1, forbid_out=1, no_spaces=1,
              enrg_out='on', J_out='on', term_out='on', conf_out='on', A_out=0, bibrefs=1,
              order_out=0, remove_js='on', page_size=15)
WINDOWS = [  # (spectrum, low A, high A, tag)
    ('Li I', 6707.0, 6709.0, 'Li_I_6708'),
    ('H I', 6562.0, 6564.0, 'H_I_Halpha'),
    ('D I', 6560.0, 6563.0, 'D_I_Halpha'),
    ('Hg I', 5460.0, 5461.5, 'Hg_I_5461'),
    ('Ba II', 4553.5, 4554.5, 'Ba_II_4554'),
    ('Ca II', 8541.5, 8542.6, 'Ca_II_8542'),
    ('Fe I', 5250.0, 5251.0, 'Fe_I_5250'),
]
out = sys.argv[1]
os.makedirs(out, exist_ok=True)
for sp, lo, hi, tag in WINDOWS:
    fn = os.path.join(out, tag + '.csv')
    if os.path.exists(fn) and os.path.getsize(fn) > 200:
        continue
    q = dict(COMMON, spectra=sp, low_w=f'{lo:.3f}', upp_w=f'{hi:.3f}')
    url = BASE + urllib.parse.urlencode(q)
    req = urllib.request.Request(url, headers={'User-Agent': 'spectral-line-wiki-research/0.1'})
    txt = urllib.request.urlopen(req, timeout=180).read().decode('utf-8', 'replace')
    open(fn, 'w').write(txt)
    open(fn + '.url', 'w').write(url + '\n')
    print(tag, len(txt), 'bytes', txt.count('\n'), 'lines')
    time.sleep(2)
