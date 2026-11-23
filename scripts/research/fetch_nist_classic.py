#!/usr/bin/env python
"""Fetch NIST ASD line data (air and vacuum) in narrow windows around the classic
Fraunhofer atomic lines and the optical-glass reference lamp lines.

Usage: python fetch_nist_classic.py <outdir>
Writes <outdir>/<tag>_air.csv (show_av=2: vac<2000 A, air 2000-20000 A, vac>20000 A)
and    <outdir>/<tag>_vac.csv (show_av=3: vacuum at all wavelengths), raw ASD CSV.
"""
import os, sys, time, urllib.parse, urllib.request

OUT = sys.argv[1]
os.makedirs(OUT, exist_ok=True)

# (spectrum, low_A, high_A) -- windows in Angstrom (air for show_av=2; the same numbers are
# used for the vacuum query but widened by +3 A to cover the air->vac shift)
WINDOWS = [
    ('H I', 6561.5, 6564.0), ('H I', 4860.5, 4862.5), ('H I', 4339.8, 4341.2), ('H I', 4101.0, 4102.5),
    ('H I', 3969.5, 3970.5),
    ('Na I', 5888.5, 5897.0), ('He I', 5875.0, 5876.5),
    ('Fe I', 5268.5, 5271.0), ('Ca I', 5269.5, 5271.0),
    ('Mg I', 5166.5, 5184.5), ('Fe I', 5166.5, 5169.5), ('Fe II', 5168.5, 5169.5),
    ('Fe I', 4956.8, 4958.2), ('Fe I', 4667.6, 4668.6), ('Fe I', 4383.0, 4384.0), ('Fe I', 4325.2, 4326.3),
    ('Fe I', 4307.3, 4308.4), ('Ca I', 4307.3, 4308.4), ('Ti II', 4307.3, 4308.4),
    ('Ca I', 4226.2, 4227.3), ('Ca II', 3968.0, 3969.0), ('Ca II', 3933.2, 3934.2),
    ('Fe I', 3819.9, 3821.0), ('Fe I', 3727.1, 3728.2), ('Fe I', 3580.7, 3581.7),
    ('Fe I', 3440.0, 3441.5), ('Ti II', 3360.7, 3361.7), ('Fe I', 3360.7, 3361.7),
    ('Fe I', 3286.2, 3287.3), ('Ca II', 3178.8, 3181.8), ('Fe I', 3099.4, 3101.2),
    ('Fe I', 3047.1, 3048.1), ('Fe I', 3020.1, 3021.6), ('Fe I', 2993.9, 2995.0), ('Ni I', 2993.9, 2995.0),
    ('Fe I', 2947.4, 2948.4), ('Fe I', 7163.9, 7165.0), ('Si I', 7164.3, 7165.0),
    # optical-glass / lamp reference lines that re-use Fraunhofer-style letters
    ('Hg I', 5460.2, 5461.3), ('Hg I', 4358.0, 4358.8), ('Hg I', 4046.2, 4047.0), ('Hg I', 3650.0, 3650.4),
    ('Hg I', 10139.0, 10140.5), ('Cd I', 4799.5, 4800.4), ('Cd I', 6438.0, 6439.0), ('He I', 7065.0, 7065.8),
    ('Cs I', 8520.5, 8521.8), ('He I', 7065.0, 7066.0),
]

BASE = 'https://physics.nist.gov/cgi-bin/ASD/lines1.pl?'
COMMON = dict(unit=0, format=2, line_out=0, output=0, show_obs_wl=1, show_calc_wl=1, unc_out=1,
              show_wn=1, intens_out='on', allowed_out=1, forbid_out=1, no_spaces=1,
              enrg_out='on', J_out='on', term_out='on', conf_out='on', A_out=0, bibrefs=1,
              order_out=0, remove_js='on', page_size=15)

def get(url):
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'spectral-line-wiki/1.0 (research script)'})
            return urllib.request.urlopen(req, timeout=120).read().decode('utf-8', 'replace')
        except Exception as e:
            print('retry', e, file=sys.stderr); time.sleep(5 * (attempt + 1))
    return None

seen = set()
for sp, lo, hi in WINDOWS:
    tag = f"{sp.replace(' ', '_')}_{lo:.1f}_{hi:.1f}"
    if tag in seen:
        continue
    seen.add(tag)
    for medium, av, dlo, dhi in (('air', 2, 0.0, 0.0), ('vac', 3, 0.0, 3.5)):
        fn = os.path.join(OUT, f'{tag}_{medium}.csv')
        if os.path.exists(fn) and os.path.getsize(fn) > 100:
            continue
        q = dict(COMMON, spectra=sp, low_w=f'{lo + dlo:.2f}', upp_w=f'{hi + dhi:.2f}', show_av=av)
        url = BASE + urllib.parse.urlencode(q)
        txt = get(url)
        if txt is None:
            print('FAILED', tag, medium); continue
        if not txt.lstrip().startswith(('obs_wl', 'ritz_wl')):
            open(fn + '.err.html', 'w').write(txt)
            print(tag, medium, 'NON-CSV', len(txt))
        else:
            open(fn, 'w').write(f'# query: {url}\n' + txt)
            print(tag, medium, txt.count('\n') - 1, 'lines')
        time.sleep(1.0)
print('DONE')
