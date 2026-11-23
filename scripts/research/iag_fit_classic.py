#!/usr/bin/env python
"""Measure solar-flux line-core positions of classic Fraunhofer lines in the IAG solar flux atlas
(Reiners et al. 2016, A&A 587, A65; VizieR J/A+A/587/A65/spvis, vacuum wavenumbers, corrected
for the Sun-observer Doppler shift, i.e. solar rest frame incl. gravitational redshift and
convective shifts).

Method (mirrors Reiners et al. 2016 Sect. 6): locate the flux minimum within +/- search_kms of
the expected position, then fit a 2nd-order polynomial to the points within +/- fit_kms of the
minimum; the vertex gives the line-core wavenumber.

Usage: python iag_fit_classic.py <nist_classic_dir> <snippet_out_dir> <result_csv>
"""
import csv, glob, math, os, sys, time, urllib.parse, urllib.request
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from airvac import vac_to_air

C = 299792.458  # km/s
NIST, SNIP, OUT = sys.argv[1:4]
os.makedirs(SNIP, exist_ok=True)

# label, NIST spectrum, lab AIR wl to pick from NIST file (+-0.002), search km/s, fit km/s
TARGETS = [
    ('C (H-alpha)', 'H I', '6562.795', 30, 10), ('F (H-beta)', 'H I', '4861.3250', 30, 10),
    ("G' (H-gamma)", 'H I', '4340.4634', 30, 10), ('h (H-delta)', 'H I', '4101.7346', 30, 10),
    ('a (Rowland line)', 'Fe I', '7164.4481', 5, 3),
    ('D2', 'Na I', '5889.95095', 8, 4), ('D1', 'Na I', '5895.92424', 8, 4),
    ('E2', 'Fe I', '5269.5370', 6, 4), ('E1 Fe', 'Fe I', '5270.3560', 4, 2.5), ('E1 Ca', 'Ca I', '5270.270', 3, 2.0),
    ('b1', 'Mg I', '5183.6042', 8, 4), ('b2', 'Mg I', '5172.6843', 8, 4), ('b4 Mg', 'Mg I', '5167.3216', 5, 3),
    ('b4 Fe', 'Fe I', '5167.4881', 5, 3), ('b3 Fe I', 'Fe I', '5168.8978', 5, 3), ('b3 Fe II', 'Fe II', '5169.0282', 5, 3),
    ('c', 'Fe I', '4957.5965', 5, 3), ('c (4957.30)', 'Fe I', '4957.2983', 5, 3),
    ('d/e (4668)', 'Fe I', '4668.1341', 5, 3), ('d/e (4383)', 'Fe I', '4383.5447', 6, 4),
    ('f (4325)', 'Fe I', '4325.7616', 6, 4), ('G Fe', 'Fe I', '4307.9020', 6, 4), ('G Ca', 'Ca I', '4307.744', 3, 2),
    ('g', 'Ca I', '4226.7276', 8, 4),
]


def clean(s):
    return s.replace('=', '').replace('"', '').strip()


def nist_lines(sp):
    rows = []
    for medium in ('air', 'vac'):
        for fn in glob.glob(os.path.join(NIST, f"{sp.replace(' ', '_')}_*_{medium}.csv")):
            lines = [l for l in open(fn).read().splitlines() if not l.startswith('# query')]
            rd = list(csv.reader(lines))
            hdr = [clean(h) for h in rd[0]]
            for r in rd[1:]:
                d = {h: clean(v) for h, v in zip(hdr, r)}
                d['_medium'] = medium
                rows.append(d)
    return rows


def lab_pair(sp, wl_air):
    rows = nist_lines(sp)
    wl_air = str(wl_air)
    air = [r for r in rows if r['_medium'] == 'air' and wl_air in (r.get('obs_wl_air(A)'), r.get('ritz_wl_air(A)'))]
    if not air:
        raise SystemExit(f'no NIST air match {sp} {wl_air}')
    a = air[0]
    wn = a['wn(cm-1)']
    # vacuum partner = same wavenumber string in the vacuum file
    keys = ('wn(cm-1)', 'Aki(s^-1)', 'Ek(cm-1)', 'J_i', 'J_k', 'term_k', 'line_ref')
    vac = [r for r in rows if r['_medium'] == 'vac' and all(r.get(k, '') == a.get(k, '') for k in keys)]
    assert len(vac) == 1, (sp, wl_air, len(vac))
    v = vac[0]
    # use the same kind (observed vs Ritz) in vacuum as was matched in air
    if a.get('obs_wl_air(A)') == wl_air and v.get('obs_wl_vac(A)'):
        lv = float(v.get('obs_wl_vac(A)'))
    else:
        lv = float(v.get('ritz_wl_vac(A)'))
    return a, v, lv


def fetch(sig_lo, sig_hi):
    url = ('https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?' +
           urllib.parse.urlencode({'-source': 'J/A+A/587/A65/spvis', '-out.max': 'unlimited',
                                   'Wave': f'{sig_lo:.5f}..{sig_hi:.5f}'}))
    for k in range(8):
        try:
            txt = urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'spectral-line-wiki/1.0'}), timeout=120).read().decode()
            n = sum(1 for l in txt.splitlines() if l[:1].isdigit() or l[:2].strip()[:1].isdigit())
            if n > 20:
                return url, txt
            print('empty/err response, retry', k)
        except Exception as e:
            print('retry', e)
        time.sleep(10 * (k + 1))
    raise SystemExit('fetch failed')


res = []
for label, sp, wl_air, search, fitw in TARGETS:
    a, v, lab_vac = lab_pair(sp, wl_air)
    sig_lab = 1e8 / lab_vac
    sig_exp = sig_lab / (1 + 0.5 / C)          # expect ~ +0.5 km/s redshift
    half = sig_exp * (search + 2) / C
    step = 0.0037665
    tag = f"{sp.replace(' ', '')}_{float(wl_air):.3f}"
    fn = os.path.join(SNIP, f'iag_spvis_{tag}.tsv')
    if not os.path.exists(fn) or os.path.getsize(fn) < 5000:
        url, txt = fetch(sig_exp - half, sig_exp + half)
        open(fn, 'w').write(f'# {url}\n' + txt)
        time.sleep(0.5)
    data = []
    for ln in open(fn):
        p = ln.split('\t')
        if len(p) >= 2:
            try:
                data.append((float(p[0]), float(p[1])))
            except ValueError:
                pass
    d = np.array(sorted(data))
    sig, f = d[:, 0], d[:, 1]
    vel = C * (sig_lab / sig - 1)            # km/s relative to lab (positive = redshift)
    m = np.abs(vel - 0.5) < search
    i0 = np.where(m)[0][np.argmin(f[m])]
    sel = np.abs(vel - vel[i0]) < fitw
    x = sig[sel] - sig[i0]
    cf = np.polyfit(x, f[sel], 2)
    sig_c = sig[i0] - cf[1] / (2 * cf[0])
    depth = 1 - np.polyval(cf, sig_c - sig[i0])
    wl_vac = 1e8 / sig_c
    v_shift = C * 1e3 * (wl_vac / lab_vac - 1)  # m/s
    res.append(dict(label=label, species=sp, lab_wl_air_A=wl_air,
                    lab_wl_vac_A=f'{lab_vac:.5f}', iag_sigma_cm1=f'{sig_c:.5f}', iag_wl_vac_A=f'{wl_vac:.4f}',
                    iag_wl_air_A=f'{float(vac_to_air(wl_vac)):.4f}', iag_minus_lab_mps=f'{v_shift:.0f}',
                    core_depth=f'{depth:.3f}', n_fit_pts=int(sel.sum()), fit_halfwidth_kms=fitw,
                    snippet=os.path.basename(fn)))
    print(res[-1])

with open(OUT, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(res[0]))
    w.writeheader(); w.writerows(res)
print('wrote', OUT)
