#!/usr/bin/env python
"""Merge the NIST Handbook strong lines with their NIST ASD matches into the compact element-line
table used by the page (data/build/elements.json).

Adopted wavelength: as decided by scripts/match_asd.py (record field 'adopt'): the ASD observed or Ritz
value for exact/close matches when ASD states an uncertainty that is eligible (ASD value with at least the
Handbook's decimals, or u <= 3 units of the Handbook's last digit) or, without a stated uncertainty, when
the ASD value has more decimals than the Handbook; the ASD value of a 'loose-unique' match (one ASD line in
the window, stated u <= one Handbook unit); the ASD blend value of an unresolved multiplet; otherwise the
Handbook value.  The adopted string is stored without ASD's '+'/'*' suffixes.

Uncertainty u (A, two significant figures, never rounded to 0) and 'us':
  us = 1  the uncertainty stated by NIST ASD for the adopted observed/Ritz value;
  us = 0  implied: 2.5 units of the last published digit of the adopted string.  NIST ASD (Lines help,
          "Uncertainties of Observed and Ritz Wavelengths"): "If such information is not available, the
          implied uncertainty is defined by the number of significant figures given: it is generally
          between 2.5 and 25 units of the least significant figure."  2.5 units is the lower end (ASD's
          "rule of 25"), so implied uncertainties are optimistic, not conservative.
Grades: A u <= 1 mA, B <= 5 mA, C <= 10 mA, D > 10 mA (1e-9 A tolerance for float noise only).

dx (mA): the other source minus the adopted value, set when they differ by more than twice their combined
uncertainty (each source's stated or implied u as above).  Computed for exact, close and loose-unique
matches; not for loose/none, and not for unresolved multiplets that keep the Handbook value (no single
ASD value represents the blend).

m: match class e exact, c close, u loose-unique, l loose, '' none.  bl: number of fine-structure
components of an unresolved multiplet (absent otherwise).

Native medium: vacuum below 2000 A, standard air above (Handbook convention).  The other-medium value is
computed with Peck & Reeder (1972) eq. 3, the dispersion formula that reproduces ASD's own air/vacuum
pairs (see notes); the vacuum wavenumber is 1e8/lambda_vac.

Output keys: 'elements' (per spectrum: ground = original text, conf/term/core/val = parsed ground state,
note, ie_cm, ie_eV), 'lines', 'stats' (integer counts, keys 'after:A', 'after_stated:A', 'after_implied:B',
'before:D', 'src:ASD Ritz', 'cls:exact', 'blend:centroid', 'disagree:2sigma', 'only_0.01A:implied' ...),
'diff_stats' (Handbook-vs-ASD differences by Handbook decimals: for adopted finer ASD values, and under
'finer_available' for every exact/close match with a finer ASD string, the definition of
scripts/research/precision_rationale.py), 'unstaged' (the strong lines printed without an ionization stage,
not in 'lines') and 'parse_notes'.  The statistics quoted in src/content/EL_METHOD.html and WHY_VERDICT.html
are printed.
"""
import json, math, os, re
from collections import Counter

import numpy as np

ROOT = os.path.join(os.path.dirname(__file__), '..')
H = json.load(open(os.path.join(ROOT, 'data', 'handbook.json')))
M = json.load(open(os.path.join(ROOT, 'data', 'handbook_asd.json')))
IMPLIED_UNITS = 2.5          # implied uncertainty, units of the last published digit
C_MS = 299792458.0
MCODE = {'exact': 'e', 'close': 'c', 'loose-unique': 'u', 'loose': 'l', 'none': ''}


def n_pr72(wl_vac):
    s2 = (1e4 / wl_vac) ** 2
    return 1 + 1e-8 * (8060.51 + 2480990 / (132.274 - s2) + 17455.7 / (39.32957 - s2))


def air_to_vac(wa):
    v = wa
    for _ in range(30):
        nv = wa * n_pr72(v)
        if abs(nv - v) < 1e-12:
            return nv
        v = nv
    return v


def dec_of(s):
    s = re.sub(r'[^\d.]', '', s or '')
    return len(s.split('.')[1]) if '.' in s else 0


def implied(dec):
    return IMPLIED_UNITS * 10.0 ** (-dec)


def sig2(x):
    """two significant figures; a positive value never becomes 0"""
    return float(f'{x:.2g}')


def grade(u):
    return 'A' if u <= 0.001 + 1e-9 else 'B' if u <= 0.005 + 1e-9 else 'C' if u <= 0.010 + 1e-9 else 'D'


def main():
    lines, mt = H['lines'], {m['i']: m for m in M['matches']}
    pers = {(p['el'], p['ion'], p['wl_str']): p for p in H['persistent']}
    out, stat = [], Counter()
    dstat = {2: [], 3: [], 4: []}          # (|adopted - Handbook| mA, |dv| m/s) for Handbook lines replaced by a finer ASD value
    for i, l in enumerate(lines):
        m = mt[i]
        adopt_asd = m.get('adopt') == 'asd'
        if adopt_asd:
            wstr = re.sub(r'[^\d.]', '', m['asd_wl_str'])
            w = float(wstr)
            unc_stated = m['asd_unc'] is not None
            unc = m['asd_unc'] if unc_stated else implied(dec_of(wstr))
            src = 'ASD ' + ('Ritz' if m['asd_kind'] == 'ritz' else 'obs')
        else:
            wstr, w, unc, src, unc_stated = l['wl_str'], l['wl'], implied(l['dec']), 'Handbook', False
        med = l['medium']
        if med == 'vacuum':
            vac, air = w, None
        else:
            air, vac = w, air_to_vac(w)
        p = pers.get((l['el'], l['ion'], l['wl_str']))
        rec = {
            'Z': l['Z'], 'el': l['el'], 'ion': l['ion'],
            'hb': l['wl_str'], 'med': 'v' if med == 'vacuum' else 'a',
            'w': wstr, 'u': sig2(unc), 'us': 1 if unc_stated else 0, 'src': src,
            'air': round(air, 5) if air else None, 'vac': round(vac, 5),
            'wn': round(1e8 / vac, 4),
            'int': l['int'], 'fl': l['flags'], 'ref': l['ref'],
            'm': MCODE[m['cls']],
            'dl': m.get('dl_mA'),
        }
        if m.get('blend'):
            rec['bl'] = m['blend']
            stat[('blend', m['blend_kind'])] += 1
        if m['cls'] not in ('none', 'loose') and not (m.get('blend') and not adopt_asd):
            # do the two sources disagree by more than twice their combined (stated or implied) uncertainty?
            if adopt_asd:
                other, u_other = l['wl'], implied(l['dec'])
            else:
                ostr = re.sub(r'[^\d.]', '', m['asd_wl_str'])
                other = float(ostr)
                u_other = m['asd_unc'] if m['asd_unc'] is not None else implied(dec_of(ostr))
            if abs(other - w) > 2 * math.hypot(unc, u_other) + 1e-9:
                rec['dx'] = round((other - w) * 1000, 3)
                stat[('disagree', '2sigma')] += 1
                stat[('disagree_cls', m['cls'])] += 1
        if m['cls'] != 'none':
            rec.update({'ai': m.get('asd_int', ''), 'lo': m.get('lower', ''), 'up': m.get('upper', ''),
                        'A': m.get('Aki'), 'acc': m.get('acc', ''), 'Ei': m.get('Ei', ''), 'Ek': m.get('Ek', ''),
                        'aw': m.get('asd_wl_str', ''), 'ak': m.get('asd_kind', ''), 'au': m.get('asd_unc'),
                        'ty': m.get('type', ''), 'lr': m.get('asd_line_ref', '')})
        if p:
            rec['P'] = 1
            if p.get('levels'):
                rec['lv'] = p['levels']
            if p.get('Aki_1e8'):
                rec['hA'] = p['Aki_1e8']
        out.append(rec)
        # ---- statistics
        g = grade(rec['u'])          # from the stored (2-significant-figure) value, as the page does
        stat[('after', g)] += 1
        stat[('after_stated' if unc_stated else 'after_implied', g)] += 1
        stat[('before', grade(implied(l['dec'])))] += 1
        stat[('src', src)] += 1
        stat[('cls', m['cls'])] += 1
        stat[('dec_hb', min(l['dec'], 5))] += 1
        stat[('dec_after', min(dec_of(wstr), 5))] += 1
        if not unc_stated and dec_of(wstr) <= 2:
            stat[('only_0.01A', 'implied')] += 1           # adopted value still has <= 2 decimals, no stated u
        if adopt_asd and l['dec'] in dstat and dec_of(wstr) > l['dec']:
            dl = (w - l['wl']) * 1e3
            dstat[l['dec']].append((abs(dl), abs(dl) / l['wl'] * C_MS / 1e3))
    diff = {}
    for d_, v in dstat.items():
        a = np.array(v)
        diff[str(d_)] = {'n': len(v), 'median_abs_mA': round(float(np.median(a[:, 0])), 3),
                         'p68_abs_mA': round(float(np.percentile(a[:, 0], 68)), 3),
                         'p95_abs_mA': round(float(np.percentile(a[:, 0], 95)), 3),
                         'median_abs_m_s': round(float(np.median(a[:, 1])), 1),
                         'p95_abs_m_s': round(float(np.percentile(a[:, 1], 95)), 1)}
    diff['note'] = ('Handbook lines with d decimals whose adopted value is an ASD value with more decimals: '
                    '|adopted - Handbook|; velocity at each line\'s own wavelength')
    # same comparison as scripts/research/precision_rationale.py (quoted in WHY_VERDICT / st_effects): every
    # exact/close match whose shown ASD row has an observed or Ritz string with more decimals than the Handbook
    # (the finer one; ties -> smaller stated u), whether adopted or not
    fin = {2: [], 3: [], 4: []}
    for i, l in enumerate(lines):
        m = mt[i]
        if l['dec'] not in fin or m.get('cls') not in ('exact', 'close'):
            continue
        c_ = []
        for key, ukey in (('asd_obs', 'asd_unc_obs'), ('asd_ritz', 'asd_unc_ritz')):
            mm = re.match(r'^\[?\(?(\d+)\.(\d+)', str(m.get(key) or ''))
            if mm and len(mm.group(2)) > l['dec']:
                c_.append((-len(mm.group(2)), m.get(ukey) if m.get(ukey) is not None else 9e9,
                           float(mm.group(1) + '.' + mm.group(2))))
        if c_:
            dl = abs(l['wl'] - sorted(c_)[0][2]) * 1e3
            fin[l['dec']].append((dl, dl / l['wl'] * C_MS / 1e3))
    diff['finer_available'] = {}
    for d_, v in fin.items():
        a = np.array(v)
        diff['finer_available'][str(d_)] = {
            'n': len(v), 'median_abs_mA': round(float(np.median(a[:, 0])), 3),
            'p68_abs_mA': round(float(np.percentile(a[:, 0], 68)), 3), 'p95_abs_mA': round(float(np.percentile(a[:, 0], 95)), 3),
            'median_abs_m_s': round(float(np.median(a[:, 1])), 1), 'p95_abs_m_s': round(float(np.percentile(a[:, 1], 95)), 1)}
    for k in sorted(stat, key=str):
        print(k, stat[k])
    print('Handbook vs finer ASD value:', json.dumps(diff))
    els = []
    for e in H['elements']:
        rec = {k: e[k] for k in ('Z', 'symbol', 'name', 'weight')}
        rec['spectra'] = {ion: {k: s[k] for k in ('ground', 'conf', 'term', 'core', 'val', 'note', 'ie_cm', 'ie_eV') if k in s}
                          for ion, s in e['spectra'].items()}
        els.append(rec)
    os.makedirs(os.path.join(ROOT, 'data', 'build'), exist_ok=True)
    json.dump({'elements': els, 'lines': out, 'stats': {f'{a}:{b}': c for (a, b), c in stat.items()},
               'diff_stats': diff,
               'unstaged': [{k: u[k] for k in ('Z', 'el', 'wl_str', 'medium', 'int', 'flags', 'ref')} for u in H.get('unstaged', [])],
               'parse_notes': H.get('parse_notes', {})},
              open(os.path.join(ROOT, 'data', 'build', 'elements.json'), 'w'), separators=(',', ':'))
    print('size kB', os.path.getsize(os.path.join(ROOT, 'data', 'build', 'elements.json')) / 1e3)


if __name__ == '__main__':
    main()
