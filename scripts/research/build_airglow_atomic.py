#!/usr/bin/env python
"""Build data/curated/airglow_atomic.csv from the ASD selection + literature metadata.

usage: python -I build_airglow_atomic.py <asd_selected.json> <scripts_dir> <out_csv> <refs_csv>

Wavelength policy
-----------------
* Every wavelength string is copied verbatim from NIST ASD v5.12 (no local rounding/padding).
* wl_vac_A : ASD vacuum value (show_av=3 query).  wl_air_A : ASD standard-air value (show_av=2
  query; ASD shows air only for 2000-20000 A, so VUV lines have no air value).
* For each medium the "adopted" value is whichever of ASD observed / Ritz has the smaller stated
  uncertainty; with no stated uncertainty, the one with more decimals (tie -> Ritz).
* air_check_mA = [lambda_vac / n_PR72(lambda_vac)] - ASD_air, in mA, using the Peck & Reeder
  (1972) 5-parameter formula that ASD itself uses (consistency check only; differences of a few
  tenths of the last ASD digit are rounding).
"""
import csv, json, re, sys

sel_json, scripts_dir, out_csv, refs_csv = sys.argv[1:5]
sys.path.insert(0, scripts_dir)
from airglow_meta import M, REFS  # noqa: E402

# ---- corrections on top of airglow_meta (verification pass 2026-11-23) ----
# PALACE Table 2 reference height is 250 km for the 630/636 nm lines (300 km only for the other O lines);
# palace_lines hlay = 250 for 6302.046 / 6365.536 / 6393.500 A (vacuum).
for _k in ('OI_6300', 'OI_6364', 'OI_6392'):
    M[_k] = dict(M[_k], altitude_km='peak ~250 (PALACE ref. height 250; 200-300 typical); day >125 (Noxon68)')
# Per-row notes for the ASD level-less multiplet entries of H-alpha / H-beta (key: group, ASD Ritz air string).
# The Fraunhofer table (Sect. 6.1) uses ASD's separate blend-centroid rows c69 / c64 instead (6562.795, 4861.3250).
ROW_NOTE = {
    ('HI_Ha', '6562.819'): ('This row is the ASD n=2-3 multiplet entry (levels given only as n=2, n=3): Ritz 6562.819 +/- 0.007 A, '
                            'observed 6562.79 +/- 0.03 A (ASD line ref. L7400c29). The Fraunhofer C entry (Sect. 6.1) uses the '
                            'separate level-less ASD blend row c69, 6562.795 +/- 0.009 A (gA-weighted centroid of all '
                            'components), 24 mA shorter. Neither is the geocoronal emission centre: for geocoronal work use the '
                            '2s 2S 1/2 - 3p 2P* 3/2 and 2s 2S 1/2 - 3p 2P* 1/2 components at 6562.72483 and 6562.77153 A (air).'),
    ('HI_Hb', '4861.333'): ('This row is the ASD n=2-4 multiplet entry (levels given only as n=2, n=4): Ritz 4861.333 +/- 0.003 A, '
                            'observed 4861.35 +/- 0.05 A (ASD line ref. L7439c30). The Fraunhofer F entry (Sect. 6.1) uses the '
                            'separate level-less ASD blend row c64, 4861.3250 +/- 0.0014 A, 8 mA shorter. For geocoronal '
                            'work use the 2s 2S 1/2 - 4p 2P* components at 4861.286949 and 4861.297761 A (air); PALACE applies the '
                            'Gardner17 Halpha fine-structure fractions (~95 % in the 2s-np pair) to the other Balmer lines.'),
}

def num(s):
    s = re.sub(r'[^0-9.eE+-]', '', s or '')
    try:
        return float(s)
    except ValueError:
        return None

def ndec(s):
    s = re.sub(r'[^0-9.]', '', s or '')
    return len(s.split('.')[1]) if '.' in s else 0

def pick(obs, uo, ritz, ur):
    c = []
    if obs:
        c.append(('observed', obs, uo))
    if ritz:
        c.append(('Ritz', ritz, ur))
    if not c:
        return None
    withu = [x for x in c if num(x[2]) is not None]
    if withu:
        return min(withu, key=lambda x: (num(x[2]), x[0] != 'Ritz'))
    return max(c, key=lambda x: (ndec(x[1]), x[0] == 'Ritz'))

def n_air(lv):
    s2 = (1e4 / lv) ** 2
    return 1 + 1e-8 * (8060.51 + 2480990 / (132.274 - s2) + 17455.7 / (39.32957 - s2))

def fmt_level(conf, term, J):
    if conf.isdigit() and not term:
        return f'n={conf} (all fine-structure levels; ASD multiplet line)'
    return f'{conf} {term} {J}'.strip()

recs = json.load(open(sel_json))
rows = []
for r in recs:
    if r['status'] != 'ok':
        continue
    t = r['target']
    meta = M.get(t['group'], {})
    med = r['medium_default']
    if med == 'air':
        air = pick(r['obs_default'], r['unc_obs_default'], r['ritz_default'], r['unc_ritz_default'])
        # use the same kind (observed / Ritz) for the vacuum value when ASD provides it
        same = {'observed': (r['obs_vac'], r['unc_obs_vac']), 'Ritz': (r['ritz_vac'], r['unc_ritz_vac'])}[air[0]] if air else ('', '')
        vac = (air[0], same[0], same[1]) if same[0] else pick(r['obs_vac'], r['unc_obs_vac'], r['ritz_vac'], r['unc_ritz_vac'])
    else:
        air = None
        vac = pick(r['obs_default'], r['unc_obs_default'], r['ritz_default'], r['unc_ritz_default'])
    lv = num(vac[1]) if vac else None
    la = num(air[1]) if air else None
    main = air if air else vac
    unc = main[2] if main else ''
    # precision statement
    if num(unc) is not None:
        prec = 'yes (stated 1-sigma <= 1 mA)' if num(unc) <= 0.001 else f'no (stated 1-sigma {unc} A)'
    else:
        d = ndec(main[1]) if main else 0
        prec = (f'implied only ({d} decimals, no ASD uncertainty; ASD rule: 2.5-25 units of last digit)' if d >= 3
                else f'no ({d} decimals in ASD; 1 mA not available from ASD)')
    chk = ''
    if lv and la:
        chk = f'{(lv / n_air(lv) - la) * 1000:+.2f}'
    typ = '+'.join(sorted({p['type'] for p in r['parts']}))
    tR = meta.get('typical_R', '')
    if isinstance(tR, dict):
        w0 = la if la else lv
        best = min(tR, key=lambda k: abs(k - w0))
        tR = tR[best] if abs(best - w0) < 0.05 else ''
    Aparts = '; '.join(f"{p['type']}:{p['Aki']}({p['Acc']})" for p in r['parts'] if p['Aki'])
    rows.append({
        'species': t['spec'],
        'wl_vac_A': vac[1] if vac else '',
        'wl_air_A': air[1] if air else '',
        'uncert_A': unc,
        'wl_source': ('NIST ASD v5.12; ' + (f'air: {air[0]} (ASD std air, Peck&Reeder 1972 index); ' if air else 'air: n/a (VUV, ASD gives vacuum only); ')
                      + (f'vac: {vac[0]}' if vac else '')),
        'lower': fmt_level(r['conf_i'], r['term_i'], r['J_i']),
        'upper': fmt_level(r['conf_k'], r['term_k'], r['J_k']),
        'Aki_s-1': r['Aki_total'],
        'category': meta.get('category', ''),
        'mechanism': meta.get('mechanism', ''),
        'altitude_km': meta.get('altitude_km', ''),
        'typical_R': tR,
        'note': ' '.join(x for x in (meta.get('note', ''), ROW_NOTE.get((t['group'], air[1] if air else ''), '')) if x),
        'refs': meta.get('refs', 'ASD'),
        # ---- extra provenance columns ----
        'group': t['group'],
        'transition_type': typ,
        'Aki_parts(acc)': Aparts,
        'asd_obs_vac_A': r['obs_vac'], 'asd_unc_obs_vac': r['unc_obs_vac'],
        'asd_ritz_vac_A': r['ritz_vac'], 'asd_unc_ritz_vac': r['unc_ritz_vac'],
        'asd_obs_air_A': r['obs_default'] if med == 'air' else '', 'asd_unc_obs_air': r['unc_obs_default'] if med == 'air' else '',
        'asd_ritz_air_A': r['ritz_default'] if med == 'air' else '', 'asd_unc_ritz_air': r['unc_ritz_default'] if med == 'air' else '',
        'asd_wavenumber_cm-1': r['wn_cm'],
        'Ei_cm-1': r['Ei'], 'Ek_cm-1': r['Ek'],
        'meets_1mA': prec,
        'air_check_PR72_mA': chk,
        'asd_line_ref': ';'.join(sorted({p['line_ref'] for p in r['parts'] if p['line_ref']})),
        'asd_tp_ref': ';'.join(sorted({p['tp_ref'] for p in r['parts'] if p['tp_ref']})),
        '_sort': lv or 0,
    })

rows.sort(key=lambda x: x['_sort'])
cols = [k for k in rows[0] if k != '_sort']
with open(out_csv, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore')
    w.writeheader()
    for x in rows:
        w.writerow(x)
with open(refs_csv, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['key', 'citation'])
    for k, v in REFS.items():
        w.writerow([k, v])
print(len(rows), 'rows ->', out_csv)
bad = [x for x in rows if x['air_check_PR72_mA'] and abs(float(x['air_check_PR72_mA'])) > 0.6 * 10 ** (3 - ndec(x['wl_air_A']))]
print('air-check outliers (> 0.6 of last ASD digit):', len(bad))
for x in bad:
    print('  ', x['species'], x['wl_vac_A'], x['wl_air_A'], x['air_check_PR72_mA'])
from collections import Counter
print(Counter(x['meets_1mA'].split(' (')[0] for x in rows))
used = {(x['group'], x['wl_air_A']) for x in rows}
missing = [k for k in ROW_NOTE if k not in used]
if missing:
    raise SystemExit(f'ROW_NOTE keys not found in the ASD selection: {missing}')
