#!/usr/bin/env python
"""Convert the curated catalogues (data/curated/*.csv) into the normalised line-set JSON blocks
used by the page (data/build/<key>.json).

Normalised row keys
  a   air wavelength (A, number) or null (lambda < 2000 A)      v   vacuum wavelength (A, number)
  n   native medium of the source value: 'a' | 'v'              s   native value, verbatim string
  sa/sv  verbatim air / vacuum strings when the SOURCE gives that medium itself
  d   decimals supported by the native value                    u   uncertainty (A) or null
  us  1 if u is stated by the source (or a source-code bound), 0 if implied by digits
  id  identification text                                       st  plot strength 0..1
  lab short plot label;  + set-specific extra columns declared in meta.extra
Derived media use Peck & Reeder (1972) eq. 3 -- the formula NIST ASD uses -- for every set.
Implied uncertainty (no stated value): 2.5 units of the last published digit (NIST ASD guideline: implied
uncertainties are generally 2.5-25 units of the least significant figure); photographic solar atlases
(MMH 1966, Babcock & Moore 1947) at least 4 mA (measured MMH - APGL98 scatter).
"""
import json, math, os, re
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
CUR = os.path.join(ROOT, 'data', 'curated')
OUT = os.path.join(ROOT, 'data', 'build')


IMPLIED = 2.5      # units of the last digit for unstated uncertainties


def n_pr72(wv):
    s2 = (1e4 / wv) ** 2
    return 1 + 1e-8 * (8060.51 + 2480990 / (132.274 - s2) + 17455.7 / (39.32957 - s2))


def vac2air(wv):
    return wv / n_pr72(wv)


def air2vac(wa):
    v = wa
    for _ in range(30):
        nv = wa * n_pr72(v)
        if abs(nv - v) < 1e-12:
            break
        v = nv
    return v


def dec(s):
    s = re.sub(r'[^\d.]', '', str(s))
    return len(s.split('.')[1]) if '.' in s else 0


def num(x):
    try:
        v = float(str(x).strip().strip('()[]+'))
        return v if math.isfinite(v) else None
    except Exception:
        return None


def rd(x, k=5):
    return None if x is None else round(x, k)


def clip01(x):
    return max(0.04, min(1.0, x))


def read(name):
    return pd.read_csv(os.path.join(CUR, name), comment='#', dtype=str, keep_default_na=False)


def base_row(native, s, other_verbatim=None, u=None, us=0):
    """native 'a'|'v', s verbatim string; returns dict with a, v, n, s, d, u, us (+sa/sv)."""
    w = float(re.sub(r'[^\d.]', '', s))
    r = {'n': native, 's': s.strip(), 'd': dec(s)}
    if native == 'v':
        r['v'] = rd(w, 6)
        r['a'] = rd(vac2air(w), 6) if w >= 2000 else None
        r['sv'] = r['s']
        if other_verbatim:
            r['sa'] = other_verbatim.strip(); r['a'] = float(re.sub(r'[^\d.]', '', other_verbatim))
    else:
        r['a'] = rd(w, 6)
        r['v'] = rd(air2vac(w), 6)
        r['sa'] = r['s']
        if other_verbatim:
            r['sv'] = other_verbatim.strip(); r['v'] = float(re.sub(r'[^\d.]', '', other_verbatim))
    if u is None:
        r['u'] = IMPLIED * 10.0 ** (-r['d']); r['us'] = 0
    else:
        r['u'] = u; r['us'] = us
    return r


def compact(rows):
    """drop fields the page can re-derive: d (from the native string), s when equal to sa/sv,
    numeric a/v when the verbatim string of that medium is present; round derived numbers."""
    for r in rows:
        r.pop('d', None)
        if r.get('n') == 'a' and r.get('sa') == r.get('s'):
            r.pop('s', None)
        if r.get('n') == 'v' and r.get('sv') == r.get('s'):
            r.pop('s', None)
        if r.get('sa'):
            r.pop('a', None)
        elif r.get('a') is not None:
            r['a'] = round(r['a'], 5)
        if r.get('sv'):
            r.pop('v', None)
        elif r.get('v') is not None:
            r['v'] = round(r['v'], 5)
        if isinstance(r.get('u'), float):
            r['u'] = float(f"{r['u']:.3g}")
        for k, v in list(r.items()):
            if isinstance(v, float) and k not in ('a', 'v', 'u'):
                r[k] = float(f'{v:.4g}') if abs(v) < 1e-2 or abs(v) >= 1e5 else round(v, 4)
            if v in ('', None):
                r.pop(k)
    return rows


def pack(key, meta, rows):
    rows = compact(rows)
    cols = []
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    table = [[r.get(c) for c in cols] for r in rows]
    dicts = {}
    for j, c in enumerate(cols):          # dictionary-encode repetitive string columns
        vals = [t[j] for t in table if t[j] is not None]
        if vals and all(isinstance(v, str) for v in vals) and len(rows) > 300:
            uniq = sorted(set(vals))
            if len(uniq) < 0.35 * len(vals) and c not in ('s', 'sa', 'sv'):
                ix = {v: i for i, v in enumerate(uniq)}
                for t in table:
                    if t[j] is not None:
                        t[j] = ix[t[j]]
                dicts[c] = uniq
    data = {'meta': meta, 'cols': cols, 'dicts': dicts, 'rows': table}
    fn = os.path.join(OUT, key + '.json')
    json.dump(data, open(fn, 'w'), separators=(',', ':'), allow_nan=False)
    print(f'{key:12s} {len(rows):6d} rows  {os.path.getsize(fn) / 1e3:8.0f} kB')


def parse_R(txt):
    """first intensity in rayleigh found in a free-text field -> R (number) or None"""
    m = re.search(r'([\d.]+)(?:\s*[-–]\s*[\d.]+)?\s*(MR|kR|R)\b', txt or '')
    if not m:
        return None
    v = float(m.group(1)) * {'R': 1, 'kR': 1e3, 'MR': 1e6}[m.group(2)]
    return v


# ------------------------------------------------------------------ airglow atomic
def ag_atomic():
    d = read('airglow_atomic.csv')
    rows = []
    for _, x in d.iterrows():
        wv = num(x['wl_vac_A'])
        has_air = bool(x['wl_air_A'].strip()) and wv >= 2000
        if has_air:
            r = base_row('a', x['wl_air_A'], other_verbatim=x['wl_vac_A'])
        else:
            r = base_row('v', x['wl_vac_A'])
        u = num(x['uncert_A'])
        if u is not None:
            r['u'], r['us'] = u, 1
        tt = x['transition_type']
        sp = x['species'] if 'E1' in tt else f"[{x['species']}]"
        R = parse_R(x['typical_R'])
        r.update({'id': f"{sp} {x['lower']} → {x['upper']}", 'lab': sp, 'sp': x['species'], 'tt': tt,
                  'A': num(x['Aki_s-1']), 'cat': x['category'].replace('|', ', '), 'mech': x['mechanism'], 'alt': x['altitude_km'],
                  'I': x['typical_R'], 'note': x['note'], 'refs': x['refs'], 'ok': x['meets_1mA'].split(' (')[0],
                  'st': clip01((math.log10(R) + 1) / 5) if R else 0.35})
        rows.append(r)
    rows.sort(key=lambda r: r['v'])
    meta = {'name': 'Atomic airglow, aurora, meteor and geocoronal lines (NIST ASD)', 'short': 'Airglow atoms',
            'color': '#2ca02c', 'mode': 'emission', 'host': '#agAtomic', 'cite': 'NIST ASD v5.12; categories and intensities from the cited aeronomy literature',
            'idLabel': 'Line (lower → upper)', 'pageSize': 250,
            'extra': [{'k': 'tt', 'label': 'Type', 'title': 'E1 allowed; M1/E2 forbidden ([ ])'},
                      {'k': 'A', 'label': 'A<sub>ki</sub> (s⁻¹)', 'num': True, 'fmt': 'sig'},
                      {'k': 'cat', 'label': 'Category'}, {'k': 'mech', 'label': 'Excitation'},
                      {'k': 'alt', 'label': 'Altitude (km)'}, {'k': 'I', 'label': 'Typical intensity'},
                      {'k': 'ok', 'label': '1 mÅ?', 'title': 'yes = ASD-stated σ ≤ 1 mÅ; implied = ≥3 decimals but no stated σ; no = below target'},
                      {'k': 'note', 'label': 'Note'}, {'k': 'refs', 'label': 'Refs'}],
            'filters': [{'k': 'ok', 'label': '1 mÅ status', 'values': sorted(set(r['ok'] for r in rows))}]}
    pack('ag_atomic', meta, rows)


# ------------------------------------------------------------------ airglow NIR (HITRAN2024 + PALACE)
def read_opt(*names):
    """concatenate the curated files that exist (later files: optional extensions)"""
    frames = [read(n) for n in names if os.path.exists(os.path.join(CUR, n))]
    return pd.concat(frames, ignore_index=True).fillna('') if len(frames) > 1 else frames[0]


def ag_nir():
    d = read_opt('airglow_nir.csv', 'airglow_nir_oh_fundamental.csv')
    rows, ids = [], set()
    for _, x in d.iterrows():
        u = num(x['uncert_A'])
        r = base_row('v', x['wl_vac_A'], u=u, us=1 if u is not None else 0)
        if u is None:
            r['u'] = None; r['us'] = 0
        I = num(x['intensity_R'])
        if x['palace_ID']:
            ids.add(x['palace_ID'])
        lab = x['branch_label'].replace('PALACE:', '')
        r.update({'ib': f"{x['isotopologue']} {x['band']}", 'lab': lab, 'mol': x['molecule'], 'I': I, 'T': num(x['palace_trans_zenith']), 'A': num(x['A_einstein_s-1']),
                  'dO': num(x['d_oliva2015_mA']),
                  'dR': num(x['d_rousselot2000_doublet_mean_mA']) if lab.startswith('Q') else num(x['d_rousselot2000_mA']),
                  'dP': num(x['d_palace2025_mA']),
                  'src': x['source_key'], 'note': x['notes'],
                  'st': clip01((math.log10(I) + 2.5) / 4) if I and I > 0 else 0.04})
        rows.append(r)
    rows.sort(key=lambda r: r['v'])
    meta = {'name': 'Near-IR airglow 0.74–4.2 µm: OH Meinel & O₂ (HITRAN2024 positions; PALACE intensities to 2.5 µm)', 'short': 'NIR OH/O₂',
            'color': '#d62728', 'mode': 'emission', 'host': '#agNir', 'cite': 'HITRAN2024; Noll et al. 2025 PALACE v1.0', 'pageSize': 200,
            'idLabel': 'Isotopologue, band, line', 'sortKey': 'v', 'idCols': ['ib', 'lab'], 'chem': True,
            'extra': [{'k': 'I', 'label': 'I (R)', 'num': True, 'fmt': 'sig', 'title': 'PALACE v1.0 reference intensity, Paranal zenith climatology, above lower-atmosphere absorption; blank for the 2336 OH Δv = 1 / Δv = 2 lines beyond 2.5 µm (outside PALACE)'},
                      {'k': 'T', 'label': 'T<sub>zen</sub>', 'num': True, 'fmt': 'fix', 'n': 3, 'title': 'PALACE zenith transmission'},
                      {'k': 'A', 'label': 'A (s⁻¹)', 'num': True, 'fmt': 'sig'},
                      {'k': 'dO', 'label': 'Oliva15 − HITRAN (mÅ)', 'num': True, 'title': 'Oliva et al. 2015 minus primary, mÅ'},
                      {'k': 'dR', 'label': 'Rousselot00 − HITRAN (mÅ)', 'num': True, 'title': 'Rousselot et al. 2000 minus primary, mÅ; Q-branch rows compare Λ-doublet means because Rousselot’s Q-branch Λ-splittings are wrong'},
                      {'k': 'mol', 'label': 'Mol.'}, {'k': 'src', 'label': 'Sources'}, {'k': 'note', 'label': 'Note'}],
            'filters': [{'k': 'mol', 'label': 'Molecule', 'values': sorted(set(r['mol'] for r in rows))},
                        {'k': 'I', 'label': 'I ≥', 'type': 'min', 'values': [0.01, 0.1, 1, 10, 100]}]}
    pack('ag_nir', meta, rows)
    return ids


# ------------------------------------------------------------------ airglow optical model (PALACE, minus NIR overlap)
def ag_palace(nir_ids):
    d = read_opt('airglow_model_palace_3000_10500.csv', 'airglow_model_palace_atomic_ir.csv')
    rows, dropped = [], 0
    for _, x in d.iterrows():
        if x['palace_id'] in nir_ids:
            dropped += 1; continue
        r = base_row('v', x['wl_vac_A'], u=None)
        r['u'] = None; r['us'] = 0
        I = num(x['I_R']); IT = num(x['I_times_trans_R'])
        r.update({'id': x['ident'], 'lab': x['ident'].split('[')[0].strip(), 'sp': x['species'], 'sys': x['system'], 'band': x['band'],
                  'I': I, 'T': num(x['trans']), 'IT': IT, 'h': num(x['hlay_km']),
                  'st': clip01((math.log10(IT) + 2.5) / 4.5) if IT and IT > 0 else 0.04})
        rows.append(r)
    rows.sort(key=lambda r: r['v'])
    print('  PALACE rows dropped as duplicates of the NIR list:', dropped)
    meta = {'name': 'Airglow model lines 3000–10 500 Å plus O I to 2.13 µm (PALACE v1.0; lab positions)', 'short': 'Airglow model',
            'color': '#ff7f0e', 'mode': 'emission', 'host': '#agPalace', 'cite': 'Noll et al. 2025', 'pageSize': 200,
            'idLabel': 'Identification', 'sortKey': 'v',
            'extra': [{'k': 'sp', 'label': 'Species'}, {'k': 'sys', 'label': 'System'},
                      {'k': 'I', 'label': 'I (R)', 'num': True, 'fmt': 'sig', 'title': 'PALACE reference intensity (Paranal zenith, 100 sfu climatology)'},
                      {'k': 'T', 'label': 'T<sub>zen</sub>', 'num': True, 'fmt': 'fix', 'n': 3},
                      {'k': 'IT', 'label': 'I·T (R)', 'num': True, 'fmt': 'sig', 'title': 'intensity reaching the ground at zenith'},
                      {'k': 'h', 'label': 'Layer (km)', 'num': True}],
            'filters': [{'k': 'sp', 'label': 'Species', 'values': sorted(set(r['sp'] for r in rows))},
                        {'k': 'IT', 'label': 'I·T ≥', 'type': 'min', 'values': [0.01, 0.1, 1, 10, 100]}]}
    pack('ag_palace', meta, rows)


# ------------------------------------------------------------------ UVES observed atlas
def ag_optical():
    d = read('airglow_optical.csv')
    d = d[d['is_primary'] == '1']
    rows = []
    for _, x in d.iterrows():
        r = base_row('a', x['wl_air_A'], u=num(x['uncert_A']), us=1)
        F = num(x['flux'])
        idt = x['ident'] or 'unidentified'
        r.update({'id': idt if len(idt) < 140 else idt[:137] + '…', 'lab': (x['species_main'] or '?'), 'F': F, 'fw': num(x['fwhm_A']),
                  'mw': x['ident_model_wl_air_A'], 'dm': (num(x['ident_dlam_A']) * 1000) if num(x['ident_dlam_A']) is not None else None,
                  'fl': x['ident_flag'], 'set': x['setting'], 'hw': x['hires_wl_air_A'], 'sp': x['species_main'] or '—', 'nt': x['notes'],
                  'st': clip01((math.log10(F) + 0.5) / 4) if F and F > 0 else 0.04})
        rows.append(r)
    rows.sort(key=lambda r: r['v'])
    meta = {'name': 'UVES/VLT observed night-sky lines (Hanuschik 2003)', 'short': 'UVES sky (obs.)',
            'color': '#9467bd', 'mode': 'emission', 'host': '#agOptical', 'cite': 'Hanuschik 2003; identifications re-derived after Cosby et al. 2006', 'pageSize': 200,
            'idLabel': 'Identification (derived)', 'sortKey': 'v', 'defaultOn': False,
            'extra': [{'k': 'mw', 'label': 'model λ<sub>air</sub> (Å)', 'num': True, 'title': 'air wavelength of the identified PALACE/NIST line (lab position)'},
                      {'k': 'dm', 'label': 'obs − model (mÅ)', 'num': True, 'fmt': 'fix', 'n': 1},
                      {'k': 'fl', 'label': 'Flag', 'title': 'Cosby et al. 2006 labels: S1 single, L2 two lines, B2 blend, Mn multiple, U unidentified; w = tentative'},
                      {'k': 'F', 'label': 'Peak flux', 'num': True, 'fmt': 'sig', 'title': 'Gaussian-fit peak flux density incl. continuum, 1e-16 erg s-1 cm-2 Å-1 arcsec-2'},
                      {'k': 'fw', 'label': 'FWHM (Å)', 'num': True, 'fmt': 'fix', 'n': 4},
                      {'k': 'set', 'label': 'UVES setting'}, {'k': 'hw', 'label': 'Keck/HIRES λ (Å)'}, {'k': 'nt', 'label': 'Note'}],
            'filters': [{'k': 'sp', 'label': 'Species', 'values': sorted(set(r['sp'] for r in rows))}]}
    pack('ag_optical', meta, rows)


# ------------------------------------------------------------------ Fraunhofer classic letters
FH_CAVEATS = {   # review findings (content-writer check of fraunhofer_classic.csv, 2026-11-23)
    'B': 'Caveat: the MMH blend pair 6869.887/6870.007 does not map cleanly onto HITRAN 6869.897/6869.974 (+32.6 mÅ for this pairing, larger than the +2–6 mÅ O₂ pressure shift).',
    'U': 'Caveat: MMH 2947.70 (2 decimals) is −0.18 Å from Fe I 2947.8758, evidently a different MMH entry; Rowland 2947.99 − 0.106 ≈ 2947.88 agrees with Fe I.',
    'P': 'Caveat: MMH 3361.193r lies 19 mÅ (−1.7 km/s) shortward of Ti II 3361.2121 — the solar feature is a blend (MMH: with 3361.287r).',
    'O': 'Caveat: MMH 3441.019r is +30 mÅ above Fe I 3440.9885; MMH values below ≈3700 Å (r flags) are of historical quality (MMH − lab −1.7 to +2.7 km/s).',
    'E1': 'Caveat: MMH 5270.383 is a photographic blend centre (+17 mÅ vs APGL98 5270.3658; Ca I 5270.27 adjacent); not usable at the mÅ level.',
    'D2': 'Caveat: solar Na D2 position is fit-dependent: MMH +22 mÅ, APGL98 flux +466 m/s, IAG core fit +816 m/s relative to lab; unresolved.',
    'y': 'Tentative identification: H₂O 8987.606 Å is 0.044 Å from the commonly quoted 898.765 nm; MMH does not cover this range.',
    'Z': 'Designation source is effectively Wikipedia; MMH supplies only the H₂O identification.',
    'h': 'Rowland 4102.000 converts to ≈4101.85 Å (ISA), 0.10 Å redward of MMH 4101.748 — a historical Rowland-scale problem.',
}


def fh_classic():
    d = read('fraunhofer_classic.csv')
    rows = []
    for _, x in d.iterrows():
        wa, wv = x['wl_air_A'].strip(), x['wl_vac_A'].strip()
        r = base_row('a', wa, other_verbatim=wv or None, u=num(x['uncert_A']), us=1 if num(x['uncert_A']) is not None else 0)
        if num(x['uncert_A']) is None:
            r['u'] = IMPLIED * 10.0 ** (-dec(wa)); r['us'] = 0
        if x['lab_value_type'].strip().lower().startswith('derived'):
            r['us'] = 0
        r.update({'id': f"{x['letter']} — {x['species']}", 'lab': x['letter'], 'L': x['letter'], 'org': x['origin'],
                  'tr': x['transition'], 'mmh': x['solar_mmh1966_air_A'], 'apg': x['solar_apgl98_flux_air_A'],
                  'iag': x['solar_iag_flux_air_A_derived'], 'diag': num(x['iag_minus_lab_mps_derived']),
                  'row': x['rowland_value_rowland_scale_A'], 'wik': x['wikipedia_value_nm'],
                  'note': ' '.join(t for t in (FH_CAVEATS.get(x['letter'].strip(), ''), x['note'],
                                               '(HITRAN2020 wavenumbers are identical in HITRAN2024.)' if 'HITRAN2020' in x['refs'] else '') if t),
                  'refs': x['refs'], 'lt': x['lab_value_type'], 'st': 1.0})
        rows.append(r)
    rows.sort(key=lambda r: r['v'])
    meta = {'name': 'Classic Fraunhofer designations', 'short': 'Fraunhofer letters', 'color': '#e377c2', 'mode': 'absorption',
            'host': '#fhClassic', 'cite': 'Fraunhofer 1817; Rowland; Moore et al. 1966; lab values NIST ASD / HITRAN2020', 'pageSize': 100,
            'idLabel': 'Letter — species', 'sortKey': 'v', 'sortDesc': True,
            'extra': [{'k': 'org', 'label': 'Origin'}, {'k': 'lt', 'label': 'Lab value'},
                      {'k': 'mmh', 'label': 'Sun: MMH 1966 (Å)', 'title': 'Moore, Minnaert & Houtgast 1966 disc-centre air wavelength (includes solar shifts)'},
                      {'k': 'apg', 'label': 'Sun: APGL98 flux (Å)'}, {'k': 'iag', 'label': 'Sun: IAG fit (Å)', 'title': 'parabola fit to the IAG 2016 flux atlas line core (derived)'},
                      {'k': 'diag', 'label': 'IAG − lab (m/s)', 'num': True},
                      {'k': 'row', 'label': 'Rowland (Rowland scale)'}, {'k': 'wik', 'label': 'Wikipedia (nm)'},
                      {'k': 'tr', 'label': 'Transition'}, {'k': 'note', 'label': 'Note'}, {'k': 'refs', 'label': 'Refs'}]}
    pack('fh_classic', meta, rows)


# ------------------------------------------------------------------ solar photospheric lines
def solar():
    d = read('solar_lines.csv')
    rows = []
    for _, x in d.iterrows():
        e = num(x['e_wl_A'])
        r = base_row('a', x['wl_air_A'], u=e, us=1 if e is not None else 0)
        if e is None:      # MMH 1966, B&M 1947 and APGL98b-only rows: no catalogue sigma
            r['u'] = max(IMPLIED * 10.0 ** (-int(x['n_decimals'])), 0.004); r['us'] = 0
        ew = num(x['ew_mA']); dp = num(x['depth_flux'])
        wt = x['wl_type']
        ident = x['ident']
        ocr = x['mmh1966_ocr_check']
        ocrlab = ('OCR ✓' if ocr == 'OCR-match' else 'OCR ?' if ocr.startswith('not in OCR; nearest') else 'no OCR' if ocr.startswith('not in OCR')
                  else 'page ✓' if ocr.startswith(('page image', 'not in MMH')) else '—') if 'MMH' in x['source_key'] else ''
        r.update({'ocr': ocrlab, 'fl': x['flag_short'], 'nt': x['notes']})
        r.update({'id': ident, 'lab': ident,
                  'typ': ('flux' if 'flux' in wt else 'disc-centre') + (' (MMH 1966)' if 'MMH' in wt else ' (Babcock & Moore 1947)' if 'Babcock' in wt else ' (APGL98b)' if 'APGL98b' in wt else ' (APGL98)'),
                  'ew': ew, 'dp': dp, 'lw': x['lab_wl_air_A'], 'dv': num(x['solar_flux_minus_lab_mps']), 'cb': num(x['conv_blueshift_iag_mps']),
                  'ep': num(x['ep_eV']), 'gf': num(x['loggf_solar']), 'src': x['source_key'],
                  'el': x['element'] or '?',
                  'st': clip01(dp) if dp else (clip01(math.log10(ew) / 2.7) if ew and ew > 0 else 0.3)})
        rows.append(r)
    rows.sort(key=lambda r: r['v'])
    els = sorted(set(r['el'] for r in rows))
    meta = {'name': 'Solar photospheric absorption lines', 'short': 'Solar lines', 'color': '#8c564b', 'mode': 'absorption',
            'host': '#fhSolar', 'cite': 'Allende Prieto & García López 1998; Moore et al. 1966; Reiners et al. 2016; Nave et al. 1994', 'pageSize': 200,
            'idLabel': 'Species', 'sortKey': 'v',
            'extra': [{'k': 'typ', 'label': 'Solar λ type'}, {'k': 'ew', 'label': 'EW (mÅ)', 'num': True},
                      {'k': 'dp', 'label': 'Depth', 'num': True, 'fmt': 'fix', 'n': 2},
                      {'k': 'lw', 'label': 'Lab λ<sub>air</sub> (Å)', 'num': True},
                      {'k': 'dv', 'label': 'Sun − lab (m/s)', 'num': True, 'title': 'solar flux minus laboratory wavelength, as velocity'},
                      {'k': 'cb', 'label': 'Conv. blueshift (m/s)', 'num': True, 'title': 'IAG 2016 convective blueshift (gravitational redshift removed)'},
                      {'k': 'ep', 'label': 'χ (eV)', 'num': True}, {'k': 'gf', 'label': 'log gf', 'num': True}, {'k': 'src', 'label': 'Sources'},
                      {'k': 'ocr', 'label': 'MMH check', 'title': 'MMH 1966 value (HelioSpectrotron digitisation) checked against the scanned book: OCR ✓ matches the OCR text; page ✓ read on the page image; OCR ? only a nearby token found (possible digitisation error)'},
                      {'k': 'fl', 'label': 'Flag'}, {'k': 'nt', 'label': 'Note'}],
            'filters': [{'k': 'el', 'label': 'Element', 'values': els}]}
    pack('solar', meta, rows)


# ------------------------------------------------------------------ telluric (HITRAN2024)
def telluric():
    d = read('telluric_lines.csv')
    rows = []
    for _, x in d.iterrows():
        u = num(x['lambda_unc_max_mA'])
        r = base_row('v', x['lambda_vac_A'], u=(u / 1000 if u is not None else None), us=1)   # air derived (P&R 1972)
        if u is None:
            r['u'] = None; r['us'] = 0
        q = f"{x['upper_local_quanta'].strip()} ← {x['lower_local_quanta'].strip()}".strip(' ←')
        dp = num(x['depth_R1e5_sl'])
        r.update({'ib': f"{x['isotopologue']} {x['band']}", 'q': q, 'mol': x['molecule'], 'reg': x['region'],
                  'nu': x['wavenumber_vac_cm-1'], 'S': num(x['S296_cm-1/(molec.cm-2)']), 'El': num(x['E_lower_cm-1']),
                  'dps': num(x['dlam_pshift_eff_mA_sl']),
                  'D': dp, 'D25': num(x['depth_R1e5_2p5km']), 'EW': num(x['EW_mA_sl'])})
        rows.append(r)
    rows.sort(key=lambda r: r['v'])
    meta = {'name': 'Telluric absorption lines (HITRAN2024, depth ≥ 1 % at R = 10⁵, sea level)', 'short': 'Telluric',
            'color': '#17becf', 'mode': 'absorption', 'host': '#fhTell', 'cite': 'HITRAN2024; depths from an AFGL US-standard model', 'pageSize': 200,
            'idLabel': 'Isotopologue, band, upper ← lower', 'sortKey': 'v', 'idCols': ['ib', 'q'], 'chem': True, 'labKey': 'mol', 'stKey': 'D', 'calc': 'telluric',
            'extra': [{'k': 'D', 'label': 'Depth (sea level)', 'num': True, 'fmt': 'fix', 'n': 3, 'title': 'central depth after R = 100 000 LSF, airmass 1, sea level (PWV 14.2 mm)'},
                      {'k': 'D25', 'label': 'Depth (2.5 km)', 'num': True, 'fmt': 'fix', 'n': 3},
                      {'k': 'EW', 'label': 'EW (mÅ)', 'num': True, 'fmt': 'fix', 'n': 2},
                      {'k': 'dps', 'label': 'P-shift (mÅ)', 'num': True, 'fmt': 'fix', 'n': 2, 'title': 'pressure shift through the sea-level atmosphere (S·u-weighted); add to the zero-pressure λ'},
                      {'k': 'aobs', 'label': 'λ<sub>air</sub> + P-shift (Å)', 'num': True, 'fmt': 'fix', 'n': 4, 'title': 'zero-pressure air λ plus the sea-level pressure shift (computed in the page)'},
                      {'k': 'nu', 'label': 'ν̃ (cm⁻¹, HITRAN)', 'num': True}, {'k': 'S', 'label': 'S<sub>296</sub>', 'num': True, 'fmt': 'sig'},
                      {'k': 'El', 'label': 'E″ (cm⁻¹)', 'num': True}, {'k': 'reg', 'label': 'Region'}],
            'filters': [{'k': 'mol', 'label': 'Molecule', 'values': sorted(set(r['mol'] for r in rows))},
                        {'k': 'D', 'label': 'Depth ≥', 'type': 'min', 'values': [0.02, 0.05, 0.1, 0.3, 0.5, 0.9]}]}
    pack('telluric', meta, rows)


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    ag_atomic()
    ids = ag_nir()
    ag_palace(ids)
    ag_optical()
    fh_classic()
    solar()
    telluric()
