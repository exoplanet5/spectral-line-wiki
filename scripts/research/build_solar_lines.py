#!/usr/bin/env python
"""Build data/curated/solar_lines.csv -- an extended solar photospheric absorption-line list.

Sources (all downloaded, see data/raw/fraunhofer_solar/):
  APGL98   Allende Prieto & Garcia Lopez 1998, A&AS 131, 431 (VizieR J/A+AS/131/431): 4947 rows,
           solar FLUX (Kurucz et al. 1984 FTS atlas) and DISC-CENTRE (Brault & Neckel 1987 FTS atlas)
           central wavelengths, air, 4 decimals, with fit errors; ident/EP/log gf from Thevenin 1989/1990.
           Ten rows occur twice with every field identical (VizieR recno pairs); the copy is dropped.
  APGL98b  Allende Prieto & Garcia Lopez 1998, A&AS 129, 41 (J/A+AS/129/41): 1446 Fe I lines with the
           same solar wavelengths plus laboratory wavelengths (Nave et al. 1994) and disc-centre EW.
  IAG16    Reiners et al. 2016, A&A 587, A65 (J/A+A/587/A65/tablea1): 1251 Fe I lines, flux depth and
           absolute convective blueshift (gravitational redshift 636 m/s removed) from the IAG atlas.
  NIST     NIST ASD v5.12 air wavelengths (lab) for the other species (nearest-line match, flagged).
  MMH66    Moore, Minnaert & Houtgast 1966 (NBS Monograph 61) as digitised by the HelioSpectrotron 5000
           project (Pietrow 2026; ~4965 deepest lines per 0.5 A, 3373.65-7328.83 A), cross-checked here
           against the Internet-Archive OCR text layer of the NBS PDF and, for flagged rows, the page images.
  BM47     Babcock & Moore 1947 (6600-13495 A) HS5000 digitisation, used above the end of the HS5000 MMH
           digitisation (lambda > 7328.83 A) unless an APGL98/APGL98b line lies within +-0.02 A.

Laboratory wavelengths are attached only when the lower-level excitation potential of the lab line
agrees with the solar catalogue's EP to within 0.02 eV (EP_lab from APGL98b, or NIST Ei/8065.544).

Digitisation precision: MMH strings are replaced by the printed MMH string when the OCR text shows the same
number with trailing zeros restored (or removed); n_decimals then follows the printed string.  B&M 1947 values
with 4 decimals keep their string but n_decimals is capped at 3 (flagged); integer values are flagged.
Identifications: LaTeX Greek is converted (H$\\alpha$ -> 'H I (Hα)'); the HS5000 codes 'T I', 'O', 'Atm', 'Fe'
were read on the MMH page images; 'c' / 'Cp' (B&M range, no page images here) become 'unidentified (HS5000 code ..)'.
Columns line_kind (absorption / emission (coronal) / telluric / chromospheric He I) and flag_short (compact
page flag) summarise wl_flag, ident_flag and notes.

Usage: python build_solar_lines.py <root>
"""
import bisect, csv, glob, os, re, statistics, sys
from collections import Counter, defaultdict

ROOT = sys.argv[1]
RAW = os.path.join(ROOT, 'data', 'raw', 'fraunhofer_solar')
OUT = os.path.join(ROOT, 'data', 'curated', 'solar_lines.csv')
C = 299792458.0
GRS = 636.0           # m/s, value adopted by Reiners et al. 2016 (and APGL98b) for the solar gravitational redshift
EV_CM = 8065.544      # cm-1 per eV (CODATA 2018: 8065.543937)
EP_TOL = 0.02         # eV, max |EP_solar - EP_lab| for a lab wavelength to be attached
MMH_END = 7328.83     # A, last line of the HS5000 MMH digitisation (7328.828)
BM_APG_TOL = 0.02     # A, a B&M line within this distance of an APGL98 line is not added as a row


def viz(fn):
    rows, hdr = [], None
    for ln in open(fn):
        if ln.startswith('#') or not ln.strip():
            continue
        p = ln.rstrip('\n').split('\t')
        if hdr is None:
            hdr = [h.strip() for h in p]; continue
        if not p[0].strip().isdigit():
            continue
        rows.append({h: v.strip() for h, v in zip(hdr, p)})
    return rows


def ndec(s):
    s = re.sub(r'[^\d.]', '', s)
    return len(s.split('.')[1]) if '.' in s else 0


ROMAN = {'1': 'I', '2': 'II', '3': 'III'}


def ion_name(code):
    code = code.strip()
    if not code:
        return ''
    return f"{code[:-1].strip().capitalize()} {ROMAN.get(code[-1], code[-1])}"


# ---------------------------------------------------------------- identifications of the HS5000 digitisation
GREEK = {'alpha': 'α', 'beta': 'β', 'gamma': 'γ', 'delta': 'δ', 'epsilon': 'ε'}
# rows whose HS5000 identification was checked on the NBS Monograph 61 page images (pdf page in brackets)
MMH_IDENT_CHECK = {
    '4821.601': ('Fe I', 'MMH page image (pdf p.228) reads Fe I (EW 4; the HS5000 EW 8 belongs to 4821.479); HS5000 "T I"'),
    '5389.847': ('unidentified', 'MMH page image (pdf p.264) gives no identification (EW 12, w); HS5000 "T I"'),
    '5585.663': ('Ti I', 'MMH page image (pdf p.276) reads Ti I; HS5000 "T I"'),
    '5739.061': ('Ti I', 'MMH page image (pdf p.284) reads 5739.061r Ti I p (EP 2.30, mult. 249); HS5000 "T I"'),
    '6017': ('unidentified', 'MMH page image (pdf p.300) reads 6017.00 m, EW 1.5, no identification; HS5000 "T I"'),
    '6996.634': ('Ti I', 'MMH page image (pdf p.330) reads Ti I (EP 2.33, mult. 256); HS5000 "T I"'),
    '7311.64': ('unidentified', 'MMH page image (pdf p.343) reads "⊙?" (probably solar, unidentified); HS5000 "O"'),
    '5624.558': ('Fe I', 'MMH page image (pdf p.278) reads 5624.558S Fe I + V I blend; HS5000 "Fe"'),
    '6397.545': ('telluric (MMH: Atm)', 'MMH page image (pdf p.312) reads 6397.545r Atm (telluric, molecule not given)'),
    '6374.56': ('Fe X', 'coronal emission line (Fe X red line), not a photospheric absorption line; not in MMH 1966 '
                '(page image pdf p.311 lists only 6374.50 a, EW 2.5, unidentified); added by HS5000'),
    '6701.47': ('Ni XV', 'coronal emission line, not a photospheric absorption line; not in MMH 1966 '
                '(page image pdf p.320 lists 6701.377r, EW 1.5, "⊙" unidentified); added by HS5000'),
}
BM_IDENT_FIX = {   # B&M 1947 pages are not available here; only notation is normalised
    'c': ('unidentified (HS5000 code c)', 'HS5000 identification "c" is not a species; B&M 1947 page not available here, not checked'),
    'Cp': ('unidentified (HS5000 code Cp)', 'HS5000 identification "Cp" is not a species; B&M 1947 page not available here, not checked'),
    'Sa II': ('Sm II', 'HS5000 "Sa II": Sa is the obsolete symbol for samarium; identification not checked'),
    'N I 1': ('N I', 'HS5000 "N I 1" (trailing digit dropped)'),
    'S I 1': ('S I', 'HS5000 "S I 1" (trailing digit dropped)'),
    '': ('unidentified', 'no identification in the HS5000 digitisation'),
}


def clean_ident(raw, wl, src):
    """-> (ident, element, ident_flag).  raw = HS5000 identification string."""
    s = raw.strip()
    flag = ''
    if src == 'MMH' and wl in MMH_IDENT_CHECK:
        s, flag = MMH_IDENT_CHECK[wl]
    elif src == 'BM' and s in BM_IDENT_FIX:
        s, flag = BM_IDENT_FIX[s]
    m = re.match(r'^H\$\\(\w+)\$$', s)                       # H$\alpha$ ...
    if m:
        s = f'H I (H{GREEK[m.group(1)]})'
    m = re.match(r'^H (\d+)$', s)                             # H 8 ... H 16 (Balmer series members)
    if m:
        s = f'H I (H{m.group(1)})'
    m = re.match(r'^Paschen(?:\$\\(\w+)\$| (\d+))$', s)       # Paschen$\beta$, Paschen 11
    if m:
        s = f'H I (Pa {GREEK[m.group(1)] if m.group(1) else m.group(2)})'
    m = re.match(r'^(Ca II|Mg I|Fe I|Na I|He I) ([HK]|b\d|D\d)$', s)   # Ca II K, Mg I b1, Na I D2 ...
    if m:
        s = f'{m.group(1)} ({m.group(2)})'
    if s.startswith('unidentified'):
        el = '?'
    elif s.startswith('telluric'):
        el = 'telluric'
    else:
        el = s.split()[0]
    assert '$' not in s and s not in ('T I', 'c', 'Cp', 'O', 'Sa II', 'Atm', 'Fe') and el not in ('T', 'Sa', 'Atm'), (raw, s)
    return s, el, flag


def line_kind(ident):
    """Physical kind of the solar-spectrum feature (all rows are absorption lines unless stated)."""
    if re.match(r'^(Fe X|Ni XV)$', ident):
        return 'emission (coronal, not a photospheric absorption line)'
    if ident.startswith('telluric'):
        return 'telluric absorption'
    if ident.startswith('He I'):
        return 'chromospheric (He I)'
    return 'absorption'


def short_flags(d):
    """Compact flag text for the page table (details stay in wl_flag / ident_flag / notes)."""
    f = []
    if d.get('line_kind', 'absorption') != 'absorption':
        f.append(d['line_kind'])
    wf, idf, nt = d.get('wl_flag', ''), d.get('ident_flag', ''), d.get('notes', '')
    if 'integer value' in wf:
        f.append('integer λ in the digitisation (decimals unknown)')
    if 'more digits than B&M' in wf:
        f.append('digits beyond B&M 1947 print precision')
    if 'printed value from the MMH page image' in wf:
        f.append('λ decimals restored from the MMH page image')
    elif 'printed value from the MMH page image' in nt:
        f.append('MMH λ decimals restored from the page image')
    if 'differs from the MMH page image' in wf + nt:
        f.append('HS5000 value differs from the MMH page')
    if 'not in MMH 1966' in wf + idf:
        f.append('not in MMH 1966 (added by HS5000)')
    if idf.startswith('MMH page image'):
        f.append('identification read from the MMH page image')
    if 'is not a species' in idf:
        f.append('HS5000 code, unidentified (not checked)')
    if 'listed twice in the APGL98' in nt:
        f.append('APGL98 duplicate removed')
    if ('not attached here' in nt or 'no lab wavelength attached' in nt) and not d.get('lab_wl_air_A'):
        f.append('lab λ withheld: ' + ('EP not testable' if 'NIST gives no level energies' in nt else 'EP mismatch'))
    return '; '.join(f)


apg_all = viz(os.path.join(RAW, 'vizier', 'apgl98_table1.tsv'))
apgb = viz(os.path.join(RAW, 'vizier', 'apgl98b_fe1_shifts.tsv'))
iag = viz(os.path.join(RAW, 'vizier', 'iag2016_tablea1.tsv'))
print('APGL98', len(apg_all), 'APGL98b', len(apgb), 'IAG', len(iag))

# --- APGL98: drop rows that are exact copies (all fields except the VizieR recno)
apg, seen, dup_note = [], {}, {}
for r in apg_all:
    key = tuple(v for k, v in r.items() if k != 'recno')
    if key in seen:
        dup_note[seen[key]] = f"listed twice in the APGL98 VizieR table (recno {seen[key]} and {r['recno']}, identical in all fields); copy dropped"
        continue
    seen[key] = r['recno']
    apg.append(r)
print('APGL98 exact duplicates dropped:', len(apg_all) - len(apg))

# APGL98b lists one line twice (recno 333/334: same flux, lab, EP and EW; disc-centre 4343.2173 vs 4343.2710,
# APGL98 has 4343.2174): keep the first copy of each (Flux, Lab, EP)
apgb_dup_note, seen_b, apgb_u = {}, {}, []
for r in apgb:
    key = (r['Flux'], r['Lab'], r['EP'])
    if key in seen_b:
        apgb_dup_note[seen_b[key]] = (f"APGL98b lists this line twice (recno {seen_b[key]} and {r['recno']}; disc-centre "
                                      f"{[x for x in apgb if x['recno'] == seen_b[key]][0]['Int.']} / {r['Int.']}); second copy dropped")
        continue
    seen_b[key] = r['recno']
    apgb_u.append(r)
print('APGL98b duplicate lines dropped:', len(apgb) - len(apgb_u))
apgb = apgb_u
# APGL98b: sorted by flux wavelength for tolerant matching (the two catalogues differ by ~1e-4 A
# for the same measurements; a few APGL98b lines are absent from APGL98)
apgb.sort(key=lambda r: float(r['Flux']))
apgb_v = [float(r['Flux']) for r in apgb]
used_apgb = set()


def find_apgb(lF, lD, ep):
    """APGL98b row with the same flux/disc-centre wavelengths (+-1.2 mA) AND the same EP (+-0.02 eV)."""
    v, vd = float(lF), float(lD) if lD else None
    i = bisect.bisect_left(apgb_v, v - 0.0012)
    best, lam_only = None, None
    while i < len(apgb) and apgb_v[i] <= v + 0.0012:
        r = apgb[i]
        if vd is None or abs(float(r['Int.']) - vd) <= 0.0012:
            if ep and abs(float(r['EP']) - float(ep)) <= EP_TOL:
                if best is None or abs(apgb_v[i] - v) < abs(apgb_v[best] - v):
                    best = i
            else:
                lam_only = i
        i += 1
    return best, lam_only


# IAG tablea1 'lambda' is the Nave et al. (1994) rest wavelength in VACUUM nm (verified here: all 1251
# values match APGL98b air lab wavelengths to <1 mA only after vacuum->air conversion).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from airvac import vac_to_air
iag_air = sorted((float(vac_to_air(float(r['lambda']) * 10)), r) for r in iag)
iag_air_v = [x[0] for x in iag_air]


def find_iag(lab_air):
    v = float(lab_air); i = bisect.bisect_left(iag_air_v, v - 0.0012)
    best = None
    while i < len(iag_air) and iag_air_v[i] <= v + 0.0012:
        if best is None or abs(iag_air_v[i] - v) < abs(iag_air_v[best] - v):
            best = i
        i += 1
    return iag_air[best][1] if best is not None else None


def clean(s):
    return s.replace('=', '').replace('"', '').strip()


def level_ev(s):
    s = re.sub(r'[\[\]()?+a-z]', '', s or '')
    try:
        return float(s) / EV_CM
    except ValueError:
        return None


# NIST lab lines per species: (air wl float, verbatim string, unc, obs/Ritz, line ref, intensity, lower-level EP eV)
NIST = {}
for fn in glob.glob(os.path.join(RAW, 'nist_asd_solar_species', '*_3900_8000_air.csv')):
    sp = os.path.basename(fn).replace('_3900_8000_air.csv', '').replace('_', ' ')
    lines = [l for l in open(fn).read().splitlines() if not l.startswith('# query')]
    rd = list(csv.reader(lines)); hdr = [clean(h) for h in rd[0]]
    L = []
    for r in rd[1:]:
        d = {h: clean(v) for h, v in zip(hdr, r)}
        w = d.get('obs_wl_air(A)') or d.get('ritz_wl_air(A)')
        if not w:
            continue
        kind = 'obs' if d.get('obs_wl_air(A)') else 'Ritz'
        unc = d.get('unc_obs_wl') if kind == 'obs' else d.get('unc_ritz_wl')
        mnum = re.match(r'\d+\.\d+', w)   # NIST sometimes appends flags such as '+' or '?'
        if not mnum:
            continue
        e = [x for x in (level_ev(d.get('Ei(cm-1)')), level_ev(d.get('Ek(cm-1)'))) if x is not None]
        L.append((float(mnum.group(0)), w, unc or '', kind, d.get('line_ref', ''), d.get('intens', ''), min(e) if e else None))
    NIST[sp] = sorted(L)

# ---------------------------------------------------------------- HS5000 MMH digitisation + OCR cross-check
ocr = open(os.path.join(RAW, 'nbs_mono61', 'nbsmonograph61.txt'), encoding='utf-8', errors='replace').read()
ocr_tok = set(m.group(1) for m in re.finditer(r'(?<![\d.])(\d{4}\.\d{1,3})(?!\d)', ocr))
ocr_sorted = sorted(ocr_tok, key=float)
ocr_vals = [float(t) for t in ocr_sorted]


def ocr_resolve(s):
    """-> (printed string, check text, wl_flag).  Restores the printed MMH string when the OCR token
    differs from the HS5000 string only by trailing zeros (dropped or added by the digitisation)."""
    if s in ocr_tok:
        return s, 'OCR-match', ''
    ip, _, fr = s.partition('.')
    alts = [f'{ip}.{(fr + "000")[:k]}' for k in (1, 2, 3) if k > len(fr)]          # zeros dropped by HS5000
    alts += [f'{ip}.{fr[:k]}' for k in (2, 1) if k < len(fr) and set(fr[k:]) == {'0'}]  # zeros added by HS5000
    for a in alts:
        if a in ocr_tok:
            how = 'drops' if ndec(a) > ndec(s) else 'adds'
            return a, 'OCR-match', f'HS5000 string "{s}" {how} trailing zero(s); MMH prints {a} (OCR text)'
    v = float(s); i = bisect.bisect_left(ocr_vals, v)
    near = [j for j in (i - 1, i) if 0 <= j < len(ocr_vals)]
    j = min(near, key=lambda j: abs(ocr_vals[j] - v)) if near else None
    if j is not None and abs(ocr_vals[j] - v) < 0.03:
        return s, f'not in OCR; nearest OCR token {ocr_sorted[j]}', ''
    return s, 'not in OCR (OCR layer incomplete)', ''


# MMH entries checked by eye against the NBS Monograph 61 page images in this work (pdf page numbers)
MMH_PAGE_CHECK = {
    '5183.600': 'page image (pdf p.253) shows 5183.619m, EW 1584 -- HS5000 value is wrong',
    '4861.35': 'page image (pdf p.231) shows 4861.342m, EW 3680 -- HS5000 value/EW differ',
    '5168.91': 'page image (pdf p.252) shows 5168.908m, EW 114 (HS5000 EW 700 is wrong)',
    '5875.62': 'not in MMH 1966 (pdf p.292 lists only Atm H2O 5875.596); He I D3 added by HS5000',
}

# HS5000 integer values: the decimals were lost in the digitisation (on most of these pages the printed
# decimal point is faint or missing, e.g. "5698 340m").  Printed value read from the page image.
MMH_INT_FIX = {
    '4898': ('4898.473', 'pdf p.234 prints "4898 473r" (Cr I, EW 2)'),
    '4944': ('4944.564', 'pdf p.237 prints "4944 564" (Cr I, EW 5)'),
    '5522': ('5522.197', 'pdf p.272 prints "5522 197" (C2, EW 5.5)'),
    '5698': ('5698.340', 'pdf p.282 prints "5698 340m" (Cr I + Fe I, EW 30)'),
    '5903': ('5903.332', 'pdf p.294 prints 5903.332 (Ti I, EW 5.5)'),
    '5597': ('5597.876', 'pdf p.277 prints "5597.876" with the digits "87" faint (Cr I, EW 3, EP 3.84, mult. 239; NIST Cr I obs 5597.87)'),
}

mmh = []
for r in csv.DictReader(open(os.path.join(RAW, 'hs5000_github', 'moore_clean02012026.csv'))):
    s = r['wavelength'].strip()
    ident, el, iflag = clean_ident(r['id'], s, 'MMH')
    if s in MMH_PAGE_CHECK:
        wl, chk = s, MMH_PAGE_CHECK[s]
        wflag = (MMH_PAGE_CHECK[s] + '; HS5000 entry kept') if MMH_PAGE_CHECK[s].startswith('not in MMH') else \
            'HS5000 value differs from the MMH page image (kept verbatim): ' + MMH_PAGE_CHECK[s]
    elif s in MMH_INT_FIX:
        wl = MMH_INT_FIX[s][0]
        chk = 'page image: ' + MMH_INT_FIX[s][1]
        wflag = f'HS5000 integer "{s}" (decimals lost in the digitisation); printed value from the MMH page image: {MMH_INT_FIX[s][1]}'
    else:
        wl, chk, wflag = ocr_resolve(s)
    if ndec(wl) == 0:
        wflag = (wflag + '; ' if wflag else '') + 'integer value in the HS5000 digitisation; printed decimals unknown'
    mmh.append(dict(wl=wl, hs=s, v=float(wl), ew=r['ew'].strip(), id=ident, el=el, raw_id=r['id'].strip(), iflag=iflag,
                    chk=chk, wflag=wflag, flag=r['flag'].strip()))
mmh.sort(key=lambda x: x['v'])
mmh_v = [m['v'] for m in mmh]
mmh_chk = Counter()
for m in mmh:
    if m['chk'] == 'OCR-match':
        mmh_chk['OCR exact' if not m['wflag'] else ('OCR after restoring dropped zero(s)' if 'drops' in m['wflag'] else
                                                     'OCR after removing added zero(s)')] += 1
    elif m['chk'].startswith('not in OCR; nearest'):
        mmh_chk['nearest OCR token differs (OCR ?)'] += 1
    elif m['chk'].startswith('not in OCR'):
        mmh_chk['OCR layer gap'] += 1
    else:
        mmh_chk['checked on page image'] += 1
print('MMH digitisation vs OCR:', dict(mmh_chk), 'integer rows left:', sum(ndec(m['wl']) == 0 for m in mmh))


def add_note(d, t):
    d['notes'] = (d['notes'] + '; ' if d.get('notes') else '') + t


def attach_mmh(d):
    """MMH (disc centre) match: same element, within 0.012 A of the row's disc-centre value."""
    if not d.get('wl_disc_centre_air_A'):
        return False
    v = float(d['wl_disc_centre_air_A']); i = bisect.bisect_left(mmh_v, v - 0.012)
    best = None
    while i < len(mmh) and mmh[i]['v'] <= v + 0.012:
        if mmh[i]['el'] == d['element']:
            if best is None or abs(mmh[i]['v'] - v) < abs(mmh[best]['v'] - v):
                best = i
        i += 1
    if best is None:
        return False
    m = mmh[best]; used_mmh.add(best)
    d.update(mmh1966_wl_air_A=m['wl'], mmh1966_ew_mA=m['ew'], mmh1966_ident=m['id'], mmh1966_ocr_check=m['chk'])
    if m['wflag']:
        add_note(d, 'MMH value: ' + m['wflag'])
    if m['iflag']:
        add_note(d, 'MMH identification: ' + m['iflag'])
    if not d.get('ew_mA') and m['ew']:
        d.update(ew_mA=m['ew'], ew_source='MMH 1966 disc-centre EW (HS5000 digitisation)')
    return True


rows = []
used_mmh = set()
nist_stats = Counter()
for r in apg:
    ion = ion_name(r['Ion'])
    lF, lD, ep = r['lambdaF'], r['lambdaD'], r['EP']
    d = dict(wl_air_A=lF, medium='air (standard air; Kitt Peak FTS flux atlas, Kurucz et al. 1984)',
             wl_type='solar flux (Sun-as-a-star) line centre, 4th-order poly fit (APGL98)', e_wl_A=r['e_lambdaF'],
             wl_disc_centre_air_A=lD, e_wl_disc_centre_A=r['e_lambdaD'], ident=ion, element=ion.split()[0] if ion else '?',
             ep_eV=ep, loggf_solar=r['loggf'])
    if r['recno'] in dup_note:
        add_note(d, dup_note[r['recno']])
    src = ['APGL98']
    bi, lam_only = find_apgb(lF, lD, ep) if (r['Ion'].strip() == 'FE1' and lF) else (None, None)
    b = apgb[bi] if bi is not None else None
    if b:
        used_apgb.add(bi)
        if b['Flux'] != lF or b['Int.'] != lD:
            add_note(d, f"APGL98b lists flux {b['Flux']} / disc-centre {b['Int.']} for this line")
        d.update(lab_wl_air_A=b['Lab'], lab_source='Nave et al. 1994 (via APGL98b)', lab_unc_A='', lab_ep_eV=b['EP'])
        if b['recno'] in apgb_dup_note:
            add_note(d, apgb_dup_note[b['recno']])
        if b['EW']:
            d.update(ew_mA=f"{float(b['EW']) * 10:.0f}", ew_source='APGL98b disc-centre EW (pm x10)')
        src.append('APGL98b')
        q = find_iag(b['Lab'])
        if q:
            d.update(depth_flux=q['relDepth'], conv_blueshift_iag_mps=q['ConvBS'], lab_wl_vac_A=f"{float(q['lambda']) * 10:.3f}")
            if not d.get('ew_mA') and q['EW']:
                d.update(ew_mA=q['EW'], ew_source='IAG16 (from Allende Prieto+ 1998, disc centre)')
            src.append('IAG16')
    else:
        if lam_only is not None:
            nist_stats['apgb_ep_mismatch'] += 1
            add_note(d, f"APGL98b assigns this solar feature to Fe I EP {apgb[lam_only]['EP']} eV (Nave lab {apgb[lam_only]['Lab']}); "
                        f"EP differs by > {EP_TOL} eV, so that lab wavelength is not attached here")
        if ion in NIST and lF:
            v = float(lF)
            win = [x for x in NIST[ion] if -0.015 < v - x[0] < 0.025]
            cands = [x for x in win if x[6] is not None and ep and abs(x[6] - float(ep)) <= EP_TOL]
            if cands:
                exp = v / (1 + 300.0 / C)  # expect ~+0.3 km/s net solar shift
                best = min(cands, key=lambda x: abs(x[0] - exp))
                flag = 'unique' if len(cands) == 1 else f'{len(cands)} NIST candidates with matching EP; nearest to +0.3 km/s chosen'
                if len(cands) > 1:
                    nist_stats['multi'] += 1
                d.update(lab_wl_air_A=best[1], lab_unc_A=best[2], lab_ep_eV=f'{best[6]:.3f}',
                         lab_source=f'NIST ASD v5.12 {ion} ({best[3]}; ref {best[4] or "n/a"}); match: {flag}; EP agrees within {EP_TOL} eV')
                src.append('NIST_ASD')
                nist_stats['matched'] += 1
            elif win:
                nist_stats['ep_reject' if any(x[6] is not None for x in win) else 'ep_untestable'] += 1
                add_note(d, f"{len(win)} NIST {ion} line(s) within -15/+25 mA, none with lower-level EP within {EP_TOL} eV of {ep} eV"
                            + (' (NIST gives no level energies)' if all(x[6] is None for x in win) else '') + '; no lab wavelength attached')
            else:
                nist_stats['nomatch'] += 1
        elif lF:
            nist_stats['no_nist_species'] += 1
    if d.get('lab_wl_air_A') and lF:
        d['solar_flux_minus_lab_mps'] = f"{C * (float(lF) / float(re.match(r'[\d.]+', d['lab_wl_air_A']).group(0)) - 1):.0f}"
    if attach_mmh(d):
        src.append('MMH66_HS5000')
    if d.get('conv_blueshift_iag_mps') and d.get('lab_wl_air_A'):
        d['iag_flux_wl_air_A_derived'] = f"{float(d['lab_wl_air_A']) * (1 + (float(d['conv_blueshift_iag_mps']) + GRS) / C):.4f}"
    d['n_decimals'] = ndec(lF)
    d['source_key'] = '+'.join(src)
    rows.append(d)

# APGL98b Fe I lines not present in APGL98 with the same wavelength AND excitation potential
n_apgb_added = Counter()
for bi, b in enumerate(apgb):
    if bi in used_apgb:
        continue
    d = dict(wl_air_A=b['Flux'], medium='air (standard air; Kitt Peak FTS flux atlas, Kurucz et al. 1984)',
             wl_type='solar flux line centre (APGL98b)', wl_disc_centre_air_A=b['Int.'], ident='Fe I', element='Fe', ep_eV=b['EP'],
             loggf_solar=b['log(gf)'], lab_wl_air_A=b['Lab'], lab_source='Nave et al. 1994 (via APGL98b)', lab_ep_eV=b['EP'],
             n_decimals=ndec(b['Flux']), source_key='APGL98b')
    if b['recno'] in apgb_dup_note:
        add_note(d, apgb_dup_note[b['recno']])
    same = [r for r in apg if r['lambdaF'] and abs(float(r['lambdaF']) - float(b['Flux'])) <= 0.0012
            and (not r['lambdaD'] or abs(float(r['lambdaD']) - float(b['Int.'])) <= 0.0012)]
    if same:
        n_apgb_added['ep'] += 1
        d['notes'] = ('APGL98 lists this solar feature only as ' +
                      ', '.join(sorted({f"{ion_name(r['Ion'])} EP {r['EP']}" for r in same})) +
                      f'; APGL98b assigns it to this Fe I line (EP {b["EP"]} eV)')
    else:
        n_apgb_added['lambda'] += 1
    if b['EW']:
        d.update(ew_mA=f"{float(b['EW']) * 10:.0f}", ew_source='APGL98b disc-centre EW (pm x10)')
    d['solar_flux_minus_lab_mps'] = f"{C * (float(b['Flux']) / float(b['Lab']) - 1):.0f}"
    q = find_iag(b['Lab'])
    if q:
        d.update(depth_flux=q['relDepth'], conv_blueshift_iag_mps=q['ConvBS'], lab_wl_vac_A=f"{float(q['lambda']) * 10:.3f}")
        if not d.get('ew_mA') and q['EW']:
            d.update(ew_mA=q['EW'], ew_source='IAG16 (from Allende Prieto+ 1998, disc centre)')
        d['iag_flux_wl_air_A_derived'] = f"{float(b['Lab']) * (1 + (float(q['ConvBS']) + GRS) / C):.4f}"
        d['source_key'] += '+IAG16'
    if attach_mmh(d):
        d['source_key'] += '+MMH66_HS5000'
    rows.append(d)
print('APGL98b rows added:', dict(n_apgb_added))

# MMH lines not represented in APGL98
for i, m in enumerate(mmh):
    if i in used_mmh:
        continue
    d = dict(wl_air_A=m['wl'], medium='air (ISA scale)', wl_type='solar disc-centre, photographic (MMH 1966; HS5000 digitisation)',
             e_wl_A='', ident=m['id'], element=m['el'], ident_flag=m['iflag'], wl_flag=m['wflag'],
             ew_mA=m['ew'], ew_source='MMH 1966' if m['ew'] else '',
             mmh1966_wl_air_A=m['wl'], mmh1966_ew_mA=m['ew'], mmh1966_ident=m['id'], mmh1966_ocr_check=m['chk'],
             n_decimals=ndec(m['wl']), source_key='MMH66_HS5000',
             notes=('HS5000 flag ' + m['flag']) if m['flag'] else '')
    if m['id'] != m['raw_id']:
        d['ident_verbatim'] = m['raw_id']
    rows.append(d)

# Babcock & Moore 1947 above the end of the HS5000 MMH digitisation, unless an APGL98 line is within 0.02 A
apg_lam = sorted({(float(x), i) for i, d in enumerate(rows) if 'APGL98' in d['source_key']
                  for x in (d['wl_air_A'], d.get('wl_disc_centre_air_A')) if x})
apg_lam_v = [x[0] for x in apg_lam]
bm_stats = Counter()
# Ca II infrared triplet: MMH 1966 disc-centre EWs read from the NBS Monograph 61 OCR text (lines 38556, 38660, 38846)
BM_MMH_EW = {'8498.062': '1470', '8542.144': '3670', '8662.170': '2600'}
for r in csv.DictReader(open(os.path.join(RAW, 'hs5000_github', 'babcock_clean_02012026.csv'))):
    s = r['wav'].strip(); v = float(s)
    if v <= MMH_END:
        bm_stats['below MMH end'] += 1
        continue
    j = bisect.bisect_left(apg_lam_v, v - BM_APG_TOL)
    near = sorted({apg_lam[k][1] for k in range(j, len(apg_lam)) if apg_lam_v[k] <= v + BM_APG_TOL})
    ident, el, iflag = clean_ident(r['id'] or '', s, 'BM')
    strength = (r['strength'] or '').strip()
    if near:
        bm_stats['skipped (APGL98 within 0.02 A)'] += 1
        for k in near:
            add_note(rows[k], f'Babcock & Moore 1947 lists {s} ({ident or "unidentified"}, strength {strength}) here; not added as a row')
        continue
    nd = ndec(s); wflag = ''
    if nd > 3:
        wflag = f'HS5000 value has {nd} decimals, more digits than B&M 1947 printed; precision taken as 3 decimals'
        t3 = s[:s.index('.') + 4]
        if v <= 8770 and t3 in ocr_tok:
            wflag += f' (MMH 1966, which reprints B&M here, prints {t3}: OCR text)'
        nd = 3
    elif nd == 0:
        wflag = 'integer value in the HS5000 digitisation; printed decimals unknown (B&M 1947 page not available here)'
    d = dict(wl_air_A=s, medium='air (ISA scale)',
             wl_type='solar disc-centre, photographic (Babcock & Moore 1947; HS5000 digitisation)', ident=ident, element=el,
             ident_flag=iflag, wl_flag=wflag, n_decimals=nd, source_key='BM47_HS5000',
             notes=f"HS5000 strength (B&M Rowland-type intensity, not EW) = {strength}" + (f"; HS5000 flag {r['flag'].strip()}" if (r['flag'] or '').strip() else ''))
    if ident != (r['id'] or '').strip():
        d['ident_verbatim'] = (r['id'] or '').strip()
    if v <= 8770:   # MMH 1966 (to 8770 A) reprints B&M wavelengths above 6600 A: look for the same string in its OCR text
        if s in ocr_tok:
            bm_stats['same value in MMH 1966 OCR'] += 1
            add_note(d, 'same wavelength printed in MMH 1966 (OCR text)')
        else:
            bm_stats['not found in MMH 1966 OCR'] += 1
    if s in BM_MMH_EW:
        d.update(ew_mA=BM_MMH_EW[s], ew_source='MMH 1966 disc-centre EW (OCR text; not in the HS5000 digitisation)')
        add_note(d, 'the HS5000 strengths 1 / 25 / 100 for Ca II 8498 / 8542 / 8662 are not in the order of the MMH 1966 EWs '
                    '(1470 / 3670 / 2600 mA); B&M 1947 page not available here, strength not checked')
    rows.append(d)
    bm_stats['added'] += 1
print('B&M 1947:', dict(bm_stats))

rows.sort(key=lambda d: (float(d['wl_air_A']), d.get('ident', '')))
# one solar feature, several catalogue entries of the same species and EP (EP test cannot separate them):
# the same lab wavelength is attached to each; say so
shared = defaultdict(list)
for d in rows:
    if d.get('lab_wl_air_A'):
        shared[(d['wl_air_A'], d['lab_wl_air_A'], d.get('lab_ep_eV'), d.get('lab_source'))].append(d)
n_shared = 0
for grp in shared.values():
    if len(grp) > 1:
        n_shared += len(grp) - 1
        gfs = ', '.join(x.get('loggf_solar') or 'n/a' for x in grp)
        eps = ', '.join(x.get('ep_eV') or 'n/a' for x in grp)
        for x in grp:
            add_note(x, f"APGL98 lists {len(grp)} {x['ident']} entries at this solar feature (EP {eps} eV; log gf {gfs}), all within "
                        f'{EP_TOL} eV of the lab line; the EP test cannot tell them apart, so the same lab wavelength is given for each')
print('rows sharing a lab wavelength with another entry of the same EP (extra rows):', n_shared)
for d in rows:
    d['line_kind'] = line_kind(d.get('ident', ''))
    d['flag_short'] = short_flags(d)
fields = ['wl_air_A', 'medium', 'wl_type', 'e_wl_A', 'n_decimals', 'wl_flag', 'ident', 'element', 'line_kind', 'flag_short',
          'ident_flag', 'ident_verbatim',
          'ew_mA', 'ew_source', 'depth_flux',
          'lab_wl_air_A', 'lab_unc_A', 'lab_source', 'lab_ep_eV', 'lab_wl_vac_A', 'solar_flux_minus_lab_mps', 'conv_blueshift_iag_mps',
          'iag_flux_wl_air_A_derived',
          'wl_disc_centre_air_A', 'e_wl_disc_centre_A', 'ep_eV', 'loggf_solar',
          'mmh1966_wl_air_A', 'mmh1966_ew_mA', 'mmh1966_ident', 'mmh1966_ocr_check', 'source_key', 'notes']
for d in rows:
    for k, v in d.items():
        assert '#' not in str(v), (k, v)      # build_sets reads the CSV with comment='#'
with open(OUT, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fields)
    w.writeheader()
    for d in rows:
        w.writerow({k: d.get(k, '') for k in fields})
print('wrote', OUT, len(rows), 'rows;', dict(nist_stats), 'MMH used', len(used_mmh), 'of', len(mmh))
