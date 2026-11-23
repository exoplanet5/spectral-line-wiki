"""Build data/curated/airglow_nir.csv: near-IR (0.9-2.5 um) OH Meinel and O2 airglow lines
(plus the O2 b-X 762/865 nm bands and other O2 bands 0.74-2.5 um).

Primary line positions: HITRAN line-by-line records (HITRANonline, fetched 2026-11-23; 16OH
Meinel positions are the HITRAN2020 PGOPHER list = Brooke+2016 fit refit with Noll+2020 UVES
airglow lines; O2 from HITRAN2020/2024).  wl_vac_A = 1e8 / nu_vac (derived; nu verbatim kept).
Intensities: PALACE v1.0 reference intensities (Noll+2025, GMD 18, 4353), rayleigh.
Cross-checks: Oliva+2015 (A&A 581, A47, table1 OH, table2 O2) and Rousselot+2000 (A&A 354,
1134; ESO ISAAC list_v2.0.dat, Tables 2 and 3 of the paper).

Usage: python -I build_airglow_nir.py <lib_dir> <raw_dir> <out_csv> <out_meta_json> [<static_meta_json>]
"""
import sys, re, json, csv, math
from collections import defaultdict, Counter
import numpy as np

sys.path.insert(0, sys.argv[1])
from airglow_nir_lib import (read_par_files, oh_quanta, o2_quanta, IERR_NU, read_vizier_tsv,
                             parse_oliva_id, read_rousselot, vac_to_air_morton2000, palace_id_from_hitran)
from astropy.io import fits

RAW, OUT_CSV, OUT_META = sys.argv[2], sys.argv[3], sys.argv[4]
OH_MIN, OH_MAX = 9000.0, 25000.0      # task range 0.9-2.5 um for OH
O2_MIN, O2_MAX = 7400.0, 25000.0      # O2: include b-X(0-0) 762 nm and (0-1) 865 nm bands
stats = {}

# ------------------------------------------------------------------ HITRAN
hit = []
for pat, mol in ((RAW + '/hitran_lbl/OH16_*.par', 'OH'), (RAW + '/hitran_lbl/O2_*.par', 'O2')):
    for r in read_par_files(pat):
        r['mol'] = mol
        r['lamA'] = 1e8 / r['nu']
        if mol == 'OH':
            r.update(oh_quanta(r)); r['isotop'] = '16OH'
        else:
            q = o2_quanta(r)
            if q is None:
                continue
            r.update(q); r['isotop'] = {'1': '16O2', '2': '16O18O', '3': '16O17O'}[r['iso']]
        hit.append(r)
# de-duplicate identical records that may appear in two overlapping windows
seen = set(); hu = []
for r in hit:
    k = (r['molec'], r['iso'], r['nu_str'], r['Vu'], r['Vl'], r['Ql'])
    if k in seen:
        continue
    seen.add(k); hu.append(r)
hit = hu
stats['hitran_records_parsed'] = Counter(r['isotop'] for r in hit)

def hit_label(r):
    if r['mol'] == 'OH':
        par = r['sym'][0]                      # verified: first sym char = upper-state parity
        if r['F_u'] == r['F_l'] and r['dN'] == r['dJ']:
            return '%s%d%s(%.1f)' % (r['dJ'], r['F_u'], par, r['J_l'])
        return '%s%s%d%d%s(%.1f)' % (r['dN'], r['dJ'], r['F_u'], r['F_l'], par, r['J_l'])
    return '%s%s(%d)' % (r['dN'], r['dJ'], r['N_l'])

def hit_band(r):
    if r['mol'] == 'OH':
        return 'OH(%d-%d)' % (r['v_u'], r['v_l'])
    return 'O2 %s-%s(%d-%d)' % (r['st_u'], r['st_l'], r['v_u'], r['v_l'])

# sorted index per isotopologue for wavelength lookups
idx = defaultdict(list)
for i, r in enumerate(hit):
    idx[r['isotop']].append(i)
srt = {}
for iso, ii in idx.items():
    ii = sorted(ii, key=lambda i: hit[i]['lamA'])
    srt[iso] = (np.array([hit[i]['lamA'] for i in ii]), np.array(ii))

def nearest_hit(iso, lamA):
    arr, ii = srt[iso]
    j = np.searchsorted(arr, lamA)
    ks = [k for k in (j - 1, j) if 0 <= k < len(arr)]
    k = min(ks, key=lambda k: abs(arr[k] - lamA))
    return ii[k], arr[k] - lamA

# ------------------------------------------------------------------ PALACE (FITS double + ASCII verbatim)
pal_fits = fits.open(RAW + '/zenodo_palace_v1.0/unzipped/PALACE/src/palace/data/palace_lines.fits')[1].data
pal_ascii = {}
with open(RAW + '/zenodo_palace_v1.0/unzipped_pmd/PMD/pmd_reflineint_mol.dat') as fh:
    for l in fh:
        if l.startswith('#') or l.strip().startswith('lam'):
            continue
        p = l.split()
        pal_ascii[(p[4], p[0])] = dict(lam_str=p[0], trans=p[1], fH2O=p[2], I_str=p[3], ID=p[4],
                                       src=p[9], vclass=p[10])
PAL_ISO = {'OH': '16OH', 'QH': '18OH', 'O2': '16O2', 'OQ': '16O18O', 'OG': '16O17O'}
pal_rows = []
for r in pal_fits:
    chem = r['chem'].strip(); ID = r['ID'].strip(); lamA = r['lam'] * 1e4
    if chem == 'OH' and not (OH_MIN <= lamA <= OH_MAX):
        continue
    if chem == 'O2' and not (O2_MIN <= lamA <= O2_MAX):
        continue
    if chem not in ('OH', 'O2'):
        continue
    lam_str = '%.7f' % r['lam']
    a = pal_ascii.get((ID, lam_str))
    if a is None:   # tolerate last-digit formatting differences
        cands = [v for (i2, ls), v in pal_ascii.items() if i2 == ID and abs(float(ls) - r['lam']) < 2e-7]
        a = cands[0] if cands else None
    pal_rows.append(dict(ID=ID, chem=chem, isotop=PAL_ISO[ID[:2]], lamA=lamA, I=float(r['I']),
                         trans=float(r['trans']), ascii=a))
stats['palace_rows_in_range'] = Counter(p['isotop'] for p in pal_rows)

def palace_vib(ID):
    def v(c):
        if c == 'X':
            return 10
        if c in 'ABCDEF':
            return 10 + 'ABCDEF'.index(c)
        return int(c)
    return v(ID[4]), v(ID[5])

pal_match = {}       # hitran index -> palace row
pal_only = []
pal_dl = []
hit_by_pid = {}
for i, r in enumerate(hit):
    if r['mol'] == 'O2' and r['st_u'] not in ('a', 'b'):
        continue
    hit_by_pid.setdefault((r['isotop'], palace_id_from_hitran(r)), []).append(i)
for p in pal_rows:
    cand = hit_by_pid.get((p['isotop'], p['ID']), [])
    cand = [i for i in cand if abs(hit[i]['lamA'] - p['lamA']) < 0.05]
    if cand:
        i = min(cand, key=lambda i: abs(hit[i]['lamA'] - p['lamA']))
        if i in pal_match:
            raise RuntimeError('double PALACE match ' + p['ID'])
        pal_match[i] = p; pal_dl.append((p['lamA'] - hit[i]['lamA'], p['ID'])); continue
    pal_only.append(p)
stats['palace_matched_to_hitran'] = len(pal_match)
dd = np.array([x[0] for x in pal_dl])
stats['palace_minus_hitran_maxabs_A_OH'] = float(max(abs(x[0]) for x in pal_dl if x[1].startswith('OH')))
stats['palace_minus_hitran_mA_O2_by_band'] = {b: dict(n=int(len(v)), median=round(float(np.median(v)), 3), max_abs=round(float(np.max(np.abs(v))), 3))
    for b, v in ((b, [x[0] * 1e3 for x in pal_dl if x[1][:6] == b]) for b in sorted(set(x[1][:6] for x in pal_dl if not x[1].startswith('OH'))))}
stats['palace_not_in_hitran'] = Counter(p['isotop'] + ' ' + p['ID'][2:6] for p in pal_only)

# ------------------------------------------------------------------ Oliva+2015 table1 (OH, by label)
oh_lab = {}
for i, r in enumerate(hit):
    if r['mol'] == 'OH' and r['F_u'] == r['F_l'] and r['dN'] == r['dJ']:
        oh_lab[(r['v_u'], r['v_l'], r['dJ'], r['F_u'], r['sym'][0], r['J_l'])] = i
oliva = {}
oliva_unmatched = []
for row in read_vizier_tsv(RAW + '/vizier_oliva2015/table1.tsv'):
    for lk, ik in (('lambda1', 'Iden1'), ('lambda2', 'Iden2')):
        p = parse_oliva_id(row[ik])
        key = None if p is None else (p['v_u'], p['v_l'], p['dJ'], p['F'], p['parity'], p['J_l'])
        if key in oh_lab:
            oliva[oh_lab[key]] = dict(lam_str=row[lk].strip(), iden=row[ik].strip(), flux=row['Flux'].strip(),
                                      recno=row['recno'].strip())
        else:
            oliva_unmatched.append(row[ik])
stats['oliva2015_table1_components'] = 2 * len(read_vizier_tsv(RAW + '/vizier_oliva2015/table1.tsv'))
stats['oliva2015_matched_by_label'] = len(oliva)
stats['oliva2015_unmatched'] = oliva_unmatched

# ------------------------------------------------------------------ Rousselot+2000 list_v2.0 (OH, by doublet)
rou = read_rousselot(RAW + '/eso_isaac_rousselot2000/list_v2.0.dat')
groups = defaultdict(list)
for k, x in enumerate(rou):
    groups[x[2]].append(k)
rpairs = []           # (k1, k2) with k1 shorter wavelength
paired = set()
for s, ks in groups.items():
    ks = sorted(ks, key=lambda k: rou[k][1])
    if len(ks) == 2:
        if rou[ks[1]][1] - rou[ks[0]][1] < 30:
            rpairs.append(tuple(ks)); paired.update(ks)
        continue
    # larger groups: greedy nearest-neighbour pairing
    free = list(ks)
    cand = sorted(((rou[b][1] - rou[a][1], a, b) for ii, a in enumerate(free) for b in free[ii + 1:]))
    for dlam, a, b in cand:
        if dlam < 30 and a not in paired and b not in paired:
            rpairs.append((a, b)); paired.update((a, b))
rou_single = [k for k in range(len(rou)) if k not in paired]
# HITRAN main-branch doublets in Rousselot's stated domain (dv=2-5, v'<=9, J<=16.5) plus J up to 17.5 margin
hd = defaultdict(list)
for i, r in enumerate(hit):
    if (r['mol'] == 'OH' and r['F_u'] == r['F_l'] and r['dN'] == r['dJ'] and 2 <= r['v_u'] - r['v_l'] <= 5
            and r['v_u'] <= 9 and r['J_l'] <= 17.5 and r['lamA'] >= 8790):
        hd[(r['v_u'], r['v_l'], r['dJ'], r['F_u'], r['J_l'])].append(i)
hdoub = [(np.mean([hit[i]['lamA'] for i in v]), sorted(v, key=lambda i: hit[i]['lamA'])) for v in hd.values() if len(v) == 2]
hdoub.sort()
hm = np.array([x[0] for x in hdoub])
rmean = np.array([(rou[a][1] + rou[b][1]) / 2 for a, b in rpairs])
ro = np.argsort(rmean); rms = rmean[ro]
def nearest_idx(arr, x):
    j = np.searchsorted(arr, x); ks = [k for k in (j - 1, j) if 0 <= k < len(arr)]
    return min(ks, key=lambda k: abs(arr[k] - x))
rous = {}
dmean = []
for pi, (a, b) in enumerate(rpairs):
    if rmean[pi] < 8790:
        continue
    hk = nearest_idx(hm, rmean[pi])
    back = ro[nearest_idx(rms, hm[hk])]
    if back != pi or abs(hm[hk] - rmean[pi]) > 0.30:
        continue
    i1, i2 = hdoub[hk][1]
    dm = (rmean[pi] - hm[hk]) * 1e3
    dmean.append((dm, hit[i1]))
    for kr, ih in ((a, i1), (b, i2)):
        rous[ih] = dict(lam_str=rou[kr][0], inten=rou[kr][2], d_mean_mA=dm,
                        mean_vac=(rou[a][1] + rou[b][1]) / 2)
stats['rousselot_lines'] = len(rou)
stats['rousselot_lines_ge_8790A'] = int(sum(1 for x in rou if x[1] >= 8790))
stats['rousselot_pairs'] = len(rpairs)
stats['rousselot_unpaired'] = len(rou_single)
stats['rousselot_components_matched'] = len(rous)

# Rousselot 2000 Table 3 (vacuum -> air for doublet means / resolved components) from paper text
t3 = []
txt = open(RAW + '/papers/rousselot2000.txt').read().split('\n')
on = False
for l in txt:
    if 'Table 3.' in l:
        on = True
    if not on:
        continue
    nums = re.findall(r'(?<![\d.])(\d{4,5}\.\d{3})(?![\d.])', l)
    for j in range(0, len(nums) - 1, 2):
        v, a = float(nums[j]), float(nums[j + 1])
        if 0 < v - a < 8 and abs((v - a) / v - 2.73e-4) < 0.15e-4:
            t3.append((nums[j], nums[j + 1]))
stats['rousselot_table3_pairs_parsed'] = len(t3)
# attach: doublet mean (or single component) of matched Rousselot entries
t3_by_vac = {round(float(v), 3): (v, a) for v, a in t3}
air_attach = {}
for ih, d in rous.items():
    for key, how in ((round(d['mean_vac'], 3), 'doublet mean'), (round(float(d['lam_str']), 3), 'component')):
        if key in t3_by_vac:
            air_attach[ih] = (t3_by_vac[key], how); break
        # Rousselot rounding of a mean may differ by 1 in the last digit
        for dd in (-0.001, 0.001):
            k2 = round(key + dd, 3)
            if how == 'doublet mean' and k2 in t3_by_vac:
                air_attach[ih] = (t3_by_vac[k2], how + ' (+/-1 in last digit)'); break
        if ih in air_attach:
            break
stats['rousselot_table3_attached_components'] = len(air_attach)

# ------------------------------------------------------------------ O2 cross-checks
# Oliva+2015 table2 (measured GIANO, 2 decimals): O2, O2?, O2_Band-head
oliva_o2 = {}
ol2 = read_vizier_tsv(RAW + '/vizier_oliva2015/table2.tsv')
o2_strong = [i for i in pal_match if hit[i]['mol'] == 'O2']
o2_arr = np.array(sorted(o2_strong, key=lambda i: hit[i]['lamA']))
o2_lam = np.array([hit[i]['lamA'] for i in o2_arr])
ol2_d = []
for row in ol2:
    iden = row['Iden'].strip()
    if not iden.startswith('O2'):
        continue
    lam = float(row['lambda'])
    w = np.where(np.abs(o2_lam - lam) < 0.25)[0]
    if len(w) == 0:
        continue
    best = max(w, key=lambda k: pal_match[o2_arr[k]]['I'])      # brightest PALACE line within 0.25 A
    ih = o2_arr[best]
    d = (lam - hit[ih]['lamA']) * 1e3
    if ih not in oliva_o2 or abs(d) < abs((float(oliva_o2[ih]['lam_str']) - hit[ih]['lamA']) * 1e3):
        oliva_o2[ih] = dict(lam_str=row['lambda'].strip(), iden=iden, flux=row['Flux'].strip(),
                            flag=row['Flag'].strip(), nblend=int(len(w)))
    ol2_d.append(d)
stats['oliva2015_table2_O2_entries'] = int(sum(1 for r in ol2 if r['Iden'].strip().startswith('O2')))
stats['oliva2015_table2_O2_matched'] = len(oliva_o2)

# Rousselot 2000 Table 2 (O2 a-X(0-0), vacuum with air in parentheses, 2 decimals)
t2 = []
on = False; prev_vac = None
for l in txt:
    if l.startswith('Table 2.'):
        on = True; continue
    if on and l.startswith('1138'):
        break
    if not on:
        continue
    vac = re.findall(r'(?<![\d.(])(1\d{4}\.\d{2})(?![\d)])', l)
    air = re.findall(r'\((1\d{4}\.\d{2})\)', l)
    if vac and not air:
        prev_vac = vac
    elif air and prev_vac is not None and len(air) == len(prev_vac):
        t2.extend(zip(prev_vac, air)); prev_vac = None
stats['rousselot_table2_O2_pairs_parsed'] = len(t2)
rous_o2 = {}
ax00 = [i for i in o2_strong if hit[i]['isotop'] == '16O2' and hit[i]['st_u'] == 'a' and hit[i]['v_u'] == 0 and hit[i]['v_l'] == 0]
ax_lam = np.array([hit[i]['lamA'] for i in ax00])
for v, a in t2:
    w = np.where(np.abs(ax_lam - float(v)) < 0.25)[0]
    if len(w) == 0:
        continue
    k = min(w, key=lambda k: abs(ax_lam[k] - float(v)))
    ih = ax00[k]
    if ih not in rous_o2 or abs(float(v) - hit[ih]['lamA']) < abs(float(rous_o2[ih]['lam_str']) - hit[ih]['lamA']):
        rous_o2[ih] = dict(lam_str=v, air_str=a)

# ------------------------------------------------------------------ assemble rows
COLS = ['wl_vac_A', 'wl_vac_basis', 'wl_air_A', 'wl_air_note', 'wl_air_A_derived',
        'uncert_A', 'uncert_A_lo', 'uncert_basis', 'nu_vac_cm-1', 'hitran_ierr',
        'intensity_R', 'palace_trans_zenith', 'palace_ID', 'palace_src',
        'intensity_oliva2015_rel', 'intensity_rousselot2000_rel', 'A_einstein_s-1', 'E_lower_cm-1',
        'molecule', 'isotopologue', 'band', 'branch_label', 'F_upper', 'F_lower', 'parity_upper',
        'N_upper', 'N_lower', 'J_upper', 'J_lower', 'transition_type', 'source_key',
        'd_palace2025_mA', 'wl_vac_A_oliva2015', 'd_oliva2015_mA', 'oliva2015_iden',
        'wl_vac_A_rousselot2000', 'd_rousselot2000_mA', 'd_rousselot2000_doublet_mean_mA', 'notes']

def fmt(x, nd=4):
    return '' if x is None else ('%.*f' % (nd, x))

noxon = [i for i, r in enumerate(hit) if r['mol'] == 'O2' and r['st_u'] == 'b' and r['st_l'] == 'a']
stats['hitran_noxon_b-a_lines'] = len(noxon)
include = set(pal_match) | set(oliva) | set(rous) | set(oliva_o2) | set(rous_o2) | set(noxon)
rows = []
for i in include:
    h = hit[i]
    lo, hi_ = IERR_NU[int(h['ierr'][0])]
    lam = h['lamA']
    if h['mol'] == 'OH' and not (OH_MIN <= lam <= OH_MAX):
        continue
    if h['mol'] == 'O2' and not (O2_MIN <= lam <= O2_MAX):
        continue
    p = pal_match.get(i); o = oliva.get(i) or oliva_o2.get(i); ro_ = rous.get(i) or rous_o2.get(i)
    notes = []
    row = dict.fromkeys(COLS, '')
    row['wl_vac_A'] = fmt(lam)
    row['wl_vac_basis'] = 'derived: 1e8/nu_vac (HITRAN wavenumber, vacuum)'
    row['nu_vac_cm-1'] = h['nu_str']
    row['hitran_ierr'] = h['ierr']
    row['uncert_A'] = '%.2g' % (lam ** 2 * hi_ / 1e8) if np.isfinite(hi_) else ''
    row['uncert_A_lo'] = '%.2g' % (lam ** 2 * lo / 1e8)
    row['uncert_basis'] = 'HITRAN nu code %s: %g <= dnu < %g cm-1 (uncert_A = upper bound, uncert_A_lo = lower bound)' % (h['ierr'][0], lo, hi_)
    row['wl_air_A_derived'] = fmt(float(vac_to_air_morton2000(lam)))
    row['A_einstein_s-1'] = '%.4g' % h['A']
    row['E_lower_cm-1'] = '%.4f' % h['Elow']
    row['molecule'] = h['mol']; row['isotopologue'] = h['isotop']
    row['band'] = hit_band(h); row['branch_label'] = hit_label(h)
    row['J_upper'] = '%g' % h['J_u']; row['J_lower'] = '%g' % h['J_l']
    if h['mol'] == 'OH':
        row['F_upper'] = h['F_u']; row['F_lower'] = h['F_l']; row['parity_upper'] = h['sym'][0]
        row['N_upper'] = '%g' % (h['J_u'] - 0.5 if h['F_u'] == 1 else h['J_u'] + 0.5)
        row['N_lower'] = '%g' % (h['J_l'] - 0.5 if h['F_l'] == 1 else h['J_l'] + 0.5)
        row['transition_type'] = 'main' if (h['F_u'] == h['F_l'] and h['dN'] == h['dJ']) else 'satellite'
    else:
        row['N_upper'] = h['N_u']; row['N_lower'] = h['N_l']
        row['transition_type'] = {'d': 'magnetic dipole', 'q': 'electric quadrupole'}.get(h['mtype'], h['mtype'])
    row['source_key'] = 'HITRAN2024_OH_Meinel' if h['mol'] == 'OH' else 'HITRAN2024_O2'
    if i in noxon:
        notes.append('Noxon band b(0)->a(0), new in HITRAN2024; dayglow/twilight emission (O3 photolysis), not in PALACE nightglow model')
    if p is not None:
        row['intensity_R'] = '%.4g' % p['I']; row['palace_trans_zenith'] = '%.4f' % p['trans']
        row['palace_ID'] = p['ID']; row['palace_src'] = p['ascii']['src'] if p['ascii'] else ''
        row['source_key'] += ';PALACE2025'
        row['d_palace2025_mA'] = '%.2f' % (round((p['lamA'] - lam) * 1e3, 2) + 0.0)
    if o is not None:
        row['wl_vac_A_oliva2015'] = o['lam_str']
        row['d_oliva2015_mA'] = '%.1f' % ((float(o['lam_str']) - lam) * 1e3)
        row['oliva2015_iden'] = o['iden']
        row['intensity_oliva2015_rel'] = o['flux']
        row['source_key'] += ';Oliva2015'
        if h['mol'] == 'O2':
            notes.append('Oliva2015 table2 measured GIANO position (0.01 A format, ~0.07 A rms calib.), matched to brightest PALACE line within 0.25 A (%d lines in window)' % o['nblend'])
        else:
            notes.append('Oliva2015 table1 computed (Bernath & Colin 2009 constants); flux is for the Lambda doublet')
    if ro_ is not None:
        row['wl_vac_A_rousselot2000'] = ro_['lam_str']
        row['d_rousselot2000_mA'] = '%.1f' % ((float(ro_['lam_str']) - lam) * 1e3)
        row['source_key'] += ';Rousselot2000'
        if h['mol'] == 'OH':
            row['intensity_rousselot2000_rel'] = ro_['inten']
            row['d_rousselot2000_doublet_mean_mA'] = '%.1f' % ro_['d_mean_mA']
            if h['dJ'] == 'Q':
                notes.append('Rousselot2000 Q-branch Lambda-doublet splitting differs from HITRAN/Oliva; compare doublet mean')
        else:
            row['wl_air_A'] = ro_['air_str']
            row['wl_air_note'] = 'Rousselot2000 Table 2 air value (2 decimals; vac %s; Herzberg & Herzberg 1947 wavenumbers, Allen 1973 conversion)' % ro_['lam_str']
    if i in air_attach:
        (v, a), how = air_attach[i]
        row['wl_air_A'] = a
        row['wl_air_note'] = 'Rousselot2000 Table 3 air value for %s vac %s (Allen 1973 formula)' % (how, v)
    if h['mol'] == 'OH' and h['v_u'] in (5, 6) and h['J_l'] >= 15.5:
        notes.append('high-J v=5/6: HITRAN2020 positions refit to Noll+2020 airglow data; older lists deviate by up to ~0.7 A')
    row['notes'] = '; '.join(notes)
    rows.append(row)

# PALACE lines not in HITRAN (18OH, minor O2 bands from Cosby+2006)
for p in pal_only:
    row = dict.fromkeys(COLS, '')
    a = p['ascii']
    row['wl_vac_A'] = fmt(p['lamA'])
    row['wl_vac_basis'] = 'PALACE v1.0 palace_lines.fits lam x1e4 (ASCII release %s um)' % (a['lam_str'] if a else '?')
    row['wl_air_A_derived'] = fmt(float(vac_to_air_morton2000(p['lamA'])))
    row['uncert_basis'] = 'not stated by source (ASCII release gives 1e-7 um = 1e-3 A format)'
    row['intensity_R'] = '%.4g' % p['I']; row['palace_trans_zenith'] = '%.4f' % p['trans']
    row['palace_ID'] = p['ID']; row['palace_src'] = a['src'] if a else ''
    ID = p['ID']; vu, vl = palace_vib(ID)
    if p['chem'] == 'OH':
        row['molecule'] = 'OH'; row['band'] = 'OH(%d-%d)' % (vu, vl)
        row['F_upper'] = ID[8]; row['F_lower'] = ID[9]; row['N_upper'] = str(int(ID[10:12])); row['parity_upper'] = ID[12]
        Nu = int(ID[10:12]); Fu = int(ID[8]); Ju = Nu + 0.5 if Fu == 1 else Nu - 0.5
        row['J_upper'] = '%g' % Ju
        row['branch_label'] = 'PALACE:%s' % ID[6:]
    else:
        row['molecule'] = 'O2'; row['band'] = 'O2 %s-%s(%d-%d)' % (ID[2], ID[3], vu, vl)
        row['N_upper'] = str(int(ID[10:12])); row['branch_label'] = 'PALACE:%s' % ID[6:]
        row['transition_type'] = {'d': 'magnetic dipole', 'q': 'electric quadrupole'}.get(ID[12], ID[12])
    row['isotopologue'] = p['isotop']
    row['source_key'] = 'PALACE2025;Cosby2006'
    row['notes'] = 'not in HITRAN; PALACE takes wavelength from Cosby et al. 2006 (UVES) model positions'
    rows.append(row)

rows.sort(key=lambda r: float(r['wl_vac_A']))
with open(OUT_CSV, 'w', newline='') as fh:
    fh.write('# NIR airglow OH + O2 line list (0.9-2.5 um OH; 0.74-2.5 um O2). Vacuum wavelengths in Angstrom.\n')
    fh.write('# Primary positions: HITRAN2024 via HITRANonline (accessed 2026-11-23; 16OH Meinel positions unchanged since HITRAN2020) -> wl_vac_A = 1e8/nu (derived, 4 decimals); nu verbatim in nu_vac_cm-1.\n')
    fh.write('# intensity_R: PALACE v1.0 reference intensity in rayleigh (Cerro Paranal, zenith, nocturnal annual mean, 100 sfu, above lower-atmosphere absorption).\n')
    fh.write('# intensity_oliva2015_rel: Oliva+2015 observed relative photon flux of the Lambda doublet (brightest = 10000; O2: table2 scale). intensity_rousselot2000_rel: theoretical relative photon flux (Trot=190 K, Tvib=9000 K).\n')
    fh.write('# d_*_mA = (other source - primary) in milli-Angstrom. wl_air_A only where a source gives an air value; wl_air_A_derived = Morton (2000) eq. 8 (= Birch & Downs 1994) applied to wl_vac_A.\n')
    fh.write('# See airglow_nir_meta.json for sources, citations and statistics.\n')
    w = csv.DictWriter(fh, fieldnames=COLS); w.writeheader(); w.writerows(rows)
stats['rows_written'] = len(rows)
stats['rows_by_isotopologue'] = Counter(r['isotopologue'] for r in rows)
stats['rows_by_molecule_band_top'] = Counter(r['band'] for r in rows).most_common(60)

# ------------------------------------------------------------------ accuracy statistics
def summ(a):
    a = np.asarray(a, float)
    if len(a) == 0:
        return {}
    return dict(n=int(len(a)), median=round(float(np.median(a)), 2), mean=round(float(a.mean()), 2),
                p68_abs=round(float(np.percentile(np.abs(a), 68)), 2), p95_abs=round(float(np.percentile(np.abs(a), 95)), 2),
                max_abs=round(float(np.abs(a).max()), 1))
oh_rows = [r for r in rows if r['molecule'] == 'OH']
d_ol = [(float(r['d_oliva2015_mA']), float(r['J_lower']), r) for r in oh_rows if r['d_oliva2015_mA']]
stats['oliva2015_minus_hitran_mA_all'] = summ([x[0] for x in d_ol])
stats['oliva2015_minus_hitran_mA_Jlow<=8.5'] = summ([x[0] for x in d_ol if x[1] <= 8.5])
stats['oliva2015_minus_hitran_mA_Jlow>=12.5'] = summ([x[0] for x in d_ol if x[1] >= 12.5])
d_ro = [(float(r['d_rousselot2000_mA']), float(r['J_lower']), r) for r in oh_rows if r['d_rousselot2000_mA']]
stats['rousselot2000_minus_hitran_mA_components_PR'] = summ([x[0] for x in d_ro if not x[2]['branch_label'].startswith('Q')])
stats['rousselot2000_minus_hitran_mA_components_Q'] = summ([x[0] for x in d_ro if x[2]['branch_label'].startswith('Q')])
stats['rousselot2000_minus_hitran_mA_doublet_mean'] = summ([float(r['d_rousselot2000_doublet_mean_mA']) for r in oh_rows if r['d_rousselot2000_doublet_mean_mA'] and r['parity_upper'] == 'e'])
stats['rousselot2000_minus_hitran_mA_PR_Jlow<=8.5'] = summ([x[0] for x in d_ro if x[1] <= 8.5 and not x[2]['branch_label'].startswith('Q')])
both = [(float(r['wl_vac_A_rousselot2000']) - float(r['wl_vac_A_oliva2015'])) * 1e3 for r in oh_rows
        if r['wl_vac_A_rousselot2000'] and r['wl_vac_A_oliva2015'] and not r['branch_label'].startswith('Q')]
stats['rousselot2000_minus_oliva2015_mA_PR'] = summ(both)
o2r = [r for r in rows if r['molecule'] == 'O2']
stats['oliva2015_O2_measured_minus_hitran_mA'] = summ([float(r['d_oliva2015_mA']) for r in o2r if r['d_oliva2015_mA']])
stats['rousselot2000_O2_table2_minus_hitran_mA'] = summ([float(r['d_rousselot2000_mA']) for r in o2r if r['d_rousselot2000_mA']])
stats['hitran_nu_code_counts_OH'] = Counter(r['hitran_ierr'][:1] for r in oh_rows if r['hitran_ierr'])
stats['hitran_nu_code_counts_O2_by_band'] = {b: dict(Counter(r['hitran_ierr'][:1] for r in o2r if r['band'] == b and r['hitran_ierr']))
                                            for b in sorted(set(r['band'] for r in o2r))}
# morton vs palace-edlen air check at a few wavelengths (documentation only)
stats['air_formula_note'] = ('PALACE uses Edlen 1966 (n-1)*1e8 = 8342.13 + 2406030/(130-s2) + 15997/(38.9-s2); wl_air_A_derived uses '
                             'Morton 2000 eq. 8 (= Birch & Downs 1994); Morton - Edlen 1966 = -1.3e-4 to -2.9e-4 A over 1.0-2.2 um '
                             '(next entry); the page air column uses Peck & Reeder 1972')
lam_t = np.array([10000., 16000., 22000.])
s2 = (1e4 / lam_t) ** 2
n_ed = 1 + 1e-8 * (8342.13 + 2406030 / (130 - s2) + 15997 / (38.9 - s2))
stats['air_morton_minus_edlen1966_A_at_1.0_1.6_2.2um'] = [round(float(x), 6) for x in (vac_to_air_morton2000(lam_t) - lam_t / n_ed)]
# Rousselot Table 3 air vs Morton derived from Rousselot vac (checks their Allen 1973 conversion)
t3d = [(float(a) - float(vac_to_air_morton2000(float(v)))) * 1e3 for v, a in t3]
stats['rousselot_table3_air_minus_morton2000_from_same_vac_mA'] = summ(t3d)

meta = json.load(open(sys.argv[5])) if len(sys.argv) > 5 else {}
# wording corrections applied on top of the static metadata (verification pass 2026-11-23)
if 'medium' in meta:
    meta['medium'] = meta['medium'].replace('with the IAU standard formula of Morton (2000, ApJS 130, 403, eq. 8).',
                                            'with Morton (2000, ApJS 130, 403) eq. 8 (= Birch & Downs 1994); the IAU-adopted standard-air formula is Edlen 1953.')
if 'columns' in meta:
    meta['columns']['wl_vac_A'] = meta['columns']['wl_vac_A'].replace(
        'PALACE-only rows: PALACE lam [um] x 1e4.',
        'PALACE-only rows: palace_lines.fits lam [um, double] x 1e4; the PMD ASCII release prints the same value rounded to 1e-7 um '
        '(0.001 A), so the 4th decimal is not a published digit.')
if 'caveats' in meta:
    meta['caveats'] = [c.replace('Their wavelengths are Cosby et al. 2006 model positions.',
                                 'Their wavelengths are Cosby et al. 2006 model positions as carried by PALACE (palace_lines.fits; the ASCII '
                                 'release rounds them to 0.001 A).') for c in meta['caveats']]
meta['stats'] = stats
json.dump(meta, open(OUT_META, 'w'), indent=1, ensure_ascii=False, default=lambda o: dict(o) if isinstance(o, Counter) else str(o))
print(json.dumps(stats, indent=1, default=lambda o: dict(o) if isinstance(o, Counter) else str(o)))
