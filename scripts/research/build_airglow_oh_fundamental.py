"""Build data/curated/airglow_nir_oh_fundamental.csv: 16OH Meinel lines between 2.5 and 4.2 um
(Delta v = 1 fundamental bands (1-0)...(9-8), plus the few Delta v = 2 P-branch lines of
(7-5)...(9-7) that fall beyond 2.5 um), in the column layout of data/curated/airglow_nir.csv.

Positions: HITRAN2024 line-by-line records (HITRANonline, global isotopologue id 48 = 16OH, accessed
2026-11-23): raw/hitran_lbl_dv1/OH16_2000-2900cm.par, OH16_2900-3800cm.par and the earlier
raw/hitran_lbl/OH16_3800-4800cm.par.  wl_vac_A = 1e8 / nu_vac (derived; nu verbatim kept).
Intensities: none.  PALACE v1.0 (Noll et al. 2025) is valid only from 0.3 to 2.5 um.
HITRAN2024 keeps the HITRAN2020 16OH wavenumbers and replaces the intensities / A values by those of
Ushakov et al. (HITRAN2024 paper, Sect. 2.13).  Consistency check: every selected record is compared with
raw/airglow_atomic/bands/hitran_oh/oh48_0-30000.par, which was fetched from the same API earlier the
same day (so it is HITRAN2024 content as well, whatever its label elsewhere); expected: identical.

Selection (see SELECTION below and the CSV header for the numbers that justify it):
  16OH, v' <= 9, Delta v = 1 or 2, 25000 < lambda_vac <= 42000 A, N'' <= 20, A >= 0.01 s-1.

Usage: python -I build_airglow_oh_fundamental.py <lib_dir> <raw_nir_dir> <oh48_par> <nir_csv> <out_csv> <out_meta_json>
  lib_dir: scripts/research (airglow_nir_lib.py, airvac_formulas.py); nir_csv: data/curated/airglow_nir.csv
  (used for the column list, the overlap check and the PALACE-based justification of the N'' cut).
"""
import sys, csv, json, glob
from collections import Counter, defaultdict
import numpy as np

sys.path.insert(0, sys.argv[1])
from airglow_nir_lib import read_par_files, oh_quanta, IERR_NU, vac_to_air_morton2000   # noqa: E402
import airvac_formulas as AV                                                             # noqa: E402

RAW, OH48, NIR_CSV, OUT_CSV, OUT_META = sys.argv[2:7]
LAM_MIN, LAM_MAX = 25000.0, 42000.0     # exclusive / inclusive, A (airglow_nir.csv ends at 25000 A inclusive)
VU_MAX, NL_MAX, A_MIN = 9, 20, 0.01
SELECTION = (f"16OH, v' <= {VU_MAX}, Delta v = 1 (plus Delta v = 2 lines in the window), {LAM_MIN:.0f} < lambda_vac <= {LAM_MAX:.0f} A, "
             f"N'' <= {NL_MAX}, A >= {A_MIN} s-1")
NOTE = 'no PALACE intensity (PALACE ends at 2.5 um)'
stats = {'selection': SELECTION}

# ------------------------------------------------------------------ column layout = airglow_nir.csv
with open(NIR_CSV) as fh:
    nir_rows = list(csv.DictReader(l for l in fh if not l.startswith('#')))
COLS = list(nir_rows[0].keys())


def n_lower(J, F):
    return J - 0.5 if F == 1 else J + 0.5        # OH X2Pi: F1 = 2Pi3/2 (J = N + 1/2), F2 = 2Pi1/2 (J = N - 1/2)


def read_oh(paths):
    seen, out = set(), []
    for p in paths:
        for r in read_par_files(p):
            if r['molec'] != 13 or r['iso'] != '1':
                continue
            k = (r['nu_str'], r['Vu'], r['Vl'], r['Ql'])
            if k in seen:
                continue
            seen.add(k)
            r.update(oh_quanta(r))
            r['lamA'] = 1e8 / r['nu']
            r['N_l'] = n_lower(r['J_l'], r['F_l']); r['N_u'] = n_lower(r['J_u'], r['F_u'])
            r['main'] = r['F_u'] == r['F_l'] and r['dN'] == r['dJ']
            out.append(r)
    return out


def hit_label(r):          # identical to build_airglow_nir.py
    par = r['sym'][0]
    if r['main']:
        return '%s%d%s(%.1f)' % (r['dJ'], r['F_u'], par, r['J_l'])
    return '%s%s%d%d%s(%.1f)' % (r['dN'], r['dJ'], r['F_u'], r['F_l'], par, r['J_l'])


# ------------------------------------------------------------------ HITRAN2024
files24 = sorted(glob.glob(RAW + '/hitran_lbl_dv1/OH16_*.par')) + [RAW + '/hitran_lbl/OH16_3800-4800cm.par']
hit = read_oh(files24)
stats['hitran2024_files'] = [f.split('/raw/')[-1] for f in files24]
stats['hitran2024_records_2000-4800cm'] = len(hit)
stats['hitran2024_records_by_dv'] = dict(sorted(Counter(r['v_u'] - r['v_l'] for r in hit).items()))
win = [r for r in hit if LAM_MIN < r['lamA'] <= LAM_MAX and 1 <= r['v_u'] - r['v_l'] <= 2 and r['v_u'] <= VU_MAX]
sel = [r for r in win if r['N_l'] <= NL_MAX and r['A'] >= A_MIN]
stats['window_records_vu<=9_dv1-2'] = len(win)
stats['dropped_N_lower>%d' % NL_MAX] = sum(1 for r in win if r['N_l'] > NL_MAX)
stats['dropped_A<%g_with_N_lower<=%d' % (A_MIN, NL_MAX)] = dict(Counter('main' if r['main'] else 'satellite'
                                                                    for r in win if r['N_l'] <= NL_MAX and r['A'] < A_MIN))
stats['selected'] = len(sel)
stats['selected_by_band'] = {f'OH({a}-{b})': n for (a, b), n in sorted(Counter((r['v_u'], r['v_l']) for r in sel).items(),
                                                                         key=lambda kv: (kv[0][0] - kv[0][1], kv[0][0]))}
stats['selected_by_type'] = dict(Counter('main' if r['main'] else 'satellite' for r in sel))
stats['selected_hitran_nu_code'] = dict(sorted(Counter(r['ierr'][0] for r in sel).items()))
stats['selected_lambda_range_A'] = [round(min(r['lamA'] for r in sel), 4), round(max(r['lamA'] for r in sel), 4)]
stats['selected_A_range_s-1'] = {'main': [min(r['A'] for r in sel if r['main']), max(r['A'] for r in sel if r['main'])],
                                 'satellite': [min(r['A'] for r in sel if not r['main']), max(r['A'] for r in sel if not r['main'])]}
# bands cut by the 4.2 um limit (lines with N'' <= 20, A >= 0.01 beyond 42000 A)
beyond = Counter(f"OH({r['v_u']}-{r['v_l']})" for r in hit if r['lamA'] > LAM_MAX and r['v_u'] - r['v_l'] == 1
                 and r['v_u'] <= VU_MAX and r['N_l'] <= NL_MAX and r['A'] >= A_MIN)
stats['same_cuts_beyond_4.2um_not_listed'] = dict(sorted(beyond.items()))
stats['note_download_lower_edge'] = '2000 cm-1 (5.0 um); the 4.2 um window edge lies well inside the downloaded range'

# ------------------------------------------------------------------ justification of the N'' cut from PALACE (Delta v = 2, 3)
oh = [x for x in nir_rows if x['molecule'] == 'OH' and x['isotopologue'] == '16OH' and x['intensity_R']]
def dv_of(b):
    u, l = b[3:-1].split('-')
    return int(u) - int(l)
oh23 = [x for x in oh if dv_of(x['band']) in (2, 3)]
Itot = sum(float(x['intensity_R']) for x in oh23)
just = {}
for N in (10, 15, 20, 25):
    keep = [x for x in oh23 if float(x['N_lower']) <= N]
    just[f"N''<={N}"] = round(sum(float(x['intensity_R']) for x in keep) / Itot, 4)
keepA = [x for x in oh23 if float(x['N_lower']) <= NL_MAX and float(x['A_einstein_s-1']) >= A_MIN]
just[f"N''<={NL_MAX} and A>={A_MIN}"] = round(sum(float(x['intensity_R']) for x in keepA) / Itot, 5)
just['brightest_dropped_by_A_cut_R'] = max(float(x['intensity_R']) for x in oh23
                                           if float(x['N_lower']) <= NL_MAX and float(x['A_einstein_s-1']) < A_MIN)
just['brightest_line_R'] = max(float(x['intensity_R']) for x in oh23)
just['brightest_beyond_N_cut_R'] = max(float(x['intensity_R']) for x in oh23 if float(x['N_lower']) > NL_MAX)
just['n_lines'] = len(oh23)
stats['palace_dv2_3_intensity_fraction_retained_by_cuts'] = just

# ------------------------------------------------------------------ download consistency (same transition, by quanta)
h20 = read_oh([OH48])
k20 = defaultdict(list)
for r in h20:
    k20[(r['v_u'], r['v_l'], r['br'], r['F_u'], r['F_l'], r['J_l'], r['sym'])].append(r)
dnu, dA, unmatched, multi = [], [], 0, 0
for r in sel:
    c = k20.get((r['v_u'], r['v_l'], r['br'], r['F_u'], r['F_l'], r['J_l'], r['sym']), [])
    c = [x for x in c if abs(x['nu'] - r['nu']) < 0.5]
    if not c:
        unmatched += 1; continue
    if len(c) > 1:
        multi += 1
    x = min(c, key=lambda x: abs(x['nu'] - r['nu']))
    r['nu20'] = x['nu']; r['A20'] = x['A']
    dnu.append(r['nu'] - x['nu']); dA.append(r['A'] / x['A'] - 1)
dnu = np.array(dnu); dA = np.array(dA)
lam_sel = np.array([r['lamA'] for r in sel if 'nu20' in r])
dlam_mA = -lam_sel ** 2 * dnu / 1e8 * 1e3
stats['consistency_vs_oh48_0-30000.par'] = dict(
    matched=int(len(dnu)), unmatched=unmatched, multiple_candidates=multi,
    max_abs_dnu_cm=float(np.max(np.abs(dnu))) if len(dnu) else None,
    max_abs_dlam_mA=float(np.max(np.abs(dlam_mA))) if len(dnu) else None,
    n_position_changed=int(np.sum(np.abs(dnu) > 0)),
    A_ratio_minus_1_median=round(float(np.median(dA)), 4), A_ratio_minus_1_p05_p95=[round(float(np.percentile(dA, 5)), 4), round(float(np.percentile(dA, 95)), 4)])

# ------------------------------------------------------------------ air-formula extrapolation check
lam_t = np.array([25000.0, 30000.0, 35000.0, 42000.0])
stats['air_morton2000_minus_peckreeder1972_5p_mA'] = dict(zip(['2.5um', '3.0um', '3.5um', '4.2um'],
    [round(float(x) * 1e3, 3) for x in (vac_to_air_morton2000(lam_t) - AV.vac_to_air_peckreeder1972_5p(lam_t))]))
stats['air_morton2000_minus_ciddor1996_mA'] = dict(zip(['2.5um', '3.0um', '3.5um', '4.2um'],
    [round(float(x) * 1e3, 3) for x in (vac_to_air_morton2000(lam_t) - AV.vac_to_air_ciddor1996(lam_t))]))

# ------------------------------------------------------------------ rows
nir_nu = {x['nu_vac_cm-1'] for x in nir_rows if x['molecule'] == 'OH' and x['nu_vac_cm-1']}
rows = []
for h in sorted(sel, key=lambda r: r['lamA']):
    lo, hi_ = IERR_NU[int(h['ierr'][0])]
    lam = h['lamA']
    row = dict.fromkeys(COLS, '')
    row['wl_vac_A'] = '%.4f' % lam
    row['wl_vac_basis'] = 'derived: 1e8/nu_vac (HITRAN wavenumber, vacuum)'
    row['nu_vac_cm-1'] = h['nu_str']
    row['hitran_ierr'] = h['ierr']
    row['uncert_A'] = '%.2g' % (lam ** 2 * hi_ / 1e8) if np.isfinite(hi_) else ''
    row['uncert_A_lo'] = '%.2g' % (lam ** 2 * lo / 1e8)
    row['uncert_basis'] = 'HITRAN nu code %s: %g <= dnu < %g cm-1 (uncert_A = upper bound, uncert_A_lo = lower bound)' % (h['ierr'][0], lo, hi_)
    row['wl_air_A_derived'] = '%.4f' % float(vac_to_air_morton2000(lam))
    row['A_einstein_s-1'] = '%.4g' % h['A']
    row['E_lower_cm-1'] = '%.4f' % h['Elow']
    row['molecule'] = 'OH'; row['isotopologue'] = '16OH'
    row['band'] = 'OH(%d-%d)' % (h['v_u'], h['v_l']); row['branch_label'] = hit_label(h)
    row['F_upper'] = h['F_u']; row['F_lower'] = h['F_l']; row['parity_upper'] = h['sym'][0]
    row['N_upper'] = '%g' % h['N_u']; row['N_lower'] = '%g' % h['N_l']
    row['J_upper'] = '%g' % h['J_u']; row['J_lower'] = '%g' % h['J_l']
    row['transition_type'] = 'main' if h['main'] else 'satellite'
    row['source_key'] = 'HITRAN2024_OH_Meinel'
    notes = [NOTE]
    if h['v_u'] - h['v_l'] == 2:
        notes.append('Delta v = 2 line beyond 2.5 um (outside airglow_nir.csv)')
    row['notes'] = '; '.join(notes)
    if h['nu_str'] in nir_nu:
        raise RuntimeError('line already in airglow_nir.csv: ' + h['nu_str'])
    rows.append(row)

cc = stats['selected_hitran_nu_code']
J = stats['palace_dv2_3_intensity_fraction_retained_by_cuts']
f15, f20 = J["N''<=15"], J["N''<=%d" % NL_MAX]
with open(OUT_CSV, 'w', newline='') as fh:
    n2 = sum(1 for r in sel if r['v_u'] - r['v_l'] == 2)
    fh.write(f"# OH Meinel airglow lines 2.5-4.2 um: Delta v = 1 fundamental bands (1-0)...(9-8) of 16OH ({len(sel) - n2} lines), plus {n2} high-N Delta v = 2 lines of (7-5), (8-6), (9-7) beyond 2.5 um. Vacuum wavelengths in Angstrom.\n")
    fh.write("# Primary positions: HITRAN2024 via HITRANonline (global iso id 48, accessed 2026-11-23; raw/airglow_nir/hitran_lbl_dv1 + hitran_lbl/OH16_3800-4800cm.par) -> wl_vac_A = 1e8/nu (derived, 4 decimals); nu verbatim in nu_vac_cm-1.\n")
    fh.write(f"# Selection: {SELECTION}. Justification from the PALACE v1.0 OH model in airglow_nir.csv (Delta v = 2-3, {J['n_lines']} lines): N'' <= 15 keeps {f15:.1%} of the intensity, "
             f"N'' <= {NL_MAX} keeps {f20:.1%} (hot rotational population; brightest line beyond N'' = {NL_MAX} is {J['brightest_beyond_N_cut_R']:.3g} R vs {J['brightest_line_R']:.4g} R); "
             f"A >= {A_MIN} s-1 removes only weak satellites and high-J Q lines (brightest PALACE line so removed: {J['brightest_dropped_by_A_cut_R']:.2g} R).\n")
    fh.write(f"# intensity_R and all palace_* columns blank: {NOTE}. A_einstein_s-1 and E_lower_cm-1: HITRAN2024 verbatim values (%.4g / %.4f formatting).\n")
    fh.write(f"# uncert_A = upper bound of the HITRAN position-code interval (lambda^2 dnu_hi / 1e8): code 3 (1e-3..1e-2 cm-1) 0.06-0.18 A, code 2 (1e-2..0.1) 0.6-1.8 A, code 1 (0.1..1) 6-18 A over 2.5-4.2 um. Rows per code: " + ', '.join(f'code {k}: {v}' for k, v in cc.items()) + '.\n')
    fh.write("# wl_air_A_derived = Morton (2000) eq. 8 (= Birch & Downs 1994) applied to wl_vac_A, i.e. extrapolated beyond the wavelengths the standard-air formulas were fitted to; "
             f"Morton - Peck & Reeder 1972 = {stats['air_morton2000_minus_peckreeder1972_5p_mA']['3.5um']:+.2f} mA, Morton - Ciddor 1996 = {stats['air_morton2000_minus_ciddor1996_mA']['3.5um']:+.2f} mA at 3.5 um. Use vacuum in the IR.\n")
    fh.write("# HITRAN2024 keeps the HITRAN2020 16OH wavenumbers (Brooke et al. 2016 line list with the Noll et al. 2020 airglow position corrections); "
             "its A values are new (Ushakov et al.; HITRAN2024 Sect. 2.13: low-J P-line intensities changed by ~2 %). "
             "Thermal-IR emitters (NO 5.3 um, CO2 4.3/15 um, O3 9.6 um) are not included. See airglow_nir_oh_fundamental_meta.json.\n")
    w = csv.DictWriter(fh, fieldnames=COLS); w.writeheader(); w.writerows(rows)
stats['rows_written'] = len(rows)

meta = {
    'file': 'data/curated/airglow_nir_oh_fundamental.csv',
    'title': 'OH Meinel airglow lines 2.5-4.2 um (Delta v = 1 fundamental plus the Delta v = 2 tail beyond 2.5 um; HITRAN2024 positions, no intensities)',
    'built_by': 'scripts/research/build_airglow_oh_fundamental.py (python -I); parsers in scripts/research/airglow_nir_lib.py',
    'columns': 'identical to data/curated/airglow_nir.csv (see airglow_nir_meta.json "columns"); intensity, PALACE and cross-check columns are blank',
    'caveats': [
        'No intensities: PALACE v1.0 covers 0.3-2.5 um only, and its Einstein-A corrections for OH were derived for Delta v >= 2 bands.',
        'HITRAN position codes are interval bounds, not measured errors; at 2.5-4.2 um code 2 means 0.6-1.8 A upper bound.',
        'Large parts of 2.5-2.9 um and of 4.2-4.4 um are opaque from the ground (H2O, CO2). The L-band sky beyond ~3.5-4 um is increasingly dominated by thermal emission of the atmosphere and telescope.',
        'Bands (8-7) and (9-8) are cut by the 4.2 um limit (their P branches extend to ~5 um).',
        'Thermal-IR emitters (NO 5.3 um, CO2 4.3 and 15 um, O3 9.6 um) are out of scope.',
    ],
    'stats': stats,
}
json.dump(meta, open(OUT_META, 'w'), indent=1, ensure_ascii=False)
print(json.dumps(stats, indent=1))
