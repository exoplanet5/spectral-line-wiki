#!/usr/bin/env python
"""Build data/curated/fraunhofer_classic.csv -- the classic Fraunhofer letter designations with
modern laboratory wavelengths (air + vacuum) and solar-atlas wavelengths.

Every number is taken from a downloaded source file (NIST ASD v5.12 CSV, HITRAN2020 .par, VizieR
tables) or transcribed verbatim from a scanned source (Moore, Minnaert & Houtgast 1966 = NBS
Monograph 61; Smithsonian Physical Tables 1920 Table 312 = Rowland's values; Wikipedia), or is
derived by a documented formula (air<->vacuum: Peck & Reeder 1972; IAG line-core fits).

Usage: python build_fraunhofer_classic.py <root>   (root = project directory)
"""
import csv, glob, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from airvac import vac_to_air
from hitran_lookup import parse_par

ROOT = sys.argv[1]
RAW = os.path.join(ROOT, 'data', 'raw', 'fraunhofer_solar')
OUT = os.path.join(ROOT, 'data', 'curated', 'fraunhofer_classic.csv')
C_MS = 299792458.0


def clean(s):
    return s.replace('=', '').replace('"', '').strip()


def nist_rows(sp, medium):
    rows = []
    for fn in glob.glob(os.path.join(RAW, 'nist_asd_classic', f"{sp.replace(' ', '_')}_*_{medium}.csv")):
        lines = [l for l in open(fn).read().splitlines() if not l.startswith('# query')]
        rd = list(csv.reader(lines))
        hdr = [clean(h) for h in rd[0]]
        for r in rd[1:]:
            d = {h: clean(v) for h, v in zip(hdr, r)}
            if d not in rows:          # windows may overlap -> de-duplicate
                rows.append(d)
    return rows


def nist(sp, wl_str):
    """Return lab values for the NIST line whose observed or Ritz AIR wavelength string == wl_str."""
    air = [r for r in nist_rows(sp, 'air') if wl_str in (r.get('obs_wl_air(A)'), r.get('ritz_wl_air(A)'))]
    if not air:
        raise SystemExit(f'NIST line not found: {sp} {wl_str}')
    a = air[0]
    kind = 'obs' if a.get('obs_wl_air(A)') == wl_str else 'Ritz'
    keys = ('wn(cm-1)', 'Aki(s^-1)', 'Ek(cm-1)', 'J_i', 'J_k', 'term_k', 'line_ref')
    vac = [r for r in nist_rows(sp, 'vac') if all(r.get(k, '') == a.get(k, '') for k in keys)]
    if len(vac) != 1:
        raise SystemExit(f'ambiguous/missing vacuum partner for {sp} {wl_str}: {len(vac)}')
    v = vac[0]
    if kind == 'obs' and v.get('obs_wl_vac(A)'):
        wv, uv = v['obs_wl_vac(A)'], v.get('unc_obs_wl', '')
        ua = a.get('unc_obs_wl', '')
    else:
        if kind == 'obs':  # air observed but vacuum observed missing -> use Ritz in both for consistency
            kind = 'Ritz'
            wl_str = a.get('ritz_wl_air(A)')
        wv, uv = v.get('ritz_wl_vac(A)'), v.get('unc_ritz_wl', '')
        ua = a.get('unc_ritz_wl', '')
    term = f"{a.get('conf_i','')} {a.get('term_i','')} {a.get('J_i','')} - {a.get('conf_k','')} {a.get('term_k','')} {a.get('J_k','')}".strip()
    if term.replace('-', '').strip() == '':
        term = 'unresolved fine-structure centroid (NIST blend row ' + a.get('line_ref', '') + ')'
    return dict(wl_air=wl_str, wl_vac=wv, unc_air=ua, unc_vac=uv, kind=kind, transition=term,
                lab_src=f"NIST ASD v5.12 {sp} ({kind}; line ref {a.get('line_ref','') or 'n/a'}; Ei={a.get('Ei(cm-1)','')} cm-1)")


HIT = []
for fn in glob.glob(os.path.join(RAW, 'hitran', '*.par')):
    HIT += parse_par(fn)
IERR = {3: 1e-2, 4: 1e-3, 5: 1e-4, 6: 1e-5, 7: 1e-6}


def hitran(nu_str):
    r = [h for h in HIT if h['nu_str'] == nu_str][0]
    lv = 1e8 / r['nu']
    la = float(vac_to_air(lv))
    dnu = IERR.get(int(r['ierr'][0]))
    unc = lv * lv * dnu / 1e8 if dnu else None
    mol = {7: 'O2', 1: 'H2O'}[r['mol']]
    return dict(wl_air=f'{la:.4f}', wl_vac=f'{lv:.4f}', unc_air=(f'{unc:.1e}' if unc else ''), unc_vac=(f'{unc:.1e}' if unc else ''),
                kind='HITRAN2020 (vac wavenumber); air derived (Peck & Reeder 1972)',
                transition=f"{mol} {r['gq_up']} <- {r['gq_lo']}; lower-state local quanta '{r['lq_lo']}'",
                lab_src=f"HITRAN2020 nu={r['nu_str']} cm-1 (ierr code {r['ierr'][0]}), S296={r['S']:.3e} cm/molec, delta_air={r['d_air']:+.6f} cm-1/atm",
                S=r['S'], d_air=r['d_air'])


# --- APGL98 (solar flux + disc centre, air) by recno
APG = {}
for ln in open(os.path.join(RAW, 'vizier', 'apgl98_table1.tsv')):
    p = ln.rstrip('\n').split('\t')
    if len(p) >= 8 and p[0].strip().isdigit():
        APG[int(p[0])] = dict(lD=p[1].strip(), eD=p[2].strip(), lF=p[3].strip(), eF=p[4].strip(), ion=p[5].strip())

# --- IAG line-core fits (derived, this work)
IAG = {}
fn = os.path.join(ROOT, 'data', 'curated', 'fraunhofer_classic_iag_linecores.csv')
if os.path.exists(fn):
    for r in csv.DictReader(open(fn)):
        IAG[r['label']] = r

R_FR = 'Fraunhofer1817'; R_ROW = 'Rowland(Smithsonian1920_T312)'; R_WP = 'Wikipedia_FraunhoferLines_rev1358424640'
R_MMH = 'MMH1966'; R_NIST = 'NIST_ASD_v5.12'; R_HIT = 'HITRAN2020'; R_APG = 'APGL1998'; R_IAG = 'IAG2016(fit:this work)'
R_SCH = 'SCHOTT2003_T9.1'

# Each row: letter, species, origin, lab (callable result), mmh (str,ew,pdfpage,ident), apgl recno, iag label,
#           rowland value (Rowland scale, verbatim), wikipedia (nm, verbatim), designation sources, note
ROWS = [
 ('A', 'O2 b1Sg+ - X3Sg- (0,0) A-band: R-branch head, line 1', 'telluric', hitran('13165.166455'),
  ('7593.695m', '193 (blend)', 350, 'Atm O2 R27+R29 (0,0)'), None, None, '7594.06* / 7621.28*', '759.370',
  f'{R_FR};{R_ROW};{R_WP}',
  'A is a BAND (O2 A-band, 7594-7684 A in MMH 1966), not a single line. Listed here: the band-head lines; MMH 7593.695 is the blend of the two HITRAN lines in rows A/A. '
  'Rowland gives two A lines: 7621.28 "beginning at the head of A, outside edge" and 7594.06 "single line beginning at the tail of A" (Rowland scale, ~+0.2 A vs ISA). '
  'Telluric: rest wavelength = lab; pressure shift delta_air ~ -0.005 to -0.010 cm-1/atm (= +3 to +6 mA at 1 atm).'),
 ('A', 'O2 b-X (0,0) A-band: R-branch head, line 2', 'telluric', hitran('13165.249023'),
  ('7593.695m', '193 (blend)', 350, 'Atm O2 R27+R29 (0,0)'), None, None, '', '', f'{R_MMH}',
  'Second component of the MMH 7593.695 band-head blend.'),
 ('a', 'H2O band group 0.72 um (telluric) -- Rowland reference line = solar Fe I', 'telluric band (Fraunhofer); solar Fe I (Rowland line)', nist('Fe I', '7164.4481'),
  ('7164.432', '153', 337, 'Fe I (EP 4.19 eV, mult. 1051)'), 4770, 'a (Rowland line)', '7164.725', '627.661 (labelled a, O2; = Rowland/Smithsonian alpha group, see row alpha)',
  f'{R_FR};{R_ROW};{R_MMH}',
  'Caveat: Fraunhofer\'s a is the telluric H2O band, and the solar values for this Fe I line are H2O-blended and disagree. MMH 7164.432 (disc centre) and the IAG core fit (-386 m/s) lie blueward of NIST Fe I 7164.4481, while the APGL98 flux value 7164.4578 is +406 m/s above it (+385 m/s above the Nave et al. 1994 value 7164.4486 used in APGL98b). Not a clean single-line reference. '
  'Fraunhofer (1817): "At a there are heaped together many lines which form a band" between A and B = telluric H2O bands (MMH 1966: H2O 400/301/202/221/103 bands, 6845-7480 A). '
  'Rowland/Smithsonian list a = 7164.725 (no element); MMH 1966 identifies the line at 7164.432 as solar Fe I (EW 153 mA). Strong telluric H2O lines lie within 0.1 A (HITRAN 7164.348, 7164.403 A air), so the solar position is blend-affected. '
  'Wikipedia/German WP "a = O2 627.661 nm" is the Greek-alpha O2 gamma-band group (Rowland alpha), a different feature.'),
 ('B', 'O2 b-X (1,0) B-band: first line of head', 'telluric', hitran('14557.976642'),
  ('6867.187m', '153 (blend with 6867.252m)', 324, 'Atm O2 R21 (1,0)'), None, None, '', '686.719',
  f'{R_FR};{R_MMH};{R_WP}', 'MMH 1966 note: "First line in head of B group: 6867.187 A" (Babcock & Herzberg 1948 O2 measurements). B is a band (6867-6967 A).'),
 ('B', "O2 b-X (1,0): Rowland's 'principal line in the head of B'", 'telluric', hitran('14552.079908'),
  ('6870.007m', '382 (blend with 6869.887m)', 324, 'Atm O2 R11 (+R13) (1,0)'), None, None, '6870.182 (dagger: principal line in head of B)', '',
  f'{R_ROW};{R_MMH}', 'Rowland value converted to ISA (~-0.21 A) matches the MMH 6869.887/6870.007 blend (HITRAN R13R13 6869.897 and R11Q12 6869.974 A air); the stronger HITRAN component is given.'),
 ('C', 'H I Balmer alpha (n=2-3), 7 fine-structure components', 'solar', nist('H I', '6562.795'),
  ('6562.808m', '4020', 317, 'H-alpha'), None, 'C (H-alpha)', '6563.045', '656.281',
  f'{R_FR};{R_ROW};{R_WP};{R_SCH}',
  'Lab value = ASD level-less blend row c69, 6562.795 +/- 0.009 A (Ritz): the gA-weighted centroid of the 7 fine-structure components (reproduced by this work as 6562.79523). '
  'ASD also lists an n=2-3 multiplet line, Ritz 6562.819 +/- 0.007 A (observed 6562.79 +/- 0.03 A, line ref. L7400c29); the airglow atomic table (Sect. 4) uses that entry. '
  'The two composite values differ by 24 mA (1.1 km/s), more than either stated uncertainty: the difference comes from how the components are combined, not from a measurement. '
  'Components span 6562.7097-6562.9094 A (air); strongest 2p2P3/2-3d2D5/2 = 6562.85175 A. At 1e-3 A a Balmer line has no single wavelength. Neither composite is the geocoronal emission centre (2s-3p pair 6562.72483 / 6562.77153 A). '
  'SCHOTT optics "C" = 656.2725 nm. Solar core is chromospheric and asymmetric.'),
 ('alpha (Greek)', 'O2 b-X (2,0) gamma-band: chief line of the alpha group', 'telluric', hitran('15924.045748'),
  ('6278.073m', '52 (blend with 6278.126m)', 308, 'Atm O2 R9 (+R11) (2,0)'), None, None, '6278.303 (double dagger: chief line in the alpha group)', '',
  f'{R_ROW};{R_MMH}', "Smithsonian/Rowland use Greek alpha for this O2 group; English & German Wikipedia list it as Latin 'a' (627.661 nm), conflicting with Fraunhofer's own a (between A and B)."),
 ('alpha (Greek) / Wikipedia a', 'O2 b-X (2,0): first line of the alpha group', 'telluric', hitran('15927.804100'),
  ('6276.590m', '24 (blend with 6276.633m)', 308, 'Atm O2 R15 (+R17) (2,0)'), None, None, '', '627.661 (as a)',
  f'{R_MMH};{R_WP}', 'MMH 1966 note: "First line of a group: 6276.590 A". Wikipedia 627.661 nm = mean of the MMH 6276.590/6276.633 blend.'),
 ('D1', 'Na I 3s 2S1/2 - 3p 2P1/2', 'solar', nist('Na I', '5895.92424'),
  ('5895.940m', '564', 293, 'Na I (D1)'), 4045, 'D1', '5896.155', '589.592', f'{R_FR};{R_ROW};{R_WP}',
  "Fraunhofer's D was the pair. Interstellar/telluric Na can add narrow components in stellar spectra."),
 ('D2', 'Na I 3s 2S1/2 - 3p 2P3/2', 'solar', nist('Na I', '5889.95095'),
  ('5889.973m', '752', 293, 'Na I (D2)'), 4039, 'D2', '5890.186', '588.995', f'{R_FR};{R_ROW};{R_WP}', ''),
 ('D', 'Na I D1+D2 doublet centre (arithmetic mean)', 'solar (doublet) / optics reference', dict(
      wl_air=f"{(float(nist('Na I', '5889.95095')['wl_air']) + float(nist('Na I', '5895.92424')['wl_air'])) / 2:.5f}",
      wl_vac=f"{(float(nist('Na I', '5889.95095')['wl_vac']) + float(nist('Na I', '5895.92424')['wl_vac'])) / 2:.5f}",
      unc_air='0.00003', unc_vac='0.000015', kind='derived: mean of NIST D1 and D2', transition='Na I 3s-3p doublet',
      lab_src='derived from NIST ASD v5.12 Na I D1, D2 (observed)'),
  ('', '', None, ''), None, None, '', '', f'{R_FR};{R_SCH}',
  "Derived row, not an observed line: the arithmetic mean of the NIST D1 and D2 values (air and vacuum); its uncertainty is the D1/D2 value 0.00003 A, taken as fully correlated (same source, L5918). "
  "Fraunhofer's original D (a pair). SCHOTT 2003 'D' = 589.2938 nm (centre of the double line); n_D is still quoted for glasses."),
 ('D3 (or d)', 'He I 2p 3P - 3d 3D (main blend of 4 components)', 'solar chromosphere (emission at limb/prominences; not a disc-integrated absorption line)', nist('He I', '5875.621'),
  ('not listed', '', 292, 'MMH has only Atm H2O 5875.596 here'), None, None, '5875.985', '587.5618 (D3 or d)', f'{R_ROW};{R_WP};{R_SCH}',
  'NIST observed value is for the unresolved J=2->1,2,3 / 1->1,2 blend (no uncertainty given); Ritz components 5875.5987-5875.6404 A plus the weak 3P0-3D1 component at 5875.9663 A. '
  'gA-weighted centroid of all six components (derived, this work) = 5875.6610 A air. Discovered in the 1868 eclipse chromosphere (Janssen, Lockyer). SCHOTT optics "d" = 587.5618 nm (He lamp).'),
 ('E1', 'Fe I a3F2 - z3D*1 (blend with Ca I 5270.270)', 'solar', nist('Fe I', '5270.3560'),
  ('5270.383', '255 (blend with Ca I 5270.269)', 258, 'Fe I (EP 1.61, mult. 37)'), 3111, 'E1 Fe', '5270.558 (Fe) / 5270.438 (Ca)', '527.039 (labelled E2)',
  f'{R_ROW};{R_WP}', "Fraunhofer: 'E itself consists of several lines'. Smithsonian/Rowland E1 = Fe 5270.558 + Ca 5270.438 (Rowland scale) = Fe I 5270.356 + Ca I 5270.270; Wikipedia calls 527.039 nm 'E2' (it is Rowland's E1 Fe line)."),
 ('E1', 'Ca I 3d4s 3D3 - 3d4p 3P*2', 'solar', nist('Ca I', '5270.270'),
  ('5270.269', '255 (blend with Fe I 5270.383)', 258, 'Ca I (EP 2.52, mult. 22)'), None, None, '5270.438 (Ca)', '', f'{R_ROW}',
  'IAG fit not usable (core merges with Fe I 5270.356).'),
 ('E2', 'Fe I a5F5 - z5D*4', 'solar', nist('Fe I', '5269.5370'),
  ('5269.550m', '478 (blend group 5269.418-5269.701)', 258, 'Fe I (EP 0.86, mult. 15)'), 3110, 'E2', '5269.723', '', f'{R_ROW}', ''),
 ('b1', 'Mg I 3s3p 3P*2 - 3s4s 3S1', 'solar', nist('Mg I', '5183.6042'),
  ('5183.619m', '1584', 253, 'Mg I (mult. 2)'), 2926, 'b1', '5183.791', '518.362', f'{R_FR};{R_ROW};{R_WP}',
  "Fraunhofer's b: '3 very strong lines'. NIST gives no uncertainty for this Mg I value (4 decimals). HS5000 digitization has 5183.600 (error)."),
 ('b2', 'Mg I 3s3p 3P*1 - 3s4s 3S1', 'solar', nist('Mg I', '5172.6843'),
  ('5172.698m', '1259', 252, 'Mg I (mult. 2)'), 2914, 'b2', '5172.856', '517.270', f'{R_FR};{R_ROW};{R_WP}', 'MMH value from OCR text layer, identical in HS5000 digitization (not checked on page image).'),
 ('b3', 'Fe I a5D3 - z7D*3', 'solar', nist('Fe I', '5168.8978'),
  ('5168.908m', '114', 252, 'Fe I (EP 0.05, mult. 1)'), 2907, 'b3 Fe I', '5169.069 (Fe)', '516.891', f'{R_ROW};{R_WP}',
  'Rowland b3 = two lines (5169.220, 5169.069 Rowland) = Fe II 5169.028 + Fe I 5168.898. HS5000 gives EW 700 for 5168.91 (MMH page shows 114).'),
 ('b3', 'Fe II a6S5/2 - z6P*7/2 (multiplet 42)', 'solar', nist('Fe II', '5169.0282'),
  ('5169.050m', '154', 252, 'Fe II (EP 2.89, mult. 42)'), None, 'b3 Fe II', '5169.220 (Fe)', '', f'{R_ROW}', ''),
 ('b4', 'Mg I 3s3p 3P*0 - 3s4s 3S1', 'solar', nist('Mg I', '5167.3216'),
  ('5167.327', '935 (blend with Fe I 5167.508m)', 252, 'Mg I (mult. 2)'), 2904, 'b4 Mg', '5167.497 (Mg)', '516.733', f'{R_ROW};{R_WP}', ''),
 ('b4', 'Fe I a3F4 - z3D*3', 'solar', nist('Fe I', '5167.4881'),
  ('5167.508m', '935 (blend with Mg I 5167.327)', 252, 'Fe I (EP 1.48, mult. 37)'), 2905, 'b4 Fe', '5167.678 (Fe)', '', f'{R_ROW}', ''),
 ('c', 'Fe I z7F*6 - e7D5 (multiplet 318)', 'solar', nist('Fe I', '4957.5965'),
  ('4957.613m', '696 (group 4957.307m-4957.697)', 238, 'Fe I (EP 2.81, mult. 318)'), 2537, 'c', '', '495.761', f'{R_WP}',
  "Not in Fraunhofer's own list nor in the Smithsonian/Rowland table; found in modern lists (Wikipedia). Companion Fe I 4957.2983 (MMH 4957.307m) is in the same blend."),
 ('F', 'H I Balmer beta (n=2-4)', 'solar', nist('H I', '4861.3250'),
  ('4861.342m', '3680', 231, 'H-beta'), None, 'F (H-beta)', '4861.527', '486.134', f'{R_FR};{R_ROW};{R_WP};{R_SCH}',
  'Lab value = ASD level-less blend row c64, 4861.3250 +/- 0.0014 A (Ritz): the gA-weighted centroid of the 7 components (reproduced as 4861.32497); components span 4861.2786-4861.3749 A. '
  'ASD\'s n=2-4 multiplet line, Ritz 4861.333 +/- 0.003 A (observed 4861.35 +/- 0.05 A, line ref. L7439c30), used in the airglow atomic table (Sect. 4), is 8 mA (0.49 km/s) longer. '
  'SCHOTT optics "F" = 486.1327 nm. HS5000 digitization 4861.35/EW 3000 disagrees with MMH page (4861.342m/3680).'),
 ('d', 'Fe I a3F4 - z5G*5 (Rowland/Smithsonian d)', 'solar', nist('Fe I', '4383.5447'),
  ('4383.557m', '1008', 192, 'Fe I (EP 1.48, mult. 41)'), 1168, 'd/e (4383)', '4383.721', '438.355 (labelled e)', f'{R_ROW};{R_WP}',
  "Designation conflict: Smithsonian/Rowland call Fe 4383.7 'd'; Wikipedia calls 438.355 nm 'e' and gives d = Fe 466.814 nm; D3 (He) is also called d."),
 ('d (Wikipedia)', 'Fe I z5D*2 - e5F3 (blend with Fe I 4668.061)', 'solar', nist('Fe I', '4668.1341'),
  ('4668.149', '112 (blend with 4668.073)', 217, 'Fe I (EP 3.26, mult. 554)'), 1948, 'd/e (4668)', '', '466.814', f'{R_WP}',
  'Only in modern lists (Wikipedia); a weak blended line (EW 112 mA) -- an odd choice for a lettered line. IAG fit affected by the blend.'),
 ('e (Wikipedia)', 'Fe I a3F4 - z5G*5', 'solar', nist('Fe I', '4383.5447'),
  ('4383.557m', '1008', 192, 'Fe I (EP 1.48, mult. 41)'), 1168, 'd/e (4383)', '4383.721 (as d)', '438.355', f'{R_WP}', 'Same line as Rowland d.'),
 ('e (optics)', 'Hg I 6s6p 3P*2 - 6s7s 3S1 (lamp line, not solar)', 'lab reference (optical glass)', nist('Hg I', '5460.750'),
  ('', '', None, ''), None, None, '', '546.073', f'{R_WP};{R_SCH}',
  'Glass-dispersion reference line borrowed the letter e. SCHOTT 2003 Table 9.1 gives 546.0740 nm; NIST ASD natural-Hg value 5460.750 A -- a 0.010 A difference (isotope composition / source).'),
 ("G'", 'H I Balmer gamma (n=2-5)', 'solar', nist('H I', '4340.4634'),
  ('4340.475m', '2855', 188, 'H-gamma'), None, "G' (H-gamma)", '4340.634', "434.047 (G')", f'{R_ROW};{R_WP}',
  "Lab value = ASD level-less blend row c61, 4340.4634 +/- 0.0011 A (Ritz; gA-weighted centroid of the 7 components, reproduced as 4340.46344). "
  "ASD's n=2-5 multiplet line is Ritz 4340.471 +/- 0.003 A, 7.6 mA longer. Some lists (e.g. SEOS) call H-gamma 'f'."),
 ('f', 'Fe I a3F2 - z3G*3', 'solar', nist('Fe I', '4325.7616'),
  ('4325.775m', '793', 186, 'Fe I (EP 1.61, mult. 42)'), 1039, 'f (4325)', '4325.939', '', f'{R_ROW}', ''),
 ('G', 'Fe I a3F3 - z3G*4 (in the CH G band)', 'solar', nist('Fe I', '4307.9020'),
  ('4307.912m', '723', 184, 'Fe I (EP 1.56, mult. 42) + Ti II + CH'), 996, 'G Fe', '4308.081 (Fe)', '430.790', f'{R_FR};{R_ROW};{R_WP}',
  "Fraunhofer: 'At G there are massed together many lines' = the CH A-X (0,0) 'G band' 4290-4315 A plus Fe I/Ca I lines."),
 ('G', 'Ca I 4s4p 3P*1 - 4p2 3P0', 'solar', nist('Ca I', '4307.744'),
  ('4307.748', '59', 184, 'Ca I (EP 1.89, mult. 5)'), 995, 'G Ca', '4307.907 (Ca)', '430.774', f'{R_ROW};{R_WP}', 'NIST observed value is 4307.74 (2 decimals); Ritz 4307.744 used.'),
 ('g', 'Ca I 4s2 1S0 - 4s4p 1P*1 (resonance line)', 'solar', nist('Ca I', '4226.7276'),
  ('4226.740m', '1476', 176, 'Ca I (EP 0.00, mult. 2)'), 737, 'g', '4226.904', '', f'{R_ROW}', 'SCHOTT optics "g" is instead Hg 435.8 nm.'),
 ('h', 'H I Balmer delta (n=2-6)', 'solar', nist('H I', '4101.7346'),
  ('4101.748m', '3133', 163, 'H-delta (+In I)'), None, 'h (H-delta)', '4102.000', '410.175', f'{R_ROW};{R_WP}',
  'Lab value = ASD level-less blend row c58, 4101.7346 +/- 0.0010 A (Ritz; gA-weighted centroid of the 7 components, reproduced as 4101.73463); ASD\'s n=2-6 multiplet line is Ritz 4101.7415 +/- 0.0024 A, 6.9 mA longer. '
  'IAG core fit pulled by Fe I 4101.68 blend. SCHOTT optics "h" = Hg 404.6561 nm.'),
 ('H', 'Ca II 4s 2S1/2 - 4p 2P*1/2', 'solar (chromospheric core)', nist('Ca II', '3968.469'),
  ('3968.492m', '15467', 149, 'Ca II (H)'), None, None, '3968.625', '396.847', f'{R_FR};{R_ROW};{R_WP}',
  "NIST ASD gives only 3 decimals (no uncertainty) for Ca II H & K; Fraunhofer's H covered both H and K ('the two bands at H'). Core shows chromospheric H2V/H3 structure; H-epsilon (3970.07) lies in its wing."),
 ('K', 'Ca II 4s 2S1/2 - 4p 2P*3/2', 'solar (chromospheric core)', nist('Ca II', '3933.663'),
  ('3933.682m', '20253', 146, 'Ca II (K)'), None, None, '3933.825', '393.366', f'{R_ROW};{R_WP}', 'Strongest line in the visible solar spectrum (EW ~20 A).'),
 ('L', 'Fe I a5F5 - y5D*4', 'solar', nist('Fe I', '3820.4249'),
  ('3820.436m', '1712', 132, 'Fe I (EP 0.86, mult. 20)'), None, None, '3820.586', '382.044', f'{R_ROW};{R_WP}', ''),
 ('M', 'Fe I a5F3 - y5F*2', 'solar', nist('Fe I', '3727.6188'),
  ('3727.634m', '632', 121, 'Fe I (EP 0.96, mult. 21)'), None, None, '3727.778', '', f'{R_ROW}', ''),
 ('N', 'Fe I a5F5 - z5G*6', 'solar', nist('Fe I', '3581.1928'),
  ('3581.209r', '2144', 104, 'Fe I (EP 0.86, mult. 23)'), None, None, '3581.349', '358.121', f'{R_ROW};{R_WP}', ''),
 ('O', 'Fe I a5D3 - z5P*2', 'solar', nist('Fe I', '3440.9885'),
  ('3441.019r', '634', 89, 'Fe I (EP 0.05, mult. 6)'), None, None, '3441.155', '', f'{R_ROW}',
  'Rowland O converts to the weaker 3441.0 line; the stronger multiplet partner Fe I 3440.6057 (MMH 3440.626r, EW 1243) is adjacent.'),
 ('P', 'Ti II a4F7/2 - z4G*9/2 (blend with Ti I/Sc II 3361.29)', 'solar', nist('Ti II', '3361.2121'),
  ('3361.193r', '939 (blend with 3361.287r)', 79, 'Ti II (EP 0.03, mult. 1)'), None, None, '3361.327 (listed as Fe)', '336.112', f'{R_ROW};{R_WP}',
  "Rowland/Smithsonian attribute P to Fe; MMH 1966 and modern lists identify Ti II. Wikipedia's 336.112 nm does not match Ti II 336.121 nm or MMH 336.119 nm (apparent digit transposition)."),
 ('Q', 'Fe I a5P3 - v5P*3', 'solar', nist('Fe I', '3286.7528'),
  ('3286.772r', '660', 70, 'Fe I (EP 2.18, mult. 91)'), None, None, '3286.898', '', f'{R_ROW}', ''),
 ('R', 'Ca II 4p 2P*3/2 - 4d 2D5/2', 'solar', nist('Ca II', '3179.331'),
  ('3179.342r', '580', 57, 'Ca II (EP 3.15, mult. 4)'), None, None, '3179.453 (Ca)', '', f'{R_ROW}', 'NIST gives 3 decimals only.'),
 ('R', 'Ca II 4p 2P*3/2 - 4d 2D3/2', 'solar', nist('Ca II', '3181.275'),
  ('3181.276r', '182', 58, 'Ca II (EP 3.15, mult. 4)'), None, None, '3181.387 (Ca)', '', f'{R_ROW}', ''),
 ('S1/S2', 'Fe I a5F3 - x5D*3', 'solar', nist('Fe I', '3100.6647'),
  ('3100.682r', '216', 48, 'Fe I (EP 0.96, mult. 28) + Ti I'), None, None, '3100.787', '', f'{R_ROW}',
  'Smithsonian brackets S1, S2 over three Fe lines 3100.787/3100.430/3100.046 (Rowland scale).'),
 ('S1/S2', 'Fe I a5F2 - x5D*2', 'solar', nist('Fe I', '3100.3033'),
  ('3100.325r', '263', 48, 'Fe I (EP 0.99, mult. 28)'), None, None, '3100.430', '', f'{R_ROW}', ''),
 ('S1/S2', 'Fe I a5F4 - x5D*4 (blend with Fe I 3099.8943)', 'solar', nist('Fe I', '3099.9679'),
  ('3099.987r', '465 (blend with 3099.896r)', 48, 'Fe I (EP 0.91, mult. 28)'), None, None, '3100.046', '', f'{R_ROW}', ''),
 ('s', 'Fe I a5D2 - y5D*3', 'solar', nist('Fe I', '3047.6043'),
  ('3047.614', '[35] (Rowland intensity, not EW)', 43, 'Fe I (EP 0.09, mult. 9)'), None, None, '3047.725', '', f'{R_ROW}', ''),
 ('T', 'Fe I a5D4 - y5D*4 (Rowland/Smithsonian T)', 'solar', nist('Fe I', '3020.6388'),
  ('3020.656', '[40] (Rowland intensity)', 40, 'Fe I (EP 0.00, mult. 9) (+Cr I)'), None, None, '3020.76', '', f'{R_ROW}',
  'Conflict: Wikipedia T = Fe 302.108 nm, i.e. the multiplet partner Fe I 3021.0725 (MMH 3021.077, [30]).'),
 ('T (Wikipedia)', 'Fe I a5D3 - y5D*3', 'solar', nist('Fe I', '3021.0725'),
  ('3021.077', '[30] (Rowland intensity)', 40, 'Fe I (EP 0.05, mult. 9)'), None, None, '', '302.108', f'{R_WP}', ''),
 ('t', 'Fe I a5D3 - y5D*2 (blend with Fe I 2994.502 and Ni I 2994.453)', 'solar', nist('Fe I', '2994.4266'),
  ('2994.436', '[40] (Rowland intensity)', 38, 'Fe I (mult. 9, 11) + Ni I (mult. 27)'), None, None, '2994.53 (Fe)', '299.444 (Ni)', f'{R_ROW};{R_WP}',
  'Rowland/Smithsonian: Fe; Wikipedia: Ni. MMH 1966 lists a blend dominated by Fe I with Ni I 2994.46 (NIST obs, 2 decimals) contributing.'),
 ('U', 'Fe I a5D3 - y5F*3', 'solar', nist('Fe I', '2947.8758'),
  ('2947.70', '[30] (Rowland intensity)', 35, 'Fe II- , Fe I (mult. 10), (V II)'), None, None, '2947.99', '', f'{R_ROW}',
  'MMH gives only 2 decimals here and its value is 0.18 A shortward of Fe I 2947.876; MMH lists the feature as a blend (Fe II, Fe I, V II).'),
 ('Z', 'H2O (2,1,1)<-(0,0,0) line (NOT O2)', 'telluric', hitran('12151.823462'),
  ('8226.962', '347', 366, 'Atm H2O Q2 (211 band)'), None, None, '', '822.696 (labelled O2)', f'{R_WP};{R_MMH}',
  'Wikipedia lists Z as O2; MMH 1966 and HITRAN2020 identify the 8226.96 A line as H2O (no O2 transition in HITRAN within 500 cm-1).'),
 ('y', 'H2O (0,0,3)<-(0,0,0) line (NOT O2)', 'telluric', hitran('11123.379544'),
  ('beyond MMH range', '', None, ''), None, None, '', '898.765 (labelled O2)', f'{R_WP}',
  'Wikipedia lists y = O2 898.765 nm; nearest strong HITRAN line is H2O at 8987.606 A (0.044 A away; S=5.8e-23). No O2 lines in HITRAN near 11123 cm-1. Identification tentative.'),
 # --- optical-glass reference lines (lamp lines that re-use Fraunhofer-style letters; not solar)
 ('i (optics)', 'Hg I 6s6p 3P*2 - 6s6d 3D3', 'lab reference (optical glass)', nist('Hg I', '3650.158'), ('', '', None, ''), None, None, '', '', R_SCH, 'SCHOTT 2003: 365.0146 nm.'),
 ('h (optics)', 'Hg I 6s6p 3P*0 - 6s7s 3S1', 'lab reference (optical glass)', nist('Hg I', '4046.565'), ('', '', None, ''), None, None, '', '', R_SCH, 'SCHOTT 2003: 404.6561 nm. (Solar h = H-delta.)'),
 ('g (optics)', 'Hg I 6s6p 3P*1 - 6s7s 3S1', 'lab reference (optical glass)', nist('Hg I', '4358.335'), ('', '', None, ''), None, None, '', '', R_SCH, 'SCHOTT 2003: 435.8343 nm. (Solar g = Ca I 4226.7.)'),
 ("F' (optics)", 'Cd I 5s5p 3P*1 - 5s6s 3S1', 'lab reference (optical glass)', nist('Cd I', '4799.9123'), ('', '', None, ''), None, None, '', '', R_SCH, 'SCHOTT 2003: 479.9914 nm.'),
 ('d (optics)', 'He I 2p 3P - 3d 3D blend (= D3)', 'lab reference (optical glass)', nist('He I', '5875.621'), ('', '', None, ''), None, None, '', '587.5618', R_SCH, 'SCHOTT 2003: 587.5618 nm; defines n_d and the Abbe number V_d.'),
 ("C' (optics)", 'Cd I 5s5p 1P*1 - 5s5d 1D2', 'lab reference (optical glass)', nist('Cd I', '6438.4695'), ('', '', None, ''), None, None, '', '', R_SCH, 'SCHOTT 2003: 643.8469 nm (also the former primary length standard, red Cd line).'),
 ('r (optics)', 'He I 2p 3P*2,1 - 3s 3S1', 'lab reference (optical glass)', nist('He I', '7065.190'), ('', '', None, ''), None, None, '', '', R_SCH, 'SCHOTT 2003: 706.5188 nm.'),
 ('s (optics)', 'Cs I 6s 2S1/2 - 6p 2P*3/2 (D2)', 'lab reference (optical glass)', nist('Cs I', '8521.13165'), ('', '', None, ''), None, None, '', '', R_SCH, 'SCHOTT 2003: 852.11 nm. (Solar s = Fe I 3047.6.)'),
 ('t (optics)', 'Hg I 6s6p 1P*1 - 6s7s 1S0', 'lab reference (optical glass)', nist('Hg I', '10139.787'), ('', '', None, ''), None, None, '', '', R_SCH,
  'SCHOTT 2003: 1013.98 nm. NIST Ritz value used (10139.787 +/- 0.007 A air, 10142.567 A vac); the NIST observed value 10139.75 +/- 0.05 A is 7x less precise. (Solar t = Fe I 2994.4.)'),
]

fields = ['letter', 'species', 'wl_air_A', 'wl_vac_A', 'uncert_A', 'origin', 'note', 'refs',
          'transition', 'lab_value_type', 'lab_source_detail', 'uncert_vac_A',
          'solar_mmh1966_air_A', 'mmh1966_ew_mA', 'mmh1966_pdf_page', 'mmh1966_ident',
          'solar_apgl98_flux_air_A', 'solar_apgl98_flux_err_A', 'solar_apgl98_disc_air_A',
          'solar_iag_flux_air_A_derived', 'solar_iag_flux_vac_A_derived', 'iag_minus_lab_mps_derived', 'iag_core_depth',
          'rowland_value_rowland_scale_A', 'wikipedia_value_nm', 'designation_sources']
out = []
for (letter, species, origin, lab, mmh, apg, iag, row_v, wp, dsrc, note) in ROWS:
    refs = [R_NIST if 'NIST' in lab['lab_src'] else R_HIT]
    if mmh[0] and mmh[2]:
        refs.append(R_MMH)
    d = dict(letter=letter, species=species, wl_air_A=lab['wl_air'], wl_vac_A=lab['wl_vac'],
             uncert_A=lab['unc_air'] if lab['unc_air'] else 'not given (precision = last decimal)',
             origin=origin, note=note, transition=lab['transition'], lab_value_type=lab['kind'],
             lab_source_detail=lab['lab_src'], uncert_vac_A=lab['unc_vac'],
             solar_mmh1966_air_A=mmh[0], mmh1966_ew_mA=mmh[1], mmh1966_pdf_page=(mmh[2] or ''), mmh1966_ident=mmh[3],
             rowland_value_rowland_scale_A=row_v, wikipedia_value_nm=wp, designation_sources=dsrc)
    if apg:
        a = APG[apg]
        d.update(solar_apgl98_flux_air_A=a['lF'], solar_apgl98_flux_err_A=a['eF'], solar_apgl98_disc_air_A=a['lD'])
        refs.append(R_APG)
    if iag and iag in IAG and iag != 'E1 Ca':
        r = IAG[iag]
        d.update(solar_iag_flux_air_A_derived=r['iag_wl_air_A'], solar_iag_flux_vac_A_derived=r['iag_wl_vac_A'],
                 iag_minus_lab_mps_derived=r['iag_minus_lab_mps'], iag_core_depth=r['core_depth'])
        refs.append(R_IAG)
    if row_v:
        refs.append(R_ROW)
    if wp:
        refs.append(R_WP)
    if 'optics' in letter or 'SCHOTT' in note:
        refs.append(R_SCH)
    if species.startswith('H I Balmer'):   # say which ASD composite entry is used
        ref = lab['lab_src'].split('line ref ')[1].split(';')[0]
        d['lab_value_type'] = f'Ritz: ASD level-less blend centroid (row {ref})'
    d['refs'] = ';'.join(dict.fromkeys(refs))
    out.append(d)

with open(OUT, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fields)
    w.writeheader(); w.writerows(out)
print('wrote', OUT, len(out), 'rows')
