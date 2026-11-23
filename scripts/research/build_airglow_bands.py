#!/usr/bin/env python
"""Build data/curated/airglow_bands.csv : molecular airglow / auroral band systems.

usage: python -I build_airglow_bands.py <oh_origins_csv> <scripts_dir> <out_csv>

Sources of band positions (all downloaded; see airglow_meta.REFS):
 * LK77  Lofthus & Krupenie 1977 (N2, N2+) - values transcribed by eye from the table page images
         of the NIST JPCRD reprint (Tables 7a, 27, 28, 29, 41, 42).  LK77 convention: air above
         2000 A, vacuum below 2000 A (Table 7 footnote); verified below against the tabulated
         wavenumbers where both are given.
 * Krupenie72 (O2) Tables 3, 4, 5, 6, 9, 11 - transcribed the same way.
 * Meier91 Tables II/III (NO gamma/delta, LBH intensities) - nominal 1-A band positions.
 * OH Meinel band origins: derived (oh_band_origins.py) from HITRAN2020 16OH level energies.
 * PALACE (Noll et al. 2025): intensities/heights for OH, O2, FeO, HO2, O2 UV continuum.
Wavelength in cm-1 -> vacuum A = 1e8/sigma; vacuum -> air with Peck & Reeder (1972) as in NIST ASD.
"""
import csv, sys
oh_csv, scripts_dir, out_csv = sys.argv[1:4]

def n_air(lv):
    s2 = (1e4 / lv) ** 2
    return 1 + 1e-8 * (8060.51 + 2480990 / (132.274 - s2) + 17455.7 / (39.32957 - s2))

ROWS = []
def add(mol, system, trans, band, wl, kind, sigma='', medium='', prec='', src='', cat='', mech='', alt='', R='', note='', refs='', sigma_kind='origin'):
    ROWS.append(dict(molecule=mol, system=system, transition=trans, band=band, wl_A=wl, wl_kind=kind,
                     sigma_cm=sigma, sigma_kind=(sigma_kind if sigma else ''), medium=medium, precision_A=prec, wl_source=src, category=cat,
                     mechanism=mech, altitude_km=alt, typical_R=R, note=note, refs=refs))

# ---------------- N2+ First Negative B2Su+ - X2Sg+  (LK77 Table 42, p.233) ----------------
S = 'LK77 Table 42 (JPCRD 6, p.233): head lambda_H of violet-degraded bands + band origin sigma_0 (heads and origins are different features)'
CAT = 'aurora|dayglow|twilight'
MECH = 'aurora: electron-impact ionization-excitation of N2; day: photoionization-excitation + resonant scattering by N2+; twilight: resonance scattering of sunlight by N2+'
for band, wl, sg, I, R in [
    ('0-0', '3914.4', '25566.0', 1500, 'day 2-5 kR zenith [Noxon68]; day 6700 R nadir 600 km [Meier91 T.II]; aurora model 0.189 kR per erg cm-2 s-1 @1 keV [Barthelemy25]'),
    ('0-1', '4278.1', '23391.2', 400, 'aurora model 0.072 kR per erg cm-2 s-1 @1 keV [Barthelemy25]; aurora peak altitude ~116.6 km (mean) [Whiter23]'),
    ('0-2', '4709.2', '21249.0', 50, ''), ('0-3', '5228.3', '19139.7', 2, ''),
    ('1-0', '3582.1', '', 300, ''), ('1-1', '3884.3', '', 250, ''), ('1-2', '4236.5', '23620.5', 250, ''),
    ('2-2', '3857.9', '25939.8', 40, ''), ('1-3', '4651.8', '21511.3', 50, ''), ('2-3', '4199.1', '23830.0', 100, '')]:
    add('N2+', 'First Negative (1N)', 'B2Sigma_u+ - X2Sigma_g+', band, wl, 'head (violet-degraded)', sg, 'air', '~0.5 (LK77: measurement uncertainty about 0.5 A)', S, CAT, MECH,
        'aurora ~100-150 (427.8 nm peak 116.55 mean, Whiter23); day 100-400 (Noxon68)', R,
        f'Herzberg lab intensity {I} (LK77 col. Ic).', 'LK77;Noxon68;Meier91;Barthelemy25;Whiter23;Hunten67')
# ---------------- N2+ Meinel A2Pi_u - X2Sg+ (LK77 Table 41, p.232) ----------------
S = 'LK77 Table 41 (JPCRD 6, p.232), data of Douglas: Q1 and R2 heads + origin sigma_0'
for band, q1, r2, sg, Ia, Ib, meas in [
    ('1-0', '9212.06', '9145.3', '10889.49', '', '100', 'long-wavelength head measured'),
    ('2-0', '7874.62', '7825.73', '12732.75', '105', '68', 'measured'),
    ('3-1', '8105.35', '8053.64', '12371.51', '53', '38', 'measured'),
    ('4-2', '8348.26', '8293.42', '12012.73', '13', '15', 'measured'),
    ('3-0', '6890.5', '6853.0', '', '7', '11', 'long-wavelength head measured'),
    ('0-0', '11133.9', '11036.2', '', '35', '', 'extrapolated by Douglas'),
    ('1-2', '15297', '15114', '', '80', '', 'extrapolated by Douglas'),
    ('0-1', '14692', '14523', '', '25', '', 'extrapolated by Douglas')]:
    add('N2+', 'Meinel', 'A2Pi_u - X2Sigma_g+', band, q1, 'Q1 head', sg, 'air', '0.01-0.1 (as printed); ' + meas, S, 'aurora|dayglow',
        'electron-impact ionization-excitation of N2 (aurora); photoionization-excitation (day)', 'aurora E region', 'aurora model (all Meinel 700-1000 nm) 1.022 kR per erg cm-2 s-1 @1 keV [Barthelemy25]',
        f'R2 head {r2} A. Auroral relative intensities Ia={Ia} (Meinel), Ib={Ib} (Fedorova) from LK77.', 'LK77;Barthelemy25;Noxon68')
# ---------------- N2 First Positive B3Pi_g - A3Su+ (LK77 Table 28, p.225) ----------------
S = 'LK77 Table 28 (JPCRD 6, p.225): P1-branch heads (red-degraded) + origin sigma_0 (Benesch et al. via LK77)'
for band, wl, sg, Ia in [('0-0', '10510.04', '', '1000'), ('1-0', '8912.39', '11257.15', '1260'), ('2-1', '8723.03', '11500.40', '862'),
                         ('3-2', '8542.54', '11742.74', '305'), ('2-0', '7753.67', '12933.31', '719'), ('3-1', '7626.76', '', '833'),
                         ('4-2', '7504.70', '13361.13', '872'), ('5-3', '7387.19', '13572.98', '590'), ('6-4', '7273.98', '13783.55', '362'),
                         ('3-0', '6875.24', '14580.89', '142'), ('4-1', '6788.60', '14766.37', '335'), ('5-2', '6704.79', '14950.40', '504'),
                         ('6-3', '6623.57', '15133.21', '548'), ('7-4', '6544.88', '15314.58', '506'), ('8-5', '6468.62', '15494.59', '408'),
                         ('9-6', '6394.66', '15673.30', '318'), ('10-7', '6322.85', '15850.77', '211')]:
    add('N2', 'First Positive (1P)', 'B3Pi_g - A3Sigma_u+', band, wl, 'P1 head', sg, 'air', '0.01 (as printed)', S, 'aurora|dayglow',
        'electron impact on N2 (direct + cascade); populates A3Sigma_u+', 'aurora E region',
        'aurora model (whole 1P 500-1000 nm) 0.303 kR per erg cm-2 s-1 @1 keV [Barthelemy25]',
        f'Lab relative intensity Ia={Ia} (Turner & Nicholls via LK77). ' + ('7-4 band near 654.5 nm contaminates H-alpha photometry (Barthelemy25). ' if band == '7-4' else '')
        + ('P1 head 0.45 A from Ca II 8542.09 (meteor/solar; derived). ' if band == '3-2' else ''), 'LK77;Barthelemy25')
# ---------------- N2 Second Positive C3Pi_u - B3Pi_g (LK77 Table 29, p.226) ----------------
S = 'LK77 Table 29 (JPCRD 6, p.226): principal (P3) heads (violet-degraded) + origin sigma_0'
for band, wl, sg, Ib in [('0-0', '3371.3', '29671.05', '100'), ('0-1', '3576.9', '', '67'), ('0-2', '3804.9', '', '22'), ('1-0', '3159.3', '', '67'),
                         ('1-2', '3536.7', '', '30'), ('1-3', '3755.4', '26636.59', '25'), ('2-1', '3136.0', '', '29'), ('2-0', '2976.8', '', '12'),
                         ('3-1', '2962.0', '', '23'), ('1-4', '3998.4', '25018.20', '12'), ('0-3', '4059.4', '24642.19', '7.0'), ('2-4', '3710.5', '26958.99', '8.6')]:
    add('N2', 'Second Positive (2P)', 'C3Pi_u - B3Pi_g', band, wl, 'P3 head (violet-degraded)', sg, 'air', '0.1 (as printed)', S, 'aurora|dayglow',
        'photoelectron / auroral electron impact on N2 (spin-forbidden excitation of C3Pi_u)', 'day obs >170 (Noxon68)',
        ('day 0.4 kR zenith (strongest band) [Noxon68]; aurora model (whole 2P 300-500 nm) 0.273 kR per erg cm-2 s-1 @1 keV [Barthelemy25]' if band == '0-0' else ''),
        f'Relative intensity Ib={Ib} (Tyte via LK77).', 'LK77;Noxon68;Barthelemy25')
# ---------------- N2 Vegard-Kaplan A3Su+ - X1Sg+ (LK77 Table 27, p.224) ----------------
S = 'LK77 Table 27 (JPCRD 6, p.224): R heads; band origins (b) from Wilkinson / Miller (unc. 0.003 cm-1 for Miller)'
for band, wl, sg, Id in [('0-3', '2332.8', '', '0.224'), ('1-4', '2377.5', '', '0.168'), ('0-4', '2461.6', '40607.246', '0.831'), ('1-5', '2509.8', '', '0.59'),
                         ('0-5', '2603.8', '38392.167', '1.40'), ('0-6', '2760.6', '36206.048', '2.33'), ('0-7', '2935.7', '', '2.39'), ('1-8', '2997.0', '', '0.453'),
                         ('0-8', '3133', '', '1.92'), ('1-9', '3197.5', '', '0.816'), ('0-9', '3351.5', '', '1.26'), ('1-10', '3424.6', '', '0.875'),
                         ('0-10', '3603.0', '', '0.641'), ('1-11', '3682.4', '', '0.664')]:
    add('N2', 'Vegard-Kaplan (VK)', 'A3Sigma_u+ - X1Sigma_g+', band, wl, 'head (red-degraded)', sg, 'air', '0.1 (as printed)', S, 'aurora|dayglow',
        'spin-forbidden emission of metastable N2(A) produced by electron impact + cascade from B, C, W; quenched by O below ~150 km (Meier91)', '>~150 (quenching below)',
        ('day 6820 R total A-state excitation, quenching not included [Meier91 T.II]; aurora model (VK 200-600 nm) 0.365 kR per erg cm-2 s-1 @1 keV [Barthelemy25]' if band == '0-5' else ''),
        f'Relative intensity Id={Id} (Chandraiah & Shepherd via LK77).', 'LK77;Meier91;Barthelemy25')
# ---------------- N2 LBH a1Pi_g - X1Sg+ (LK77 Table 7a cont., p.204) ----------------
S = 'LK77 Table 7(a) (JPCRD 6, p.204): R heads (vacuum below 2000 A) + origin sigma_0'
MR = {'3-0': '222', '4-0': '200', '2-0': '267 (with 5-2)', '1-0': '207 (with 4-2)', '1-1': '144', '1-2': '127 (with 4-4)', '0-2': '87 (with 4-5)',
      '2-4': '73', '0-4': '98 (with 4-7)', '1-5': '116 (with 5-8)', '1-6': '113 (with 5-9)', '2-7': '127 (with 6-10)', '2-8': '100 (with 6-11)',
      '5-0': '(111)', '6-0': '<70'}
for band, wl, sg, I in [('6-0', '1273.2', '', '2.27'), ('5-0', '1298.5', '', '3.14'), ('4-0', '1325.3', '', ''), ('3-0', '1353.6', '', '17.2'),
                        ('2-0', '1383.8', '72256.08', '40.4'), ('1-0', '1415.9', '70617.58', '36.8'), ('0-0', '1450.1', '68951.24', '23.7'),
                        ('1-1', '1464.2', '', '62.6'), ('0-1', '1500.8', '', '72.2'), ('1-2', '1515.3', '', '23.7'), ('2-3', '1530.0', '', '19.4'),
                        ('0-2', '1554.5', '64319.98', '100'), ('2-4', '1584.4', '', '13.7'), ('0-3', '1611.4', '62047.63', '92.6'), ('1-4', '1626.6', '61470.02', '20.7'),
                        ('0-4', '1671.9', '59803.71', '52.9'), ('1-5', '1687.4', '59254.97', '34.2'), ('2-6', '1703.1', '58707.35', '6.93'),
                        ('0-5', '1736.1', '57588.64', '23.7'), ('1-6', '1752.0', '57068.85', '28.7'), ('2-7', '1768.1', '', '15.3'),
                        ('1-7', '1820.8', '', '13.1'), ('2-8', '1837.2', '', '12.7')]:
    add('N2', 'Lyman-Birge-Hopfield (LBH)', 'a1Pi_g - X1Sigma_g+', band, wl, 'R head', sg, 'vacuum', '0.1 (McEwen calc. vs obs. within 0.1 A)', S, 'dayglow|aurora',
        'photoelectron / auroral electron impact on N2', 'thermosphere (~150 dayglow peak region not quantified here)',
        (f'day {MR[band]} R nadir 600 km [Meier91 T.II]' if band in MR else '') + ('; day LBH total 4560 R [Meier91 T.II]' if band == '2-0' else ''),
        f'Relative intensity {I} (McEwen via LK77).', 'LK77;Meier91')
# ---------------- O2 Atmospheric b1Sg+ - X3Sg- (Krupenie72 Tables 4, 5) ----------------
S4 = 'Krupenie72 Table 4 (JPCRD 1, p.459): band origin lambda_0 / sigma_0 (Babcock & Herzberg; refit Albritton et al.)'
S5 = 'Krupenie72 Table 5 (JPCRD 1, p.460): R heads'
MECH = 'night: O+O+M -> O2*, collisional cascade to b (Barth-type) [PALACE]; day/twilight: resonance absorption of sunlight + O(1D) energy transfer'
for band, wl, sg, kind, src, R in [
    ('0-0', '7619.33', '13120.9085', 'origin', S4, 'night: strong but self-absorbed, not observable from ground [PALACE]; day ~300 kR (7600+8640) [Noxon68]'),
    ('0-0', '7593.73', '13165.14', 'head (red-degraded)', S5, ''),
    ('1-0', '6882.47', '14525.6609', 'origin', S4, ''), ('1-0', '6867.2', '14558', 'head (red-degraded)', S5, ''),
    ('1-1', '7708.41', '12969.2821', 'origin', S4, ''), ('1-1', '7683.85', '13010.73', 'head (red-degraded)', S5, ''),
    ('2-0', '6286.61', '15902.4174', 'origin', S4, ''), ('2-1', '6968.63', '14346.076', 'origin', S4, ''),
    ('3-0', '5795.13', '17251.0922', 'origin', S4, ''),
    ('0-1', '(8623)', '11594', 'head (red-degraded, uncertain)', S5, 'night 340 R (Paranal) [PALACE]'),
    ('0-1', '~8645', '', 'origin (Kaplan; Meinel)', S5 + ' text', 'night 340 R (Paranal) [PALACE]')]:
    add('O2', 'Atmospheric (A band; b-X)', 'b1Sigma_g+ - X3Sigma_g-', band, wl, kind, sg, 'air', '0.01 (origins, Babcock & Herzberg); heads 0.01-several A' if 'origin' == kind else 'as printed', src,
        'nightglow|dayglow|twilight', MECH, '~94 peak (Yee et al. 1997 via PALACE); day 40-130 (Noxon68)', R,
        'Total b-X model (all bands incl. self-absorbed 0-0) 6.9 kR at night [PALACE].' if (band, kind) == ('0-0', 'origin') else '', 'Krupenie72;PALACE;Noxon68',
        sigma_kind=('head' if 'head' in kind else 'origin'))
# ---------------- O2 IR Atmospheric a1Dg - X3Sg- (Krupenie72 Table 3) ----------------
S3 = 'Krupenie72 Table 3 (JPCRD 1, p.459): band origins lambda_0 / sigma_0 (Herzberg & Herzberg); 0-1 is a twilight head'
for band, wl, sg, kind, R in [('0-0', '12683.0', '7882.39', 'origin', 'night 155 kR without self-absorption (~17 kR at ground) [PALACE]; day ~25 MR [Noxon68]'),
                              ('1-0', '10674.1', '9365.89', 'origin', 'night (1-0) 0.012 R model [PALACE]'),
                              ('0-1', '(15800)', '6327', 'head (twilight, Vallance Jones & Harrison)', 'night (0-1) 2.0 kR [PALACE]')]:
    # the (0-1) head is printed in parentheses (approximate) and lies 0.96 A from the air wavelength of its own sigma (medium_check)
    pr = '>~1 (approximate head, printed in parentheses; 0.96 A from its own sigma)' if band == '0-1' else '0.1 (origins as printed)'
    add('O2', 'Infrared Atmospheric (a-X)', 'a1Delta_g - X3Sigma_g-', band, wl, kind, sg, 'air', pr, S3, 'nightglow|dayglow|twilight',
        'night: O recombination + other pathways; day/early night: O3 photolysis -> O2(a), lifetime ~1 h [PALACE]', '~90 night (PALACE); day 40-90 (Noxon68)', R, '', 'Krupenie72;PALACE;Noxon68;Hunten67',
        sigma_kind=('head' if 'head' in kind else 'origin'))
# ---------------- O2 Herzberg I A3Su+ - X3Sg- (Krupenie72 Table 11, emission heads) ----------------
S11 = 'Krupenie72 Table 11 (JPCRD 1, p.462): R heads observed in emission (Broida & Gaydon; Barth & Kaplan); uncertain by 2-5 A'
for band, wl, sg, I in [('4-4', '3142', '31818', 7), ('3-4', '3211', '31134', 10), ('2-4', '3285', '30433', 7), ('3-5', '3370', '29665', 10),
                        ('2-5', '3453', '28952', 8), ('3-6/1-5', '3542', '28225', 8), ('2-6', '3633', '27518', 8), ('1-6', '3734', '26773', 8),
                        ('2-7', '3829', '26109', 8), ('1-7', '3938', '25386', 7), ('0-8', '4309', '23201', 7)]:
    add('O2', 'Herzberg I', 'A3Sigma_u+ - X3Sigma_g-', band, wl, 'head (red-degraded)', sg, 'air', '2-5 (stated)', S11, 'nightglow',
        'three-body recombination O+O+M -> O2(A)+M; quenched at low altitude (Meier91)', '~94 (PALACE, O2 non a-X systems)',
        ('night total Herzberg I 225 R nadir 600 km [Meier91 T.III]; 600 R zenith (Hennes 1966, via Meier91)' if band == '3-5' else ''),
        f'Lab afterglow intensity {I} (Broida & Gaydon via Krupenie72). Herzberg I : II = 4 : 1 and Herzberg I : Chamberlain = 2 : 1 (Murtagh et al. 1986 via Meier91).', 'Krupenie72;Meier91;PALACE',
        sigma_kind='head')
# ---------------- O2 Herzberg II c1Su- - X3Sg- (Krupenie72 Table 6a) ----------------
for band, wl, sg in [('0-7', '4491.5', '22258.3'), ('0-8', '4791.4', '20864.7')]:
    add('O2', 'Herzberg II', 'c1Sigma_u- - X3Sigma_g-', band, wl, 'origin (emission)', sg, 'air', 'as printed (0.1)', 'Krupenie72 Table 6(a) (JPCRD 1, p.460), data of Degen', 'nightglow',
        'O+O+M recombination', '~94 (PALACE)', '', 'Calculated origins for other bands in Krupenie72 Table 7.', 'Krupenie72;Meier91;PALACE')
# ---------------- O2 Chamberlain A'3Du - a1Dg (Krupenie72 Table 9, tentative) ----------------
add('O2', 'Chamberlain', "A'3Delta_u - a1Delta_g", 'various (tentative)', '3698-4378', 'head range', '', 'air', 'tentative heads (to 1 A)', 'Krupenie72 Table 9 (JPCRD 1, p.461), data of Chamberlain; quantum numbering tentative',
    'nightglow', 'O+O+M recombination', '~94 (PALACE)', '', 'Krupenie72 labels upper state C 3Delta_u.', 'Krupenie72;Meier91')
# ---------------- NO gamma / delta (Meier91 Tables II, III) ----------------
for band, wl, Rd, Rn in [('1-0', '2149', '752', ''), ('2-0', '2047', '197', ''), ('3-1', '2030', '29', ''), ('2-2', '2216', '82', ''), ('1-1', '2239', '206', ''),
                         ('0-0', '2262', '265', '1.8'), ('0-1', '2363', '369', '2.4'), ('0-2', '2471', '', '3.8'), ('0-3', '2587', '', '1.5'), ('0-4', '2713', '', '0.5')]:
    add('NO', 'gamma', 'A2Sigma+ - X2Pi', band, wl, 'nominal band position', '', 'unspecified', '~1-10 (15 A resolution spectra)', 'Meier91 Tables II/III', 'dayglow|nightglow',
        'day: solar-stimulated fluorescence; night: N(4S)+O(3P) radiative association', '~80-140 (day, Noxon68)',
        '; '.join(x for x in [f'day {Rd} R nadir [Meier91 T.II]' if Rd else '', f'night {Rn} R nadir [Meier91 T.III]' if Rn else '',
                               'day strongest NO-gamma band ~1 kR zenith [Noxon68]' if band == '1-0' else ''] if x), '', 'Meier91;Noxon68')
for band, wl, Rn in [('0-0', '1900', '5.5'), ('0-1', '1980', '7.9'), ('0-2', '2055', '6.5'), ('0-3', '2134', '4.0'), ('0-4', '2219', '2.1'),
                     ('0-5', '2310', '0.93'), ('0-6', '2406', '0.38'), ('0-7', '2509', '0.16'), ('0-8', '2620', '0.06')]:
    add('NO', 'delta', 'C2Pi - X2Pi', band, wl, 'nominal band position', '', 'unspecified', '~1-10', 'Meier91 Table III', 'nightglow',
        'N(4S)+O(3P) radiative association (night)', '~115-150 (rocket view, Meier91 Fig.11)', f'night {Rn} R nadir [Meier91 T.III]', '', 'Meier91')
# ---------------- continua ----------------
add('FeO', 'orange arc', "D5Delta_i, D'5Delta_i - X5Delta_i", 'unresolved', '5000-7200 (main peak ~5950)', 'pseudo-continuum range / peak', '', 'n/a', 'peak to ~50 A', 'PALACE Sect. 5.2', 'nightglow',
    'Fe + O3 -> FeO* + O2', '~88 (WACCM centroid, PALACE)', '', 'Discovered with Odin/OSIRIS (Evans10). Correlated shallow maximum near 8000 A may include other molecules (NiO, OFeOH).', 'PALACE;Evans10')
add('NO2', 'air afterglow continuum', 'NO + O -> NO2* (chemiluminescence)', 'continuum', 'visible-NIR', 'continuum', '', 'n/a', 'n/a', 'PALACE Sect. 5.2 (citing Gattinger10)', 'nightglow',
    'NO + O -> NO2 + hv', '', '', 'Named as a minor possible contributor to the visible nightglow continuum; no band positions (true continuum).', 'PALACE;Gattinger10')
add('HO2', 'NIR bands', "X2A'' (200-000), A2A'(001)-X2A''(000)", 'unresolved', '~15100 peak; features ~12700 and ~16200', 'feature peaks', '', 'n/a', '~100', 'PALACE Sect. 5.1', 'nightglow',
    'H + O2 + M -> HO2* + M', '~81 (PALACE)', 'night 13.1 kR (0.3-2.5 um model continuum) [PALACE]', '', 'PALACE')
add('O2', 'UV-blue band continuum (c-b, A\'-a, A-X unresolved)', 'c1Sigma_u-, A\'3Delta_u, A3Sigma_u+ upper states', 'unresolved', '3000-6800', 'range', '', 'n/a', 'n/a', 'PALACE Sect. 5.3', 'nightglow',
    'O+O+M recombination', '~94 (PALACE)', 'night 405 R (339 R excluding 300-320 nm) [PALACE]', '', 'PALACE')

# ---------------- OH Meinel (derived origins) ----------------
src = 'derived: T_v from fit of 2Pi (T,A,B,D,gamma) Hamiltonian to HITRAN2020 16OH X-state level energies (J<=12.5), nu0 = T_v\' - T_v\'\'; scripts/research/oh_band_origins.py'
pal = {'4-2': 'night 92 kR (strongest band) [PALACE]', '8-5': 'night 24 kR (strongest dv=3) [PALACE]', '9-5': 'night 3.8 kR [PALACE]',
       '9-4': 'night 570 R [PALACE]', '9-3': 'night 61 R [PALACE]', '9-2': 'night 5.1 R [PALACE]'}
for r in csv.DictReader(open(oh_csv)):
    b = r['band']
    note = (f"vacuum {r['wl_vac_A']} A. Fit unc (rms-based) {r['nu0_fit_unc_cm-1']} cm-1. Huber-Herzberg/Maillard G(v) polynomial gives {r['nu0_HH_poly_cm-1']} cm-1 "
            f"(diff {r['HH_minus_fit_cm-1']} cm-1; polynomial valid only v<=5). A band origin is a model quantity: no line lies at it; individual lines are listed in HITRAN / sky atlases.")
    add('OH', f"Meinel dv={r['dv']}", 'X2Pi (v\') - X2Pi (v\'\')', b, r['wl_air_A'] or r['wl_vac_A'], 'band origin (derived)', r['nu0_cm-1'],
        'air' if r['wl_air_A'] else 'vacuum', f"~{max(0.01, round(float(r['nu0_fit_unc_cm-1']) * (float(r['wl_vac_A']) / 1e4) ** 2 * 1e0, 3))} (from fit; definition-dependent)",
        src, 'nightglow|dayglow', 'H + O3 -> OH(v<=9) + O2 (PALACE R1)', '~87, FWHM ~8 (Baker & Stair 1988 via PALACE); higher v emit a few km higher',
        pal.get(b, ''), note + (' Total OH 0.3-2.5 um 715 kR; 642-858 nm 3.1 kR [PALACE].' if b == '4-2' else ''), 'HITRAN2020;HH79-WebBook;PALACE')

# ---- origin wavelengths from sigma_0 and medium check (only where lambda and sigma refer to the same feature) ----
out = []
for r in ROWS:
    chk = ''
    r['origin_wl_vac_A'] = r['origin_wl_air_A'] = ''
    try:
        sg = float(r['sigma_cm']); lv = 1e8 / sg
    except (ValueError, ZeroDivisionError):
        out.append(r); r['medium_check'] = ''; continue
    if r['sigma_kind'] == 'origin':
        r['origin_wl_vac_A'] = f'{lv:.3f}'
        r['origin_wl_air_A'] = f'{lv / n_air(lv):.3f}' if 2000 < lv < 20000 else ''
    same = (r['sigma_kind'] == 'head') == ('head' in r['wl_kind'])
    if same and r['molecule'] != 'OH':
        try:
            wl = float(r['wl_A'].strip('()~'))
            chk = f'lambda - vac(sigma) = {wl - lv:+.2f} A' + (f'; lambda - air(sigma) = {wl - lv / n_air(lv):+.2f} A' if lv > 2000 else '')
        except ValueError:
            pass
    r['medium_check'] = chk
    out.append(r)
cols = ['molecule', 'system', 'transition', 'band', 'wl_A', 'wl_kind', 'medium', 'precision_A', 'sigma_cm', 'sigma_kind', 'origin_wl_vac_A', 'origin_wl_air_A',
        'medium_check', 'wl_source', 'category', 'mechanism', 'altitude_km', 'typical_R', 'note', 'refs']
with open(out_csv, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore')
    w.writeheader()
    for r in out:
        w.writerow({c: r.get(c, '') for c in cols})
print(len(out), 'band rows ->', out_csv)
for r in out:
    if r['medium_check']:
        print(f"{r['molecule']:4s} {r['system'][:22]:22s} {r['band']:6s} {r['wl_kind'][:12]:12s} {r['wl_A']:>10s} {r['medium']:7s} {r['medium_check']}")
