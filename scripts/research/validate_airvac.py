#!/usr/bin/env python
"""Validate airvac_formulas.py against published test values.

Usage (run isolated, downloaded files are passed as paths):
  python -I validate_airvac.py <specutils_wcs_utils.py> <asd_check_dir> <out_csv>

Every reference value below is transcribed from the cited source (file in data/raw/airvac/...).
"""
import os, sys, csv, ast, re, glob, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # our own scripts dir only
import numpy as np
import airvac_formulas as F

SPECUTILS_SRC, ASD_DIR, OUT_CSV = sys.argv[1], sys.argv[2], sys.argv[3]
rows = []


def rec(test, source, inputs, quantity, published, computed, unit, tol, note='', expect='match'):
    """expect='match': formula should reproduce the value; expect='differ': discriminator test that a
    rival formula should NOT reproduce (PASS means it indeed differs by more than tol)."""
    diff = computed - published
    ok = abs(diff) <= tol
    if expect == 'match':
        status = 'PASS' if ok else 'FAIL'
    elif expect == 'differ':
        status = 'PASS' if not ok else 'FAIL'
    else:  # 'info' or 'misprint?': reported, not scored
        status = 'NOTE'
    rows.append(dict(test=test, source=source, inputs=inputs, quantity=quantity, published=published,
                     computed=round(float(computed), 9), diff=round(float(diff), 9), unit=unit, tol=tol,
                     expect=expect, status=status, note=note))
    print(f"{status} [{expect:6s}] {test:38s} {inputs:44s} pub={published!s:>14} calc={computed:.6f} d={diff:+.6f} {unit}")


# ---------------- 1. Ciddor (1996) Table 1: dry air, 450 ppm, lambda = 633 nm --------------------
for t, pk, val in [(20, 80, 21458.0), (20, 100, 26824.4), (20, 120, 32191.6), (10, 100, 27774.7), (30, 100, 25937.2)]:
    n = F.n_ciddor1996(6330.0, t_C=t, p_Pa=pk * 1e3, rh=0, xc=450)
    rec('Ciddor1996 Table1 (Eq.5)', 'Ciddor 1996 AO 35,1566 Table 1', f't={t}C p={pk}kPa dry 450ppm 633nm',
        '1e8(n-1)', val, (n - 1) * 1e8, '1e-8', 0.06)

# ---------------- 2. Ciddor (1996) Table 2 (moist air, 633 nm), column 'Calculated, Eq. (5)' ----------
for t, p, pw, xc, val in [(19.526, 102094.8, 1065, 510, 27392.9), (19.517, 102096.8, 1065, 510, 27394.0),
                          (19.173, 102993.0, 641, 450, 27682.4), (19.173, 103006.0, 642, 440, 27685.8),
                          (19.188, 102918.8, 706, 450, 27658.7), (19.189, 102927.8, 708, 440, 27660.8),
                          (19.532, 103603.2, 986, 600, 27802.0), (19.534, 103596.2, 962, 600, 27800.8),
                          (19.534, 103599.2, 951, 610, 27802.2)]:
    n = F.n_ciddor1996(6330.0, t_C=t, p_Pa=p, pw_Pa=pw, xc=xc)
    mis = (t == 19.517)
    rec('Ciddor1996 Table2 (Eq.5)', 'Ciddor 1996 Table 2', f't={t} p={p} pw={pw} xc={xc}', '1e8(n-1)', val,
        (n - 1) * 1e8, '1e-8', 0.06, expect='misprint?' if mis else 'match',
        note=('printed 27394.0; row-to-row change of Ciddor\'s own Ref.4 column (+1.4) and of B&D 1993 Table 4 '
              '(+1.4) implies ~27394.3 = our value') if mis else '')

# ---------------- 3. Ciddor (1996) Table 3 extreme conditions (Eq.5 column; B&D column) --------------
# B&D values: Ciddor footnote d: B&D 1993 eq.(5) with denominator constant 96094.69, 1994 dispersion.
for t, pk, rh, pw, ciddor, bd in [(20, 80, 75, 1754, 21394.0, 21394.8), (20, 120, 75, 1754, 32127.8, 32129.2),
                                  (40, 80, 75, 5539, 19996.5, 19883.6), (40, 120, 75, 5539, 29941.8, 29928.6),
                                  (50, 80, 100, 12352, 19058.4, 19012.9), (50, 120, 100, 12352, 28792.4, 28745.2)]:
    n = F.n_ciddor1996(6330.0, t_C=t, p_Pa=pk * 1e3, pw_Pa=pw, xc=450)
    mis = (t == 40 and pk == 80)
    rec('Ciddor1996 Table3 (Eq.5)', 'Ciddor 1996 Table 3', f't={t} p={pk}kPa RH={rh}% pw={pw}', '1e8(n-1)', ciddor,
        (n - 1) * 1e8, '1e-8', 0.06, expect='misprint?' if mis else 'match',
        note='printed 19996.5 conflicts with printed B&D 19883.6 and difference -12.9 -> 19896.5 implied (= our 19896.49)' if mis else '')
    # Ciddor's B&D column is reproduced by Edlen-1966 dispersion x 1.0000799 (CO2 factor only) and D = 96094.69
    s2_ = F._s2(6330.0); p_ = pk * 1e3
    X_ = (1 + 1e-8 * (0.601 - 0.00972 * t) * p_) / (1 + 0.0036610 * t)
    nbh = p_ * F.nm1_edlen1966(s2_) * 1.0000799 * X_ / 96094.69 - 1e-10 * (3.7345 - 0.0401 * s2_) * pw
    rec('Ciddor Tab3 B&D column (reconstructed)', 'Ciddor 1996 Table 3, col. Refs. 3 and 4', f't={t} p={pk}kPa pw={pw}',
        '1e8(n-1)', bd, nbh * 1e8, '1e-8', 0.15,
        note='reproduced with Edlen1966 x 1.0000799 dispersion, D=96094.69 (footnote d), B&D water term')
    nb = F.n_birchdowns(6330.0, t_C=t, p_Pa=pk * 1e3, pv_Pa=pw, version=1994, D=96095.43)
    rec('BirchDowns1994 (as published) vs Tab3', 'Ciddor 1996 Table 3, col. Refs. 3 and 4', f't={t} p={pk}kPa pw={pw}',
        '1e8(n-1)', bd, (nb - 1) * 1e8, '1e-8', 1.2, expect='info',
        note='B&D 1994 constants give values 0.7-1.2e-8 lower (dispersion factor 1.0000486 vs 1.0000799)')

# ---------------- 4. NIST EMToolBox Table 1 (Ciddor & modified Edlen, 450 ppm) ------------------------
emt = [(20, 0, 101.325, 633, 1.000271800, 1.000271799), (20, 0, 60, 633, 1.000160924, 1.000160920),
       (20, 0, 120, 633, 1.000321916, 1.000321918), (50, 0, 100, 633, 1.000243285, 1.000243270),
       (5, 0, 100, 633, 1.000282756, 1.000282750), (-40, 0, 100, 633, 1.000337580, 1.000337471),
       (50, 100, 120, 633, 1.000287924, 1.000287864), (40, 75, 120, 633, 1.000299418, 1.000299406),
       (20, 100, 100, 633, 1.000267394, 1.000267394), (40, 100, 110, 1700, 1.000270247, 1.000270237),
       (20, 0, 101.325, 1700, 1.000268479, 1.000268483), (40, 100, 110, 300, 1.000289000, 1.000288922),
       (20, 0, 101.325, 300, 1.000286581, 1.000286579), (-40, 0, 120, 300, 1.000427233, 1.000427072)]
for t, rh, pk, lam, nc, ne in emt:
    T = t + 273.15
    svp = F.svp_iapws(T) if t >= 0 else F.svp_iapws_ice(T)
    n1 = F.n_ciddor1996(lam * 10.0, t_C=t, p_Pa=pk * 1e3, rh=rh / 100, xc=450, variant='nist')
    rec('NIST EMToolBox Ciddor', 'emtoolbox.nist.gov Wavelength/Documentation.asp Table 1',
        f't={t} RH={rh} p={pk}kPa lam={lam}nm', '1e8(n-1)', round((nc - 1) * 1e8, 1), (n1 - 1) * 1e8, '1e-8', 0.1)
    n2 = F.n_birchdowns(lam * 10.0, t_C=t, p_Pa=pk * 1e3, pv_Pa=rh / 100 * svp, version=1994, nist_water=True)
    rec('NIST EMToolBox mod. Edlen', 'emtoolbox Table 1', f't={t} RH={rh} p={pk}kPa lam={lam}nm', '1e8(n-1)',
        round((ne - 1) * 1e8, 1), (n2 - 1) * 1e8, '1e-8', 0.1)
for t, val in [(100, 101418), (50, 12351), (40, 7384), (20, 2339), (-10, 260)]:
    T = t + 273.15
    v = F.svp_iapws(T) if t >= 0 else F.svp_iapws_ice(T)
    rec('NIST EMToolBox svp (IAPWS)', 'emtoolbox Table 2', f't={t}C', 'svp Pa', val, float(v), 'Pa', 0.6)

# ---------------- 5. Birch & Downs (1993) Table 3, 633 nm (632.991 nm vac He-Ne) --------------------
for pk, t, orig, rev in [(80, 20, 21456.8, 21459.0), (100, 20, 26823.2, 26826.0), (120, 20, 32190.5, 32193.8),
                         (100, 10, 27773.3, 27776.1), (100, 30, 25935.7, 25938.5)]:
    ne = F.n_edlen1966(6329.91, t_C=t, p_Pa=pk * 1e3, f_Pa=0, x_co2=0.0003)
    rec('Edlen1966 T,p (B&D93 Tab3 orig.)', 'Birch & Downs 1993 Table 3', f'p={pk}kPa t={t} dry 632.991nm',
        '1e8(n-1)', orig, (ne - 1) * 1e8, '1e-8', 0.15)
    nb = F.n_birchdowns(6329.91, t_C=t, p_Pa=pk * 1e3, pv_Pa=0, version=1993)
    rec('BirchDowns1993 (Tab3 revised)', 'Birch & Downs 1993 Table 3', f'p={pk}kPa t={t} dry 632.991nm',
        '1e8(n-1)', rev, (nb - 1) * 1e8, '1e-8', 0.15)
rec('Edlen1966 (n-1)s 633nm', 'Birch & Downs 1993 App. A.1', '632.991 nm, 300 ppm', '1e6(n-1)', 276.5175,
    F.nm1_edlen1966(F._s2(6329.91)) * 1e6, '1e-6', 0.0003)
rec('Edlen1966+CO2 eq17 633nm', 'Birch & Downs 1993 App. A.1', '632.991 nm, 450 ppm', '1e6(n-1)', 276.5399,
    (F.n_edlen1966(6329.91, x_co2=0.00045) - 1) * 1e6, '1e-6', 0.0003)

# ---------------- 6. Peck & Reeder (1972) Tables III and V ------------------------------------------
# Table III: vacuum wavelength (um), Eq.(2) prediction, and [Eq.(1)-Eq.(2)] (Eq.(1) = Edlen 1966)
tab3 = [(1.6945208, 27314.19, 0.52), (1.5300153, 27326.37, 0.50), (1.3722327, 27342.42, 0.47),
        (1.0142573, 27410.90, 0.36), (0.7247163, 27557.37, 0.18), (0.6440249, 27638.092, 0.102),
        (0.5462271, 27789.843, 0.009), (0.5017074, 27891.525, -0.026), (0.4923304, 27916.714, -0.032),
        (0.4472732, 28062.082, -0.048), (0.4359562, 28106.335, -0.048), (0.4047714, 28249.611, -0.036),
        (0.3889751, 28336.843, -0.023), (0.3391680, 28705.90, 0.03), (0.2926300, 29264.37, -0.02),
        (0.2675750, 29719.74, -0.23), (0.2464820, 30246.41, -0.36), (0.2302890, 30787.68, 0.45)]
for lam, eq2, d12 in tab3:
    dig = 3 if abs(eq2 * 1000 - round(eq2 * 1000)) < 1e-6 and len(f'{eq2}'.split('.')[1]) == 3 else 2
    tol = 0.0015 if dig == 3 else 0.015
    rec('PeckReeder1972 Eq.(2)', 'Peck & Reeder 1972 Table III', f'lam_vac={lam}um', '1e8(n-1)', eq2,
        F.nm1_peckreeder1972(1 / lam ** 2) * 1e8, '1e-8', tol)
    rec('Edlen1966 = P&R Eq.(1)', 'Peck & Reeder 1972 Table III (Eq2 + [Eq1-Eq2])', f'lam_vac={lam}um',
        '1e8(n-1)', round(eq2 + d12, 3), F.nm1_edlen1966(1 / lam ** 2) * 1e8, '1e-8', 2 * tol)
# Table V: far UV, Eq.(3) column; Eq.(2) and Eq.(1) from data minus residual columns
tab5 = [(0.214506, 31496.8, 31496.6, 5.5, -0.3), (0.202605, 32214.8, 32217.4, 18.8, -3.6),
        (0.199052, 32479.5, 32478.5, 32.7, -0.4), (0.193585, 32939.7, 32937.5, 60.4, -0.1),
        (0.186277, 33697.4, 33702.7, 129, -11.1), (0.185473, 33805.5, 33801.4, 152, -2.4)]
for lam, data, eq3, dd2, dd1 in tab5:
    rec('PeckReeder1972 Eq.(3) 5-param', 'Peck & Reeder 1972 Table V', f'lam={lam}um', '1e8(n-1)', eq3,
        F.nm1_peckreeder1972_5p(1 / lam ** 2) * 1e8, '1e-8', 0.06)
    rec('PeckReeder1972 Eq.(2) far-UV', 'P&R 1972 Table V (Data-[Data-Eq2])', f'lam={lam}um', '1e8(n-1)',
        round(data - dd2, 1), F.nm1_peckreeder1972(1 / lam ** 2) * 1e8, '1e-8', 0.6)
    rec('Edlen1966 far-UV', 'P&R 1972 Table V (Data-[Data-Eq1])', f'lam={lam}um', '1e8(n-1)', round(data - dd1, 1),
        F.nm1_edlen1966(1 / lam ** 2) * 1e8, '1e-8', 0.06)

# ---------------- 7. Edlen (1953) Table I 'formula (1)' = CO2-free air? ------------------------------
# Edlen 1966 states Edlen 1953 used factor 1.000162 to pass from CO2-free to standard air.
for s2, val in [(2.4110, 27633.6), (3.3516, 27785.6), (5.2616, 28102.4), (10.1905, 28977.3), (16.4601, 30240.1),
                (21.7330, 31491.0), (28.8192, 33697.3)]:
    calc = F.nm1_edlen1953(s2) * 1e8 / 1.000162
    rec('Edlen1953 / 1.000162 vs Tab.I', 'Edlen 1953 JOSA 43,339 Table I', f'sigma2={s2}', '1e8(n-1)', val, calc,
        '1e-8', 0.1, 'Table I compares with CO2-free air: standard-air formula / 1.000162 (Edlen 1966 p.79)')

# ---------------- 8. Single-wavelength published examples (A) ----------------------------------------
rec('Edlen1953 (Morton 1991)', 'Morton 1991 ApJS 77,119 p.~121', 'lam_air=2000 A', 'lam_vac A', 2000.648,
    float(F.air_to_vac_edlen1953(2000.0)), 'A', 0.0005)
rec('NIST ASD (P&R)', 'ASD lineshelp: 200.0648 nm vac = 200.0000 nm air', 'lam_vac=2000.648 A', 'lam_air A',
    2000.000, float(F.vac_to_air_peckreeder1972_5p(2000.648)), 'A', 0.0005)
rec('NIST ASD (P&R)', 'ASD lineshelp: 199.9352 nm air = 200.0000 nm vac', 'lam_vac=2000.000 A', 'lam_air A',
    1999.352, float(F.vac_to_air_peckreeder1972_5p(2000.0)), 'A', 0.0005)
rec('NIST ASD (P&R eq2, 4-param)', 'same, testing 4-param alternative', 'lam_vac=2000.000 A', 'lam_air A',
    1999.352, float(F.vac_to_air_peckreeder1972(2000.0)), 'A', 0.0005, expect='differ')
rec('IDL vactoair (Ciddor)', 'IDLAstro vactoair.pro header example', 'lam_vac=2000 A', 'lam_air A', 1999.353,
    float(F.vac_to_air_idl_astrolib_literal(2000.0)), 'A', 0.0005)
rec('IDL airtovac (Ciddor)', 'IDLAstro airtovac.pro header example (Kr line)', 'lam_air=6056.125 A',
    'lam_vac A', 6057.8019, float(F.air_to_vac_idl_astrolib_literal(6056.125)), 'A', 0.00005)
for nm, nm1e8, corr in [(200, 32409, 0.06480), (300, 29157, 0.08745), (400, 28276, 0.11307), (500, 27897, 0.13945),
                        (600, 27698, 0.16614), (700, 27581, 0.19301), (800, 27505, 0.21998), (900, 27454, 0.24701),
                        (1000, 27417, 0.2741), (1500, 27331, 0.4099), (2000, 27301, 0.5459)]:
    lv = nm * 10.0
    rec('CRC87 table (B&D 1994)', 'CRC Handbook 87th ed. p.10-252', f'lam_vac={nm}nm', '1e8(n-1)', nm1e8,
        F.nm1_birchdowns1994(F._s2(lv)) * 1e8, '1e-8', 0.6)
    rec('CRC87 table (B&D 1994)', 'CRC Handbook 87th ed. p.10-252', f'lam_vac={nm}nm', 'lam_vac-lam_air nm', corr,
        float(lv - F.vac_to_air_birchdowns1994(lv)) / 10.0, 'nm', 0.000006 if nm < 1000 else 0.00006)

# SDSS line tables (vac -> air), DR9+ = Ciddor; DR7 = polynomial
sdss17 = [('Hbeta', 4861.325, 4862.683), ('[OIII]', 4958.911, 4960.295), ('[OIII]', 5006.843, 5008.239),
          ('[NII]', 6548.05, 6549.86), ('Halpha', 6562.801, 6564.614), ('[NII]', 6583.45, 6585.27),
          ('[SII]', 6716.44, 6718.29), ('[SII]', 6730.82, 6732.68)]
for name, air, vac in sdss17:
    nd = len(f'{air}'.split('.')[1])
    tol = 0.5 * 10 ** -nd + 0.5 * 10 ** -len(f'{vac}'.split('.')[1])
    rec('SDSS DR17 table (Ciddor)', 'sdss4.org/dr17/spectro/spectro_basics', f'{name} vac={vac}', 'lam_air A', air,
        float(F.vac_to_air_ciddor1996(vac)), 'A', tol)
for name, air, vac in [('Hbeta', 4861.363, 4862.721), ('[OIII]', 5006.843, 5008.239), ('Halpha', 6562.801, 6564.614)]:
    rec('SDSS DR7 table vs DR7 poly', 'classic.sdss.org/dr7/products/spectra/vacwavelength', f'{name} vac={vac}',
        'lam_air A', air, float(F.vac_to_air_sdss_dr7(vac)), 'A', 0.0005, expect='differ',
        note='the DR7 table pairs are NOT generated by the DR7 polynomial (poly is 1-3 mA low)')
    rec('SDSS DR7 table vs Ciddor', 'classic.sdss.org/dr7/products/spectra/vacwavelength', f'{name} vac={vac}',
        'lam_air A', air, float(F.vac_to_air_ciddor1996(vac)), 'A', 0.001 if name != 'Hbeta' else 0.05,
        note='DR7 Hbeta air 4861.363 is the older (pre-NIST-2009) air value; vac 4862.721 not 4862.683')

# ---------------- 9. Na D lines in NIST ASD (vacuum vs air output of the same lines) ------------------
for name, vac, air in [('Na D2 obs', 5891.583264, 5889.95095), ('Na D2 Ritz', 5891.583253, 5889.95094),
                       ('Na D1 obs', 5897.558147, 5895.92424)]:
    for fname, fn in [('P&R eq3 (5p)', F.vac_to_air_peckreeder1972_5p), ('P&R eq2 (4p)', F.vac_to_air_peckreeder1972),
                      ('Ciddor', F.vac_to_air_ciddor1996), ('Morton2000', F.vac_to_air_morton2000),
                      ('Edlen1953', F.vac_to_air_edlen1953)]:
        rec(f'ASD Na D: {fname}', 'NIST ASD lines1.pl show_av=3 vs 2 (2026-11-23)', f'{name} vac={vac}', 'lam_air A',
            air, float(fn(vac)), 'A', 0.0000101, expect='match' if fname.startswith('P&R eq3') else 'info',
            note='tolerance = half ULP of vac (1e-6) + air (1e-5) rounding')

# ---------------- 10. astropy.wcs (WCSLIB) AWAV -> WAVE ----------------------------------------------
try:
    from astropy import wcs as awcs
    import astropy
    for lam_air in [2000.0, 3000.0, 5000.0, 6562.8, 10000.0, 20000.0, 25000.0]:
        w = awcs.WCS(naxis=1)
        w.wcs.ctype = ['AWAV']
        w.wcs.cunit = ['m']
        w.wcs.crval = [lam_air * 1e-10]
        w.wcs.cdelt = [1e-13]
        w.wcs.crpix = [1]
        w.wcs.set()
        w2 = w.deepcopy()
        w2.wcs.sptr('WAVE-A2W')
        w2.wcs.set()
        lv = w2.wcs.crval[0] * 1e10
        rec('astropy.wcs AWAV->WAVE vs wcslib', f'astropy {astropy.__version__} wcs.sptr(WAVE-A2W)',
            f'lam_air={lam_air}', 'lam_vac A', round(float(F.air_to_vac_wcslib(lam_air)), 7), lv, 'A', 2e-6,
            'live call; reference column is our re-implementation of spx.c')
        rec('astropy.wcs AWAV->WAVE vs Edlen1953 exact', 'Edlen 1953 (vacuum-argument)', f'lam_air={lam_air}',
            'lam_vac A', round(float(F.air_to_vac_edlen1953(lam_air)), 7), lv, 'A', 2e-4,
            note='difference = effect of WCSLIB evaluating Edlen 1953 at the air wavenumber')
except Exception as e:  # noqa
    print('astropy.wcs test skipped:', e)

# ---------------- 11. specutils 2.4.0 functions, extracted from the wheel and executed ---------------
src = open(SPECUTILS_SRC).read()
tree = ast.parse(src)
want = {'refraction_index', 'vac_to_air', 'air_to_vac'}
mod = ast.Module(body=[n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in want], type_ignores=[])
import astropy.units as u
from astropy.utils.exceptions import AstropyDeprecationWarning
ns = dict(u=u, warnings=warnings, AstropyDeprecationWarning=AstropyDeprecationWarning, np=np)
exec(compile(mod, SPECUTILS_SRC, 'exec'), ns)
grid = np.array([2000.0, 3000.0, 5000.0, 10000.0, 20000.0])
for meth, mine in [('Morton2000', F.vac_to_air_morton2000), ('Edlen1953', F.vac_to_air_edlen1953),
                   ('Edlen1966', F.vac_to_air_edlen1966), ('PeckReeder1972', F.vac_to_air_peckreeder1972),
                   ('Ciddor1996', F.vac_to_air_ciddor1996), ('Greisen2006', F.vac_to_air_greisen2006_literal)]:
    su = ns['vac_to_air'](grid * u.AA, method=meth).to_value(u.AA)
    for g, s, m in zip(grid, su, mine(grid)):
        rec(f'specutils vac_to_air {meth}', 'specutils 2.4.0 wcs_utils (executed)', f'lam_vac={g}', 'lam_air A',
            round(float(s), 7), float(m), 'A', 1e-6)
    sinv = ns['air_to_vac'](grid * u.AA, scheme='inversion', method=meth).to_value(u.AA)
    exact = {'Morton2000': F.air_to_vac_morton2000, 'Edlen1953': F.air_to_vac_edlen1953,
             'Edlen1966': F.air_to_vac_edlen1966, 'PeckReeder1972': F.air_to_vac_peckreeder1972,
             'Ciddor1996': F.air_to_vac_ciddor1996, 'Greisen2006': F.air_to_vac_greisen2006}[meth](grid)
    for g, s, m in zip(grid, sinv, exact):
        rec(f'specutils air_to_vac inversion {meth}', 'specutils 2.4.0 (executed) vs exact inverse',
            f'lam_air={g}', 'lam_vac A', round(float(m), 7), float(s), 'A', 2e-4,
            'difference = error of the "inversion" scheme (n evaluated at the air wavelength)')
    sit = ns['air_to_vac'](grid * u.AA, scheme='iteration', method=meth).to_value(u.AA)
    exact_it = exact if meth != 'Greisen2006' else np.array([
        F._invert_vacarg(lambda s2: F.nm1_iag1999_phase(1.0 / np.sqrt(s2)), g) for g in grid])
    for g, s, m in zip(grid, sit, exact_it):
        rec(f'specutils air_to_vac iteration {meth}', 'specutils 2.4.0 (executed)', f'lam_air={g}', 'lam_vac A',
            round(float(m), 7), float(s), 'A', 1e-6)
sp = ns['air_to_vac'](grid * u.AA, scheme='Piskunov').to_value(u.AA)
for g, s, m in zip(grid, sp, F.air_to_vac_vald(grid)):
    rec('specutils air_to_vac Piskunov', 'specutils 2.4.0 (executed)', f'lam_air={g}', 'lam_vac A',
        round(float(m), 7), float(s), 'A', 1e-6)


# ---------------- 12. NIST ASD: which formula reproduces its air wavelengths? ------------------------
def _clean(v):
    v = v.strip().replace('=', '').replace('"', '')
    v = re.sub(r'[^0-9.eE+-]', '', v.rstrip('+*?'))
    return v


def load_asd(fn):
    """Parse ASD CSV; a new header line appears whenever the wavelength medium switches.
    Returns list of dicts with keys: obs, ritz (strings), medium ('vac'/'air'), key tuple."""
    out = []
    hdr = None
    with open(fn) as fh:
        for row in csv.reader(fh):
            if not row:
                continue
            if row[0].startswith('obs_wl') or row[0].startswith('ritz_wl'):
                hdr = row
                continue
            if hdr is None or len(row) < 5:
                continue
            d = dict(zip(hdr, row))
            med = 'air' if any('_air(' in h for h in hdr) else 'vac'
            obs = next((d[h] for h in hdr if h.startswith('obs_wl')), '')
            ritz = next((d[h] for h in hdr if h.startswith('ritz_wl')), '')
            key = (d.get('Ei(cm-1)'), d.get('Ek(cm-1)'), d.get('J_i'), d.get('J_k'), d.get('wn(cm-1)'),
                   d.get('conf_i'), d.get('conf_k'))
            out.append(dict(obs=_clean(obs), ritz=_clean(ritz), medium=med, key=key))
    return out


pairs = []
switch_report = []
for f3 in sorted(glob.glob(os.path.join(ASD_DIR, '*_av3.csv'))):
    spec = os.path.basename(f3)[:-8]
    r3 = load_asd(f3)
    for av in (4, 2):
        fx = f3.replace('_av3', f'_av{av}')
        if not os.path.exists(fx):
            continue
        rx = load_asd(fx)
        idx = {}
        for d in rx:
            idx.setdefault(d['key'], []).append(d)
        for d in r3:
            if d['medium'] != 'vac' or d['key'] not in idx or len(idx[d['key']]) != 1:
                continue
            e = idx[d['key']][0]
            for col in ('ritz', 'obs'):
                sv, sx = d[col], e[col]
                if not sv or not sx:
                    continue
                try:
                    lv, lx = float(sv), float(sx)
                except ValueError:
                    continue
                ndv = len(sv.split('.')[1]) if '.' in sv else 0
                ndx = len(sx.split('.')[1]) if '.' in sx else 0
                pairs.append((spec, av, col, lv, lx, ndv, ndx, e['medium']))
# medium switch boundaries actually used by ASD
for av in (4, 2):
    for med in ('vac', 'air'):
        lv = [p[3] for p in pairs if p[1] == av and p[7] == med]
        if lv:
            print(f'ASD show_av={av}: medium={med}: lam_vac range of lines {min(lv):.3f} - {max(lv):.3f} A (n={len(lv)})')
            switch_report.append((av, med, min(lv), max(lv), len(lv)))
asd_rows = []
P = [p for p in pairs if p[1] == 4 and p[7] == 'air' and p[5] >= 4 and p[6] >= 4]
print(f'ASD pairs (show_av=4 air vs show_av=3 vac, >=4 decimals both): {len(P)}')
if P:
    lv = np.array([p[3] for p in P]); la = np.array([p[4] for p in P])
    ndv = np.array([p[5] for p in P]); nda = np.array([p[6] for p in P])
    quant = 0.5 * 10.0 ** -nda + 0.5 * 10.0 ** -ndv * la / lv
    cands = {'peckreeder1972_5p (eq.3)': F.vac_to_air_peckreeder1972_5p,
             'peckreeder1972 (eq.2)': F.vac_to_air_peckreeder1972,
             'ciddor1996': F.vac_to_air_ciddor1996, 'birchdowns1994/morton2000': F.vac_to_air_birchdowns1994,
             'edlen1966': F.vac_to_air_edlen1966, 'edlen1953': F.vac_to_air_edlen1953}
    bins = [(1850, 2000), (2000, 3000), (3000, 5000), (5000, 10000), (10000, 20000), (20000, 26000)]
    for name, fn in cands.items():
        res = la - fn(lv)
        print(f'  {name:28s} all n={len(res)}: mean={np.mean(res)*1e3:+.4f} mA rms={np.std(res)*1e3:.4f} mA '
              f'within-rounding={np.mean(np.abs(res) <= quant + 1e-9)*100:.2f}%')
        for lo, hi in bins:
            m = (lv >= lo) & (lv < hi)
            if m.sum() < 3:
                continue
            row = dict(formula=name, lam_lo_A=lo, lam_hi_A=hi, n_lines=int(m.sum()),
                       mean_resid_mA=round(float(np.mean(res[m])) * 1e3, 5),
                       rms_resid_mA=round(float(np.std(res[m])) * 1e3, 5),
                       max_abs_resid_mA=round(float(np.max(np.abs(res[m]))) * 1e3, 5),
                       pct_within_rounding=round(float(np.mean(np.abs(res[m]) <= quant[m] + 1e-9)) * 100, 2))
            asd_rows.append(row)
            print(f'      {lo:6d}-{hi:6d} A n={m.sum():5d} mean={row["mean_resid_mA"]:+.4f} mA rms={row["rms_resid_mA"]:.4f} '
                  f'max={row["max_abs_resid_mA"]:.4f} within={row["pct_within_rounding"]:.1f}%')
for av, med, lo, hi, n in switch_report:
    asd_rows.append(dict(formula=f'MEDIUM show_av={av} {med}', lam_lo_A=round(lo, 3), lam_hi_A=round(hi, 3),
                         n_lines=n, mean_resid_mA='', rms_resid_mA='', max_abs_resid_mA='', pct_within_rounding=''))

with open(OUT_CSV, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
with open(OUT_CSV.replace('.csv', '_asd.csv'), 'w', newline='') as fh:
    if asd_rows:
        w = csv.DictWriter(fh, fieldnames=list(asd_rows[0].keys()))
        w.writeheader()
        w.writerows(asd_rows)
npass = sum(r['status'] == 'PASS' for r in rows)
print(f'\n{npass}/{len(rows)} checks pass')
