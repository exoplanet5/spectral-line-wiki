#!/usr/bin/env python
"""Compare air<->vacuum formulas over 2000-25000 A and compute Ciddor sensitivities at 5000 A.

Usage: python compare_airvac.py <out_comparison_csv> <out_sensitivity_csv>
All numbers are DERIVED from the formulas in airvac_formulas.py (no external data read).
"""
import os, sys, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import airvac_formulas as F

OUT, OUT_SENS = sys.argv[1], sys.argv[2]
C_KMS = 299792.458
wl = np.arange(2000.0, 25000.0 + 1, 100.0)
s2 = F._s2(wl)
ref_air = F.vac_to_air_ciddor1996(wl)          # reference: Ciddor 1996 standard air (15 C, 101325 Pa, dry, 450 ppm)

cols = {}
cols['wl_vac_A'] = wl
cols['sigma2_um-2'] = s2
cols['nm1_ciddor1996_x1e8'] = F.nm1_ciddor1996(s2) * 1e8
cols['wl_air_ciddor1996_A'] = ref_air
cols['vac_minus_air_ciddor1996_A'] = wl - ref_air
cols['vac_minus_air_ciddor1996_kms'] = (wl - ref_air) / wl * C_KMS

# --- standard-air formulas: lambda_air(formula) - lambda_air(Ciddor) for the same vacuum wavelength, mA
std = [('edlen1953', F.vac_to_air_edlen1953), ('edlen1966', F.vac_to_air_edlen1966),
       ('peckreeder1972_eq2', F.vac_to_air_peckreeder1972), ('peckreeder1972_eq3_NIST_ASD', F.vac_to_air_peckreeder1972_5p),
       ('birchdowns1993', F.vac_to_air_birchdowns1993), ('birchdowns1994_morton2000_VALD', F.vac_to_air_birchdowns1994),
       ('ciddor1996_xc420ppm', lambda w: F.vac_to_air_ciddor1996_co2(w, 420.0)),
       ('ciddor1996_xc300ppm', lambda w: F.vac_to_air_ciddor1996_co2(w, 300.0)),
       ('wcslib_astropy_wcs_AWAV', F.vac_to_air_wcslib),
       ('greisen2006_iag1999_0C', F.vac_to_air_greisen2006),
       ('sdss_dr7_polynomial', F.vac_to_air_sdss_dr7),
       ('iraf_disptrans_default_t15_p760_f4mmHg', F.vac_to_air_iraf_literal),
       ('iraf_disptrans_dry_t15_p760_f0', lambda w: F.vac_to_air_iraf_literal(w, 15.0, 760.0, 0.0))]
for name, fn in std:
    cols[f'dair_{name}_mA'] = (fn(wl) - ref_air) * 1e3

# --- inversion / round-trip errors of software implementations, in lambda_vac (mA)
lam_air_m = F.vac_to_air_morton2000(wl)
cols['err_specutils_inversion_morton2000_mA'] = (F.air_to_vac_specutils_inversion_literal(lam_air_m, 'Morton2000') - wl) * 1e3
cols['err_vald_piskunov_inverse_mA'] = (F.air_to_vac_vald(lam_air_m) - wl) * 1e3
lam_air_c = F.vac_to_air_ciddor1996(wl)
cols['err_idl_airtovac_2iter_mA'] = (F.air_to_vac_idl_astrolib_literal(lam_air_c) - wl) * 1e3
lam_air_g = F.vac_to_air_greisen2006(wl)
cols['err_greisen2006_vac2air_shortcut_mA'] = (F.vac_to_air_greisen2006_literal(wl) - lam_air_g) * 1e3
lam_air_i = F.vac_to_air_iraf_literal(wl)
cols['err_iraf_air2vac_roundtrip_mA'] = (F.air_to_vac_iraf_literal(lam_air_i) - wl) * 1e3
cols['err_wcslib_literal_4iter_mA'] = (F.vac_to_air_wcslib_literal(wl) - F.vac_to_air_wcslib(wl)) * 1e3

# --- non-standard air (Ciddor full model) relative to standard air, mA (illustrative conditions, DERIVED)
cols['dair_lab_20C_101325Pa_RH50_420ppm_mA'] = (F.vac_to_air_ciddor1996_full(wl, t_C=20.0, p_Pa=101325.0, rh=0.5, xc=420.0) - ref_air) * 1e3
p2600 = float(F.isa_pressure_Pa(2600.0))
cols['dair_site_2600m_ISA_0C_RH30_420ppm_mA'] = (F.vac_to_air_ciddor1996_full(wl, t_C=0.0, p_Pa=p2600, rh=0.3, xc=420.0) - ref_air) * 1e3

names = list(cols)
with open(OUT, 'w', newline='') as fh:
    fh.write('# Air/vacuum formula comparison. DERIVED by scripts/research/compare_airvac.py from airvac_formulas.py.\n')
    fh.write('# Reference = Ciddor (1996) standard air (15 C, 101325 Pa, 0 %RH, 450 ppm CO2). dair_* = lambda_air(formula)'
             ' - lambda_air(Ciddor) at fixed lambda_vac, in milli-Angstrom (mA). err_* = error of a software'
             ' implementation in lambda_vac (or lambda_air for err_greisen/err_wcslib) w.r.t. the exact inverse of its'
             ' own formula, mA. Wavelengths in Angstrom; lambda_vac grid step 100 A.\n')
    fh.write(f'# ISA pressure at 2600 m used for site column: {p2600:.1f} Pa (derived from barometric formula).\n')
    w = csv.writer(fh)
    w.writerow(names)
    for i in range(len(wl)):
        out = []
        for k in names:
            v = cols[k][i]
            out.append(f'{v:.1f}' if k == 'wl_vac_A' else (f'{v:.6f}' if k.endswith('_A') else f'{v:.5f}'))
        w.writerow(out)
print('wrote', OUT, len(wl), 'rows')

# --- summary at key wavelengths
key = [2000.0, 3000.0, 4000.0, 5000.0, 6564.6, 8000.0, 10000.0, 15000.0, 20000.0, 25000.0]
kk = np.array(key)
print('\nlambda_vac  vac-air(A)  ' + '  '.join(n for n, _ in std))
for x in kk:
    ra = F.vac_to_air_ciddor1996(x)
    print(f'{x:8.1f} {x - ra:9.5f} ' + ' '.join(f'{(fn(np.array([x]))[0] - ra) * 1e3:+8.4f}' for _, fn in std))
print('\nmax |d| over 2000-25000 A (mA):')
for name, _ in std:
    v = cols[f'dair_{name}_mA']
    print(f'  {name:42s} min={v.min():+9.4f} max={v.max():+9.4f}')
for k in names:
    if k.startswith('err_') or k.startswith('dair_lab') or k.startswith('dair_site'):
        v = cols[k]
        print(f'  {k:42s} min={v.min():+10.5f} max={v.max():+10.5f}  @2000={v[0]:+.5f} @5000={v[wl == 5000.0][0]:+.5f}')

# --- Ciddor sensitivities at 5000 A (vacuum), central differences
lam = 5000.0
conds = {'standard_air_15C_101325Pa_RH0_450ppm': dict(t_C=15.0, p_Pa=101325.0, rh=0.0, xc=450.0),
         'lab_20C_101325Pa_RH50_420ppm': dict(t_C=20.0, p_Pa=101325.0, rh=0.5, xc=420.0),
         'site_2600m_ISA_0C_RH30_420ppm': dict(t_C=0.0, p_Pa=p2600, rh=0.3, xc=420.0),
         'site_4200m_ISA_0C_RH20_420ppm': dict(t_C=0.0, p_Pa=float(F.isa_pressure_Pa(4200.0)), rh=0.2, xc=420.0)}
srows = []
for cname, c in conds.items():
    n0 = F.n_ciddor1996(lam, **c)
    la0 = lam / n0

    def lair(**kw):
        cc = dict(c); cc.update(kw)
        return lam / F.n_ciddor1996(lam, **cc)
    dT = (lair(t_C=c['t_C'] + 0.5) - lair(t_C=c['t_C'] - 0.5))           # per K
    dP = (lair(p_Pa=c['p_Pa'] + 50.0) - lair(p_Pa=c['p_Pa'] - 50.0))     # per hPa
    rh_lo = max(c['rh'] - 0.005, 0.0)
    dRH = (lair(rh=rh_lo + 0.01) - lair(rh=rh_lo))                       # per 1 % RH
    dC = (lair(xc=c['xc'] + 50.0) - lair(xc=c['xc'] - 50.0))             # per 100 ppm
    dla_std = la0 - F.vac_to_air_ciddor1996(lam)
    row = dict(condition=cname, t_C=c['t_C'], p_Pa=round(c['p_Pa'], 1), RH_pct=c['rh'] * 100, xc_ppm=c['xc'],
               nm1_x1e8=round((n0 - 1) * 1e8, 3), lam_air_A=round(la0, 6),
               lam_air_minus_standard_air_mA=round(dla_std * 1e3, 3),
               dlam_dT_mA_per_K=round(dT * 1e3, 4), dlam_dP_mA_per_hPa=round(dP * 1e3, 4),
               dlam_dRH_mA_per_pctRH=round(dRH * 1e3, 5), dlam_dCO2_mA_per_100ppm=round(dC * 1e3, 5),
               dv_dT_ms_per_K=round(dT / lam * C_KMS * 1e3, 2), dv_dP_ms_per_hPa=round(dP / lam * C_KMS * 1e3, 2),
               dv_dRH_ms_per_pctRH=round(dRH / lam * C_KMS * 1e3, 3), dv_dCO2_ms_per_100ppm=round(dC / lam * C_KMS * 1e3, 3))
    srows.append(row)
    print('\n', row)
with open(OUT_SENS, 'w', newline='') as fh:
    fh.write('# Ciddor (1996) full-model sensitivities of the air wavelength at lambda_vac = 5000 A. DERIVED '
             '(central differences: dT=1 K, dp=100 Pa, dRH=1 %, dCO2=100 ppm). Velocity = c*dlambda/lambda.\n')
    w = csv.DictWriter(fh, fieldnames=list(srows[0].keys()))
    w.writeheader()
    w.writerows(srows)
print('wrote', OUT_SENS)
