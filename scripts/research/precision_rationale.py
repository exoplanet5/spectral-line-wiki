#!/usr/bin/env python
"""
precision_rationale.py -- the numbers behind data/notes/precision_rationale.md:
WHY the Spectral Line Wiki quotes wavelengths to 1e-3 Angstrom (1 mA), where that is not enough,
and where it is overkill.

Usage (from anywhere; ROOT = project root):
    ~/.venvs/astro313/bin/python -I scripts/research/precision_rationale.py ROOT

Reads (read-only):
    ROOT/data/handbook.json                                  NIST Handbook lines (decimal counts)
    ROOT/data/handbook_asd.json                              Handbook <-> NIST ASD cross-match (snapshot)
    ROOT/data/notes/murphy_berengut_2014_mm_transitions.json q coefficients (Murphy & Berengut 2014)
    ROOT/data/raw/precision_rationale/nist_asd/*.csv         NIST ASD windows (downloaded; parsed as text only)
    scripts/research/airvac_formulas.py                      project implementation of the air/vacuum formulas
Writes:
    ROOT/data/curated/spectrograph_precision.csv
    ROOT/data/notes/precision_rationale_numbers.json
and prints the markdown tables used in the report.

Provenance flags used in the output:
    lit      value typed verbatim from the cited source (string kept as printed)
    derived  computed here from the stated formula / constants
    assumed  an explicit assumption (not from a source); always says so
Constants: astropy.constants (CODATA 2018; IAU 2015 nominal solar values; IAU 2012 au).
"""
import csv, importlib.util, io, json, math, os, re, sys
from collections import Counter

import numpy as np
from scipy.special import erf
from astropy import constants as const
from astropy import units as u
from astropy.cosmology import Planck18

ROOT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), '..', '..'))
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('airvac_formulas', os.path.join(HERE, 'airvac_formulas.py'))
AV = importlib.util.module_from_spec(spec)
spec.loader.exec_module(AV)

C = const.c.value                      # 299792458 m/s exactly
H_PLANCK = const.h.value
E_CHARGE = const.e.si.value
M_E = const.m_e.value
GM_SUN = const.GM_sun.value            # IAU 2015 B3 nominal 1.3271244e20 m^3 s^-2
R_SUN = const.R_sun.value              # IAU 2015 B3 nominal 6.957e8 m
AU = const.au.value                    # IAU 2012 B2
OMEGA_E = 7.292115e-5                  # rad/s, IERS Conventions 2010 (TN36) Table 1.2
A_E = 6378137.0                        # m, GRS80/WGS84 equatorial radius

OUT = {'generated_by': 'scripts/research/precision_rationale.py', 'constants': {
    'c_m_s': C, 'h_J_s': H_PLANCK, 'e_C': E_CHARGE, 'm_e_kg': M_E, 'GM_sun': GM_SUN, 'R_sun_m': R_SUN,
    'au_m': AU, 'omega_E_rad_s': OMEGA_E, 'a_E_m': A_E}}


def mA_of_v(v_ms, lam_A):
    """Doppler-equivalent wavelength step in mA: dlam = lam * v / c (first order)."""
    return lam_A * v_ms / C * 1e3


def v_of_mA(dl_mA, lam_A):
    return C * dl_mA * 1e-3 / lam_A


def md_table(header, rows):
    s = '| ' + ' | '.join(header) + ' |\n|' + '---|' * len(header) + '\n'
    for r in rows:
        s += '| ' + ' | '.join(str(x) for x in r) + ' |\n'
    return s


def g(x, nd=3):
    """Compact formatting with nd significant digits."""
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return ''
    if x == 0:
        return '0'
    e = math.floor(math.log10(abs(x)))
    if -3 <= e < 5:
        return f'{x:.{max(nd - 1 - e, 0)}f}'
    return f'{x:.{nd - 1}e}'


PRINT = []


def show(title, text):
    PRINT.append(f'\n### {title}\n\n{text}')


# =====================================================================================
# (1) Velocity, frequency and energy equivalents of 1 mA
# =====================================================================================
LAMS = [1200, 3000, 4000, 5000, 6563, 8000, 10000, 16000, 22000]
rows, sec1 = [], []
for lam in LAMS:
    dv1 = v_of_mA(1.0, lam)                         # m/s per mA
    lam_m = lam * 1e-10
    dnu = C * 1e-13 / lam_m ** 2                    # Hz for 1 mA (vacuum wavelength)
    dsig = 1e-11 / (lam * 1e-8) ** 2                # cm^-1: dsigma = dlam/lam^2 with lam in cm, dlam = 1e-11 cm
    dE = H_PLANCK * dnu / E_CHARGE                  # eV
    # medium check: is c*dlam/lam the same for air and vacuum wavelengths?  lam_air = lam_vac/n(lam_vac)
    # => dlam_air/lam_air = (dlam_vac/lam_vac) * (1 - lam_vac n'(lam_vac)/n).  Evaluated with Ciddor standard air.
    if lam >= 2000:
        h_ = 0.01
        lv = np.array([lam - h_, lam + h_])
        la = AV.vac_to_air_ciddor1996(lv)
        ratio = ((la[1] - la[0]) / AV.vac_to_air_ciddor1996(lam)) / ((lv[1] - lv[0]) / lam)
        medium_factor = float(ratio)
    else:
        medium_factor = None
    rec = dict(lam_A=lam, dv_per_mA_m_s=dv1, dv_per_10mA_m_s=10 * dv1, dv_per_0p1mA_m_s=0.1 * dv1,
               dlam_for_1m_s_A=lam / C, dlam_for_10cm_s_A=0.1 * lam / C, dlam_for_1cm_s_A=0.01 * lam / C,
               dnu_per_mA_MHz=dnu / 1e6, dsigma_per_mA_cm1=dsig, dE_per_mA_eV=dE,
               airwl_over_vacwl_fractional_step_ratio_Ciddor=medium_factor, flag='derived')
    sec1.append(rec)
    rows.append([lam, f'{dv1:.2f}', f'{10 * dv1:.1f}', f'{lam / C:.3e}', f'{0.1 * lam / C:.2e}',
                 f'{dnu / 1e6:.2f}', f'{dsig:.3e}', f'{dE:.3e}',
                 '' if medium_factor is None else f'{medium_factor - 1:+.1e}'])
OUT['sec1_velocity_equivalents'] = sec1
show('1. Velocity / frequency / energy equivalent of 1 mA (derived, exact c)',
     md_table(['lambda (A)', '1 mA = m/s', '10 mA = m/s', '1 m/s = A', '10 cm/s = A', '1 mA = MHz',
               '1 mA = cm^-1', '1 mA = eV', '(dlam/lam)_air/(dlam/lam)_vac - 1'], rows))

# =====================================================================================
# (2) Spectrographs: resolution element, sampling, centroid precision
# =====================================================================================
SQ2PI = math.sqrt(2 * math.pi)
FW2S = 2 * math.sqrt(2 * math.log(2))  # FWHM / sigma


def _Phi(z):
    return 0.5 * (1 + erf(z / math.sqrt(2)))


def _phi(z):
    return np.exp(-0.5 * z * z) / SQ2PI


def crb_centroid(fwhm, npix_per_fwhm, snr=100.0, kind='absorption', depth=0.5, nphase=64):
    """Photon-noise Cramer-Rao bound on the centroid (same units as fwhm) of an unresolved Gaussian line
    (intrinsic width << instrumental, so the observed profile is the Gaussian LSF of FWHM 'fwhm'),
    integrated over pixels of width p = fwhm/npix_per_fwhm, averaged over sub-pixel phase.
    This is the single-line case of Bouchy, Pepe & Queloz (2001) eqs. 8-12 with sigma_D = 0:
        dlam_RMS = 1/sqrt(sum_i (dA0_i/dlam)^2 / A0_i).
    kind='absorption': continuum level SNR^2 photo-electrons per pixel, line depth 'depth'.
    kind='emission'  : no continuum/background; peak (continuous profile, per pixel width) = SNR^2 e-."""
    s = fwhm / FW2S
    p = fwhm / npix_per_fwhm
    Cc = snr ** 2
    halfw = 8 * s + 2 * p
    infos = []
    for ph in np.arange(nphase) / nphase:
        edges = np.arange(-halfw - ph * p, halfw + p, p)
        a, b = edges[:-1] / s, edges[1:] / s           # (x - x0)/s at the pixel edges, x0 = 0
        prof = (s * SQ2PI / p) * (_Phi(b) - _Phi(a))   # pixel-averaged exp(-x^2/2s^2)
        dprof = (SQ2PI / p) * (_phi(a) - _phi(b))      # d prof / d x0
        if kind == 'absorption':
            F = Cc * (1 - depth * prof)
            dF = -Cc * depth * dprof
        else:
            F = Cc * prof
            dF = Cc * dprof
        ok = F > 0
        infos.append(np.sum(dF[ok] ** 2 / F[ok]))
    return 1 / math.sqrt(np.mean(infos))


def rule_centroid(fwhm, npix_per_fwhm, snr=100.0):
    """Rule of thumb sigma_lam ~ FWHM / (SNR sqrt(N_pix)), N_pix = pixels per FWHM."""
    return fwhm / (snr * math.sqrt(npix_per_fwhm))


# Analytic limits (continuous sampling) for checking the numerics:
#   emission:   sigma = s / sqrt(N_tot), N_tot = SNR^2 sqrt(2 pi) s/p  ->  sigma = 0.4116 FWHM/(SNR sqrt(Npix))
#   absorption (D -> 0): sigma = sqrt(2 s p / sqrt(pi)) / (SNR D)
_chk_e = crb_centroid(1.0, 20.0, 100.0, 'emission')
_chk_e_an = (1 / FW2S) / (100 * math.sqrt(SQ2PI * (20.0 / FW2S)))
_chk_a = crb_centroid(1.0, 20.0, 100.0, 'absorption', depth=0.01)
_chk_a_an = math.sqrt(2 * (1 / FW2S) * (1 / 20.0) / math.sqrt(math.pi)) / (100 * 0.01)
OUT['crb_numeric_checks'] = dict(emission_numeric=_chk_e, emission_analytic=_chk_e_an,
                                 absorption_D0p01_numeric=_chk_a, absorption_D0p01_analytic=_chk_a_an,
                                 rule_over_crb_emission_continuous=1 / (1 / (FW2S * math.sqrt(SQ2PI / FW2S))))

# Instrument parameters.  'src_*' strings are verbatim from the cited page/paper; 'derived' entries say how.
REF = {
    'harps': ('Mayor M. et al. 2003, The Messenger 114, 20 "Setting New Standards with HARPS" (2003Msngr.114...20M); '
              'ESO HARPS overview page', 'https://www.eso.org/sci/publications/messenger/archive/no.114-dec03/messenger-no114-20-24.pdf ; '
              'https://www.eso.org/sci/facilities/lasilla/instruments/harps/overview.html'),
    'espresso': ('Pepe F. et al. 2021, A&A 645, A96 (doi:10.1051/0004-6361/202038306), Table 3; ESO ESPRESSO overview page',
                 'https://arxiv.org/abs/2010.00316 ; https://www.eso.org/sci/facilities/paranal/instruments/espresso/overview.html'),
    'expres': ('Blackman R. T. et al. 2020, AJ 159, 238 (doi:10.3847/1538-3881/ab811d); Jurgenson C. et al. 2016, Proc. SPIE 9908, 99086T '
               '(doi:10.1117/12.2233002)', 'https://arxiv.org/abs/2003.08852 ; https://arxiv.org/abs/1606.04413'),
    'neid': ('WIYN NEID characteristics page; NEID Observers and Proposers FAQ (Penn State); Schwab C. et al. 2016, Proc. SPIE 9908, '
             '99087H (doi:10.1117/12.2234411)', 'https://www.wiyn.org/Instruments/wiynneid.html ; https://neid.psu.edu/observers-and-proposers-faq'),
    'pepsi': ('Strassmeier K. G. et al. 2015, Astron. Nachr. 336, 324 (doi:10.1002/asna.201512172), Table 1',
              'https://arxiv.org/abs/1505.06492'),
    'uves': ('ESO UVES instrument description and overview pages; Dekker H. et al. 2000, Proc. SPIE 4008, 534 (doi:10.1117/12.395512)',
             'https://www.eso.org/sci/facilities/paranal/instruments/uves/inst.html ; https://www.eso.org/sci/facilities/paranal/instruments/uves/overview.html'),
    'hires': ('Keck HIRES "Available Slits and Resolutions" and "Binning Options" pages; Vogt S. S. et al. 1994, Proc. SPIE 2198, 362 '
              '(doi:10.1117/12.176725)', 'https://www2.keck.hawaii.edu/inst/hires/slitres.html ; https://www2.keck.hawaii.edu/inst/hires/binning.html'),
    'hds': ('Subaru HDS instrument pages; Noguchi K. et al. 2002, PASJ 54, 855 (doi:10.1093/pasj/54.6.855)',
            'https://subarutelescope.org/Instruments/HDS/ ; https://subarutelescope.org/en/about/instrument/hds/index.html'),
    'crires': ('Dorn R. J. et al. 2023, A&A 671, A24 (doi:10.1051/0004-6361/202245217), Sect. 4.9; ESO CRIRES overview page',
               'https://arxiv.org/abs/2301.08048 ; https://www.eso.org/sci/facilities/paranal/instruments/crires/overview.html'),
    'igrins': ('Gemini IGRINS page; McDonald Obs. IGRINS page (Park C. et al. 2014, Proc. SPIE 9147, 91471D, doi:10.1117/12.2056431); '
               'Le H. A. N. et al. 2015, Adv. Space Res. 55, 2509 (doi:10.1016/j.asr.2015.03.007)',
               'https://www.gemini.edu/instrumentation/igrins ; https://arxiv.org/abs/1501.03249'),
    'apogee': ('Wilson J. C. et al. 2019, PASP 131, 055001 (doi:10.1088/1538-3873/ab0075), Table 1', 'https://arxiv.org/abs/1902.00928'),
    'gaia': ('Cropper M. et al. 2018, A&A 616, A5 (doi:10.1051/0004-6361/201832763), Sect. 4.2 & Table 1', 'https://arxiv.org/abs/1804.09369'),
    'desi': ('DESI Collaboration (Abareshi B. et al.) 2022, AJ 164, 207 (doi:10.3847/1538-3881/ac882b), Sect. 7', 'https://arxiv.org/abs/2205.10939'),
    'boss': ('Smee S. A. et al. 2013, AJ 146, 32 (doi:10.1088/0004-6256/146/2/32), abstract, Sect. 2.1.5, Table 3', 'https://arxiv.org/abs/1208.2233'),
    'lamost': ('LAMOST DR11 v1.1 LRS data-production description; Cui X.-Q. et al. 2012, RAA 12, 1197 (doi:10.1088/1674-4527/12/9/003)',
               'https://www.lamost.org/dr11/v1.1/doc/lr-data-production-description'),
    'lamost_mrs': ('Liu C. et al. 2020, arXiv:2005.07210 (LAMOST-MRS survey plan)', 'https://arxiv.org/abs/2005.07210'),
    '4most': ('de Jong R. S. et al. 2019, The Messenger 175, 3 (doi:10.18727/0722-6691/5117), Table 1 and text', 'https://arxiv.org/abs/1903.02464'),
    'pfs': ('Subaru PFS "Instrument parameters" page (as of July 2024, page updated 2026-10-21); Tamura N. et al. 2016, Proc. SPIE 9908, '
            '99081M (doi:10.1117/12.2232103)', 'https://pfs.ipmu.jp/research/parameters.html'),
    'weave': ('Jin S. et al. 2024, MNRAS 530, 2688 (doi:10.1093/mnras/stad557), Table 3', 'https://arxiv.org/abs/2212.03981'),
}

# fields: name, R, lam_ref_A, sampling_px (px per FWHM), fwhm_A (if the source quotes it), source strings, notes, ref key
INSTR = [
    dict(name='HARPS (ESO 3.6 m)', R=115000, lam=5500, samp=3.2, ref='harps',
         src='R=115,000; "each spectral element is still sampled by 3.2 CCD pixels"; 378-691 nm; 15 um pixels',
         notes='Lovis & Pepe (2007) measured R=110,000, 0.015 A/pixel, ~3.3 pixels per resolution element.'),
    dict(name='ESPRESSO HR 1-UT', R=140000, lam=5500, samp=4.5, ref='espresso',
         src='HR 1-UT: median R 140,000 (ESO page; 138,000 in Pepe+2021 Table 3), spectral sampling (average) 4.5 px; 380-788 nm', notes=''),
    dict(name='ESPRESSO UHR 1-UT', R=190000, lam=5500, samp=2.5, ref='espresso',
         src='UHR 1-UT: median R 190,000 (">190 000" in Pepe+2021), spectral sampling 2.5 px', notes=''),
    dict(name='EXPRES (LDT)', R=137500, lam=5500, samp=4.0, ref='expres',
         src='median R = 137,500 +- 6,100 (measured with LFC), "sampling of four pixels"; 3800-7800 A', notes='design R=150,000 (Jurgenson+2016)'),
    dict(name='NEID HR (WIYN)', R=110000, lam=5500, samp=5.0, ref='neid',
         src='HR mode R>110,000 (WIYN page); "spectral resolution element in HR mode is 5 pixels wide" (NEID FAQ); 380-930 nm', notes='R taken as the quoted lower bound 110,000'),
    dict(name='PEPSI 270k (LBT)', R=270000, lam=5500, samp=2.0, ref='pepsi',
         src='100 um fibre, 7 slices: R=270 000, sampling 2.0 pix; 383-907 nm', notes=''),
    dict(name='PEPSI 120k (LBT)', R=120000, lam=5500, samp=4.2, ref='pepsi',
         src='200 um fibre, 5 slices: R=120 000, sampling 4.2 pix', notes=''),
    dict(name='PEPSI 43k (LBT)', R=43000, lam=5500, samp=12.1, ref='pepsi',
         src='300 um fibre, 3 slices: R=43 000, sampling 12.1 pix', notes=''),
    dict(name='UVES red, 1" slit (VLT)', R=38700, lam=6000, samp=1.0 / 0.16, ref='uves',
         src='red-arm resolution-slit product 38,700 (R~40,000 for 1" slit); red CCD 15 um pixels, 0.16 arcsec/pix',
         notes='sampling derived = slit width / pixel scale = 1.0/0.16 = 6.25 px (assumes equal spectral and spatial scale)'),
    dict(name='UVES red, 0.3" slit (VLT)', R=110000, lam=6000, samp=2.0, ref='uves',
         src='"maximum (two-pixel) resolution is ... 110,000 in the Red Arm" (~110,000 with 0.3" slit)', notes=''),
    dict(name='Keck HIRES 0.861" (C2/B5)', R=50000, lam=5240, samp=4.8, ref='hires',
         src='0.861" deckers: R calculated 48,000, measured 49,000 (UV XD, 4100 A), 50,000 (red XD, 5240 A); binning page: 0.86" slit projects to 2.4 super-pixels when binned x2',
         notes='sampling derived = 2 x 2.4 = 4.8 unbinned 15 um pixels'),
    dict(name='Keck HIRES 0.574" (B1-B4)', R=66400, lam=5240, samp=4.8 * 0.574 / 0.861, ref='hires',
         src='0.574" deckers: R measured 66,400 (red XD, 5240 A)', notes='sampling derived by scaling 4.8 px x 0.574/0.861 = 3.2 px'),
    dict(name='Subaru HDS 0.2" slit', R=160000, lam=5500, samp=2.0, ref='hds',
         src='"maximum (2-pixel FWHM) resolution of 160 000 (1.8 km/s) ... with a 0.2\'\' slit or Image Slicer#3"; 13.5 um pixels, 0.138"/pixel', notes=''),
    dict(name='CRIRES+ 0.2" slit (VLT)', R=104800, lam=15820, samp=0.2 / 0.056, ref='crires',
         src='H1582 setting, 0.2" slit: average R = 104,800 +- 2600 (minimum 100,000); science detector 0.056 arcsec/pixel',
         notes='sampling derived = 0.2"/0.056" = 3.6 px (slit projection; AO-sharpened PSF can be narrower)'),
    dict(name='IGRINS (H band)', R=45000, lam=16000, samp=3.66, ref='igrins',
         src='R~45,000 (Gemini page; McDonald page); "pslit is 3.66 pixels in the spatial and spectral resolutions" (Le+2015, who adopt R=40,000)',
         notes='sampling from the ETC paper'),
    dict(name='APOGEE (green chip)', R=23700, lam=16100, samp=None, disp=0.282, ref='apogee',
         src='median native (lambda/FWHM) resolution 22,800; 23,700; 22,000 (blue; green; red); dispersion 0.326; 0.282; 0.235 A/pixel at 1.54; 1.61; 1.66 um; requirement R>=22,500',
         notes='sampling derived = (16100/23700)/0.282 = 2.4 px; blue end undersampled, remedied by 0.5-pixel dithering'),
    dict(name='Gaia RVS (HR mode)', R=11500, lam=8600, samp=None, disp=0.245, ref='gaia',
         src='R = 11500 adopted (average 10 500-12 500 required); dispersion linear term ~0.0245 nm/pixel; 847-874 nm',
         notes='sampling derived = (8600/11500)/0.245 = 3.05 px'),
    dict(name='DESI blue arm', R=2500, lam=4750, samp=None, disp=(5930 - 3600) / 4096, ref='desi',
         src='blue 360-593 nm, R ~ 2000-3000; 4096x4096 15 um CCD; fibre (107 um) imaged with magnification 0.48',
         notes='R = mid-point of quoted 2000-3000 (assumed); dispersion derived = (5930-3600 A)/4096 px; fibre image 0.48x107/15 = 3.4 px'),
    dict(name='DESI red arm', R=4000, lam=6660, samp=None, disp=(7720 - 5600) / 4114, ref='desi',
         src='red 560-772 nm, R ~ 3500-4500; 4114x4128 15 um CCD', notes='R = mid-point (assumed); dispersion derived = (7720-5600)/4114'),
    dict(name='DESI NIR arm', R=4750, lam=8635, samp=None, disp=(9800 - 7470) / 4114, ref='desi',
         src='NIR 747-980 nm, R ~ 4000-5500', notes='R = mid-point (assumed); dispersion derived = (9800-7470)/4114'),
    dict(name='SDSS/BOSS', R=2000, lam=6000, samp=3.0, ref='boss',
         src='"resolving power R = lambda/FWHM ~ 2000"; "3 pixels per resolution element on the CCD"; measured BOSS R 1550-2650 for lambda>4900 A (Table 3)',
         notes='pipeline output grid 1e-4 dex (69 km/s) per pixel'),
    dict(name='LAMOST LRS', R=1800, lam=5500, samp=None, disp=5500 * math.log(10) * 1e-4, ref='lamost',
         src='"resolution of 1800 at the 5500 A"; slit mode x2/3; COEFF1 = 0.0001 (log10 dispersion per pixel)',
         notes='sampling derived on the 1e-4 dex output grid (1.266 A/px at 5500 A) -> 2.4 output px; native CCD sampling not given in source'),
    dict(name='LAMOST MRS', R=7500, lam=5150, samp=3.0, ref='lamost_mrs',
         src='R ~ 7500; blue 4950-5350 A, red 6300-6800 A; RV precision ~1 km/s for late-type stars',
         notes='sampling NOT in cited source: 3 px ASSUMED for the centroid estimate'),
    dict(name='4MOST HRS', R=20000, lam=5450, samp=3.0, ref='4most',
         src='<R> = 20 000; passbands 3926-4355, 5160-5730, 6100-679x A; "sampled with about three pixels per resolution element"; velocity accuracy < 1 km/s', notes=''),
    dict(name='4MOST LRS', R=6500, lam=6600, samp=3.0, ref='4most',
         src='<R> = 6500; 3700-9500 A; about three pixels per resolution element', notes=''),
    dict(name='Subaru PFS blue', R=2500, lam=5000, samp=None, fwhm=2.1, disp=0.7, ref='pfs',
         src='380-650 nm; dispersion ~0.7 A/pix; spectral resolution ~2.1 A; R ~2500 (@500nm)', notes='sampling derived = 2.1/0.7 = 3.0 px'),
    dict(name='Subaru PFS red LR', R=3000, lam=8000, samp=None, fwhm=2.7, disp=0.9, ref='pfs',
         src='630-970 nm; ~0.9 A/pix; ~2.7 A; R ~3000 (@800nm)', notes='sampling derived = 3.0 px'),
    dict(name='Subaru PFS red MR', R=5500, lam=8000, samp=None, fwhm=1.6, disp=0.4, ref='pfs',
         src='710-885 nm; ~0.4 A/pix; ~1.6 A; R ~5500 (@800nm)', notes='sampling derived = 4.0 px'),
    dict(name='Subaru PFS NIR', R=4500, lam=11000, samp=None, fwhm=2.4, disp=0.8, ref='pfs',
         src='940-1260 nm; ~0.8 A/pix; ~2.4 A; R ~4500 (@1100nm)', notes='sampling derived = 3.0 px'),
    dict(name='WEAVE HR', R=20000, lam=5000, samp=3.0, ref='weave',
         src='HR modes R = 20 000 (Blue-1 4040-4650, Blue-2 4730-5450, Red 5950-6850 A); 15 um pixels',
         notes='sampling NOT in cited source: 3 px ASSUMED'),
    dict(name='WEAVE LR', R=5000, lam=6600, samp=3.0, ref='weave',
         src='LR R = 5 000 (3660-9590 A)', notes='sampling NOT in cited source: 3 px ASSUMED'),
]

SNR0 = 100.0
inst_rows, csv_rows = [], []
for it in INSTR:
    lam, R = it['lam'], it['R']
    fwhm = it.get('fwhm') or lam / R
    fwhm_flag = 'lit (quoted FWHM)' if it.get('fwhm') else 'derived lambda/R'
    samp = it['samp'] if it['samp'] is not None else fwhm / it['disp']
    pix = fwhm / samp
    rot = rule_centroid(fwhm, samp, SNR0)
    crb_abs = crb_centroid(fwhm, samp, SNR0, 'absorption', 0.5)
    crb_em = crb_centroid(fwhm, samp, SNR0, 'emission')
    fwhm_kms = C * fwhm / lam / 1e3
    rec = dict(instrument=it['name'], R=R, ref_wavelength_A=lam, FWHM_A=round(fwhm, 5), FWHM_flag=fwhm_flag,
               FWHM_kms=round(fwhm_kms, 3), pixel_sampling_px_per_FWHM=round(samp, 2), pixel_size_A=round(pix, 5),
               pixel_size_kms=round(fwhm_kms / samp, 3),
               centroid_mA_SNR100_rule=round(rot * 1e3, 3), centroid_ms_SNR100_rule=round(v_of_mA(rot * 1e3, lam), 1),
               centroid_mA_SNR100_CRB_abs_depth0p5=round(crb_abs * 1e3, 3),
               centroid_mA_SNR100_CRB_emission=round(crb_em * 1e3, 3),
               source_values=it['src'], derivation_notes=it['notes'],
               citation=REF[it['ref']][0], url=REF[it['ref']][1])
    csv_rows.append(rec)
    inst_rows.append([it['name'], f'{R:,}', lam, g(fwhm, 3), g(fwhm_kms, 3), g(samp, 2), g(pix * 1e3, 3),
                      g(rot * 1e3, 2), g(v_of_mA(rot * 1e3, lam), 2), g(crb_abs * 1e3, 2), g(crb_em * 1e3, 2)])

os.makedirs(os.path.join(ROOT, 'data', 'curated'), exist_ok=True)
csv_path = os.path.join(ROOT, 'data', 'curated', 'spectrograph_precision.csv')
cols = ['instrument', 'R', 'ref_wavelength_A', 'FWHM_A', 'FWHM_flag', 'FWHM_kms', 'pixel_sampling_px_per_FWHM', 'pixel_size_A',
        'pixel_size_kms', 'centroid_mA_SNR100_rule', 'centroid_ms_SNR100_rule', 'centroid_mA_SNR100_CRB_abs_depth0p5',
        'centroid_mA_SNR100_CRB_emission', 'source_values', 'derivation_notes', 'citation', 'url']
with open(csv_path, 'w', newline='') as f:
    f.write('# Spectrograph resolution elements and single-line centroid precision at SNR = 100 per pixel. Built by '
            'scripts/research/precision_rationale.py.\n')
    f.write('# R, sampling, dispersion and quoted FWHM are from the cited source (source_values, verbatim where quoted); '
            'FWHM_A = ref_wavelength/R unless FWHM_flag says it is quoted. Vacuum/air immaterial at this level (3e-4).\n')
    f.write('# centroid_mA_SNR100_rule = FWHM/(SNR*sqrt(N_pix)), N_pix = pixels per FWHM (rule of thumb, task definition). '
            'CRB columns = photon-noise Cramer-Rao bound for one unresolved Gaussian line (Bouchy, Pepe & Queloz 2001, '
            'A&A 374, 733, eqs. 8-12 with sigma_D=0): absorption depth 0.5 in a continuum of SNR 100/pixel, or emission '
            'with peak 1e4 e-/pixel and no background; detector noise, blends, calibration and LSF errors excluded.\n')
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    for r in csv_rows:
        w.writerow(r)
OUT['sec2_spectrographs'] = csv_rows
show('2a. Spectrographs at SNR = 100 per pixel (FWHM from R unless quoted; CRB = Cramer-Rao bound)',
     md_table(['instrument', 'R', 'lambda_ref (A)', 'FWHM (A)', 'FWHM (km/s)', 'px/FWHM', 'pixel (mA)',
               'rule FWHM/(SNR sqrt N) (mA)', 'rule (m/s)', 'CRB abs. D=0.5 (mA)', 'CRB emission (mA)'], inst_rows))

# generic grid: centroid precision vs R and SNR at 5000 A, 3 px sampling
grid_R = [1000, 2000, 5000, 10000, 20000, 50000, 100000, 140000, 200000]
grid_S = [10, 30, 100, 300, 1000]
grid, grow = {}, []
for R in grid_R:
    fw = 5000.0 / R
    row = [f'{R:,}', g(fw, 3)]
    for S in grid_S:
        r_ = rule_centroid(fw, 3.0, S) * 1e3
        a_ = crb_centroid(fw, 3.0, S, 'absorption', 0.5) * 1e3
        grid[f'R{R}_SNR{S}'] = dict(rule_mA=r_, crb_abs_mA=a_)
        row.append(f'{g(r_, 2)} / {g(a_, 2)}')
    grow.append(row)
OUT['sec2_grid_5000A_3px'] = grid
show('2b. Single-line centroid precision at 5000 A, 3 px/FWHM: rule-of-thumb / CRB (absorption D=0.5), in mA',
     md_table(['R', 'FWHM (A)'] + [f'SNR {S}' for S in grid_S], grow))
# Where is 1 mA reached? minimal R for rule <= 1 mA at SNR 100 and 3 px:  FWHM <= 0.1*sqrt(3) A
OUT['sec2_threshold'] = dict(fwhm_for_1mA_rule_SNR100_3px_A=1e-3 * 100 * math.sqrt(3),
                             R_for_1mA_rule_SNR100_3px_at5000=5000 / (1e-3 * 100 * math.sqrt(3)))

# =====================================================================================
# (3) Physical / instrumental effects at the mA level (5000 A unless stated)
# =====================================================================================
L0 = 5000.0
eff = []


def add(effect, v_ms=None, dl_mA=None, lam=L0, flag='derived', source='', note=''):
    if dl_mA is None:
        dl_mA = mA_of_v(v_ms, lam)
    if v_ms is None:
        v_ms = v_of_mA(dl_mA, lam)
    eff.append(dict(effect=effect, lam_A=lam, v_m_s=v_ms, dlam_mA=dl_mA, flag=flag, source=source, note=note))


v_rot = OMEGA_E * A_E
add('Earth rotation, equator (omega_E a_E)', v_rot, source='IERS TN36 omega_E; GRS80 a_E')
v_orb = math.sqrt(GM_SUN / AU)
add('Earth orbital speed, circular sqrt(GM_sun/au)', v_orb, source='IAU nominal GM_sun, au',
    note='DE440s 2026: 29.294-30.278 km/s (velocity_frames.md)')
add('Earth orbital speed, DE440s minimum 2026', 29294.0, flag='derived (DE440s)', source='data/notes/velocity_frames.md')
add('Earth orbital speed, DE440s maximum 2026', 30278.0, flag='derived (DE440s)', source='data/notes/velocity_frames.md')
v_gr = GM_SUN / (R_SUN * C)
v_gr1 = GM_SUN / C * (1 / R_SUN - 1 / AU)
add('Solar gravitational redshift, surface to infinity GM/(Rc)', v_gr, source='IAU 2015 B3 nominal GM_sun, R_sun')
add('Solar gravitational redshift seen at 1 au GM/c (1/R - 1/au)', v_gr1, source='idem; Gonzalez Hernandez et al. 2020 quote 633.35 m/s')
add('Solar GRS measured with HARPS-LFC (15 Fe I lines, 150<EW<550 mA)', 639.0, flag='lit', source='Gonzalez Hernandez et al. 2020, A&A 643, A146',
    note='639 +- 14 m/s')
# Reiners et al. 2016, A&A 587, A65, Table 5: absolute convective blueshift vs relative line depth (Fe I, IAG atlas)
REINERS_T5 = [('0.10', '-547', '303'), ('0.20', '-513', '199'), ('0.30', '-515', '232'), ('0.40', '-491', '155'),
              ('0.50', '-450', '284'), ('0.60', '-419', '224'), ('0.70', '-334', '215'), ('0.80', '-207', '193'),
              ('0.90', '2', '128')]
for d, v, sg in REINERS_T5:
    add(f'Solar convective blueshift, Fe I lines of relative depth {d} (mean)', float(v), flag='lit',
        source='Reiners et al. 2016, A&A 587, A65, Table 5', note=f'sigma {sg} m/s')
add('Convective blueshift, Fe I 5250.2084 A bisector maximum, disk centre (LARS)', -290.0, lam=5250.2084, flag='lit',
    source='Lohner-Bottcher et al. 2019, A&A 624, A57, Sect. 4', note='-100 m/s near continuum to -290 m/s at I=0.76')
add('Sun about the SSB (max, 1990-2040)', 16.08, flag='derived (DE440s)', source='data/notes/velocity_frames.md')
add('Earth about the Earth-Moon barycentre (mean)', 12.43, flag='derived (DE440s)', source='data/notes/velocity_frames.md')
add('Earth-induced reflex of the Sun', 0.09, flag='lit', source='Fischer et al. 2016, PASP 128, 066001 ("9 cm/s")')
add('HARPS ThAr: strongest lines, photon noise 0.002 px', 1.6, lam=5500, flag='lit', source='Lovis & Pepe 2007, A&A 468, 1115')
add('HARPS ThAr: CCD systematics floor 0.01 px', 8.0, lam=5500, flag='lit', source='Lovis & Pepe 2007')
add('HARPS ThAr: typical line-list uncertainty', 10.0, lam=5500, flag='lit', source='Lovis & Pepe 2007 ("~10 m/s ... 0.18 mA")')
add('HARPS ThAr: faintest lines 0.2 px', 160.0, lam=5500, flag='lit', source='Lovis & Pepe 2007')
add('Ar lines in ThAr lamps: lamp-to-lamp (pressure) shifts, "a few tens of m/s" (10-50 m/s shown)', 10.0, lam=5500, flag='lit',
    source='Lovis & Pepe 2007, Sect. 3 (citing Whaling et al. 2002)')
add('  (upper end of "a few tens")', 50.0, lam=5500, flag='lit', source='Lovis & Pepe 2007')
add('Palmer & Engleman 1983 Th atlas quoted accuracy 0.001-0.005 cm^-1', 16.0, lam=5500, flag='lit',
    source='Lovis & Pepe 2007 ("16-82 m/s ... 0.29-1.50 mA")')
add('  (faint end)', 82.0, lam=5500, flag='lit', source='Lovis & Pepe 2007')
add('IAG solar flux atlas absolute accuracy', 10.0, flag='lit', source='Reiners et al. 2016 ("+-10 m/s in 405-1065 nm")')
add('Kitt Peak solar atlases: wavelength-scale deviations up to', 100.0, flag='lit', source='Reiners et al. 2016 ("several ten to 100 m/s")')
OUT['sec3_effects'] = eff

# Zeeman: dlam = g e B lam^2/(4 pi m_e c) (normal Zeeman pattern, sigma-component offset)
kZ = E_CHARGE / (4 * math.pi * M_E * C)            # m^-1 T^-1
kZ_AG = kZ * 1e-10 * 1e-4                          # A^-1 G^-1 for lam in A: dlam[A] = kZ_AG * g * B[G] * lam[A]^2
zee = []
for B in [1, 100, 1000, 3000]:
    for gL, lam in [(1.0, 5000.0), (3.0, 5250.2084), (1.0, 16000.0)]:
        dl = kZ_AG * gL * B * lam ** 2 * 1e3
        zee.append(dict(B_G=B, g_eff=gL, lam_A=lam, dlam_mA=dl, v_m_s=v_of_mA(dl, lam)))
OUT['sec3_zeeman'] = dict(k_A_per_G_per_A2=kZ_AG, rows=zee,
                          note='g_eff = 3.00 for Fe I 5250.2084 A from Lohner-Bottcher et al. 2019 Table 1 (lit)')

# Air / vacuum: formula differences (fixed vacuum wavelength -> air wavelength), mA, relative to Ciddor 1996
AV_FORMS = [('Edlen 1953', AV.vac_to_air_edlen1953), ('Edlen 1966', AV.vac_to_air_edlen1966),
            ('Peck & Reeder 1972 eq.2', AV.vac_to_air_peckreeder1972),
            ('Peck & Reeder 1972 eq.3 (NIST ASD)', AV.vac_to_air_peckreeder1972_5p),
            ('Birch & Downs 1994 = Morton 2000 (VALD)', AV.vac_to_air_birchdowns1994),
            ('Ciddor 1996 (SDSS DR9+, IDL)', AV.vac_to_air_ciddor1996)]
AVL = [2000, 3000, 4000, 5000, 6563, 8000, 10000, 16000, 22000]
avrows, avrec = [], []
for lam in AVL:
    ref = float(AV.vac_to_air_ciddor1996(lam))
    vals = {nm: float(fn(lam)) for nm, fn in AV_FORMS}
    d = {nm: (v - ref) * 1e3 for nm, v in vals.items()}
    spread = max(vals.values()) - min(vals.values())
    avrec.append(dict(lam_vac_A=lam, lam_air_ciddor_A=ref, vac_minus_air_A=lam - ref,
                      diff_vs_ciddor_mA=d, spread_mA=spread * 1e3, spread_m_s=v_of_mA(spread * 1e3, lam),
                      note=('outside the 0.23-1.69 um fit range of Peck & Reeder eq.2 / Ciddor (extrapolation)'
                            if lam < 2300 or lam > 16900 else '')))
    avrows.append([lam, f'{ref:.4f}', f'{lam - ref:.4f}'] + [f'{d[nm]:+.3f}' for nm, _ in AV_FORMS[:-1]]
                  + [f'{spread * 1e3:.3f}', f'{v_of_mA(spread * 1e3, lam):.1f}'])
OUT['sec3_airvac_formulas'] = avrec
show('3a. Air wavelength from a fixed vacuum wavelength: formula minus Ciddor (1996), mA (standard air of each formula)',
     md_table(['lambda_vac (A)', 'lambda_air Ciddor (A)', 'vac-air (A)'] + [nm for nm, _ in AV_FORMS[:-1]]
              + ['spread (mA)', 'spread (m/s)'], avrows))

# Ciddor full-model sensitivities at 5000 A (standard air 15 C, 101325 Pa, dry, 450 ppm) and lab air
sens = {}
for tag, kw in [('standard_15C_101325Pa_dry_450ppm', dict(t_C=15.0, p_Pa=101325.0, rh=0.0, xc=450.0)),
                ('lab_20C_101325Pa_RH50_420ppm', dict(t_C=20.0, p_Pa=101325.0, rh=0.5, xc=420.0))]:
    la = float(AV.vac_to_air_ciddor1996_full(L0, **kw))
    def dl(**dk):
        k2 = dict(kw); k2.update(dk)
        return float(AV.vac_to_air_ciddor1996_full(L0, **k2))
    sens[tag] = dict(lam_air_A=la,
                     dmA_per_K=(dl(t_C=kw['t_C'] + 0.5) - dl(t_C=kw['t_C'] - 0.5)) * 1e3,
                     dmA_per_hPa=(dl(p_Pa=kw['p_Pa'] + 50) - dl(p_Pa=kw['p_Pa'] - 50)) * 1e3,
                     dmA_per_10pctRH=(dl(rh=kw['rh'] + 0.05) - dl(rh=max(kw['rh'] - 0.05, 0.0))) * 1e3 * (0.1 / (min(kw['rh'] + 0.05, 1) - max(kw['rh'] - 0.05, 0))),
                     dmA_per_100ppm_CO2=(dl(xc=kw['xc'] + 50) - dl(xc=kw['xc'] - 50)) * 1e3)
sens['lab_minus_standard_mA'] = (sens['lab_20C_101325Pa_RH50_420ppm']['lam_air_A'] - sens['standard_15C_101325Pa_dry_450ppm']['lam_air_A']) * 1e3
OUT['sec3_airvac_sensitivity_5000A'] = sens

# Isotope shifts, hyperfine and fine structure: verbatim literature / NIST ASD strings
iso = []


def add_split(what, s1, s2, medium, source, note=''):
    a, b = float(re.sub(r'\(.*\)', '', s1)), float(re.sub(r'\(.*\)', '', s2))
    lam = 0.5 * (a + b)
    d = (b - a) * 1e3
    iso.append(dict(what=what, wl1=s1, wl2=s2, medium=medium, split_mA=d, v_m_s=v_of_mA(d, lam), source=source, note=note))


MB = 'Murphy & Berengut 2014, MNRAS 438, 388 (arXiv:1311.2949v3)'
add_split('Mg II 2796: 24Mg -> 26Mg isotope shift', '2796.35509903(42)', '2796.34704565(42)', 'vacuum', MB + ' Table 2')
add_split('Mg II 2796: 25Mg hyperfine components F=2 -> F=3', '2796.353449(50)', '2796.349030(50)', 'vacuum', MB + ' Table 2')
add_split('Ca II 3934 (K): 40Ca -> 48Ca isotope shift', '3934.7747160(26)', '3934.76575(65)', 'vacuum', MB + ' Table 5')
add_split('Ca II 3934 (K): 43Ca hyperfine F=3 vs F=4', '3934.77998(11)', '3934.76438(11)', 'vacuum', MB + ' Table 5')
add_split('Na I 5891 (D2): 23Na ground-state hfs F=1 vs F=2', '5891.570875', '5891.590672', 'vacuum', MB + ' Table 1')
add_split('Mn II 2576: 55Mn hfs, extreme components', '2576.890898', '2576.856181', 'vacuum', MB + ' Table 9 (14-16 components)')

# Li 6708 D2: centre-of-gravity frequencies (Brown et al. 2013, PRA 87, 032504, Table I)
NU6 = '446799648.870(15)'
NU7 = '446810183.163(16)'
nu6 = float(NU6.split('(')[0]) * 1e6
nu7 = float(NU7.split('(')[0]) * 1e6
l6v, l7v = C / nu6 * 1e10, C / nu7 * 1e10
l6a = float(AV.vac_to_air_peckreeder1972_5p(l6v))
l7a = float(AV.vac_to_air_peckreeder1972_5p(l7v))
li = dict(nu_6Li_D2_cog_MHz=NU6, nu_7Li_D2_cog_MHz=NU7, source='Brown et al. 2013, PRA 87, 032504, Table I',
          isotope_shift_MHz_lit='10534.293(22) (Brown+2013 Table IV); 10534.357(29) Sansonetti+2011 corrected (PRL 109, 259901)',
          lam_vac_6Li_A=l6v, lam_vac_7Li_A=l7v, lam_air_6Li_A_PR72eq3=l6a, lam_air_7Li_A_PR72eq3=l7a,
          shift_mA=(l6v - l7v) * 1e3, shift_m_s=C * (nu7 - nu6) / nu7)
# 7Li D2 hyperfine: ground-state F=2 vs F=1 (components 2->3 and 1->2 of Table I)
nu_a, nu_b = 446809874.895e6, 446810687.873e6
li['7Li_D2_ground_hfs_components_MHz'] = ['446809874.895(20) [F=2->3]', '446810687.873(25) [F=1->2]']
li['7Li_D2_ground_hfs_split_mA'] = (C / nu_a - C / nu_b) * 1e13
iso.append(dict(what='Li I 6708 D2: 7Li -> 6Li isotope shift (cog)', wl1=f'{l7v:.6f}', wl2=f'{l6v:.6f}', medium='vacuum (derived from frequencies)',
                split_mA=li['shift_mA'], v_m_s=li['shift_m_s'], source=li['source'], note='derived lambda = c/nu'))
iso.append(dict(what='Li I 6708 D2: 7Li ground-state hfs (F=2 vs F=1)', wl1='', wl2='', medium='vacuum (derived)',
                split_mA=li['7Li_D2_ground_hfs_split_mA'], v_m_s=v_of_mA(li['7Li_D2_ground_hfs_split_mA'], 6709.7),
                source=li['source'], note='components F=2->3 and F=1->2'))
OUT['sec3_lithium'] = li


def read_asd(fn):
    rows_ = []
    txt = open(fn, encoding='utf-8', errors='replace').read()
    rd = csv.reader(io.StringIO(txt))
    hdr = next(rd)
    for r in rd:
        if not r:
            continue
        rows_.append({h: re.sub(r'^="(.*)"$', r'\1', v).strip() for h, v in zip(hdr, r)})
    return hdr, rows_


ASD = {}
asd_dir = os.path.join(ROOT, 'data', 'raw', 'precision_rationale', 'nist_asd')
for fn in sorted(os.listdir(asd_dir)):
    if fn.endswith('.csv'):
        hdr, rr = read_asd(os.path.join(asd_dir, fn))
        ASD[fn[:-4]] = [{k: v for k, v in r.items() if k and v} for r in rr]
OUT['sec3_nist_asd_windows'] = ASD
# H I / D I Halpha fine-structure components (air) and composite
def _num(s):
    m = re.match(r'^[\[(]?([-+]?\d+\.?\d*)', s or '')
    return float(m.group(1)) if m else None


hcomp = [(_num(r.get('ritz_wl_air(A)')), r.get('ritz_wl_air(A)'), r.get('unc_ritz_wl'), r.get('J_i', ''), r.get('J_k', ''),
          r.get('conf_i', ''), r.get('conf_k', ''), r.get('obs_wl_air(A)', '')) for r in ASD['H_I_Halpha'] if r.get('Aki(s^-1)')]
fs = [h for h in hcomp if h[5] in ('2s', '2p')]
comp = [h for h in hcomp if h[5] == '2']
dcomp = [(_num(r.get('ritz_wl_air(A)')), r.get('ritz_wl_air(A)'), r.get('conf_i', '')) for r in ASD['D_I_Halpha'] if r.get('Aki(s^-1)')]
d_comp = [h for h in dcomp if h[2] == '2']
halpha = dict(fine_structure_components_air=[h[1] for h in fs],
              fine_structure_span_mA=(max(h[0] for h in fs) - min(h[0] for h in fs)) * 1e3,
              composite_ritz_air=comp[0][1] if comp else None, composite_obs_air=comp[0][7] if comp else None,
              D_composite_ritz_air=d_comp[0][1] if d_comp else None)
if comp and d_comp:
    halpha['H_minus_D_composite_mA'] = (comp[0][0] - d_comp[0][0]) * 1e3
    halpha['H_minus_D_v_m_s'] = v_of_mA(halpha['H_minus_D_composite_mA'], comp[0][0])
OUT['sec3_halpha'] = halpha
iso.append(dict(what='H I Halpha: fine-structure components span (2s/2p - 3s/3p/3d)', wl1=min(fs)[1], wl2=max(fs)[1], medium='air',
                split_mA=halpha['fine_structure_span_mA'], v_m_s=v_of_mA(halpha['fine_structure_span_mA'], 6562.8),
                source='NIST ASD (H I), downloaded 2026-11-23', note=f"composite Ritz {halpha['composite_ritz_air']}, observed {halpha['composite_obs_air']}"))
if 'H_minus_D_composite_mA' in halpha:
    iso.append(dict(what='Halpha: H -> D isotope shift (composites)', wl1=comp[0][1], wl2=d_comp[0][1], medium='air',
                    split_mA=-halpha['H_minus_D_composite_mA'], v_m_s=-halpha['H_minus_D_v_m_s'], source='NIST ASD (H I, D I)', note=''))
OUT['sec3_isotope_hfs'] = iso
show('3b. Isotope, hyperfine and fine-structure splittings (verbatim wavelength strings; split = wl2 - wl1)',
     md_table(['structure', 'wl1 (A)', 'wl2 (A)', 'medium', 'split (mA)', 'split (m/s)', 'source'],
              [[r['what'], r['wl1'], r['wl2'], r['medium'], f"{r['split_mA']:+.3f}", f"{r['v_m_s']:+.0f}", r['source']] for r in iso]))

effrows = [[e['effect'], g(e['lam_A'], 6), g(e['v_m_s'], 4), g(e['dlam_mA'], 3), e['flag'], e['source'] + (('; ' + e['note']) if e['note'] else '')]
           for e in eff]
show('3c. Effects on line positions, expressed at lambda (mA)', md_table(['effect', 'lambda (A)', 'v (m/s)', 'dlambda (mA)', 'flag', 'source / note'], effrows))
show('3d. Zeeman sigma-component offset dlam = g e B lam^2/(4 pi m_e c)',
     md_table(['B (G)', 'g_eff', 'lambda (A)', 'dlam (mA)', 'v (m/s)'],
              [[z['B_G'], z['g_eff'], z['lam_A'], g(z['dlam_mA'], 3), g(z['v_m_s'], 3)] for z in zee]))

# =====================================================================================
# (4) Catalogue rounding: NIST Handbook decimals, and Handbook vs ASD
# =====================================================================================
HB = json.load(open(os.path.join(ROOT, 'data', 'handbook.json')))
lines = HB['lines']
dec = Counter(l['dec'] for l in lines)
decm = Counter((l['medium'], l['dec']) for l in lines)
cat = dict(n_lines=len(lines), by_decimals={str(k): dec[k] for k in sorted(dec)},
           by_medium_decimals={f'{m}_{d}': n for (m, d), n in sorted(decm.items())},
           n_le_1_decimal=dec[0] + dec[1], n_2_decimals=dec[2], n_3_decimals=dec[3], n_ge_4_decimals=sum(v for k, v in dec.items() if k >= 4))
cat['frac_2_decimals'] = dec[2] / len(lines)
for d in (1, 2, 3, 4):
    step = 10.0 ** (-d) * 1e3
    cat[f'rounding_{d}dec'] = dict(step_mA=step, max_err_mA=step / 2, rms_err_mA=step / math.sqrt(12),
                                   step_m_s_5000A=v_of_mA(step, 5000), rms_m_s_5000A=v_of_mA(step / math.sqrt(12), 5000))
try:
    M = json.load(open(os.path.join(ROOT, 'data', 'handbook_asd.json')))['matches']
    stats = {}
    for d in (2, 3, 4):
        dl, un, v, n_all, n_finer, n_not_finer = [], [], [], 0, 0, 0
        for m in M:
            l = lines[m['i']]
            if l['dec'] != d or m.get('cls') not in ('exact', 'close') or m.get('dl_mA') is None:
                continue
            n_all += 1
            if m.get('asd_unc') is not None:
                un.append(m['asd_unc'] * 1e3)
            # the rounding error of the Handbook string is only visible against an ASD value (observed or
            # Ritz) that is quoted with MORE decimals; take the finer one, ties broken by smaller stated unc.
            cands = []
            for key, ukey in (('asd_obs', 'asd_unc_obs'), ('asd_ritz', 'asd_unc_ritz')):
                sv = str(m.get(key) or '').strip()
                mm = re.match(r'^\[?\(?(\d+)\.(\d+)', sv)
                if mm and len(mm.group(2)) > d:
                    cands.append((-len(mm.group(2)), m.get(ukey) if m.get(ukey) is not None else 9e9,
                                  float(mm.group(1) + '.' + mm.group(2))))
            if not cands:
                n_not_finer += 1
                continue
            n_finer += 1
            best = sorted(cands)[0][2]
            dlm = (l['wl'] - best) * 1e3
            dl.append(abs(dlm))
            v.append(abs(v_of_mA(dlm, l['wl'])))
        dl, un, v = np.array(dl), np.array(un), np.array(v)
        q = (lambda a, p: float(np.percentile(a, p)) if a.size else None)
        stats[str(d)] = dict(n_exact_or_close=n_all, n_asd_more_decimals=n_finer, n_asd_not_more_decimals=n_not_finer,
                             median_abs_dl_mA=q(dl, 50), p68_abs_dl_mA=q(dl, 68), p95_abs_dl_mA=q(dl, 95),
                             rms_dl_mA=float(np.sqrt(np.mean(dl ** 2))) if dl.size else None,
                             median_abs_dv_m_s=q(v, 50), p95_abs_dv_m_s=q(v, 95),
                             n_with_asd_unc=int(un.size), median_asd_unc_mA=q(un, 50),
                             frac_asd_unc_le_1mA=float(np.mean(un <= 1.0)) if un.size else None,
                             note='dl statistics only over matches where the ASD wavelength string has more decimals than the Handbook')
    cls = Counter(m.get('cls') for m in M)
    cat['handbook_vs_asd'] = dict(classes=dict(cls), by_handbook_decimals=stats,
                                  note='snapshot of data/handbook_asd.json (built by scripts/match_asd.py); dl = Handbook - ASD')
except Exception as ex:  # the cross-match is produced by another process; tolerate its absence
    cat['handbook_vs_asd'] = dict(error=str(ex))
OUT['sec4_catalogue'] = cat
# examples: same line at different precisions (NIST ASD windows)
fe = ASD.get('Fe_I_5250', [])
OUT['sec4_examples'] = dict(
    Fe_I_5250p2=dict(asd_obs_air=fe[0].get('obs_wl_air(A)') if fe else None, asd_unc_obs=fe[0].get('unc_obs_wl') if fe else None,
                     asd_ritz_air=fe[0].get('ritz_wl_air(A)') if fe else None, asd_unc_ritz=fe[0].get('unc_ritz_wl') if fe else None,
                     lars_lab_air='5250.2084', lars_source='Lohner-Bottcher et al. 2019 Table 1 (LARS hollow-cathode, LFC-calibrated)'),
    Li_I_6708=ASD.get('Li_I_6708'), Ca_II_8542=ASD.get('Ca_II_8542'), Hg_I_5461=ASD.get('Hg_I_5461'), Ba_II_4554=ASD.get('Ba_II_4554'))

# =====================================================================================
# (5) Where 1 mA is not enough, and where it is overkill
# =====================================================================================
need = []
def addn(what, v_ms, lam, source):
    need.append(dict(what=what, v_m_s=v_ms, lam_A=lam, dlam_A=lam * v_ms / C, dlam_mA=lam * v_ms / C * 1e3, source=source))


addn('EPRV goal 10 cm/s', 0.10, 5000, 'Fischer et al. 2016, PASP 128, 066001')
addn('Earth-twin reflex of the Sun 9 cm/s', 0.09, 5000, 'Fischer et al. 2016')
addn('ideal calibrator single-measurement precision ~1 cm/s', 0.01, 5000, 'Fischer et al. 2016')
addn('LFC-vs-LFC agreement ~1 cm/s (HARPS)', 0.01, 5000, 'Probst et al. 2020, Nature Astronomy 4, 603 (doi:10.1038/s41550-020-1010-x)')
addn('1 m/s exoplanet RV', 1.0, 5000, 'definition')
addn('Lab rest-frame requirement for Delta-alpha/alpha work, dv <~ 20 m/s, at Fe II 2382', 20.0, 2382.7642, MB + ' Sect. 1')
addn('Redshift drift ~6 cm/s per decade at z=4 (Lya forest observed near 6080 A)', 0.06, 1215.67 * 5, 'Liske et al. 2008, MNRAS 386, 1192')
# alpha: Delta v for Delta-alpha/alpha = 1 ppm from q (Murphy & Berengut 2014 via data/notes json)
try:
    MBJ = json.load(open(os.path.join(ROOT, 'data', 'notes', 'murphy_berengut_2014_mm_transitions.json')))
    dvs = []
    for t in MBJ:
        w0 = float(t['omega0_cm-1'].split('(')[0])
        dvs.append((-2 * C * t['q_cm-1'] * 1e-6 / w0, float(t['lambda0_vac_A'].split('(')[0]), t['ion'], t['transition']))
    dv_arr = np.array([d[0] for d in dvs])
    OUT['sec5_alpha'] = dict(n=len(dvs), dv_1ppm_min_m_s=float(dv_arr.min()), dv_1ppm_max_m_s=float(dv_arr.max()),
                             dlam_1ppm_min_mA=float(min(d[0] * d[1] / C for d in dvs) * 1e3),
                             dlam_1ppm_max_mA=float(max(d[0] * d[1] / C for d in dvs) * 1e3),
                             formula='dv = -2 c q (Delta-alpha/alpha) / omega0', source=MB)
except Exception as ex:
    OUT['sec5_alpha'] = dict(error=str(ex))
# redshift drift with Planck18: vdot = c zdot/(1+z), zdot = (1+z) H0 - H(z)
drift = []
for z in [1, 2, 3, 4, 5]:
    H0 = Planck18.H0.to(1 / u.yr).value
    Hz = Planck18.H(z).to(1 / u.yr).value
    zdot = (1 + z) * H0 - Hz
    vdot = C * zdot / (1 + z)                     # m/s per yr
    lam_obs = 1215.67 * (1 + z)
    drift.append(dict(z=z, vdot_cm_s_per_decade=vdot * 10 * 100, lam_obs_Lya_A=lam_obs,
                      dlam_per_decade_A=lam_obs * vdot * 10 / C, dlam_per_20yr_A=lam_obs * vdot * 20 / C))
OUT['sec5_redshift_drift_planck18'] = drift
# overkill: survey requirements in mA
over = []
def addo(what, v_ms, lam, source):
    over.append(dict(what=what, v_m_s=v_ms, lam_A=lam, dlam_mA=lam * v_ms / C * 1e3, ratio_to_1mA=lam * v_ms / C * 1e3, source=source))


addo('DESI LRG/ELG random redshift error < 0.0005(1+z) = 150 km/s rms', 150e3, 7000, 'DESI Collaboration 2022, AJ 164, 207, req. L2.2.2/L2.3.2')
addo('DESI systematic mean-redshift error < 0.0002(1+z) = 60 km/s', 60e3, 7000, 'DESI Collaboration 2022, req. L2.2.3')
addo('DESI QSO redshift error 0.0025(1+z) = 750 km/s rms', 750e3, 7000, 'DESI Collaboration 2022, req. L2.4.2')
addo('4MOST HRS/LRS velocity accuracy < 1 km/s', 1e3, 5450, 'de Jong et al. 2019, Messenger 175, 3, Table 1')
addo('LAMOST-MRS RV precision ~1 km/s (late-type stars)', 1e3, 5150, 'Liu et al. 2020, arXiv:2005.07210')
addo('Gaia RVS RV systematic error after calibration <= 300 m/s', 300.0, 8600, 'Cropper et al. 2018, A&A 616, A5, Table 1')
addo('Galaxy velocity dispersion 100-200 km/s (BOSS resolution driver), lower end', 100e3, 6000, 'Smee et al. 2013, Sect. 2.1.5')
OUT['sec5_not_enough'] = need
OUT['sec5_overkill'] = over
show('5a. Requirements finer than 1 mA', md_table(['requirement', 'v (m/s)', 'lambda (A)', 'dlambda (A)', 'source'],
     [[n['what'], g(n['v_m_s'], 3), g(n['lam_A'], 6), f"{n['dlam_A']:.2e}", n['source']] for n in need]))
show('5b. Redshift drift (Planck18, derived): vdot = c zdot/(1+z), zdot = (1+z)H0 - H(z)',
     md_table(['z', 'vdot (cm/s per decade)', 'Lya observed (A)', 'dlam per decade (A)', 'dlam per 20 yr (A)'],
              [[d['z'], f"{d['vdot_cm_s_per_decade']:+.2f}", f"{d['lam_obs_Lya_A']:.1f}", f"{d['dlam_per_decade_A']:+.2e}",
                f"{d['dlam_per_20yr_A']:+.2e}"] for d in drift]))
show('5c. Survey requirements coarser than 1 mA', md_table(['requirement', 'v (m/s)', 'lambda (A)', 'dlambda (mA)', 'source'],
     [[o['what'], g(o['v_m_s'], 3), g(o['lam_A'], 5), g(o['dlam_mA'], 3), o['source']] for o in over]))

json.dump(OUT, open(os.path.join(ROOT, 'data', 'notes', 'precision_rationale_numbers.json'), 'w'), indent=1, default=float)
print('\n'.join(PRINT))
print('\n### extra numbers\n')
print(json.dumps({k: OUT[k] for k in ['crb_numeric_checks', 'sec2_threshold', 'sec3_airvac_sensitivity_5000A', 'sec3_lithium',
                                      'sec3_halpha', 'sec4_catalogue', 'sec5_alpha']}, indent=1, default=float))
print('\nZeeman k =', kZ_AG, 'A/(G A^2)')
print('wrote', csv_path)
