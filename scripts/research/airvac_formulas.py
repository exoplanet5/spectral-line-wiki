"""
airvac_formulas.py -- air <-> vacuum wavelength conversion formulas, researched to 1e-4 A.

All wavelengths are in Angstrom (A).  sigma is the wavenumber in micrometre^-1, sigma^2 in um^-2.
n is the phase refractive index of air; lambda_air = lambda_vac / n.

Every coefficient below is copied from the primary paper or the cited implementation
(see data/notes/airvac.md for verbatim sources, URLs, DOIs and bibcodes).  Nothing is fitted here.

Naming convention
-----------------
  vac_to_air_<name>(wl_vac_A)  -> air wavelength (A)
  air_to_vac_<name>(wl_air_A)  -> vacuum wavelength (A)
Formulas whose argument is the VACUUM wavenumber (Edlen, Peck & Reeder, Birch & Downs, Ciddor,
Morton) are explicit for vac->air; the inverse is solved exactly by fixed-point iteration
(convergence factor ~ lambda*dn/dlambda ~ 1e-6 per step) to |delta| < 1e-12 A.
Formulas whose argument is the AIR wavelength (Greisen et al. 2006 eq. 65, WCSLIB AWAV, VALD
Piskunov inverse) are explicit for air->vac; their vac->air is solved by iteration.
Functions named *_literal reproduce a published software implementation exactly as coded
(including its approximations), so that their errors can be quantified.
"""
import numpy as np

TORR_PA = 101325.0 / 760.0          # 1 Torr (mmHg, conventional) in Pa
_ITER_TOL = 1e-12                   # A
_ITER_MAX = 60


def _a(x):
    return np.asarray(x, dtype=float)


def _s2(wl_A):
    """sigma^2 in um^-2 for a wavelength in A."""
    return (1.0e4 / _a(wl_A)) ** 2


def _invert_vacarg(nm1_of_s2, wl_air_A):
    """Solve lambda_vac = lambda_air * (1 + nm1(sigma_vac^2)) exactly (fixed point)."""
    wa = _a(wl_air_A)
    wv = wa * (1.0 + nm1_of_s2(_s2(wa)))
    for _ in range(_ITER_MAX):
        new = wa * (1.0 + nm1_of_s2(_s2(wv)))
        if np.all(np.abs(new - wv) < _ITER_TOL):
            return new
        wv = new
    return wv


def _invert_airarg(nm1_of_s2air, wl_vac_A):
    """Solve lambda_air = lambda_vac / (1 + nm1(sigma_air^2)) exactly (fixed point)."""
    wv = _a(wl_vac_A)
    wa = wv / (1.0 + nm1_of_s2air(_s2(wv)))
    for _ in range(_ITER_MAX):
        new = wv / (1.0 + nm1_of_s2air(_s2(wa)))
        if np.all(np.abs(new - wa) < _ITER_TOL):
            return new
        wa = new
    return wa


# =====================================================================================
# 1. Standard-air dispersion formulas, argument = VACUUM wavenumber squared (um^-2)
# =====================================================================================
def nm1_edlen1953(s2):
    """Edlen (1953) JOSA 43, 339, abstract:  (n-1)1e8 = 6432.8 + 2949810/(146-s^2) + 25540/(41-s^2).
    Standard air: dry, 15 C, 760 Torr, 0.03 % CO2.  IAU standard (Oosterhoff 1957, Trans. IAU 9, 69 & 202;
    re-affirmed IAU 1991 Comm. 44 resolution C15); used by Morton (1991) eq. 3."""
    return 1e-8 * (6432.8 + 2949810.0 / (146.0 - s2) + 25540.0 / (41.0 - s2))


def nm1_edlen1966(s2):
    """Edlen (1966) Metrologia 2, 71, eq. (1): (n-1)_s 1e8 = 8342.13 + 2406030/(130-s^2) + 15997/(38.9-s^2).
    Standard air: dry, 15 C, 760 Torr, 0.0003 (300 ppm) CO2 by volume."""
    return 1e-8 * (8342.13 + 2406030.0 / (130.0 - s2) + 15997.0 / (38.9 - s2))


def nm1_peckreeder1972(s2):
    """Peck & Reeder (1972) JOSA 62, 958, eq. (2) (two-term Sellmeier, 4 parameters, 0.23-1.69 um):
    (n-1)1e8 = 5791817/(238.0185-s^2) + 167909/(57.362-s^2).  15 C, 101325 Pa, dry, 0.03% CO2
    (their far-UV table is 'recalculated for 15 C and 0.033% CO2')."""
    return 1e-8 * (5791817.0 / (238.0185 - s2) + 167909.0 / (57.362 - s2))


def nm1_peckreeder1972_5p(s2):
    """Peck & Reeder (1972) eq. (3) (5-parameter fit extended to 0.185 um), the formula NIST ASD
    says it uses ('five-parameter formula'):  (n-1)1e8 = 8060.51 + 2480990/(132.274-s^2)
    + 17455.7/(39.32957-s^2).  Coefficients as given by refractiveindex.info (Peck page); verified here
    against the Eq.(3) column of P&R Table V (see validate_airvac.py)."""
    return 1e-8 * (8060.51 + 2480990.0 / (132.274 - s2) + 17455.7 / (39.32957 - s2))


def nm1_birchdowns1993(s2):
    """Birch & Downs (1993) Metrologia 30, 155, eq. (8) (Edlen 1966 x 1.0001099; ITS-90 sign error):
    (n-1)_s 1e8 = 8343.05 + 2406294/(130-s^2) + 15999/(38.9-s^2).  15 C, 101325 Pa, dry, 450 ppm CO2.
    Superseded by the 1994 correction."""
    return 1e-8 * (8343.05 + 2406294.0 / (130.0 - s2) + 15999.0 / (38.9 - s2))


def nm1_birchdowns1994(s2):
    """Birch & Downs (1994) Metrologia 31, 315 (corrected; overall factor 1.0000486 on Edlen 1966):
    (n-1)_s 1e8 = 8342.54 + 2406147/(130-s^2) + 15998/(38.9-s^2).  15 C, 101325 Pa, dry, 450 ppm CO2.
    Identical to Morton (2000) ApJS 130, 403, eq. 8 and to VALD, CRC Handbook."""
    return 1e-8 * (8342.54 + 2406147.0 / (130.0 - s2) + 15998.0 / (38.9 - s2))


def nm1_morton2000(s2):
    """Morton (2000) eq. 8 as written on the VALD wiki:
    n = 1 + 0.0000834254 + 0.02406147/(130 - s^2) + 0.00015998/(38.9 - s^2), s = 1e4/lambda_vac(A)."""
    return 0.0000834254 + 0.02406147 / (130.0 - s2) + 0.00015998 / (38.9 - s2)


def nm1_ciddor1996(s2, xc=450.0):
    """Ciddor (1996) Appl. Opt. 35, 1566, eqs. (1)-(2), standard air 15 C, 101325 Pa, 0 % RH:
    (n_as-1)1e8 = k1/(k0-s^2) + k3/(k2-s^2), k0=238.0185, k1=5792105, k2=57.362, k3=167917 (um^-2);
    (n_axs-1) = (n_as-1)[1 + 0.534e-6 (xc-450)], xc in ppm CO2."""
    return 1e-8 * (5792105.0 / (238.0185 - s2) + 167917.0 / (57.362 - s2)) * (1.0 + 0.534e-6 * (xc - 450.0))


def nm1_sdss_dr7_poly(wl_vac_A):
    """SDSS DR7 'vacwavelength' page (attributed to Morton 1991 IAU standard; same Cauchy form as the
    pre-2011 IDL Astronomy Library VACTOAIR):  AIR = VAC/(1.0 + 2.735182E-4 + 131.4182/VAC^2
    + 2.76249E8/VAC^4), VAC in A."""
    v = _a(wl_vac_A)
    return 2.735182e-4 + 131.4182 / v ** 2 + 2.76249e8 / v ** 4


# =====================================================================================
# 2. Formulas whose argument is the AIR wavelength
# =====================================================================================
def nm1_iag1999_phase(lam_um):
    """IUGG/IAG (1999) Resolution 3 phase refractivity of 'standard air' with 0.0375 % CO2 at
    T = 273.15 K (0 C), p = 1013.25 hPa, e = 0:  N_ph = (n_ph-1)1e6 = 287.6155 + 1.62887/l^2 + 0.01360/l^4.
    Greisen et al. (2006) A&A 446, 747, eq. (65) adopt it with l = air wavelength (um)."""
    l = _a(lam_um)
    return 1e-6 * (287.6155 + 1.62887 / l ** 2 + 0.01360 / l ** 4)


def nm1_iag1999_group(lam_um):
    """IAG (1999) group refractivity N_g = 287.6155 + 4.88660/l^2 + 0.06800/l^4 (0 C, 1013.25 hPa, 375 ppm)."""
    l = _a(lam_um)
    return 1e-6 * (287.6155 + 4.88660 / l ** 2 + 0.06800 / l ** 4)


def nm1_wcslib_airarg(s2air):
    """WCSLIB 8.x spx.c (astropy.wcs, AWAV): n = 1.000064328 + 294.981e8/(1.46e14 - s) + 2.554e8/(0.41e14 - s),
    s = (1/lambda_air[m])^2  -- i.e. Edlen 1953 coefficients evaluated at the AIR wavenumber."""
    return 6.4328e-5 + 2.94981e-2 / (146.0 - s2air) + 2.5540e-4 / (41.0 - s2air)


def nm1_vald_piskunov(s2air):
    """VALD3 air->vac (N. Piskunov): n = 1 + 0.00008336624212083 + 0.02408926869968/(130.1065924522 - s^2)
    + 0.0001599740894897/(38.92568793293 - s^2), s = 1e4/lambda_air(A)."""
    return (0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s2air)
            + 0.0001599740894897 / (38.92568793293 - s2air))


# =====================================================================================
# 3. Public vac_to_air_* / air_to_vac_* pairs (exact inverses)
# =====================================================================================
def _make_vacarg_pair(nm1):
    def v2a(wl_vac_A):
        wv = _a(wl_vac_A)
        return wv / (1.0 + nm1(_s2(wv)))

    def a2v(wl_air_A):
        return _invert_vacarg(nm1, wl_air_A)
    return v2a, a2v


vac_to_air_edlen1953, air_to_vac_edlen1953 = _make_vacarg_pair(nm1_edlen1953)
vac_to_air_edlen1966, air_to_vac_edlen1966 = _make_vacarg_pair(nm1_edlen1966)
vac_to_air_peckreeder1972, air_to_vac_peckreeder1972 = _make_vacarg_pair(nm1_peckreeder1972)
vac_to_air_peckreeder1972_5p, air_to_vac_peckreeder1972_5p = _make_vacarg_pair(nm1_peckreeder1972_5p)
vac_to_air_birchdowns1993, air_to_vac_birchdowns1993 = _make_vacarg_pair(nm1_birchdowns1993)
vac_to_air_birchdowns1994, air_to_vac_birchdowns1994 = _make_vacarg_pair(nm1_birchdowns1994)
vac_to_air_morton2000, air_to_vac_morton2000 = _make_vacarg_pair(nm1_morton2000)
vac_to_air_ciddor1996, air_to_vac_ciddor1996 = _make_vacarg_pair(nm1_ciddor1996)
# NIST ASD: Peck & Reeder five-parameter formula (eq. 3); SDSS DR9+: Ciddor standard; IDL vactoair: Ciddor
vac_to_air_nist_asd, air_to_vac_nist_asd = vac_to_air_peckreeder1972_5p, air_to_vac_peckreeder1972_5p
vac_to_air_sdss_dr9, air_to_vac_sdss_dr9 = vac_to_air_ciddor1996, air_to_vac_ciddor1996


def vac_to_air_ciddor1996_co2(wl_vac_A, xc):
    wv = _a(wl_vac_A)
    return wv / (1.0 + nm1_ciddor1996(_s2(wv), xc))


def air_to_vac_ciddor1996_co2(wl_air_A, xc):
    return _invert_vacarg(lambda s2: nm1_ciddor1996(s2, xc), wl_air_A)


def vac_to_air_sdss_dr7(wl_vac_A):
    wv = _a(wl_vac_A)
    return wv / (1.0 + nm1_sdss_dr7_poly(wv))


def air_to_vac_sdss_dr7(wl_air_A):
    wa = _a(wl_air_A)
    wv = wa * (1.0 + nm1_sdss_dr7_poly(wa))
    for _ in range(_ITER_MAX):
        new = wa * (1.0 + nm1_sdss_dr7_poly(wv))
        if np.all(np.abs(new - wv) < _ITER_TOL):
            return new
        wv = new
    return wv


# Greisen et al. 2006 (FITS WCS Paper III) eq. 64-65: lambda = n(lambda_a) * lambda_a, lambda_a in um
def air_to_vac_greisen2006(wl_air_A):
    wa = _a(wl_air_A)
    return wa * (1.0 + nm1_iag1999_phase(wa * 1e-4))


def vac_to_air_greisen2006(wl_vac_A):
    """Exact inverse of eq. 65.  (Greisen et al. themselves recommend lambda_a = lambda/n(lambda),
    i.e. evaluating n at the vacuum wavelength; see vac_to_air_greisen2006_literal.)"""
    return _invert_airarg(lambda s2a: nm1_iag1999_phase(1.0 / np.sqrt(s2a)), wl_vac_A)


def vac_to_air_greisen2006_literal(wl_vac_A):
    wv = _a(wl_vac_A)
    return wv / (1.0 + nm1_iag1999_phase(wv * 1e-4))


def dlamvac_dlamair_greisen2006(wl_air_A):
    """Greisen et al. (2006) eq. (66): d lambda/d lambda_a = 1 + 1e-6 (287.6155 - 1.62887/l^2 - 0.04080/l^4)."""
    l = _a(wl_air_A) * 1e-4
    return 1.0 + 1e-6 * (287.6155 - 1.62887 / l ** 2 - 0.04080 / l ** 4)


# WCSLIB / astropy.wcs AWAV <-> WAVE (Edlen 1953 coefficients with air-wavenumber argument)
def air_to_vac_wcslib(wl_air_A):
    wa = _a(wl_air_A)
    return wa * (1.0 + nm1_wcslib_airarg(_s2(wa)))


def vac_to_air_wcslib(wl_vac_A):
    return _invert_airarg(nm1_wcslib_airarg, wl_vac_A)


def vac_to_air_wcslib_literal(wl_vac_A):
    """spx.c waveawav(): 4 fixed iterations n_k from s = (n_{k-1}/lambda_vac)^2, n_0 = 1."""
    wv = _a(wl_vac_A)
    n = np.ones_like(wv)
    for _ in range(4):
        n = 1.0 + nm1_wcslib_airarg((n * 1e4 / wv) ** 2)
    return wv / n


# VALD: vac->air = Morton 2000 (vacuum argument); air->vac = Piskunov closed form (air argument)
vac_to_air_vald = vac_to_air_morton2000


def air_to_vac_vald(wl_air_A):
    wa = _a(wl_air_A)
    return wa * (1.0 + nm1_vald_piskunov(_s2(wa)))


# IAG 1999 phase formula used at its stated conditions with the VACUUM ('carrier') wavelength
vac_to_air_iag1999, air_to_vac_iag1999 = _make_vacarg_pair(lambda s2: nm1_iag1999_phase(1.0 / np.sqrt(s2)))


# =====================================================================================
# 4. Literal re-implementations of software packages (to quantify their approximations)
# =====================================================================================
def air_to_vac_idl_astrolib_literal(wl_air_A):
    """IDL Astronomy Library AIRTOVAC (Landsman/Schlegel 2011): Ciddor standard air, 2 iterations,
    starting from wave_vac = wave_air; wavelengths < 2000 A unchanged."""
    wa = _a(wl_air_A)
    wv = wa.copy()
    for _ in range(2):
        wv = wa * (1.0 + nm1_ciddor1996(_s2(wv)))
    return np.where(wa >= 2000.0, wv, wa)


def vac_to_air_idl_astrolib_literal(wl_vac_A):
    wv = _a(wl_vac_A)
    return np.where(wv >= 2000.0, wv / (1.0 + nm1_ciddor1996(_s2(wv))), wv)


def _specutils_nm1(method):
    m = method.lower()
    return {'edlen1953': nm1_edlen1953, 'edlen1966': nm1_edlen1966, 'morton2000': nm1_morton2000,
            'peckreeder1972': nm1_peckreeder1972, 'ciddor1996': nm1_ciddor1996,
            'greisen2006': lambda s2: nm1_iag1999_phase(1.0 / np.sqrt(s2))}[m]


def air_to_vac_specutils_inversion_literal(wl_air_A, method='Morton2000'):
    """specutils 2.x wcs_utils.air_to_vac(scheme='inversion'): lambda_vac = lambda_air * n(lambda_air)
    (n evaluated at the AIR wavelength although the formulas are defined for the vacuum one)."""
    wa = _a(wl_air_A)
    return wa * (1.0 + _specutils_nm1(method)(_s2(wa)))


def vac_to_air_specutils_literal(wl_vac_A, method='Morton2000'):
    wv = _a(wl_vac_A)
    return wv / (1.0 + _specutils_nm1(method)(_s2(wv)))


def nm1_iraf_disptrans(wl_A, t=15.0, p=760.0, f=4.0):
    """IRAF noao.onedspec.disptrans (dispcor/t_disptrans.x dt_airvac) and noao/lib/linelists/airvac.cl,
    'formulas in Allen's Astrophysical Quantities (1973)': x2 = 1/lambda[um]^2,
    a = 64.328 + 29498.1/(146-x2) + 255.4/(41-x2)
    a = a*p*(1+(1.049-0.0157 t)1e-6 p)/(720.883 (1+0.003661 t))
    a = a - (0.0624 - 0.000680 x2)/(1+0.003661 t) * f ;  n = 1 + a/1e6
    t in C, p and f (water-vapour pressure) in mmHg.  Defaults t=15, p=760, f=4 (NOT dry air)."""
    x2 = _s2(wl_A)
    a = 64.328 + 29498.1 / (146.0 - x2) + 255.4 / (41.0 - x2)
    a = a * p * (1.0 + (1.049 - 0.0157 * t) * 1e-6 * p) / (720.883 * (1.0 + 0.003661 * t))
    a = a - (0.0624 - 0.000680 * x2) / (1.0 + 0.003661 * t) * f
    return a / 1e6


def vac_to_air_iraf_literal(wl_vac_A, t=15.0, p=760.0, f=4.0):
    wv = _a(wl_vac_A)
    return wv / (1.0 + nm1_iraf_disptrans(wv, t, p, f))


def air_to_vac_iraf_literal(wl_air_A, t=15.0, p=760.0, f=4.0):
    """IRAF evaluates n at the input (air) wavelength -- not the exact inverse of its vac2air."""
    wa = _a(wl_air_A)
    return wa * (1.0 + nm1_iraf_disptrans(wa, t, p, f))


def air_to_vac_iraf_astlvac_literal(wl_air_A):
    """IRAF noao/astutil/asttools/astlvac.x: Edlen 1966 with s2 = 1e8/lambda_air^2 (air argument)."""
    wa = _a(wl_air_A)
    return wa * (1.0 + nm1_edlen1966(_s2(wa)))


# =====================================================================================
# 5. Full T / p / humidity / CO2 models
# =====================================================================================
def svp_davis1992(T_K):
    """Saturation vapour pressure over liquid water used by Ciddor (1996) App. A (BIPM 1981/91, Davis 1992):
    svp = exp(A T^2 + B T + C + D/T) Pa, A=1.2378847e-5 K^-2, B=-1.9121316e-2 K^-1, C=33.93711047, D=-6.3431645e3 K."""
    T = _a(T_K)
    return np.exp(1.2378847e-5 * T ** 2 - 1.9121316e-2 * T + 33.93711047 - 6.3431645e3 / T)


def svp_ice_marti1993(T_K):
    """Ciddor (1996) App. C, over ice (Marti & Mauersberger 1993): log10(svp) = -2663.5/T + 12.537 (Pa)."""
    return 10.0 ** (-2663.5 / _a(T_K) + 12.537)


_IAPWS_K = (1.16705214528e3, -7.24213167032e5, -1.70738469401e1, 1.20208247025e4, -3.23255503223e6,
            1.49151086135e1, -4.82326573616e3, 4.05113405421e5, -2.38555575678e-1, 6.50175348448e2)


def svp_iapws(T_K):
    """NIST EMToolBox eqs. A1-A7 (IAPWS-IF97 saturation line, Huang 1998), over water, Pa."""
    K1, K2, K3, K4, K5, K6, K7, K8, K9, K10 = _IAPWS_K
    T = _a(T_K)
    Om = T + K9 / (T - K10)
    A = Om ** 2 + K1 * Om + K2
    B = K3 * Om ** 2 + K4 * Om + K5
    C = K6 * Om ** 2 + K7 * Om + K8
    X = -B + np.sqrt(B ** 2 - 4 * A * C)
    return 1e6 * (2 * C / X) ** 4


def svp_iapws_ice(T_K):
    """NIST EMToolBox eqs. A8-A13 (over ice): Theta = T/273.16, Y = A1(1-Theta^-1.5) + A2(1-Theta^-1.25),
    psv = 611.657 e^Y, A1=-13.928169, A2=34.7078238."""
    Th = _a(T_K) / 273.16
    Y = -13.928169 * (1 - Th ** -1.5) + 34.7078238 * (1 - Th ** -1.25)
    return 611.657 * np.exp(Y)


def _Z_bipm(p, T, xw):
    """Compressibility of moist air, BIPM 1981/91 (Ciddor 1996 eq. 12)."""
    t = T - 273.15
    a0, a1, a2 = 1.58123e-6, -2.9331e-8, 1.1043e-10
    b0, b1 = 5.707e-6, -2.051e-8
    c0, c1 = 1.9898e-4, -2.376e-6
    d, e = 1.83e-11, -0.765e-8
    return (1 - (p / T) * (a0 + a1 * t + a2 * t ** 2 + (b0 + b1 * t) * xw + (c0 + c1 * t) * xw ** 2)
            + (p / T) ** 2 * (d + e * xw ** 2))


def n_ciddor1996(wl_vac_A, t_C=15.0, p_Pa=101325.0, rh=0.0, xc=450.0, pw_Pa=None, variant='ciddor'):
    """Full Ciddor (1996) phase refractive index, 'n_prop' (eq. 5), App. A/B procedure.
    rh: fractional humidity 0..1 (ignored if pw_Pa given).  variant='ciddor': R = 8.314510, Davis svp
    (ice formula below 0 C); variant='nist': NIST EMToolBox (R = 8.314472, IAPWS svp, ice < 0 C)."""
    s2 = _s2(wl_vac_A)
    T = t_C + 273.15
    if variant == 'nist':
        R = 8.314472
        svp = svp_iapws(T) if t_C >= 0 else svp_iapws_ice(T)
    else:
        R = 8.314510
        svp = svp_davis1992(T) if t_C >= 0 else svp_ice_marti1993(T)
    f = 1.00062 + 3.14e-8 * p_Pa + 5.6e-7 * t_C ** 2
    pw = rh * svp if pw_Pa is None else pw_Pa
    xw = f * pw / p_Pa
    nas = 1e-8 * (5792105.0 / (238.0185 - s2) + 167917.0 / (57.362 - s2))
    naxs = nas * (1 + 0.534e-6 * (xc - 450.0))
    w0, w1, w2, w3 = 295.235, 2.6422, -0.032380, 0.004028
    nws = 1.022e-8 * (w0 + w1 * s2 + w2 * s2 ** 2 + w3 * s2 ** 3)
    Ma = 1e-3 * (28.9635 + 12.011e-6 * (xc - 400.0))
    Mw = 0.018015
    Za = _Z_bipm(101325.0, 288.15, 0.0)
    Zw = _Z_bipm(1333.0, 293.15, 1.0)
    rho_axs = 101325.0 * Ma / (Za * R * 288.15)
    rho_ws = 1333.0 * Mw / (Zw * R * 293.15)
    Z = _Z_bipm(p_Pa, T, xw)
    rho_a = p_Pa * Ma * (1 - xw) / (Z * R * T)
    rho_w = p_Pa * Mw * xw / (Z * R * T)
    return 1.0 + (rho_a / rho_axs) * naxs + (rho_w / rho_ws) * nws


def n_edlen1966(wl_vac_A, t_C=15.0, p_Pa=101325.0, f_Pa=0.0, x_co2=0.0003):
    """Edlen (1966): eq.(12) (n-1)_tp = (n-1)_s p[1 + p(0.817-0.0133 t)1e-6]/(720.775(1+0.0036610 t)) (p in Torr);
    eq.(17) CO2: (n-1)_x = [1 + 0.540(x-0.0003)](n-1)_s ; eq.(22) water: n_tpf - n_tp = -f(5.722-0.0457 s^2)1e-8 (f in Torr)."""
    s2 = _s2(wl_vac_A)
    p = p_Pa / TORR_PA
    f = f_Pa / TORR_PA
    ns = nm1_edlen1966(s2) * (1 + 0.540 * (x_co2 - 0.0003))
    ntp = ns * p * (1 + p * (0.817 - 0.0133 * t_C) * 1e-6) / (720.775 * (1 + 0.0036610 * t_C))
    return 1.0 + ntp - f * (5.722 - 0.0457 * s2) * 1e-8


def n_birchdowns(wl_vac_A, t_C=15.0, p_Pa=101325.0, pv_Pa=0.0, version=1994, nist_water=False, D=96095.43):
    """Birch & Downs (1993, corrected 1994): (n-1)_tp = p (n-1)_s [1 + 1e-8 (0.601 - 0.00972 t) p] / (D (1 + 0.0036610 t)),
    D = 96095.43 (Pa); n_tpf - n_tp = -f (3.7345 - 0.0401 s^2) 1e-10 (f in Pa).  450 ppm CO2 fixed.
    nist_water=True multiplies the water term by 292.75/(t+273.15) (NIST EMToolBox eq. A49)."""
    s2 = _s2(wl_vac_A)
    ns = nm1_birchdowns1994(s2) if version == 1994 else nm1_birchdowns1993(s2)
    X = (1 + 1e-8 * (0.601 - 0.00972 * t_C) * p_Pa) / (1 + 0.0036610 * t_C)
    ntp = p_Pa * ns * X / D
    wt = 1e-10 * (3.7345 - 0.0401 * s2) * pv_Pa
    if nist_water:
        wt = wt * 292.75 / (t_C + 273.15)
    return 1.0 + ntp - wt


def n_nist_shopfloor(t_C, p_kPa, RH_percent):
    """NIST EMToolBox App. B shop-floor formula (633 nm): n = 1 + 7.86e-4 p/(273+t) - 1.5e-11 RH (t^2+160)."""
    return 1 + 7.86e-4 * p_kPa / (273 + t_C) - 1.5e-11 * RH_percent * (t_C ** 2 + 160)


def vac_to_air_ciddor1996_full(wl_vac_A, **kw):
    wv = _a(wl_vac_A)
    return wv / n_ciddor1996(wv, **kw)


def air_to_vac_ciddor1996_full(wl_air_A, **kw):
    wa = _a(wl_air_A)
    wv = wa * n_ciddor1996(wa, **kw)
    for _ in range(_ITER_MAX):
        new = wa * n_ciddor1996(wv, **kw)
        if np.all(np.abs(new - wv) < _ITER_TOL):
            return new
        wv = new
    return wv


def isa_pressure_Pa(h_m):
    """ICAO/ISA standard-atmosphere barometric formula (troposphere): p = 101325 (1 - 2.25577e-5 h)^5.25588 Pa.
    Used only to illustrate observatory-altitude pressures (derived, not measured)."""
    return 101325.0 * (1 - 2.25577e-5 * _a(h_m)) ** 5.25588


ALL_STANDARD = {
    'edlen1953': (vac_to_air_edlen1953, air_to_vac_edlen1953),
    'edlen1966': (vac_to_air_edlen1966, air_to_vac_edlen1966),
    'peckreeder1972': (vac_to_air_peckreeder1972, air_to_vac_peckreeder1972),
    'peckreeder1972_5p': (vac_to_air_peckreeder1972_5p, air_to_vac_peckreeder1972_5p),
    'birchdowns1993': (vac_to_air_birchdowns1993, air_to_vac_birchdowns1993),
    'birchdowns1994': (vac_to_air_birchdowns1994, air_to_vac_birchdowns1994),
    'morton2000': (vac_to_air_morton2000, air_to_vac_morton2000),
    'ciddor1996': (vac_to_air_ciddor1996, air_to_vac_ciddor1996),
    'greisen2006': (vac_to_air_greisen2006, air_to_vac_greisen2006),
    'wcslib': (vac_to_air_wcslib, air_to_vac_wcslib),
    'sdss_dr7': (vac_to_air_sdss_dr7, air_to_vac_sdss_dr7),
    'vald': (vac_to_air_vald, air_to_vac_vald),
}
