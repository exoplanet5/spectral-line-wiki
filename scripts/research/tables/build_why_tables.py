#!/usr/bin/env python
"""Build data/build/st_spectrographs.json and data/build/st_effects.json for section 3 of the Spectral Line Wiki.

Inputs: data/curated/spectrograph_precision.csv, data/notes/precision_rationale_numbers.json (numbers are
copied or converted from these and from the notes cited row by row; nothing is typed that is not in the notes).
"""
import csv, html, json, math, os, sys

ROOT = sys.argv[1]
C = 299792458.0
OUT = os.path.join(ROOT, 'data', 'build')


def sig(x, n=3):
    if x == 0:
        return 0
    d = n - 1 - int(math.floor(math.log10(abs(x))))
    v = round(x, d)
    return int(v) if d <= 0 else v


def cite(key, text):
    return f'<a class="cite" href="#ref-{key}">{text}</a>'


def ext(url, text):
    return f'<a href="{html.escape(url, quote=True)}" target="_blank" rel="noopener">{text}</a>'


# ------------------------------------------------------------------ spectrographs
rows_in = []
with open(os.path.join(ROOT, 'data', 'curated', 'spectrograph_precision.csv'), encoding='utf-8') as f:
    lines = [l for l in f if not l.startswith('#')]
for r in csv.DictReader(lines):
    rows_in.append(r)

# citation cell per instrument (paper key + the URL(s) of the CSV row)
CITE = [
    ('HARPS', [('c', 'mayor2003', 'Mayor et al. 2003'), ('u', 1, 'ESO')]),
    ('ESPRESSO', [('c', 'pepe2021', 'Pepe et al. 2021'), ('u', 1, 'ESO')]),
    ('EXPRES', [('c', 'blackman2020', 'Blackman et al. 2020'), ('u', 0, 'arXiv')]),
    ('NEID', [('u', 0, 'WIYN'), ('u', 1, 'NEID FAQ'), ('c', 'schwab2016', 'Schwab et al. 2016')]),
    ('PEPSI', [('c', 'strassmeier2015', 'Strassmeier et al. 2015'), ('u', 0, 'arXiv')]),
    ('UVES', [('u', 0, 'ESO UVES'), ('c', 'dekker2000', 'Dekker et al. 2000')]),
    ('Keck HIRES', [('u', 0, 'Keck slits'), ('u', 1, 'binning'), ('c', 'vogt1994', 'Vogt et al. 1994')]),
    ('Subaru HDS', [('u', 0, 'Subaru HDS'), ('c', 'noguchi2002', 'Noguchi et al. 2002')]),
    ('CRIRES+', [('c', 'dorn2023', 'Dorn et al. 2023'), ('u', 1, 'ESO')]),
    ('IGRINS', [('u', 0, 'Gemini'), ('c', 'le2015', 'Le et al. 2015'), ('c', 'park2014', 'Park et al. 2014')]),
    ('APOGEE', [('c', 'wilson2019', 'Wilson et al. 2019'), ('u', 0, 'arXiv')]),
    ('Gaia RVS', [('c', 'cropper2018', 'Cropper et al. 2018'), ('u', 0, 'arXiv')]),
    ('DESI', [('c', 'desi2022', 'DESI Collaboration 2022'), ('u', 0, 'arXiv')]),
    ('SDSS/BOSS', [('c', 'smee2013', 'Smee et al. 2013'), ('u', 0, 'arXiv')]),
    ('LAMOST LRS', [('u', 0, 'LAMOST DR11'), ('c', 'cui2012', 'Cui et al. 2012')]),
    ('LAMOST MRS', [('c', 'liu2020', 'Liu et al. 2020'), ('u', 0, 'arXiv')]),
    ('4MOST', [('c', 'dejong2019', 'de Jong et al. 2019'), ('u', 0, 'arXiv')]),
    ('Subaru PFS', [('u', 0, 'PFS parameters'), ('c', 'tamura2016', 'Tamura et al. 2016')]),
    ('WEAVE', [('c', 'jin2024', 'Jin et al. 2024'), ('u', 0, 'arXiv')]),
]


def cite_cell(r):
    urls = [u.strip() for u in r['url'].split(';') if u.strip()]
    for pre, parts in CITE:
        if r['instrument'].startswith(pre):
            out = []
            for kind, a, text in parts:
                if kind == 'c':
                    out.append(cite(a, text))
                elif a < len(urls):
                    out.append(ext(urls[a], text + ' ↗'))
            return ' '.join(out)
    raise SystemExit('no citation mapping for ' + r['instrument'])


def basis_cell(r):
    """Sampling / R provenance flag, with the CSV derivation note as tooltip."""
    dn = (r['derivation_notes'] or '').strip()
    low = dn.lower()
    parts = []
    if 'mid-point' in low:
        parts.append('R mid-range')
    if r['FWHM_flag'].startswith('lit'):
        parts.append('FWHM quoted')
    if 'assumed' in low and 'sampling' in low:
        parts.append('3 px assumed')
    elif 'sampling derived' in low or 'sampling from' in low and 'etc' not in low:
        parts.append('sampling derived')
    elif 'sampling from the etc' in low:
        parts.append('sampling quoted (ETC)')
    else:
        parts.append('sampling quoted')
    if 'lower bound' in low:
        parts.insert(0, 'R lower bound')
    tip = 'source: ' + r['source_values'] + (' | derivation: ' + dn if dn else '')
    return f'<span class="small" title="{html.escape(tip, quote=True)}">{", ".join(parts)}</span>'


spec_rows = []
for r in rows_in:
    spec_rows.append([
        r['instrument'],
        int(r['R']),
        int(float(r['ref_wavelength_A'])),
        sig(float(r['FWHM_A']), 3),
        sig(float(r['FWHM_kms']), 3),
        float(r['pixel_sampling_px_per_FWHM']),
        basis_cell(r),
        sig(float(r['centroid_mA_SNR100_rule']), 2 if float(r['centroid_mA_SNR100_rule']) < 1 else 3),
        sig(float(r['centroid_ms_SNR100_rule']), 3),
        sig(float(r['centroid_mA_SNR100_CRB_abs_depth0p5']), 2 if float(r['centroid_mA_SNR100_CRB_abs_depth0p5']) < 1 else 3),
        cite_cell(r),
    ])

spec_caption = (
    'Resolution element and single-line centroid precision of 31 spectrograph modes at SNR = 100 per pixel. '
    'FWHM = λ<sub>ref</sub>/R unless quoted by the source (PFS); FWHM in km/s = c/R. '
    'σ<sub>rule</sub> = FWHM/(SNR·√N<sub>pix</sub>), N<sub>pix</sub> = pixels per FWHM; '
    'σ<sub>CRB</sub> = photon-noise Cramér–Rao bound for one unresolved Gaussian absorption line of depth 0.5, '
    'integrated over pixels (' + cite('bouchy2001', 'Bouchy et al. 2001') + ', eqs. 8–12, σ<sub>D</sub> = 0). '
    'For emission lines (peak 10⁴ e⁻/pixel, no background) the bound is ≈ 2.4× smaller than the rule. '
    'Both scale as λ/(R·SNR·√N<sub>pix</sub>): multiply by 100/SNR for other signal levels. '
    '<b>Photon noise only:</b> read noise, blends, intrinsic width and asymmetry, LSF mismatch, flat-field and pixel-size errors '
    '(HARPS floor ≈ 0.01 px = 0.15 mÅ, ' + cite('lovis2007', 'Lovis &amp; Pepe 2007') + ') and wavelength-solution errors are excluded, '
    'so real per-line accuracies are worse. R conventions differ (design vs measured, e.g. HARPS 115 000 vs 110 000); '
    'DESI R is the mid-point of each quoted range; LAMOST-MRS and WEAVE samplings (3 px) are assumptions; '
    'UVES, HIRES and CRIRES+ samplings are slit-width/pixel-scale estimates and LAMOST-LRS is sampled on its 10⁻⁴-dex output grid '
    '(hover the basis column for each derivation). Verbatim source values: <code>data/curated/spectrograph_precision.csv</code>.'
)

spec = {
    'host': '#specTable',
    'caption': spec_caption,
    'tall': False,
    'search': False,
    'cols': [
        {'k': 'inst', 'label': 'Instrument / mode', 'num': False, 'html': False},
        {'k': 'R', 'label': 'R', 'num': True, 'html': False, 'title': 'resolving power λ/FWHM as quoted by the source'},
        {'k': 'lam', 'label': 'λ<sub>ref</sub> (Å)', 'num': True, 'html': False, 'title': 'reference wavelength of the row'},
        {'k': 'fw', 'label': 'FWHM (Å)', 'num': True, 'html': False, 'title': 'λ_ref/R unless quoted'},
        {'k': 'fwv', 'label': 'FWHM (km/s)', 'num': True, 'html': False, 'title': 'c/R (or c·FWHM/λ_ref)'},
        {'k': 'spx', 'label': 'px / FWHM', 'num': True, 'html': False, 'title': 'pixels per resolution element'},
        {'k': 'basis', 'label': 'basis', 'num': False, 'html': True, 'title': 'provenance of R, FWHM and sampling; hover for the derivation'},
        {'k': 'cmA', 'label': 'σ<sub>rule</sub> (mÅ)', 'num': True, 'html': False, 'title': 'FWHM/(SNR·√N_pix) at SNR 100 per pixel'},
        {'k': 'cms', 'label': 'σ<sub>rule</sub> (m/s)', 'num': True, 'html': False, 'title': 'velocity equivalent of σ_rule at λ_ref'},
        {'k': 'crb', 'label': 'σ<sub>CRB</sub> abs. (mÅ)', 'num': True, 'html': False, 'title': 'photon-noise Cramér–Rao bound, 50 %-deep unresolved Gaussian absorption line, continuum SNR 100 per pixel'},
        {'k': 'cite', 'label': 'Source', 'num': False, 'html': True},
    ],
    'rows': spec_rows,
}

# ------------------------------------------------------------------ effects
N = json.load(open(os.path.join(ROOT, 'data', 'notes', 'precision_rationale_numbers.json')))


def ma(v, lam=5000.0):          # m/s -> mA at lam
    return v / C * lam * 1e3


def ms(dl_mA, lam=5000.0):      # mA at lam -> m/s
    return dl_mA * 1e-3 / lam * C


eff = {e['effect']: e for e in N['sec3_effects']}
iso = {e['what']: e for e in N['sec3_isotope_hfs']}
av = N['sec3_airvac_sensitivity_5000A']
cat = N['sec4_catalogue']
hb = cat['handbook_vs_asd']['by_handbook_decimals']
avf = {r['lam_vac_A']: r for r in N['sec3_airvac_formulas']}

D = 'derived'
park = cite('park2021', 'DE440s')
rows = []


def rnd(x):
    return int(round(x)) if abs(x) >= 1000 else sig(x, 3)


def add(effect, lam, mA, v, sign, note, source):
    rows.append([effect, lam, rnd(mA), rnd(v), sign, note, source])


# medium
vma = avf[5000]['vac_minus_air_A']
add('Air vs vacuum wavelength (medium confusion)', 5000, vma * 1e3, ms(vma * 1e3),
    'vac − air > 0; fixed convention',
    f'0.647 Å at 2000 Å … 6.00 Å at 2.2 µm in standard air, i.e. {ms(avf[22000]["vac_minus_air_A"]*1e3, 22000)/1e3:.0f}–{ms(avf[2000]["vac_minus_air_A"]*1e3, 2000)/1e3:.0f} km/s; 1000× the 1 mÅ target',
    D + ' (' + cite('ciddor1996', 'Ciddor 1996') + ')')
# observer motion
add('Earth orbital velocity', 5000, ma(29784.69), 29784.69,
    '± annual; projected ≈ v cos β (ecliptic latitude β)',
    'circular √(GM☉/au); DE440s 2026: 29.294–30.278 km/s = 489–505 mÅ',
    D + ' (' + park + ')')
add('Earth rotation (equator, ω⊕a⊕)', 5000, 7.757, 465.10,
    '± diurnal; ∝ cos φ′ cos δ sin H',
    '438 m/s at Mauna Kea; changes by up to 3.4 cm/s per second, so the photon-weighted exposure mid-time matters',
    D + ' (' + cite('petit2010', 'IERS TN36') + '; ' + cite('wright2014', 'Wright &amp; Eastman 2014') + ')')
add('Earth about the Earth–Moon barycentre (lunar reflex)', 5000, ma(12.43), 12.43,
    '± monthly; 11.7–13.4 m/s',
    'lost if the EMB orbit is used for the observer velocity',
    D + ' (' + park + ')')
add('Heliocentric vs barycentric frame (Sun about the SSB)', 5000, ma(16.08), 16.08,
    '≤ 16.1 m/s; 8.5–16.1 m/s over 1990–2040; years (giant-planet reflex, Jupiter alone 12.4 m/s)',
    'IAU 2000 Resolution C1 specifies the barycentre, not the Sun',
    D + ' (' + park + '; ' + cite('lindegren2003', 'Lindegren &amp; Dravins 2003') + ')')
add('Earth-induced reflex of the Sun (Earth-twin signal)', 5000, ma(0.09), 0.09,
    '± yearly', 'the EPRV benchmark', cite('fischer2016', 'Fischer et al. 2016'))
# relativity
add('Solar gravitational redshift, surface → ∞', 5000, 10.6125, 636.31,
    '+ (red); constant', 'GM☉/(R☉c) with IAU 2015 nominal values; 636.486 ± 0.024 m/s with the older R☉',
    D + ' (' + cite('lindegren2003', 'Lindegren &amp; Dravins 2003') + ')')
add('Solar gravitational redshift, received at 1 au', 5000, 10.563, 633.35,
    '+ (red); constant',
    'GM☉/c·(1/R☉ − 1/au); measured 639 ± 14 m/s (15 strong Fe I lines, HARPS-LFC Moon spectra)',
    D + '; ' + cite('gonzalezhernandez2020', 'González Hernández et al. 2020'))
add('Solar convective blueshift (Fe I, depth-dependent)', 5000, 9.123, 547,
    '− (blue): −547 m/s at relative depth 0.1 → +2 m/s at 0.9; σ 128–303 m/s per bin',
    'partly cancels the gravitational redshift; Fe I 5250.2 disk-centre bisector −100 to −290 m/s (−1.8 to −5.1 mÅ)',
    cite('reiners2016', 'Reiners et al. 2016') + '; ' + cite('lohnerbottcher2019', 'Löhner-Böttcher et al. 2019'))
add('Relativistic observer-clock term', 5000, ma(4.648), 4.648,
    '+ constant offset; 0.20 m/s peak-to-peak',
    '= L_B·c; 3.0 (solar potential) + 0.2 (Earth potential) + 1.5 (orbital speed) m/s; included in the §7.4 calculator',
    D + ' (' + park + '); ' + cite('lindegren2003', 'Lindegren &amp; Dravins 2003'))
add('Additive instead of multiplicative barycentric correction', 5000, ma(3.0), 3.0,
    '∝ v_meas·v_B/c', 'cross term; needed at the 1 m/s level',
    cite('wright2014', 'Wright &amp; Eastman 2014'))
add('Transverse Doppler, star with v_t = 100 km/s', 5000, ma(16.7), 16.7,
    '+ (red); v_t²/2c', 'Greisen et al. quote 15 m/s for the same case',
    D + '; ' + cite('greisen2006', 'Greisen et al. 2006'))
add('Solar Shapiro delay (night-time)', 5000, ma(0.0067), 0.0067,
    'varies with solar elongation: 0.59 mm/s at 90°, 6.7 mm/s at 10°',
    '≤ 1.5 mm/s for τ Cet; omitted by astropy and by the §7.4 calculator',
    D + '; ' + cite('wright2014', 'Wright &amp; Eastman 2014'))
add('Radio vs optical velocity convention at z = 0.001', 5000, ma(299), 299,
    'grows as cz²/c: 29.7 km/s at z = 0.01', 'quote VOPT, VRAD or VELO with every velocity',
    D + ' (' + cite('greisen2006', 'Greisen et al. 2006') + ')')
# air
spread5000 = max(avf[5000]['diff_vs_ciddor_mA'].values()) - min(avf[5000]['diff_vs_ciddor_mA'].values())
add('Choice of standard-air formula', 5000, spread5000, ms(spread5000),
    'fixed per formula; 0.04–0.21 mÅ over 3000 Å–1 µm, ≈ 0.6 mÅ at 2000 Å and 2.2 µm',
    'spread of Edlén 1953/1966, Peck & Reeder eq. 2/3, Birch & Downs 1994, Ciddor 1996; NIST ASD = Peck & Reeder eq. 3 (0.02 mÅ rms)',
    D + ' (' + cite('peck1972', 'Peck &amp; Reeder 1972') + '; ' + cite('ciddor1996', 'Ciddor 1996') + ')')
add('Non-standard “standard air” in software', 5000, 76.85, ms(76.85),
    '−77 mÅ: FITS-WCS eq. 65 (0 °C air); −1.2 mÅ: SDSS DR7 web polynomial; +1.2 mÅ: IRAF disptrans default (4 mmHg H₂O)',
    'check which air a pipeline assumes before comparing at 1 mÅ',
    D + ' (' + cite('greisen2006', 'Greisen et al. 2006') + ')')
s = av['standard_15C_101325Pa_dry_450ppm']
add('Air temperature (per K)', 5000, s['dmA_per_K'], 290.97,
    '+ for warmer air', 'Ciddor full model at standard air; 0.2 K error ≈ 1 mÅ',
    D + ' (' + cite('ciddor1996', 'Ciddor 1996') + ')')
add('Air pressure (per hPa)', 5000, abs(s['dmA_per_hPa']), 82.53,
    '− for higher pressure', 'Ciddor full model at standard air',
    D + ' (' + cite('ciddor1996', 'Ciddor 1996') + ')')
add('Humidity (per 10 % RH)', 5000, s['dmA_per_10pctRH'], 18.75,
    '+ for moister air', 'Ciddor full model at standard air (+0.031 mÅ per % RH)',
    D + ' (' + cite('ciddor1996', 'Ciddor 1996') + ')')
add('CO₂ (per 100 ppm)', 5000, abs(s['dmA_per_100ppm_CO2']), 4.464,
    '− for more CO₂', 'standard air is defined at 450 ppm (Ciddor) or ≈ 300 ppm (Edlén)',
    D + ' (' + cite('ciddor1996', 'Ciddor 1996') + ')')
add('Laboratory air vs standard air', 5000, av['lab_minus_standard_mA'], ms(av['lab_minus_standard_mA']),
    '+; 20 °C, 101 325 Pa, 50 % RH, 420 ppm vs 15 °C, dry, 450 ppm',
    'an “air wavelength” is defined to 1 mÅ only for standard air',
    D + ' (' + cite('ciddor1996', 'Ciddor 1996') + ')')
add('Mountain-site air vs standard air (2600 m, 0 °C)', 5000, 324.06, ms(324.06),
    '+; 30 % RH, ISA pressure',
    'absorbed by a wavelength solution fitted to lamp lines in the same air (≤ 0.02 mÅ residual, 3800–6800 Å); not absorbed: air changes between arc and science exposure',
    D + ' (' + cite('ciddor1996', 'Ciddor 1996') + ')')
# atomic structure
def isorow(name, key, lam, sign, note, src):
    e = iso[key]
    add(name, lam, abs(e['split_mA']), abs(e['v_m_s']), sign, note, src)


mb = cite('murphy2014', 'Murphy &amp; Berengut 2014')
isorow('Isotope shift Ca II K: ⁴⁰Ca → ⁴⁸Ca', 'Ca II 3934 (K): 40Ca -> 48Ca isotope shift', 3934.77,
       '− (heavier isotope bluer)', 'vacuum; a single Ca II K wavelength assumes the terrestrial isotopic mix', mb)
isorow('Isotope shift Mg II 2796: ²⁴Mg → ²⁶Mg', 'Mg II 2796: 24Mg -> 26Mg isotope shift', 2796.35,
       '− (heavier isotope bluer)', 'vacuum; component wavelengths known to 4 × 10⁻⁷ Å', mb)
isorow('Isotope shift Li I 6708 D2: ⁷Li → ⁶Li', 'Li I 6708 D2: 7Li -> 6Li isotope shift (cog)', 6709.62,
       '+ (⁶Li redder)', '≈ the D1–D2 splitting (≈ 151 mÅ): ⁶Li D2 falls within ≈ 10 mÅ of ⁷Li D1; ⁷Li ground hfs 12.2 mÅ; NIST ASD gives 6707.76 ± 0.06 Å (air)',
       D + ' (' + cite('brown2013', 'Brown et al. 2013') + '); ' + cite('kramida2024', 'NIST ASD'))
isorow('Hyperfine structure Na I D2 (ground state)', 'Na I 5891 (D2): 23Na ground-state hfs F=1 vs F=2', 5891.58,
       '± two components', 'vacuum; unresolved in most astrophysical spectra, so the composite depends on the weighting', mb)
isorow('Hyperfine structure Mn II 2576 (extreme components)', 'Mn II 2576: 55Mn hfs, extreme components', 2576.87,
       '± 14–16 components', 'vacuum; hfs-sensitive lines need component lists, not one wavelength', mb)
isorow('Hyperfine structure ⁴³Ca II K (F = 3 vs 4)', 'Ca II 3934 (K): 43Ca hyperfine F=3 vs F=4', 3934.77,
       '± two components', 'vacuum; odd isotope', mb)
isorow('Hyperfine structure ²⁵Mg II 2796 (F = 2 vs 3)', 'Mg II 2796: 25Mg hyperfine components F=2 -> F=3', 2796.35,
       '± two components', 'vacuum', mb)
add('Hα fine structure (7 components)', 6562.8, 199.74, 9124,
    'span of components; composite depends on weighting',
    'composite Ritz 6562.819 ± 0.007 Å vs observed 6562.79 ± 0.03 Å (air); H → D shift −1.807 Å',
    cite('kramida2024', 'NIST ASD'))
add('Zeeman σ components, B = 1 kG, g = 1', 5000, 11.672, 699.8,
    '± symmetric; ∝ g B λ²',
    '1.17 mÅ at 100 G; Fe I 5250.2 (g = 3) 38.6 mÅ at 1 kG; 120 mÅ at 1.6 µm',
    D + '; g from ' + cite('lohnerbottcher2019', 'Löhner-Böttcher et al. 2019'))
# telluric
hit = cite('gordon2026', 'HITRAN2024')
add('Telluric O₂ pressure shift (sea-level path)', 7600, 2.09, 81,
    '+ median; scales with effective pressure',
    'zero-pressure HITRAN positions must be shifted by δ_air·p before use as a 1 mÅ reference',
    D + ' from δ_air (' + hit + ')')
add('Telluric CO₂ pressure shift (sea-level path)', 20000, 11.0, 177,
    '+ median', 'CO +13.1 mÅ (168 m/s), N₂O +9.9 mÅ (135 m/s)',
    D + ' from δ_air (' + hit + ')')
add('Telluric H₂O pressure shift (sea-level path)', 13500, 12.4, 331,
    '+ median; effective pressure ≈ 0.81 atm', 'at 1 atm, all molecules: median +17.8 mÅ, range −62 to +231 mÅ',
    D + ' from δ_air (' + hit + ')')
add('Telluric CH₄ pressure shift (sea-level path)', 22900, 22.5, 303,
    '+ median', 'largest median shift of the selected molecules',
    D + ' from δ_air (' + hit + ')')
add('Atmospheric winds (telluric line shifts)', 5000, ma(40), 40,
    '± line-of-sight projection changes with azimuth and elevation; fitted 15–49 m/s per night',
    'high-altitude wind ≈ 40 m/s over La Silla; telluric O₂ stable to 10 m/s rms over 6 yr, 1–2 m/s after modelling',
    cite('figueira2010', 'Figueira et al. 2010'))
# laboratory and instrumental
add('Lab pressure shifts: Ar lines in Th–Ar lamps', 5500, ma(50, 5500), 50,
    'lamp-to-lamp, “a few tens of m/s” (10–50 m/s = 0.18–0.92 mÅ)',
    'Th lines stable within errors; avoid Ar lines for high precision',
    cite('lovis2007', 'Lovis &amp; Pepe 2007'))
add('Lab-to-lab: Fe I 5250.2084 (LARS vs NIST ASD)', 5250.2, 0.5, 29,
    '0.2–0.5 mÅ (11–29 m/s)', 'LFC-calibrated hollow cathode vs ASD observed/Ritz values',
    cite('lohnerbottcher2019', 'Löhner-Böttcher et al. 2019') + '; ' + cite('kramida2024', 'NIST ASD'))
add('HARPS Th–Ar line list (per line)', 5500, ma(10, 5500), 10,
    'random', 'CCD systematics floor 0.01 px = 8 m/s (0.15 mÅ); strongest lines 1.6 m/s photon noise',
    cite('lovis2007', 'Lovis &amp; Pepe 2007'))
add('Kitt Peak FTS solar atlases: wavelength-scale errors', 5000, ma(100), 100,
    'up to; smooth in λ', 'IAG flux atlas accurate to ±10 m/s (0.17 mÅ) over 405–1065 nm',
    cite('reiners2016', 'Reiners et al. 2016'))
add('UVES night-sky atlas: observed line positions', 7000, 17, ms(17, 7000),
    'single-line 1σ; red settings offset +8 to +15 mÅ vs lab/HITRAN',
    'printed to 0.1 mÅ; UVES − NIST = +13 ± 2 mÅ for 13 atomic lines',
    cite('hanuschik2003', 'Hanuschik 2003') + '; ' + cite('cosby2006', 'Cosby et al. 2006'))
# catalogue rounding
r2 = cat['rounding_2dec']; r3 = cat['rounding_3dec']
add('Rounding to 0.01 Å (uniform, rms)', 5000, r2['rms_err_mA'], r2['rms_m_s_5000A'],
    'max ±5 mÅ (±300 m/s); rms = step/√12', '43 % of the NIST Handbook lines are given to 2 decimals', D)
add('NIST Handbook 2-decimal values vs NIST ASD', None, hb['2']['median_abs_dl_mA'], hb['2']['median_abs_dv_m_s'],
    f'median |Δ|; 68 %: {hb["2"]["p68_abs_dl_mA"]:.1f} mÅ; 95 %: {hb["2"]["p95_abs_dl_mA"]:.0f} mÅ ({hb["2"]["p95_abs_dv_m_s"]/1e3:.2f} km/s)',
    f'{hb["2"]["n_asd_more_decimals"]:,} lines where ASD gives more decimals'.replace(',', ' ') + '; size (m/s) is the median over each line’s own λ',
    D + ' (' + cite('sansonetti2005', 'Sansonetti &amp; Martin 2005') + '; ' + cite('kramida2024', 'NIST ASD') + ')')
add('Rounding to 0.001 Å (uniform, rms)', 5000, r3['rms_err_mA'], r3['rms_m_s_5000A'],
    'max ±0.5 mÅ (±30 m/s)', 'the display precision of this page', D)
add('NIST Handbook 3-decimal values vs NIST ASD', None, hb['3']['median_abs_dl_mA'], hb['3']['median_abs_dv_m_s'],
    f'median |Δ|; 68 %: {hb["3"]["p68_abs_dl_mA"]:.1f} mÅ; 95 %: {hb["3"]["p95_abs_dl_mA"]:.1f} mÅ',
    f'{hb["3"]["n_asd_more_decimals"]:,} lines'.replace(',', ' ') + '; size (m/s) is the median over each line’s own λ',
    D + ' (' + cite('sansonetti2005', 'Sansonetti &amp; Martin 2005') + '; ' + cite('kramida2024', 'NIST ASD') + ')')

eff_caption = (
    'Size of physical and instrumental effects, as the magnitude |Δλ| (sign and range in the third column). '
    'Velocity-like effects are converted at λ = 5000 Å (Δλ = λ·v/c, 1 mÅ ≙ 59.96 m/s); effects tied to a line or band '
    '(isotope shifts, hyperfine and fine structure, telluric pressure shifts, UVES) are given at that wavelength, stated in the λ column '
    '(for telluric molecules the median λ of the selected lines: O₂ 0.58–1.30 µm, CO₂ 1.20–2.17 µm, H₂O 0.44–2.50 µm, CH₄ 1.13–2.50 µm). '
    'Telluric pressure shifts are medians of δ<sub>air</sub>·p<sub>eff</sub> over the selected lines, with p<sub>eff</sub> the S·u-weighted '
    'pressure along a vertical sea-level path in the US Standard atmosphere; they neglect line mixing, speed dependence and the temperature '
    'dependence of δ. Air-state sensitivities use the full Ciddor (1996) model at λ<sub>vac</sub> = 5000 Å. '
    '“derived” marks values computed for this page (DE440s, IAU 2015 nominal constants, HITRAN2024 parameters, NIST Handbook/ASD comparison); '
    'the derivations are in <code>data/notes/precision_rationale.md</code> and <code>relativity_cosmology.md</code>. '
    'Rows below ≈ 1 mÅ are kept because they decide what the barycentric correction (§7.4) must contain.'
)

effects = {
    'host': '#effTable',
    'caption': eff_caption,
    'tall': True,
    'search': True,
    'cols': [
        {'k': 'eff', 'label': 'Effect', 'num': False, 'html': False},
        {'k': 'lam', 'label': 'λ (Å)', 'num': True, 'html': False, 'title': 'wavelength at which the mÅ size is quoted'},
        {'k': 'mA', 'label': 'size (mÅ)', 'num': True, 'html': False, 'title': 'magnitude |Δλ| at the λ of the row'},
        {'k': 'ms', 'label': 'size (m/s)', 'num': True, 'html': False, 'title': 'velocity equivalent c·Δλ/λ'},
        {'k': 'sign', 'label': 'sign / variability', 'num': False, 'html': False},
        {'k': 'note', 'label': 'note', 'num': False, 'html': False},
        {'k': 'src', 'label': 'source', 'num': False, 'html': True},
    ],
    'rows': rows,
}

for name, obj in (('st_spectrographs.json', spec), ('st_effects.json', effects)):
    p = os.path.join(OUT, name)
    with open(p, 'w', encoding='utf-8') as f:
        json.dump(obj, f, ensure_ascii=False, indent=1, allow_nan=False)
    print('wrote', p, len(obj['rows']), 'rows')
