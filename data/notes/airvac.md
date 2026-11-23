# Air ↔ vacuum wavelength conversion — research notes (target: 1e-4 Å)

Compiled 2026-11-23 for the Spectral Line Wiki. Everything numeric here was copied from a primary source
(file kept under `data/raw/airvac/`) or computed by the scripts in `scripts/research/` (marked **derived**).

* Code: `scripts/research/airvac_formulas.py` (all formulas; exact iterative inverses),
  `validate_airvac.py` (304 checks vs. published values), `compare_airvac.py` (comparison table + sensitivities),
  `fetch_asd_airvac_check.py` (NIST ASD downloads used to reverse-engineer ASD's conversion).
* Outputs: `data/curated/airvac_comparison.csv` (2000–25000 Å, step 100 Å, differences in mÅ),
  `data/curated/airvac_validation.csv` (+ `_asd.csv`), `data/curated/airvac_sensitivity_ciddor_5000A.csv`.
* Re-run: `python -I scripts/research/validate_airvac.py data/raw/airvac/specutils_pkg/extracted/specutils/utils/wcs_utils.py data/raw/airvac/nist_asd_checks data/curated/airvac_validation.csv`
  and `python -I scripts/research/compare_airvac.py data/curated/airvac_comparison.csv data/curated/airvac_sensitivity_ciddor_5000A.csv`.

---------------------------------------------------------------------------------------------------

## 0. Bottom line

1. λ_air = λ_vac / n, with n the phase index of *standard air* evaluated (in every modern formula) at the
   **vacuum** wavenumber σ = 1/λ_vac [µm⁻¹]. λ_vac − λ_air ≈ 0.647 Å at 2000 Å, 1.394 Å at 5000 Å,
   1.813 Å at Hα, 2.741 Å at 1 µm, 5.459 Å at 2 µm (Ciddor 1996, derived) — i.e. (n−1)·c ≈ 83.6 km s⁻¹ at 5000 Å.
2. The modern standard-air formulas (Edlén 1953/1966, Peck & Reeder 1972, Birch & Downs 1994 = Morton 2000 =
   VALD = CRC, Ciddor 1996 = SDSS DR9+ = IDL astrolib) agree with each other to **≤ 0.65 mÅ** (max pairwise
   spread) over 2000–25000 Å, ≤ 0.15 mÅ over 3000–8000 Å and ≤ 0.21 mÅ over 3000–10000 Å (derived,
   `airvac_comparison.csv`). So at the 1e-3 Å
   level any of them is acceptable; **at the 1e-4 Å level the formula must be named**.
3. Formulas/implementations that are *not* acceptable at 1e-3 Å:
   * Greisen et al. 2006 eq. 65 (FITS WCS Paper III, = IUGG/IAG 1999 phase refractivity) is for **0 °C**
     (273.15 K, 1013.25 hPa, 375 ppm CO₂): air wavelengths −26 mÅ (2000 Å), **−77 mÅ (5000 Å)**, −375 mÅ (2.5 µm)
     vs. 15 °C standard air (≈ −4.6 km s⁻¹).
   * SDSS DR7/DR8/DR12 web polynomial `AIR = VAC/(1.0 + 2.735182E-4 + 131.4182/VAC^2 + 2.76249E8/VAC^4)` (a truncated
     expansion formerly in IDL VACTOAIR): −1.2 mÅ at 5000 Å, −3.1 mÅ at Hα, −21.6 mÅ at 2.5 µm.
   * IRAF `onedspec.disptrans` / `airvac.cl` defaults (t = 15 °C, p = 760 mmHg, **f = 4 mmHg water vapour**):
     +1.2 mÅ at 5000 Å, +6.4 mÅ at 2.5 µm (with f = 0 it reduces to Edlén 1953 to < 0.01 mÅ).
4. **NIST ASD** air wavelengths are reproduced by the **Peck & Reeder (1972) five-parameter formula (their eq. 3)**:
   35 988 ASD line pairs (vacuum vs. air output of the same lines, 1850–26000 Å) agree to within the ASD rounding
   in 99.99 % of cases, mean residual 0.0001 mÅ (derived, `airvac_validation_asd.csv`). Ciddor would be off by
   +0.05…+0.29 mÅ and the 4-parameter P&R eq. 2 by up to −2.9 mÅ below 2000 Å.
5. **NIST ASD `show_av` conventions (verified from the CSV headers)**: `show_av=2` (default) → vacuum < 2000 Å
   (σ ≥ 50 000 cm⁻¹), air 2000–20000 Å, vacuum > 20000 Å (σ ≤ 5000 cm⁻¹); `show_av=3` → vacuum everywhere;
   **`show_av=4` → vacuum < 1850 Å, AIR ≥ 1850 Å (including > 20000 Å)** — not "< 2000 Å vacuum" as assumed in
   `scripts/fetch_asd.py`. Lines between 1850 and 2000 Å fetched with `show_av=4` are air wavelengths (≈ 0.6 Å
   shorter than vacuum) and lines > 20000 Å are air (P&R extrapolated beyond its 1.69 µm validity).
6. Recommended wiki practice: store **vacuum** wavelengths as the primary quantity; display air wavelengths
   computed with one declared formula (Ciddor 1996 standard air for new work; P&R eq. 3 when the goal is to
   reproduce NIST-ASD air values exactly); never mix formulas inside one column; keep the source's own
   air/vac string verbatim alongside.

---------------------------------------------------------------------------------------------------

## 1. Why 1e-3 Å (and when 1e-4 Å) — numbers relevant to the wiki's precision target (derived)

* Velocity equivalent of δλ: δv = c·δλ/λ. **1 mÅ = 150 m s⁻¹ at 2000 Å, 60 m s⁻¹ at 5000 Å, 15 m s⁻¹ at 2 µm.**
* The air–vacuum shift itself is 0.65–6.8 Å (2000–25000 Å); a relative error of 7e-4 in (n−1) (e.g. 0.2 K in
  the assumed temperature, or the wrong standard-air definition) already produces 1 mÅ at 5000 Å.
* Laboratory wavelengths of good lines are given to 1e-5–1e-6 Å in vacuum: e.g. NIST ASD Na I D2
  `obs_wl_vac 5891.583264 ± 0.000015 Å`, Ritz `5891.583253 ± 0.000017 Å`; in air ASD prints `5889.95095 ± 0.00003 Å`
  (the air value carries the 5e-9 relative uncertainty of the air formula, added in quadrature).
* Spread among acceptable standard-air formulas: ≤ 0.65 mÅ (2000–25000 Å); NIST ASD's own stated relative
  uncertainty of P&R is 5×10⁻⁹ for λ > 400 nm (→ 0.025 mÅ at 5000 Å), and Birch & Downs (1993) quote 3σ = 3×10⁻⁸
  for ambient conditions (→ 0.15 mÅ at 5000 Å).
* Hence: a 1e-3 Å catalogue precision is meaningful only if the medium and the conversion formula are stated;
  1e-4 Å requires naming the exact formula *and* standard-air definition (CO₂, temperature scale).

---------------------------------------------------------------------------------------------------

## 2. "Standard air" definitions (they are not all the same!)

| Source | T | p | humidity | CO₂ | notes |
|---|---|---|---|---|---|
| Edlén 1953 (IAU 1957, Morton 1991) | 15 °C | 760 Torr | dry | 0.03 % (300 ppm) | Edlén 1966: "factor 1.000162 used by Edlén (1953) for transferring (n−1) for CO₂-free air to standard air" |
| Edlén 1966 | 15 °C (IPTS-48) | 760 Torr | dry | "Standard air is defined as containing CO₂ to 0.0003 parts by volume"; composition 78.09 % N₂, 20.95 % O₂, 0.93 % Ar, 0.03 % CO₂ (as quoted by Birch & Downs 1993) | |
| Peck & Reeder 1972 | 15 °C | 101325 Pa | dry | 300 ppm per Ciddor 1996; their Table V says data "recalculated for 15 °C and 0.033 % CO₂" | |
| NIST ASD (P&R) | 15 °C | 101 325 Pa | dry | "0.033 % CO₂" (lineshelp) | |
| Birch & Downs 1993/1994, CRC, VALD, Morton 2000 | 15 °C (ITS-90) | 101 325 Pa | dry | **450 ppm** ("0.045 % by volume") | B&D: mean of 6500 lab measurements was 444 ± 57 ppm |
| Ciddor 1996 (IAG 1999 procedure), SDSS DR9+, IDL astrolib | 15 °C (ITS-90) | 101 325 Pa | 0 % | **450 ppm** (`x_c` adjustable) | |
| IUGG/IAG 1999 closed formulae, Greisen 2006 eq. 65 | **0 °C** (273.15 K) | 1013.25 hPa | e = 0 | **375 ppm** | geodesy (EDM) convention |

Ambient CO₂ today is ~425 ppm (NOAA: Mauna Loa annual mean 2024 = 424.61 ppm; monthly mean Sept 2026 = 426.07 ppm);
in Ciddor's formula a change of 100 ppm moves λ_air at 5000 Å by 0.074 mÅ (derived), i.e. irrelevant at 1e-3 Å
but visible at 1e-4 Å. Changing the CO₂ convention 300→450 ppm moves λ_air by −0.05 (2000 Å) to −0.55 mÅ (2.5 µm).

---------------------------------------------------------------------------------------------------

## 3. The formulas (verbatim coefficients; σ = vacuum wavenumber in µm⁻¹ unless stated)

### 3.1 Edlén (1953) — the IAU standard (Oosterhoff 1957; Morton 1991)
Edlén, B. 1953, "The Dispersion of Standard Air", J. Opt. Soc. Am. 43, 339–344, doi:10.1364/JOSA.43.000339,
bibcode 1953JOSA...43..339E. Abstract (verbatim): "(n−1)10⁸ = 6432.8 + 2 949 810(146−σ²)⁻¹ + 25 540(41−σ²)⁻¹,
σ being vacuum wave number in µ⁻¹".
* Morton (1991, ApJS 77, 119, eq. 3; doi:10.1086/191601; 1991ApJS...77..119M): "The IAU standard for conversion
  between air and vacuum wavelengths is, according to Oosterhoff (1957) and Edlén (1953),
  λ_vac/λ_air = n = 1 + 6.4328×10⁻⁵ + 2.94981×10⁻²/(146 − σ²) + 2.5540×10⁻⁴/(41 − σ²), where σ = 10⁴/λ, with λ in
  angstroms. … Longward of λ_air = 2000 Å, corresponding to λ_vac = 2000.648 Å, the tables list air wavelengths".
  Oosterhoff, P. T. 1957, Trans. IAU 9, 69 & 202 (adoption; not accessed directly). Edlén (1966) records that the
  1953 formula was adopted by the Joint Commission for Spectroscopy, the Advisory Committee for the Definition of
  the Metre and IAU Commission 14. The IAU 1991 re-affirmation ("resolution No. C15, Commission 44, XXI GA") is
  cited by Allende Prieto (2011, APOGEE technical note) and by specutils; the IAU Transactions XXI B text was
  not accessed.
* Validation: Edlén's Table I "formula" column equals (this formula)/1.000162 (CO₂-free air) to ≤ 0.065×10⁻⁸
  (7 points); Morton's 2000 Å ↔ 2000.648 Å reproduced (2000.64806). Edlén (1966) states the 1953→1965 change
  |Δλ| does not exceed 0.00003 Å below 7700 Å — consistent with our 0.02–0.04 mÅ (derived).

### 3.2 Edlén (1966)
Edlén, B. 1966, "The Refractive Index of Air", Metrologia 2, 71–80, doi:10.1088/0026-1394/2/2/002,
1966Metro...2...71E (PDF in raw/airvac/papers/edlen1966.pdf).
* Dispersion (eq. 1): **(n − 1)ₛ × 10⁸ = 8342.13 + 2406030 (130 − σ²)⁻¹ + 15997 (38.9 − σ²)⁻¹.**
* T, p (eq. 12, as quoted by Birch & Downs 1993 eq. 1): (n − 1)ₜₚ = (n − 1)ₛ · p [1 + p (0.817 − 0.0133 t) × 10⁻⁶] /
  [720.775 (1 + 0.0036610 t)], p in Torr, t in °C. Simplified (eq. 15): (n−1)ₜₚ = (n−1)ₛ × 0.00138823 p/(1 + 0.003671 t),
  good to 1.0×10⁻⁸ for 15–30 °C, 700–800 Torr.
* CO₂ (eq. 17): (n − 1)ₓ = [1 + 0.540 (x − 0.0003)] (n − 1)ₛ, x = volume fraction.
* Water vapour (eq. 22): n_tpf − n_tp = −f (5.722 − 0.0457 σ²) × 10⁻⁸, f in Torr (valid near 20 °C, 760 Torr, f = 10 Torr).
* Stated uncertainty of absolute values ≈ ±5×10⁻⁸.
* Validation: P&R (1972) Table III Eq.(1) column (18 λ, ≤ 0.009×10⁻⁸); Birch & Downs 1993 Table 3 "original Edlén"
  (5 T,p cases, ≤ 0.10×10⁻⁸); B&D App. A: (n−1)ₛ(632.991 nm) = 276.5175×10⁻⁶ (300 ppm), 276.5399×10⁻⁶ (450 ppm).

### 3.3 Peck & Reeder (1972)
Peck, E. R. & Reeder, K. 1972, "Dispersion of Air", J. Opt. Soc. Am. 62, 958–962, doi:10.1364/JOSA.62.000958,
1972JOSA...62..958P. Eq. (1) in their paper = Edlén 1966; their new fits:
* **Eq. (2)** (two-term Sellmeier, 4 parameters, 0.23–1.69 µm, accuracy "a little better than 2 parts in 10⁹" per
  A. T. Young): **(n − 1) × 10⁸ = 5 791 817 / (238.0185 − σ²) + 167 909 / (57.362 − σ²)**.
* **Eq. (3)** (5 parameters, extended to 0.185 µm): **(n − 1) × 10⁸ = 8060.51 + 2 480 990 / (132.274 − σ²) +
  17 455.7 / (39.32957 − σ²)**. Coefficients taken from refractiveindex.info ("Peck" page, "n−1 = 8.06051×10⁻⁵ +
  2.480990×10⁻²/(132.274−λ⁻²) + 1.74557×10⁻⁴/(39.32957−λ⁻²)") — the original equation text was not accessible, but
  the coefficients reproduce all six Eq.(3) values of P&R Table V (0.1855–0.2145 µm) to ≤ 0.045×10⁻⁸, and they
  reproduce NIST ASD (§4.1).
* Validation of eq. (2): 18 values of Table III (≤ 0.0044×10⁻⁸) and six far-UV values of Table V.

### 3.4 Birch & Downs (1993) and correction (1994) — "modified Edlén"
Birch, K. P. & Downs, M. J. 1993, "An Updated Edlén Equation for the Refractive Index of Air", Metrologia 30,
155–162, doi:10.1088/0026-1394/30/3/004, 1993Metro..30..155B (PDF in raw/airvac/papers).
Birch, K. P. & Downs, M. J. 1994, "Correction to the Updated Edlén Equation for the Refractive Index of Air",
Metrologia 31, 315–316, doi:10.1088/0026-1394/31/4/006, 1994Metro..31..315B (abstract only: sign error in the
IPTS-48→ITS-90 factor; true factor 0.9999687, overall factor 1.0000486).
* 1993 eq. (8) (factor 1.0001099 on Edlén 1966; superseded): (n − 1)ₛ × 10⁸ = 8343.05 + 2 406 294 [130 − σ²]⁻¹ +
  15 999 [38.9 − σ²]⁻¹.
* **1994 (current)**: **(n − 1)ₛ × 10⁸ = 8342.54 + 2 406 147 (130 − σ²)⁻¹ + 15 998 (38.9 − σ²)⁻¹** (15 °C,
  101 325 Pa, dry, 450 ppm). Same constants in NIST EMToolBox (A = 8342.54, B = 2406147, C = 15998), CRC Handbook
  87th ed. p. 10-252, VALD, Morton (2000).
* T, p (1993 eq. 12, unchanged in 1994): (n − 1)ₜₚ = p (n − 1)ₛ [1 + 10⁻⁸ (0.601 − 0.00972 t) p] / [96095.43 (1 +
  0.0036610 t)], p in Pa, t in °C (ITS-90). (Ciddor 1996 notes the denominator constant should be 96094.69 for exact
  reduction at 15 °C/101325 Pa; effect 2×10⁻⁹.)
* Water (1993 eq. 14): n_tpf − n_tp = −f [3.7345 − 0.0401 σ²] × 10⁻¹⁰, f in Pa. NIST EMToolBox multiplies this by
  292.75/(t + 273.15) (their eq. A49).
* CO₂: no explicit term; "+1 ppm CO₂ increases n by 1.45×10⁻¹⁰" (at 633 nm).
* Stated: "3σ uncertainty ±3×10⁻⁸ … over the range of wavelengths from 350 nm to 650 nm".
* Validation: B&D Table 3 revised-Edlén values (5 cases, ≤ 0.05×10⁻⁸, 1993 constants); NIST EMToolBox Table 1
  modified-Edlén (14 cases, ≤ 0.04×10⁻⁸); CRC table (11 λ, ≤ 0.44×10⁻⁸ — the CRC table is rounded to 1×10⁻⁸ —
  and λ_vac−λ_air to ≤ 0.006 pm).

### 3.5 Ciddor (1996) — IAG 1999 recommended procedure
Ciddor, P. E. 1996, "Refractive index of air: new equations for the visible and near infrared", Appl. Opt. 35,
1566–1573, doi:10.1364/AO.35.001566, 1996ApOpt..35.1566C. (Group index: Ciddor & Hill 1999, Appl. Opt. 38,
1663–1667, doi:10.1364/AO.38.001663; sign error in their eq. B2 corrected by Pollinger 2020, Appl. Opt. 59, 9771.)
* Standard air (15 °C, 101 325 Pa, 0 % RH, 450 ppm): **10⁸ (n_as − 1) = k₁/(k₀ − σ²) + k₃/(k₂ − σ²)**, k₀ = 238.0185 µm⁻²,
  k₁ = 5 792 105 µm⁻², k₂ = 57.362 µm⁻², k₃ = 167 917 µm⁻² (P&R eq. 2 × ≈1.0000497: ITS-90 + 300→450 ppm).
* CO₂: **(n_axs − 1) = (n_as − 1)[1 + 0.534×10⁻⁶ (x_c − 450)]**, x_c in ppm.
* Water vapour (20 °C, 1333 Pa): 10⁸ (n_ws − 1) = c_f (w₀ + w₁σ² + w₂σ⁴ + w₃σ⁶), c_f = 1.022, w₀ = 295.235, w₁ = 2.6422,
  w₂ = −0.032380, w₃ = 0.004028 (µm⁻² powers).
* Density (BIPM 1981/91; Giacomo 1982, Davis 1992): ρ = (p M_a / Z R T)[1 − x_w (1 − M_w/M_a)];
  M_a = 10⁻³[28.9635 + 12.011×10⁻⁶ (x_c − 400)] kg/mol, M_w = 0.018015 kg/mol, R = 8.314510 J mol⁻¹ K⁻¹;
  svp = exp(A T² + B T + C + D/T) Pa with A = 1.2378847×10⁻⁵ K⁻², B = −1.9121316×10⁻² K⁻¹, C = 33.93711047,
  D = −6.3431645×10³ K (over ice: log₁₀ svp = −2663.5/T + 12.537); enhancement f = α + β p + γ t², α = 1.00062,
  β = 3.14×10⁻⁸ Pa⁻¹, γ = 5.6×10⁻⁷ °C⁻²; x_w = f h svp/p;
  Z = 1 − (p/T)[a₀ + a₁t + a₂t² + (b₀ + b₁t)x_w + (c₀ + c₁t)x_w²] + (p/T)²(d + e x_w²), a₀ = 1.58123×10⁻⁶ K Pa⁻¹,
  a₁ = −2.9331×10⁻⁸ Pa⁻¹, a₂ = 1.1043×10⁻¹⁰ K⁻¹ Pa⁻¹, b₀ = 5.707×10⁻⁶ K Pa⁻¹, b₁ = −2.051×10⁻⁸ Pa⁻¹, c₀ = 1.9898×10⁻⁴ K Pa⁻¹,
  c₁ = −2.376×10⁻⁶ Pa⁻¹, d = 1.83×10⁻¹¹ K² Pa⁻², e = −0.765×10⁻⁸ K² Pa⁻².
* **n_prop − 1 = (ρ_a/ρ_axs)(n_axs − 1) + (ρ_w/ρ_ws)(n_ws − 1)** (eq. 5); ρ_axs at 15 °C, 101325 Pa, x_w = 0; ρ_ws at 20 °C,
  1333 Pa, x_w = 1; ρ_a = p M_a (1 − x_w)/(Z R T), ρ_w = p M_w x_w/(Z R T).
* Stated validity: −40 to +100 °C, 80–120 kPa, 0–100 % RH, ≈ 300–1690 nm; overall uncertainty 2–5×10⁻⁸ (1σ).
  NIST ASD warns Ciddor is valid 350–1200 nm and errs by −1.0×10⁻⁸ at 228 nm and −1.6×10⁻⁶ at 185 nm.
* NIST EMToolBox variant (Stone & Zimmerman): R = 8.314472, IAPWS saturation vapour pressure (K₁…K₁₀ constants,
  eqs. A1–A13), Z_a = 0.9995922115, ρ_vs = 0.00985938.
* Validation (derived): Ciddor Table 1 (5 dry cases) ≤ 0.042×10⁻⁸; Table 2 (8 of 9 moist cases) ≤ 0.06×10⁻⁸;
  Table 3 (5 of 6) ≤ 0.06×10⁻⁸; NIST EMToolBox Table 1 (14 cases, −40…50 °C, 60–120 kPa, 300–1700 nm) ≤ 0.061×10⁻⁸;
  IAPWS svp vs EMToolBox Table 2 ≤ 0.43 Pa. **Two probable misprints in Ciddor 1996**: Table 3 (40 °C, 80 kPa) prints
  19996.5 but the printed B&D value 19883.6 and difference −12.9 imply 19896.5 (we compute 19896.49); Table 2 row 2
  prints 27394.0, the row-to-row change of both B&D columns implies ≈ 27394.3 (we compute 27394.33). Ciddor's Table 3
  "Refs. 3 and 4" column is reproduced (≤ 0.11×10⁻⁸) with Edlén-1966 dispersion × 1.0000799 and D = 96094.69 — i.e.
  without the ITS-90 dispersion factor; with the published 1994 constants it is 0.7–1.1×10⁻⁸ lower.

### 3.6 IUGG/IAG (1999) and Greisen et al. (2006) — FITS WCS Paper III
IUGG 1999, Comptes Rendus of the XXII General Assembly, Birmingham, pp. 110–111, IAG Resolution 3
(https://iugg.org/wp-content/uploads/2022/10/IUGG_GA_1999_Comptes_Rendus.pdf): recommends Ciddor (1996) and
Ciddor & Hill (1999) for < 1 ppm; closed formulae for "standard air with 0.0375 % CO₂ content at T = 273.15 K (0°C),
p = 1013.25 hPa, e = 0.0 hPa": group N_g = (n_g − 1)10⁶ = 287.6155 + 4.88660/λ² + 0.06800/λ⁴; **phase
N_ph = (n_ph − 1)10⁶ = 287.6155 + 1.62887/λ² + 0.01360/λ⁴** (λ "carrier wavelength", µm); ambient group:
N_L = (273.15/1013.25)(N_g p/T) − 11.27 e/T (p, e in hPa).
Greisen, E. W., Calabretta, M. R., Valdes, F. G. & Allen, S. L. 2006, "Representations of spectral coordinates in
FITS", A&A 446, 747–771, doi:10.1051/0004-6361:20053818, 2006A&A...446..747G (arXiv:astro-ph/0507293):
* eq. (64) λ = n(λ_a) λ_a; **eq. (65) n(λ_a) = 1 + 10⁻⁶ (287.6155 + 1.62887/λ_a² + 0.01360/λ_a⁴)**, λ_a = air wavelength in µm,
  "This formula suffices for conversions from air to vacuum when no more than 0.25 parts per million accuracy is
  required"; eq. (66) dλ/dλ_a = 1 + 10⁻⁶ (287.6155 − 1.62887/λ_a² − 0.04080/λ_a⁴); "the vacuum wavelength may be used
  to evaluate Eq. (65) … error in the index of refraction amounts to 1:10⁹". Text calls n the index "of dry air at
  standard temperature and pressure" — but the IUGG source formula is for **0 °C**. Result: Greisen λ_air is
  −77 mÅ at 5000 Å relative to 15 °C standard air (specutils docs reach the same conclusion: "Greisen (2006) assumes
  an air temperature of 0C"). Ratio IUGG/Edlén ≈ 1.000015 in n (WCSLIB note).
* Greisen footnote: "The quoted formulæ apply only at normal optical wavelengths. In the UV and IR … atmospheric
  absorption lines cause refractivity to be a strong, weather dependent function of frequency. See … Mathar (2004)"
  (Mathar, R. J. 2004, Appl. Opt. 43, 928).

### 3.7 Morton (2000) / VALD
Morton, D. C. 2000, ApJS 130, 403, doi:10.1086/317349, 2000ApJS..130..403M — eq. 8 (full text not accessible from
here; quoted via the VALD wiki and specutils): **n = 1 + 0.0000834254 + 0.02406147/(130 − s²) + 0.00015998/(38.9 − s²),
s = 10⁴/λ_vac (Å)** = Birch & Downs (1994). VALD wiki ("Air-to-vacuum conversion", N. Piskunov, last edited
2018-12-12, https://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion): vac→air uses this formula ("applies to dry
air at 1 atm pressure and 15ºC with 0.045% CO2 by volume"); **air→vac uses Piskunov's fitted inverse
n = 1 + 0.00008336624212083 + 0.02408926869968/(130.1065924522 − s²) + 0.0001599740894897/(38.92568793293 − s²),
s = 10⁴/λ_air**; air wavelengths only offered > 200 nm. Derived: Piskunov inverse reproduces the exact inverse of
Morton 2000 to ≤ 1×10⁻⁵ mÅ over 2000–25000 Å.

### 3.8 Other implementations (all checked)
* **SDSS**: DR7 page (2003) / DR8–DR12: `AIR = VAC / (1.0 + 2.735182E-4 + 131.4182 / VAC^2 + 2.76249E8 / VAC^4)`
  attributed to Morton (1991) — but Morton 1991 prints only the Edlén 3-term formula; Allende Prieto (2011, APOGEE
  technical note) traces this to IDL VACTOAIR/AIRTOVAC "approximating Eq. 1 by a truncated expansion", "off by 2 parts
  in a million in the infrared" (fixed by W. Landsman after a bug report by D. Schlegel). **Since DR9 (2012)**:
  `AIR = VAC / ( 1.0 + 5.792105E-2/(238.0185 - (1E4/VAC)^2) + 1.67917E-3/( 57.362 - (1E4/VAC)^2))` (Ciddor). SDSS
  spectra are stored in vacuum, barycentric (DR17/DR20 pages). DR7 table pairs (Hα 6562.801/6564.614 etc.) are
  consistent with Ciddor, not with the DR7 polynomial (−1.1 to −3.0 mÅ).
* **IDL Astronomy Library** (Landsman): `vactoair.pro` Ciddor standard air for λ ≥ 2000 Å; `airtovac.pro` Ciddor with
  two fixed-point iterations ("Iterate for better precision W.L./D. Schlegel Mar 2011"), λ < 2000 Å unchanged
  ("Take care within 1 A of 2000 A"). Header examples reproduced: vac 2000 → air 1999.353 (we: 1999.35262);
  air 6056.125 (Kr) → vac 6057.8019 (we: 6057.80193).
* **IRAF** `noao.onedspec.disptrans` (dispcor/t_disptrans.x `dt_airvac`) and `noao/lib/linelists/airvac.cl`:
  "formulas in Allen's Astrophysical Quantities (1973)": x2 = 1/λ[µm]², a = 64.328 + 29498.1/(146 − x2) + 255.4/(41 − x2);
  a = a·p·(1 + (1.049 − 0.0157 t)10⁻⁶ p)/(720.883 (1 + 0.003661 t)); a = a − (0.0624 − 0.000680 x2)/(1 + 0.003661 t)·f;
  n = 1 + a/10⁶; defaults t = 15 °C, p = 760 mmHg, **f = 4 mmHg**; n evaluated at the *input* wavelength (so air2vac is
  not the exact inverse of vac2air; error ≤ 0.1 mÅ). `noao/astutil/asttools/astlvac.x`: Edlén 1966 with σ from λ_air.
* **astropy**: no air/vac function in astropy core units/SpectralCoord. **astropy.wcs (WCSLIB, `spx.c`) AWAV ↔ WAVE
  uses Edlén 1953 (via Cox 2000, Allen's AQ 4th ed.) evaluated at the AIR wavenumber**:
  n = 1.000064328 + 294.981e8/(1.46e14 − s) + 2.554e8/(0.41e14 − s), s = (1/λ_air[m])² — *not* Greisen eq. 65 (the
  WCSLIB header documents this deliberate choice: "This is the IAU standard, adopted in 1957 and again in 1991 …
  it is the one used by WCSLIB"). Live astropy 7.2.0 `wcs.sptr('WAVE-A2W')` reproduces our re-implementation to
  4×10⁻⁸ Å; it differs from Edlén 1953 with vacuum argument by +0.10 mÅ at 2000 Å, +0.018 mÅ at 5000 Å.
* **specutils 2.4.0** `specutils.utils.wcs_utils` (executed from the downloaded wheel): `refraction_index(method=…)`
  with 'Morton2000' (default since 1.17), 'Greisen2006', 'Edlen1953', 'Edlen1966', 'PeckReeder1972' (= eq. 2),
  'Ciddor1996' (+co2); raises below 200 nm except Greisen. `air_to_vac(scheme='inversion')` (default) evaluates n at the
  air wavelength: error +0.098 mÅ at 2000 Å, +0.018 mÅ at 5000 Å, +0.004 mÅ at 2 µm; `scheme='iteration'` is exact;
  `scheme='Piskunov'` = VALD inverse. All specutils outputs reproduced to ≤ 5×10⁻⁸ Å.
* **PyAstronomy 0.25** `pyasl.airtovac2/vactoair2`: modes 'edlen53' (> 2000 Å), 'peckReeder' (= P&R eq. 2, 2300–16900 Å),
  'ciddor' (default, 3000–16900 Å, optional xc); air→vac iterated to 1e-12.

---------------------------------------------------------------------------------------------------

## 4. NIST ASD (what it does — documentation + reverse engineering)

### 4.1 Documentation (https://physics.nist.gov/PhysRefData/ASD/Html/lineshelp.html, accessed 2026-11-23)
* "In ASD, the index of refraction of air is derived from the five-parameter formula given by E.R. Peck and
  K. Reeder, J. Opt. Soc. Am. 62, 958 (1972). These authors fitted data between 185 nm and 1700 nm."
* Ambiguity example: "a wavelength of 200.0648 nm in vacuum corresponds to 200.0000 nm in "standard air" (i.e., 15 °C,
  101 325 Pa pressure, with 0.033 % CO₂). Conversely, an air wavelength of 199.9352 nm corresponds to 200.0000 nm in
  vacuum." Default convention by wavenumber: σ ≥ 50 000 cm⁻¹ → vacuum; 5000 < σ < 50 000 cm⁻¹ → air;
  σ ≤ 5000 cm⁻¹ → vacuum.
* Uncertainty: "The relative uncertainty δλ/λ of the formula of Peck & Reeder is 5×10⁻⁹ for wavelengths λ > 400 nm and
  increases according to the approximate formula δλ/λ ≈ (0.35734 + 38.24/(λ − 180.29) + 0.000023λ)×10⁻⁸ for shorter
  wavelengths (λ is in nm). The maximum uncertainty is about 9×10⁻⁸ at 185 nm." Air-wavelength uncertainties are the
  quadrature sum with this term (so air values often carry one digit fewer than vacuum).
* Validity: "experimentally verified in the wavelength range between 185 nm and 1690 nm … wavelengths shorter than
  185 nm … cannot be displayed in air. For wavelengths >1690 nm, output of wavelengths in air is allowed, because the
  formula … behaves smoothly in this range. However, the users should be cautious …".
* ASD on alternatives: Birch & Downs (1994) "have a five times lower uncertainty … but only in a restricted wavelength
  range from 350 nm to 650 nm"; Ciddor "valid only in the wavelength range from 350 nm to 1200 nm and give
  increasingly large errors for shorter wavelengths: about −1.0×10⁻⁸ at 228 nm and −1.6×10⁻⁶ at 185 nm."
* Form options (`show_av`): 0 = vac < 2000 Å, air 2000–10 000 Å, wavenumber > 10 000 Å; 1 = vac < 10 000 Å, wavenumber
  above; **2 = vac < 2000, air 2000–20 000, vac > 20 000 (default)**; 3 = vacuum all; **4 = vac < 1850 Å, air > 1850 Å**;
  5 (commented out) = wavenumber all.

### 4.2 Verification (derived; `fetch_asd_airvac_check.py`, `validate_airvac.py`)
Fetched Fe I, Ne I, Ar I, Kr I, Hg I, Na I, H I, Th I (1800–26000 Å) with show_av = 3, 2, 4 (2026-11-23).
Matched 35 988 line pairs (≥ 4 decimals in both media; obs and Ritz) between the vacuum (av=3) and air (av=4) output:

| formula | 1850–2000 Å | 2000–3000 | 3000–5000 | 5000–10000 | 10000–20000 | 20000–26000 |
|---|---|---|---|---|---|---|
| **P&R eq. 3 (5-param)**: mean / rms resid (mÅ), % within rounding | −0.002 / 0.017, 100 % | −0.002 / 0.022, 99.97 % | 0.000 / 0.023, 99.99 % | +0.000 / 0.026, 99.98 % | +0.001 / 0.029, 100 % | +0.004 / 0.038, 100 % |
| P&R eq. 2 | −1.63 mÅ, 0 % | −0.022, 80 % | 95 % | 92 % | 95 % | 91 % |
| Ciddor 1996 | −1.60, 0 % | +0.017, 44 % | +0.054, 46 % | +0.097, 29 % | +0.183, 2 % | +0.288, 0 % |
| B&D 1994 / Morton 2000 | +0.098 | +0.037 | +0.053 | +0.105 | +0.248 | +0.414 |
| Edlén 1966 | +0.067 | −0.001 | −0.003 | +0.014 | +0.060 | +0.111 |
| Edlén 1953 | +0.021 | −0.012 | +0.023 | +0.001 | −0.088 | −0.174 |

(residual = ASD air − formula(ASD vacuum)). Na D example: vac 5891.583264 → P&R eq. 3 gives 5889.950954, ASD air
5889.95095. Medium boundaries actually returned: show_av=4 air from 1851.25 Å (Fe I 1850.622 Å air) upward incl.
> 20 000 Å; show_av=2 air 2001.2–19986.9 Å (vacuum above, e.g. Fe I 20011.18 Å vac).

---------------------------------------------------------------------------------------------------

## 5. Comparison across 2000–25000 Å (derived; `airvac_comparison.csv`)
λ_air(formula) − λ_air(Ciddor 1996, 450 ppm) at fixed λ_vac, in mÅ:

| λ_vac (Å) | vac−air (Å) | Edlén53 | Edlén66 | P&R eq2 | P&R eq3 (ASD) | B&D93 | B&D94=Morton00=VALD | Ciddor 420 ppm | WCSLIB/astropy.wcs | Greisen06 (0 °C) | SDSS DR7 poly | IRAF default (f=4 mmHg) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2000 | 0.64738 | −0.566 | −0.564 | +0.032 | −0.538 | −0.636 | −0.596 | +0.010 | −0.664 | −26.15 | +0.31 | −0.22 |
| 3000 | 0.87445 | +0.034 | +0.043 | +0.043 | +0.038 | −0.053 | +0.000 | +0.014 | −0.003 | −47.48 | +0.11 | +0.66 |
| 5000 | 1.39448 | +0.050 | +0.070 | +0.069 | +0.073 | −0.083 | +0.002 | +0.022 | +0.032 | −76.85 | −1.21 | +1.18 |
| 6564.6 | 1.81297 | +0.093 | +0.082 | +0.090 | +0.096 | −0.117 | −0.007 | +0.029 | +0.079 | −99.86 | −3.06 | +1.61 |
| 10000 | 2.74091 | +0.180 | +0.100 | +0.136 | +0.137 | −0.202 | −0.034 | +0.044 | +0.171 | −150.84 | −6.93 | +2.52 |
| 20000 | 5.45859 | +0.405 | +0.160 | +0.270 | +0.253 | −0.441 | −0.108 | +0.087 | +0.400 | −300.23 | −16.88 | +5.12 |
| 25000 | 6.81978 | +0.513 | +0.192 | +0.338 | +0.312 | −0.558 | −0.142 | +0.109 | +0.509 | −375.06 | −21.58 | +6.42 |

Inversion/round-trip errors of software (mÅ, in λ_vac): specutils 'inversion' +0.098 (2000 Å) → +0.018 (5000 Å) →
+0.003 (25000 Å); IRAF air2vac(vac2air) identical magnitude; Greisen's recommended shortcut (n at λ_vac) +0.078 (2000 Å);
VALD Piskunov ≤ 1e-5; IDL airtovac 2-iteration ≤ 1e-5 for λ_air ≥ 2000 Å (but −0.647 Å discontinuity for
1999.35 ≤ λ_air < 2000: not converted); WCSLIB 4 fixed iterations converge exactly.

Below 2000 Å (not in the table): the formulas diverge (P&R eq. 2 and Ciddor are −1.6…−2.9 mÅ from ASD at 1850–2000 Å;
poles at λ = 1/√k: Ciddor/P&R eq. 2 1320 Å, Edlén 1966 1603 Å, Edlén 1953 1562 Å, P&R eq. 3 1595 Å). Air wavelengths
are conventionally not used below 2000 Å (ASD option down to 1850 Å).

---------------------------------------------------------------------------------------------------

## 6. Lab / observatory air vs. standard air (Ciddor full model, derived; `airvac_sensitivity_ciddor_5000A.csv`)

At λ_vac = 5000 Å:

| condition | 10⁸(n−1) | λ_air − λ_air(std) | dλ/dT | dλ/dp | dλ/dRH | dλ/dCO₂ |
|---|---|---|---|---|---|---|
| standard 15 °C, 101325 Pa, 0 %, 450 ppm | 27897.38 | 0 | **+4.85 mÅ/K** (291 m s⁻¹/K) | **−1.38 mÅ/hPa** (−82.5 m s⁻¹/hPa) | +0.031 mÅ/%RH (1.9 m s⁻¹) | −0.074 mÅ/100 ppm (−4.5 m s⁻¹) |
| lab 20 °C, 101325 Pa, 50 %, 420 ppm | 27377.71 | +25.97 mÅ | +4.81 | −1.35 | +0.042 | −0.072 |
| site 2600 m (ISA p = 73 749 Pa), 0 °C, 30 %, 420 ppm | 21412.96 | **+324.1 mÅ** | +3.96 | −1.45 | +0.012 | −0.057 |
| site 4200 m (ISA p = 60 051 Pa), 0 °C, 20 %, 420 ppm | 17435.29 | **+522.9 mÅ** | +3.22 | −1.45 | +0.012 | −0.046 |

(ISA pressures from p = 101325 (1 − 2.25577×10⁻⁵ h)^5.25588 Pa — illustrative, not site measurements.)
Edlén (1966) T,p factor and Ciddor agree for ordinary lab air; Birch & Downs vs Ciddor differ by up to 5×10⁻⁷
(Ciddor: "at high temperatures and humidities"); NIST EMToolBox: a 1 °C temperature error, a 0.4 kPa pressure error, or a
50 % humidity error each yields about 1×10⁻⁶ in n.

**When it matters**
* *Catalogue conversion* (vacuum ↔ standard air) is a convention: ambient conditions are irrelevant; only the formula
  and its standard-air definition matter (§5).
* *Spectrograph in ambient air* (not vacuum-enclosed): the dispersion happens in site air, which differs from standard
  air by ~0.3–0.5 Å at 5000 Å on a mountain — but if the wavelength solution is fitted to lamp lines through the same
  air, the difference is almost a pure scale factor and is absorbed: after a cubic fit the residual is ≤ 0.02 mÅ
  (3800–6800 Å) / ≤ 0.17 mÅ (3600–10000 Å) for the 2600 m case (derived). What is **not** absorbed is the change of air
  between calibration and science exposure: 4.9 mÅ K⁻¹ and 1.4 mÅ hPa⁻¹ at 5000 Å (≈ 290 m s⁻¹ K⁻¹, 83 m s⁻¹ hPa⁻¹) —
  the reason precision spectrographs are operated in vacuum (e.g. APOGEE, per Allende Prieto 2011; the HARPS/ESPRESSO
  class of instruments likewise uses vacuum vessels).
* *Instruments that measure wavelength absolutely* (FTS with a He-Ne reference, wavemeters, laser frequency combs used
  in air): the actual n(T, p, RH, CO₂) along the path must be computed (Ciddor/EMToolBox), not standard air.
* *Lamp reference lists* (Th–Ar, U–Ne, Fe) are mostly vacuum wavenumbers from vacuum FTS (e.g. Kerber et al. 2008;
  Lovis & Pepe 2007, as listed by Allende Prieto 2011): convert with the same formula the pipeline uses for science
  lines.
* Beyond ~1.7 µm all dispersion formulas are extrapolations; water-vapour and CO₂ resonances make n humidity-
  dependent in the IR (Greisen 2006 footnote; Mathar 2004).

---------------------------------------------------------------------------------------------------

## 7. Recommendations for the wiki page (calibration section)
1. Primary column = vacuum wavelength (Å, string as given by the source + its stated uncertainty).
2. Air column = λ_vac / n_Ciddor1996(σ_vac) for λ ≥ 2000 Å (declare: 15 °C, 101 325 Pa, dry, 450 ppm). When showing
   NIST-ASD air values, keep ASD's own string (P&R eq. 3); they differ from Ciddor by +0.05…+0.29 mÅ (3000–26000 Å).
3. For inverse (air→vac), iterate (3–4 fixed-point steps reach 1e-12 Å) — do not evaluate n at the air wavelength if
   1e-4 Å is wanted.
4. Flag, in any line table, the medium boundary (2000 Å; 1850 Å for ASD show_av=4; > 20000 Å vacuum in ASD default).
5. Interactive calculator: expose T, p, RH, CO₂ (Ciddor) and show the shift relative to standard air and the
   velocity equivalent.

## 8. Caveats / not verified
* Morton (2000) and Birch & Downs (1994) full texts, IAU Trans. 9 (Oosterhoff 1957) and the IAU 1991 resolution C15 text
  were not accessible; their content is taken from VALD wiki, specutils docs, CRC 87th ed., the B&D 1994 abstract and
  the APOGEE note (all saved/cited).
* P&R eq. (3) coefficients come from refractiveindex.info (secondary), but they reproduce P&R Table V and NIST ASD.
* The CRC table excerpt was retrieved via a web-extraction tool (direct download failed); saved as text.
* Ciddor PDF obtained from a personal-web mirror (walter.bislins.ch); text checked against rendered page images.
* "lab/site" columns use illustrative conditions and ISA pressures, not measured site data.
