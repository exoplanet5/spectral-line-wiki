# Why 1e-3 Å? A numbers-first rationale for the wiki's wavelength precision

Research notes for the Spectral Line Wiki. Written 2026-11-23.

**What produced these numbers**

| file | role |
|---|---|
| `scripts/research/precision_rationale.py` | computes every **[derived]** number below. Run it as `~/.venvs/astro313/bin/python -I scripts/research/precision_rationale.py <project root>` |
| `data/notes/precision_rationale_numbers.json` | machine-readable output of the script |
| `data/curated/spectrograph_precision.csv` | spectrograph table: 31 instrument modes, each with a source string and a URL |
| `scripts/research/airvac_formulas.py` | the project's air/vacuum formulas, reused here. NIST ASD is reproduced to 0.02 mÅ rms with Peck & Reeder eq. 3 (`data/curated/airvac_validation_asd.csv`) |
| `scripts/research/prec_fetch_asd.py`, `prec_crossref.py` | NIST ASD windows and Crossref metadata. Raw downloads are under `data/raw/precision_rationale/` |

**Provenance flags**

* **[lit]** marks a value copied from the cited source. Wavelength strings are kept exactly as printed.
* **[derived]** marks a value computed by the script.
* **[assumed]** marks an explicit assumption that no source supplies. Three spectrograph samplings are flagged this way.

**Constants.** c = 299 792 458 m/s (exact). The other constants come from `astropy.constants`: CODATA 2018 h, e and mₑ; IAU 2015 B3 nominal GM☉ and R☉; IAU 2012 au.

---

## 0. Verdict

1. **What 1 mÅ means in velocity.** 1 mÅ is 59.96 m/s at 5000 Å. Across 1200 Å–2.2 µm it is 250–14 m/s. In frequency it is 119.9 MHz (4.0 × 10⁻³ cm⁻¹) at 5000 Å.
2. **1 mÅ is the precision real spectrographs can use.** A single unblended line observed at R ≥ 50 000 and SNR ≈ 100 per pixel centroids to 0.15–0.6 mÅ. This follows from the photon-noise bound of Bouchy, Pepe & Queloz (2001). The rule of thumb FWHM/(SNR√N_pix) gives 0.14–0.62 mÅ for the 13 optical echelle modes in the table. At 3 px/FWHM and 5000 Å, the rule reaches 1 mÅ at SNR 100 once R ≳ 29 000. A catalogue wavelength rounded to 0.01 Å has an rms rounding error of 2.9 mÅ, about 10 times larger than these centroid errors. A value rounded to 0.001 Å has an rms rounding error of 0.29 mÅ, comparable to them. **Finer than 1 mÅ is better. Coarser is not good enough for R ≳ 30 000 work.**
3. **Real physical shifts are 1–10 mÅ at 5000 Å.** A 10 mÅ grid blurs all of them; a 1 mÅ grid resolves them:
   * Earth rotation: 7.8 mÅ.
   * Solar gravitational redshift: 10.6 mÅ.
   * Solar convective blueshift: 0 to −9.1 mÅ.
   * Isotope shifts: Ca II K ⁴⁰Ca→⁴⁸Ca −9.0 mÅ; Mg II 2796 −8.1 mÅ; Li 6708 +158 mÅ.
   * Hyperfine splittings: 4–35 mÅ.
   * Zeeman shift: 11.7 mÅ for g = 1 and B = 1 kG.
   * Air model: lab air versus standard air shifts the air wavelength by 26 mÅ.
4. **Laboratory data can deliver 1 mÅ.** Where NIST ASD quotes an uncertainty, the median is:
   * 0.13 mÅ for lines that the Handbook gives to 4 decimals;
   * 0.6 mÅ for lines the Handbook gives to 3 decimals;
   * 3 mÅ for lines the Handbook gives to only 2 decimals.

   HARPS Th–Ar lines are known to about 0.18 mÅ (Lovis & Pepe 2007).
5. **1 mÅ is not enough for:**
   * Extreme-precision radial velocities (EPRV): 10 cm/s is 1.7 × 10⁻⁶ Å.
   * Laser-frequency-comb calibration: 1 cm/s is 1.7 × 10⁻⁷ Å.
   * Δα/α studies: they need rest wavelengths to ≲ 20 m/s, i.e. 1.6 × 10⁻⁴ Å at 2383 Å. A 1 ppm change in α moves lines by 0.12–0.20 mÅ.
   * Redshift drift: about 1 × 10⁻⁶ Å per decade.

   These measurements are differential and rely on a comb or etalon, not on catalogue rest wavelengths.
6. **1 mÅ is overkill for low-resolution survey redshifts.**
   * The R ≈ 2000–5000 surveys (DESI, BOSS, LAMOST-LRS, PFS, WEAVE-LR, 4MOST-LRS) centroid single lines to 6–20 mÅ at SNR 100.
   * DESI's own redshift requirements correspond to 1.4–3.5 Å at 7000 Å.
   * For these surveys, 0.01 Å rest wavelengths are sufficient. The air/vacuum convention (≈ 1.4–6 Å) still matters.
7. **Display rule for the wiki.**
   * Store the source string verbatim, with its medium and its stated uncertainty.
   * Display 3 decimals only when the source has them, and never pad.
   * Flag values with fewer decimals; 43 % of NIST Handbook lines (5 122 of 11 980) have only 2.
   * Keep 4–9 decimals wherever the source gives them. Examples are Ritz wavelengths, Th/U atlases, and the Murphy & Berengut (2014) table.

---

## 1. Velocity, frequency and energy equivalents of 1 mÅ [derived]

The step is converted with δv = c δλ/λ, δν = c δλ/λ², δσ = δλ/λ² and δE = h δν. The frequency and wavenumber columns use the vacuum wavelength.

| λ (Å) | 1 mÅ = m/s | 10 mÅ = m/s | 1 m/s = Å | 10 cm/s = Å | 1 mÅ = MHz | 1 mÅ = cm⁻¹ | 1 mÅ = eV |
|---|---|---|---|---|---|---|---|
| 1200 | 249.83 | 2498.3 | 4.003e-06 | 4.00e-07 | 2081.89 | 6.944e-02 | 8.610e-06 |
| 3000 | 99.93 | 999.3 | 1.001e-05 | 1.00e-06 | 333.10 | 1.111e-02 | 1.378e-06 |
| 4000 | 74.95 | 749.5 | 1.334e-05 | 1.33e-06 | 187.37 | 6.250e-03 | 7.749e-07 |
| 5000 | 59.96 | 599.6 | 1.668e-05 | 1.67e-06 | 119.92 | 4.000e-03 | 4.959e-07 |
| 6563 | 45.68 | 456.8 | 2.189e-05 | 2.19e-06 | 69.60 | 2.322e-03 | 2.878e-07 |
| 8000 | 37.47 | 374.7 | 2.669e-05 | 2.67e-06 | 46.84 | 1.562e-03 | 1.937e-07 |
| 10000 | 29.98 | 299.8 | 3.336e-05 | 3.34e-06 | 29.98 | 1.000e-03 | 1.240e-07 |
| 16000 | 18.74 | 187.4 | 5.337e-05 | 5.34e-06 | 11.71 | 3.906e-04 | 4.843e-08 |
| 22000 | 13.63 | 136.3 | 7.338e-05 | 7.34e-06 | 6.19 | 2.066e-04 | 2.562e-08 |

* **Air or vacuum makes no difference to δv.** λ_air = λ_vac/n(λ_vac), so (δλ/λ)_air = (δλ/λ)_vac (1 − λ n′/n). With Ciddor standard air the correction is +4.2e-5 at 3000 Å, +1.3e-5 at 5000 Å and +6e-7 at 2.2 µm [derived]. The velocity equivalent of 1 mÅ is therefore the same for air and vacuum wavelengths to < 0.005 %.
* **Fixed decimals mean unequal velocity precision.** 1 mÅ is 250 m/s in the far-UV but 14 m/s at 2.2 µm. For a uniform 50 m/s target, the step is 2.0e-4 Å at 1200 Å, 8.3e-4 Å at 5000 Å and 3.7e-3 Å at 22 000 Å. So far-UV lines need 4 decimals, and 3 decimals are generous in the K band. The wiki should show the source's own decimals rather than a fixed format.

---

## 2. What spectrographs can actually measure

### 2.1 Method

The photon-noise limit on a velocity shift is given by Bouchy, Pepe & Queloz (2001, A&A 374, 733) [lit]:

* eq. (8): W(i) = λ²(i) (∂A₀(i)/∂λ(i))² / (A₀(i) + σ_D²)
* eq. (10): δV_RMS/c = 1/√ΣW(i)
* eq. (12): δV_RMS = c/(Q√N_e⁻)

For a single line, eq. (10) reduces to δλ_RMS = [Σ (∂A₀/∂λ)²/A₀]^(-1/2). That is the photon-noise Cramér–Rao bound (CRB) on the centroid. The script evaluates it numerically under these assumptions:

* The line is unresolved: its observed profile is the instrumental Gaussian of FWHM λ/R.
* The profile is integrated over pixels and averaged over sub-pixel phase.
* σ_D = 0.

Two cases are computed:

* **Absorption.** Line depth 0.5 in a continuum of SNR 100 per pixel, i.e. 10⁴ e⁻.
* **Emission.** Peak 10⁴ e⁻ per pixel and no background.

The numerical code agrees with the analytic limits. For emission the limit is σ = s/√N_tot; the numerical result is 0.9209 mÅ against 0.9204 mÅ analytic. For a weak absorption line the limit is σ = √(2sp/√π)/(SNR·D); the results are 0.1545 and 0.1548.

The rule of thumb **σ ≈ FWHM/(SNR√N_pix)**, with N_pix the number of pixels per FWHM, is tabulated for comparison:

* For a 50 %-deep absorption line it is within 20 % of the CRB.
* For emission it is a factor 2.43 more pessimistic than the CRB (continuous limit).

None of these numbers include read noise, blends, errors in the line-spread function (LSF), or calibration.

### 2.2 Instruments

All values are at SNR = 100 per pixel. FWHM = λ_ref/R unless the source quotes it. The column "CRB abs." is a 50 %-deep absorption line.

| instrument | R | λ_ref (Å) | FWHM (Å) | FWHM (km/s) | px/FWHM | pixel (mÅ) | rule (mÅ) | rule (m/s) | CRB abs. (mÅ) | CRB emission (mÅ) |
|---|---|---|---|---|---|---|---|---|---|---|
| HARPS | 115,000 | 5500 | 0.0478 | 2.61 | 3.2 | 14.9 | 0.27 | 15 | 0.32 | 0.11 |
| ESPRESSO HR | 140,000 | 5500 | 0.0393 | 2.14 | 4.5 | 8.73 | 0.19 | 10 | 0.22 | 0.077 |
| ESPRESSO UHR | 190,000 | 5500 | 0.0289 | 1.58 | 2.5 | 11.6 | 0.18 | 10 | 0.23 | 0.078 |
| EXPRES | 137,500 | 5500 | 0.0400 | 2.18 | 4.0 | 10.0 | 0.20 | 11 | 0.24 | 0.084 |
| NEID HR | 110,000 | 5500 | 0.0500 | 2.73 | 5.0 | 10.0 | 0.22 | 12 | 0.26 | 0.093 |
| PEPSI 270k | 270,000 | 5500 | 0.0204 | 1.11 | 2.0 | 10.2 | 0.14 | 7.9 | 0.18 | 0.063 |
| PEPSI 120k | 120,000 | 5500 | 0.0458 | 2.50 | 4.2 | 10.9 | 0.22 | 12 | 0.27 | 0.093 |
| PEPSI 43k | 43,000 | 5500 | 0.128 | 6.97 | 12.1 | 10.6 | 0.37 | 20 | 0.43 | 0.15 |
| UVES red 1″ | 38,700 | 6000 | 0.155 | 7.75 | 6.25ᵈ | 24.8 | 0.62 | 31 | 0.73 | 0.26 |
| UVES red 0.3″ | 110,000 | 6000 | 0.0545 | 2.73 | 2.0 | 27.3 | 0.39 | 19 | 0.49 | 0.17 |
| Keck HIRES 0.861″ | 50,000 | 5240 | 0.105 | 6.00 | 4.8ᵈ | 21.8 | 0.48 | 27 | 0.57 | 0.20 |
| Keck HIRES 0.574″ | 66,400 | 5240 | 0.0789 | 4.51 | 3.2ᵈ | 24.7 | 0.44 | 25 | 0.53 | 0.19 |
| Subaru HDS 0.2″ | 160,000 | 5500 | 0.0344 | 1.87 | 2.0 | 17.2 | 0.24 | 13 | 0.31 | 0.11 |
| CRIRES+ 0.2″ (H1582) | 104,800 | 15820 | 0.151 | 2.86 | 3.6ᵈ | 42.3 | 0.80 | 15 | 0.96 | 0.33 |
| IGRINS (H) | 45,000 | 16000 | 0.356 | 6.66 | 3.66 | 97.1 | 1.9 | 35 | 2.2 | 0.78 |
| APOGEE (green) | 23,700 | 16100 | 0.679 | 12.6 | 2.4ᵈ | 282 | 4.4 | 82 | 5.5 | 1.9 |
| Gaia RVS | 11,500 | 8600 | 0.748 | 26.1 | 3.05ᵈ | 245 | 4.3 | 149 | 5.2 | 1.8 |
| DESI blue | 2,500ᵃ | 4750 | 1.90 | 120 | 3.3ᵈ | 569 | 10 | 656 | 13 | 4.4 |
| DESI red | 4,000ᵃ | 6660 | 1.67 | 74.9 | 3.2ᵈ | 515 | 9.3 | 417 | 11 | 3.9 |
| DESI NIR | 4,750ᵃ | 8635 | 1.82 | 63.1 | 3.2ᵈ | 566 | 10 | 352 | 12 | 4.3 |
| SDSS/BOSS | 2,000 | 6000 | 3.00 | 150 | 3.0 | 1000 | 17 | 865 | 21 | 7.3 |
| LAMOST LRS | 1,800 | 5500 | 3.06 | 167 | 2.4ᵈ | 1266 | 20 | 1072 | 25 | 8.4 |
| LAMOST MRS | 7,500 | 5150 | 0.687 | 40.0 | 3ᵃ | 229 | 4.0 | 231 | 4.8 | 1.7 |
| 4MOST HRS | 20,000 | 5450 | 0.273 | 15.0 | 3.0 | 90.8 | 1.6 | 87 | 1.9 | 0.66 |
| 4MOST LRS | 6,500 | 6600 | 1.02 | 46.1 | 3.0 | 338 | 5.9 | 266 | 7.1 | 2.5 |
| PFS blue | 2,500 | 5000 | 2.1 (quoted) | 126 | 3.0ᵈ | 700 | 12 | 727 | 15 | 5.1 |
| PFS red LR | 3,000 | 8000 | 2.7 (quoted) | 101 | 3.0ᵈ | 900 | 16 | 584 | 19 | 6.6 |
| PFS red MR | 5,500 | 8000 | 1.6 (quoted) | 60.0 | 4.0ᵈ | 400 | 8.0 | 300 | 9.6 | 3.3 |
| PFS NIR | 4,500 | 11000 | 2.4 (quoted) | 65.4 | 3.0ᵈ | 800 | 14 | 378 | 17 | 5.8 |
| WEAVE HR | 20,000 | 5000 | 0.250 | 15.0 | 3ᵃ | 83.3 | 1.4 | 87 | 1.8 | 0.61 |
| WEAVE LR | 5,000 | 6600 | 1.32 | 60.0 | 3ᵃ | 440 | 7.6 | 346 | 9.3 | 3.2 |

Footnotes:

* **ᵈ derived sampling:**
  * UVES: slit width/pixel scale = 1.0″/0.16″.
  * HIRES: the binning page says the 0.86″ slit spans 2.4 binned (×2) pixels, i.e. 4.8 unbinned. The 0.574″ value is scaled from that.
  * CRIRES+: 0.2″/0.056″.
  * APOGEE, Gaia, DESI, LAMOST and PFS: FWHM divided by the quoted dispersion. DESI's dispersion is the channel range divided by the CCD columns. LAMOST uses its 10⁻⁴-dex output grid.
* **ᵃ assumed:**
  * DESI R is the mid-point of the quoted ranges (2000–3000, 3500–4500, 4000–5500).
  * The sampling of LAMOST-MRS and WEAVE is not stated in the cited sources; 3 px is assumed.
* Every row's verbatim source values, notes, citation and URL are in `data/curated/spectrograph_precision.csv`.

### 2.3 Scaling grid (5000 Å, 3 px per FWHM) [derived]

Each cell gives the rule-of-thumb value and the CRB for a 50 %-deep absorption line, in mÅ.

| R | FWHM (Å) | SNR 10 | SNR 30 | SNR 100 | SNR 300 | SNR 1000 |
|---|---|---|---|---|---|---|
| 1,000 | 5.00 | 289 / 351 | 96 / 117 | 29 / 35 | 9.6 / 12 | 2.9 / 3.5 |
| 2,000 | 2.50 | 144 / 176 | 48 / 59 | 14 / 18 | 4.8 / 5.9 | 1.4 / 1.8 |
| 5,000 | 1.00 | 58 / 70 | 19 / 23 | 5.8 / 7.0 | 1.9 / 2.3 | 0.58 / 0.70 |
| 10,000 | 0.500 | 29 / 35 | 9.6 / 12 | 2.9 / 3.5 | 0.96 / 1.2 | 0.29 / 0.35 |
| 20,000 | 0.250 | 14 / 18 | 4.8 / 5.9 | 1.4 / 1.8 | 0.48 / 0.59 | 0.14 / 0.18 |
| 50,000 | 0.100 | 5.8 / 7.0 | 1.9 / 2.3 | 0.58 / 0.70 | 0.19 / 0.23 | 0.058 / 0.070 |
| 100,000 | 0.0500 | 2.9 / 3.5 | 0.96 / 1.2 | 0.29 / 0.35 | 0.096 / 0.12 | 0.029 / 0.035 |
| 140,000 | 0.0357 | 2.1 / 2.5 | 0.69 / 0.84 | 0.21 / 0.25 | 0.069 / 0.084 | 0.021 / 0.025 |
| 200,000 | 0.0250 | 1.4 / 1.8 | 0.48 / 0.59 | 0.14 / 0.18 | 0.048 / 0.059 | 0.014 / 0.018 |

* **R ≥ 50 000, SNR ≈ 100:** one line centroids to ≤ 0.6 mÅ (rule) or ≤ 0.7 mÅ (CRB).
* **R ≈ 2000 surveys:** 14–18 mÅ at SNR 100 and 48–59 mÅ at SNR 30. This is the "10–50 mÅ suffices" regime.
* **Threshold:** the 1 mÅ rule-of-thumb level at SNR 100 and 3 px needs FWHM ≤ 0.173 Å, i.e. R ≥ 28 900 at 5000 Å [derived].

### 2.4 Check against measured HARPS Th–Ar performance [lit]

Lovis & Pepe (2007, A&A 468, 1115) measured HARPS directly:

* R = 110 000, 0.015 Å per pixel (820 m/s per pixel), about 3.3 pixels per resolution element. The derived HARPS pixel above, 14.9 mÅ, agrees.
* Strongest unsaturated lines: photon noise of 0.002 px = 1.6 m/s = 0.029 mÅ.
* Faintest detectable lines: 0.2 px = 160 m/s ≈ 2.9 mÅ.
* CCD inhomogeneities set a floor of about 0.01 px = 8 m/s = 0.15 mÅ.
* Line positions in the new list are good to about 10 m/s = 0.18 mÅ.

The emission CRB at peak SNR 100 is 0.11 mÅ ≈ 0.0075 px, and it scales as 1/SNR. It therefore falls inside the measured range (between the brightest and faintest lines). The CCD floor (0.15 mÅ) shows that **sub-mÅ per-line accuracy is the practical ceiling for a single line, even on HARPS.** RV precision below that comes only from averaging thousands of lines.

---

## 3. Physical and instrumental effects at the mÅ level

### 3.1 Velocity ladder, expressed at 5000 Å unless stated

| effect | v | Δλ (mÅ) | flag / source |
|---|---|---|---|
| Earth rotation at the equator, ω_E a_E | 465.1 m/s | 7.76 | [derived] IERS TN36 ω = 7.292115e-5 rad/s, a = 6 378 137 m |
| Earth orbit, √(GM☉/au) | 29.785 km/s | 497 | [derived]. DE440s range for 2026 is 29.294–30.278 km/s, i.e. 489–505 mÅ (`velocity_frames.md`) |
| Solar gravitational redshift, surface to ∞, GM☉/(R☉c) | 636.3 m/s | 10.6 | [derived] IAU 2015 B3 nominal values |
| Solar gravitational redshift seen at 1 au, GM☉/c (1/R☉ − 1/au) | 633.4 m/s | 10.6 | [derived]. González Hernández et al. (2020) quote 633.35 m/s |
| Solar gravitational redshift measured with HARPS-LFC Moon spectra | 639 ± 14 m/s | 10.7 | [lit] González Hernández et al. 2020, A&A 643, A146 |
| Convective blueshift of solar Fe I lines, relative depth 0.1 → 0.9 | −547 → +2 m/s | −9.12 → +0.03 | [lit] Reiners et al. 2016, Table 5. Per-bin σ is 128–303 m/s |
| Fe I 5250.2084 Å bisector, disk centre | −100 → −290 m/s | −1.8 → −5.1 (at 5250 Å) | [lit] Löhner-Böttcher et al. 2019 |
| Sun about the solar-system barycentre (max, 1990–2040) | 16.1 m/s | 0.27 | [derived] DE440s (`velocity_frames.md`) |
| Earth about the Earth–Moon barycentre | 12.4 m/s | 0.21 | [derived] DE440s |
| Kitt Peak FTS solar atlases: wavelength-scale errors | up to ~100 m/s | 1.7 | [lit] Reiners et al. 2016 |
| IAG solar flux atlas accuracy | ±10 m/s | 0.17 | [lit] Reiners et al. 2016 |
| Palmer & Engleman (1983) Th atlas | 16–82 m/s | 0.29–1.50 (at 5500 Å) | [lit] as quoted by Lovis & Pepe 2007 |
| Ar lines in Th–Ar lamps, lamp-to-lamp shifts | "a few tens of m/s"; 10–50 m/s shown | 0.18–0.92 (at 5500 Å) | [lit] Lovis & Pepe 2007, citing Whaling et al. 2002 |
| Earth-induced reflex of the Sun | 9 cm/s | 0.0015 | [lit] Fischer et al. 2016 |

**Consequences.**

* **Earth rotation (7.8 mÅ).** At 1 mÅ a topocentric-to-barycentric correction that drops the diurnal term is visibly wrong.
* **Gravitational redshift (10.6 mÅ) and convective blueshift (up to 9 mÅ).** They are the same size and partly cancel. A 2-decimal catalogue (rms rounding error 2.9 mÅ, steps of 10 mÅ) cannot separate them.
* **Barycentric versus heliocentric frame and the lunar reflex (0.2–0.3 mÅ).** These sit just below the last digit.

### 3.2 Air ↔ vacuum: the choice of formula and the state of the air

The calculation fixes λ_vac, computes λ_air with each formula for that formula's own standard air, and subtracts Ciddor (1996) [derived]. All differences are in mÅ.

| λ_vac (Å) | λ_air Ciddor (Å) | vac − air (Å) | Edlén 1953 | Edlén 1966 | Peck & Reeder 1972 eq. 2 | Peck & Reeder 1972 eq. 3 (= NIST ASD) | Birch & Downs 1994 = Morton 2000 (VALD) | spread (mÅ) | spread (m/s) |
|---|---|---|---|---|---|---|---|---|---|
| 2000 | 1999.3526 | 0.6474 | −0.566 | −0.564 | +0.032 | −0.538 | −0.596 | 0.628 | 94.2 |
| 3000 | 2999.1255 | 0.8745 | +0.034 | +0.043 | +0.043 | +0.038 | +0.000 | 0.043 | 4.3 |
| 4000 | 3998.8693 | 1.1307 | +0.026 | +0.057 | +0.056 | +0.053 | +0.002 | 0.057 | 4.3 |
| 5000 | 4998.6055 | 1.3945 | +0.050 | +0.070 | +0.069 | +0.073 | +0.002 | 0.073 | 4.4 |
| 6563 | 6561.1875 | 1.8125 | +0.093 | +0.082 | +0.090 | +0.096 | −0.007 | 0.102 | 4.7 |
| 8000 | 7997.8002 | 2.1998 | +0.130 | +0.090 | +0.109 | +0.114 | −0.018 | 0.148 | 5.5 |
| 10000 | 9997.2591 | 2.7409 | +0.179 | +0.100 | +0.136 | +0.137 | −0.034 | 0.213 | 6.4 |
| 16000 | 15995.6297 | 4.3703 | +0.317 | +0.134 | +0.216 | +0.206 | −0.080 | 0.396 | 7.4 |
| 22000 | 21993.9970 | 6.0030 | +0.448 | +0.173 | +0.297 | +0.276 | −0.122 | 0.569 | 7.8 |

**Validity ranges.** Ciddor and Peck & Reeder eq. 2 are fitted over 0.23–1.69 µm, so the 2000 Å and 22 000 Å rows are extrapolations for those formulas. Peck & Reeder eq. 3 extends to 0.185 µm.

**Consequences.**

* **The air–vacuum difference itself** is 0.65–6.0 Å, i.e. 90–97 km/s. A wiki without an explicit medium flag is wrong by 1000 × the target precision.
* **The choice among the standard formulas** matters at 0.04–0.2 mÅ (4–6 m/s) between 3000 Å and 1 µm. It is about 0.6 mÅ at 2000 Å and at 2.2 µm. This is below 1 mÅ but not negligible, so **record which formula converted any air value.** NIST ASD reproduces Peck & Reeder eq. 3 to 0.02 mÅ rms (`data/curated/airvac_validation_asd.csv`).
* **The state of the air** dominates over the formula choice. With the full Ciddor model at 5000 Å [derived]:

  | change | Δλ_air |
  |---|---|
  | temperature | +4.85 mÅ per K |
  | pressure | −1.38 mÅ per hPa |
  | humidity (dry standard air) | +0.31 mÅ per 10 % RH |
  | CO₂ | −0.074 mÅ per 100 ppm |
  | lab air (20 °C, 101 325 Pa, 50 % RH, 420 ppm) vs standard air (15 °C, dry, 450 ppm) | +26.0 mÅ |

  An "air wavelength" is defined to 1 mÅ only for standard air.

### 3.3 Isotope shifts, hyperfine structure and fine structure

Split = wl2 − wl1. Wavelength strings are verbatim.

| structure | wl1 (Å) | wl2 (Å) | medium | split (mÅ) | split (m/s) | source |
|---|---|---|---|---|---|---|
| Mg II 2796: ²⁴Mg → ²⁶Mg | 2796.35509903(42) | 2796.34704565(42) | vac | −8.053 | −863 | [lit] Murphy & Berengut 2014, Table 2 |
| Mg II 2796: ²⁵Mg hfs, F=2 vs F=3 | 2796.353449(50) | 2796.349030(50) | vac | −4.419 | −474 | [lit] idem |
| Ca II 3934 (K): ⁴⁰Ca → ⁴⁸Ca | 3934.7747160(26) | 3934.76575(65) | vac | −8.966 | −683 | [lit] idem, Table 5 |
| Ca II 3934: ⁴³Ca hfs, F=3 vs F=4 | 3934.77998(11) | 3934.76438(11) | vac | −15.600 | −1189 | [lit] idem |
| Na I 5891 (D2): ²³Na ground-state hfs | 5891.570875 | 5891.590672 | vac | +19.797 | +1007 | [lit] idem, Table 1 |
| Mn II 2576: ⁵⁵Mn hfs, extreme components (14–16 in total) | 2576.890898 | 2576.856181 | vac | −34.717 | −4039 | [lit] idem, Table 9 |
| Li I 6708 D2: ⁷Li → ⁶Li, centres of gravity | 6709.615611 | 6709.773805 | vac, from ν | +158.194 | +7068 | [derived] from Brown et al. 2013, Table I |
| Li I 6708 D2: ⁷Li ground hfs (F=2→3 vs F=1→2) | | | vac, from ν | +12.208 | +545 | [derived] from Brown et al. 2013, Table I |
| H I Hα: span of the 7 fine-structure components | 6562.70970 | 6562.90944 | air | +199.740 | +9124 | [lit] NIST ASD |
| Hα: H → D, composite Ritz values | 6562.819 | 6561.012 | air | −1807.0 | −82545 | [lit] NIST ASD (H I, D I) |

**Lithium 6708 Å.**

* Brown et al. (2013, PRA 87, 032504) give the D2 centres of gravity: ⁶Li 446 799 648.870(15) MHz and ⁷Li 446 810 183.163(16) MHz [lit]. The isotope shift is 10 534.293(22) MHz; Sansonetti et al. (2011), as corrected in their 2012 erratum, give 10 534.357(29) MHz.
* Vacuum wavelengths from these frequencies are 6709.6156 Å (⁷Li) and 6709.7738 Å (⁶Li). Converted to air with Peck & Reeder eq. 3 they are 6707.7638 Å and 6707.9219 Å [derived].
* NIST ASD lists natural Li I as **"6707.76" (D2) and "6707.91" (D1), each ± 0.06 Å, with only 2 decimals** [lit, ASD 2026-11-23].
* The isotope shift (0.158 Å) is almost the size of the 2p fine-structure splitting (≈ 0.151 Å). So ⁶Li D2 falls within ≈ 0.01 Å of ⁷Li D1. A ⁶Li/⁷Li measurement is a 10-mÅ-scale problem, and a single "6707.8" or even "6707.76" hides it.

**Hα.**

* NIST ASD lists the fine-structure components at 6562.70970, .72483, .75181, .77153, .85177, .86734 and .90944 Å (air). Each Ritz value has an uncertainty of ≈ 3e-5 Å.
* The composite is given as 6562.819 ± 0.007 Å (Ritz) and 6562.79 ± 0.03 Å (observed).
* A "1 mÅ" Hα wavelength is therefore meaningful only for a stated component, or for a stated weighting such as LS-coupling intensities or a source temperature.

**Mercury, barium and calcium in NIST ASD [lit].**

* Natural-Hg 5460.750 ± 0.001 Å (observed) and 5460.7498 ± 0.0009 Å (Ritz).
* Ba II 4554.033 Å (observed) and "4554.0332+" (Ritz).
* Ca II 8542.09 Å (2 decimals).

These are composites over isotopes and hyperfine structure. Hfs-sensitive optical lines of Sc, V, Mn and Co are a well-known abundance trap (Prochaska & McWilliam 2000, ApJ 537, L57). The 4–35 mÅ splittings above show why wiki entries for such lines must say "composite (natural isotopic mix, unresolved hfs)".

### 3.4 Zeeman effect [derived]

The σ-component offset is Δλ = g e B λ²/(4π mₑ c) = 4.669 × 10⁻¹³ g B[G] λ²[Å] Å.

| B | g_eff | λ (Å) | Δλ (mÅ) | v (m/s) |
|---|---|---|---|---|
| 1 G | 1 | 5000 | 0.0117 | 0.70 |
| 100 G | 1 | 5000 | 1.17 | 70 |
| 1 kG | 1 | 5000 | 11.7 | 700 |
| 1 kG | 3.00 (Fe I 5250.2084; [lit] Löhner-Böttcher+2019 Table 1) | 5250.2 | 38.6 | 2204 |
| 3 kG | 3.00 | 5250.2 | 116 | 6613 |
| 1 kG | 1 | 16000 | 120 | 2239 |

A 100 G field already moves the σ components by 1 mÅ. Wavelengths in the wiki are field-free, and the Landé g-factor is the quantity to store alongside them.

### 3.5 Pressure shifts in laboratory sources [lit]

Th lines are "individually stable in wavelength at least within the error bars". Ar lines shift from lamp to lamp by "a few tens of m/s", which Lovis & Pepe (2007) attribute to pressure shifts, citing Whaling et al. (2002). They recommend not using Ar lines when high precision is needed. At 5500 Å this is 0.2–0.9 mÅ. It is a concrete example of a source-condition effect just below the 1 mÅ display digit.

The laboratory reference for Fe I 5250.2084 Å, measured with the LARS hollow cathode and LFC calibration (Löhner-Böttcher et al. 2019), differs from the NIST ASD values (5250.2086 ± 0.0008 Å observed, 5250.2089 ± 0.0004 Å Ritz) by 0.2–0.5 mÅ, i.e. 11–29 m/s [lit, derived difference].

---

## 4. What a 2-decimal catalogue costs

**Rounding alone [derived].** A uniform rounding error has rms = step/√12. At 5000 Å:

| decimals | step | max error | rms error |
|---|---|---|---|
| 1 | 100 mÅ = 6.0 km/s | 50 mÅ | 28.9 mÅ = 1.73 km/s |
| 2 | 10 mÅ = 600 m/s | 5 mÅ = 300 m/s | 2.9 mÅ = 173 m/s |
| 3 | 1 mÅ = 60 m/s | 0.5 mÅ | 0.29 mÅ = 17 m/s |
| 4 | 0.1 mÅ = 6 m/s | 0.05 mÅ | 0.029 mÅ = 1.7 m/s |

**NIST Handbook census (`data/handbook.json`, 11 980 strong lines; Sansonetti & Martin 2005) [derived].**

| decimals | lines |
|---|---|
| 0 | 4 |
| 1 | 43 |
| 2 | **5 122 (42.8 %)** |
| 3 | 5 136 |
| 4 | 1 638 |
| 5 | 37 |

* That is 1 675 lines (14.0 %) at 0.0001 Å or better.
* By medium: 4 877 of the 2-decimal lines are air wavelengths and 245 are vacuum.

**Comparison with NIST ASD [derived].** This uses a snapshot of `data/handbook_asd.json` built by `scripts/match_asd.py`, file time 2026-11-23 17:48. Only exact and close matches are used. The ASD observed or Ritz value is used only when it is quoted with *more* decimals than the Handbook.

| Handbook decimals | matches | ASD finer | median \|Δ\| | 68 % | 95 % | median \|Δv\| | 95 % \|Δv\| | median ASD unc. | ASD unc. ≤ 1 mÅ |
|---|---|---|---|---|---|---|---|---|---|
| 2 | 4895 | 2567 | 4.0 mÅ | 6.0 mÅ | 22 mÅ | 303 m/s | 1.86 km/s | 3.0 mÅ (n = 1475) | 33 % |
| 3 | 4517 | 2014 | 0.5 mÅ | 0.9 mÅ | 4.2 mÅ | 45 m/s | 426 m/s | 0.6 mÅ (n = 1593) | 66 % |
| 4 | 1568 | 929 | 0.17 mÅ | 0.26 mÅ | 0.88 mÅ | 19 m/s | 135 m/s | 0.13 mÅ (n = 1289) | 95 % |

* **2-decimal lines.** Replacing a 2-decimal Handbook value with the best ASD value moves it by a median of 4 mÅ, or about 300 m/s, and by up to about 2 km/s at the 95th percentile.
* **Coarse sources.** For 2328 of the 4895 matched 2-decimal lines, ASD has nothing finer either. For these the wiki must show the coarse value with a "< 3 decimals in source" flag, not a padded "x.xx0".
* **3-decimal lines.** Their median change of 0.5 mÅ confirms that 1 mÅ is the precision level these data actually support.

---

## 5. Where 1 mÅ is not enough, and where it is overkill

### 5.1 Requirements finer than 1 mÅ

| requirement | v | λ (Å) | Δλ (Å) | source |
|---|---|---|---|---|
| EPRV goal | 10 cm/s | 5000 | 1.67e-06 | [lit] Fischer et al. 2016 |
| Earth-twin reflex of the Sun | 9 cm/s | 5000 | 1.50e-06 | [lit] Fischer et al. 2016 |
| Ideal calibrator, single measurement | ~1 cm/s | 5000 | 1.67e-07 | [lit] Fischer et al. 2016 |
| Agreement of two independent LFCs on HARPS | ~1 cm/s | 5000 | 1.67e-07 | [lit] Probst et al. 2020 |
| HARPS Th–Ar line list | ~10 m/s | 5500 | 1.8e-04 | [lit] Lovis & Pepe 2007 |
| Lab rest frequencies for Δα/α ("δv ≲ 20 m/s") | 20 m/s | 2382.76 (Fe II) | 1.59e-04 | [lit] Murphy & Berengut 2014 |
| Line shifts for Δα/α = 1 ppm, 43 MM transitions | −30.0 to +19.4 m/s | 1394–5898 | −2.0e-04 to +1.2e-04 | [derived] Δv = −2cq(Δα/α)/ω₀ with q and ω₀ from Murphy & Berengut 2014 |
| Redshift drift at z = 4 | ~6 cm/s per decade | 6078 (Lyα) | 1.2e-06 per decade | [lit] Liske et al. 2008 |

**Redshift drift with Planck18 [derived].** Here v̇ = c ż/(1+z) and ż = (1+z)H₀ − H(z).

| z | v̇ (cm/s per decade) | λ(Lyα) (Å) | Δλ per decade (Å) | Δλ per 20 yr (Å) |
|---|---|---|---|---|
| 1 | +2.25 | 2431.3 | +1.8e-07 | +3.7e-07 |
| 2 | −0.11 | 3647.0 | −1.3e-08 | −2.6e-08 |
| 3 | −2.80 | 4862.7 | −4.6e-07 | −9.1e-07 |
| 4 | −5.37 | 6078.4 | −1.09e-06 | −2.18e-06 |
| 5 | −7.76 | 7294.0 | −1.89e-06 | −3.78e-06 |

These measurements need 10⁻⁵–10⁻⁷ Å. That is 100–10 000 × finer than the wiki's display precision. They are **differential**: a laser frequency comb or a comb-referenced etalon defines the wavelength scale, and the absolute catalogue wavelength only fixes the zero point.

* Murphy & Berengut (2014) tabulate 6–9 decimals for exactly this reason. Their Mg II isotopic components carry uncertainties of 4.2e-7 Å.
* For such data the wiki must keep the full source string; displaying it rounded to 3 decimals is fine, storing it that way is not.

### 5.2 Requirements coarser than 1 mÅ (1 mÅ is overkill)

| requirement | v | λ (Å) | Δλ (mÅ) | source |
|---|---|---|---|---|
| DESI LRG/ELG random redshift error, σ_z < 0.0005(1+z) | 150 km/s rms | 7000 | 3502 | [lit] DESI Collaboration 2022, L2.2.2 |
| DESI systematic mean-redshift error, < 0.0002(1+z) | 60 km/s | 7000 | 1401 | [lit] idem, L2.2.3 |
| DESI QSO redshift error, 0.0025(1+z) | 750 km/s | 7000 | 17 512 | [lit] idem, L2.4.2 |
| Galaxy velocity dispersions (the driver for BOSS resolution) | 100–200 km/s | 6000 | 2000–4000 | [lit] Smee et al. 2013 |
| 4MOST velocity accuracy | < 1 km/s | 5450 | 18 | [lit] de Jong et al. 2019 |
| LAMOST-MRS RV precision | ~1 km/s | 5150 | 17 | [lit] Liu et al. 2020 |
| Gaia RVS systematic RV error after calibration | ≤ 300 m/s | 8600 | 8.6 | [lit] Cropper et al. 2018 |

* **Galaxy and QSO redshift surveys at R ≈ 2000–5000.** Rest wavelengths good to 0.01 Å, about 0.5 km/s, are 100–1000 × better than needed. 1 mÅ adds nothing.
* **Medium-resolution stellar surveys.** 4MOST, LAMOST-MRS, APOGEE and Gaia RVS need 0.1–1 km/s, i.e. 2–20 mÅ, so 2–3 decimals suffice.
* **Some things still matter for these surveys:**
  * the air/vacuum convention (1.4–6 Å);
  * the velocity convention: radio and optical velocities already differ by 0.3 km/s at z = 0.001 (`relativity_cosmology.md`);
  * multiplicative (1+z) composition (`relativity_cosmology.md`).

---

## 6. Recommendations for the wiki

1. **Store, do not reformat.** Keep the source wavelength string verbatim. Record its medium (vacuum, or standard air with the formula used), its stated uncertainty, the number of decimals, and whether it is a composite (natural isotopic mix, unresolved hyperfine or fine structure).
2. **Display at 3 decimals (1 mÅ) when the source supports it.** Show 4 or more decimals in the detail view when the source gives them. Never pad; flag sources with ≤ 2 decimals.
3. **Show the precision in velocity next to it:**
   * δv = c·(last-digit step)/λ;
   * for example "± 60 m/s at 5000 Å" for 1 mÅ, and "± 600 m/s" for 0.01 Å.

   A fixed number of decimals means very different velocity precision across the spectrum (250 m/s at 1200 Å, 14 m/s at 2.2 µm).
4. **Calibration section.** Once the medium is right, the dominant terms at 1 mÅ are:
   * the state of the air (lab air versus standard air: 26 mÅ);
   * Earth rotation (7.8 mÅ);
   * solar gravitational redshift (10.6 mÅ);
   * convective blueshift (≤ 9 mÅ).

   The formula choice (≤ 0.2 mÅ in the optical) is below the display digit but must be recorded. Below the display digit are:
   * the Sun's reflex about the solar-system barycentre (0.27 mÅ);
   * the Earth–Moon reflex (0.21 mÅ);
   * the relativistic and clock terms (≈ 0.08 mÅ, `velocity_frames.md`).

   Those terms are mandatory only for EPRV.

---

## 7. Sources

Every source below was read from a downloaded copy, or from the page text via the extraction tools where downloading was blocked. Journal, volume and page were checked against Crossref (`data/raw/precision_rationale/crossref/parsed.json`).

**Method**

* Bouchy F., Pepe F., Queloz D. 2001, "Fundamental photon noise limit to radial velocity measurements", A&A 374, 733–739. doi:10.1051/0004-6361:20010730, bibcode 2001A&A...374..733B, https://www.aanda.org/articles/aa/pdf/2001/29/aa1316.pdf. Eqs. 8 and 10–12 were read via the exa page fetch; aanda.org returned HTTP 403 to curl.

**Instruments**

* Mayor M. et al. 2003, "Setting New Standards with HARPS", The Messenger 114, 20. Bibcode 2003Msngr.114...20M, https://www.eso.org/sci/publications/messenger/archive/no.114-dec03/messenger-no114-20-24.pdf. Local copy: `web/www.eso.org_sci_publications_messenger_archive_no.114-dec03_messenger-no114-20-24.pdf`. ESO HARPS overview: https://www.eso.org/sci/facilities/lasilla/instruments/harps/overview.html
* Pepe F. et al. 2021, "ESPRESSO at VLT. On-sky performance and first results", A&A 645, A96. doi:10.1051/0004-6361/202038306, arXiv:2010.00316. ESO ESPRESSO overview: https://www.eso.org/sci/facilities/paranal/instruments/espresso/overview.html
* Blackman R. T. et al. 2020, "Performance Verification of the EXtreme PREcision Spectrograph", AJ 159, 238. doi:10.3847/1538-3881/ab811d, arXiv:2003.08852.
* Jurgenson C. et al. 2016, Proc. SPIE 9908, 99086T. doi:10.1117/12.2233002, arXiv:1606.04413.
* NEID: WIYN NEID characteristics, https://www.wiyn.org/Instruments/wiynneid.html. NEID Observers' and Proposers' FAQ, https://neid.psu.edu/observers-and-proposers-faq. Schwab C. et al. 2016, "Design of NEID, an extreme precision Doppler spectrograph for WIYN", Proc. SPIE 9908, 99087H, doi:10.1117/12.2234411.
* Strassmeier K. G. et al. 2015, "PEPSI: The high-resolution échelle spectrograph and polarimeter for the Large Binocular Telescope", Astron. Nachr. 336, 324–361. doi:10.1002/asna.201512172, arXiv:1505.06492.
* UVES: ESO pages https://www.eso.org/sci/facilities/paranal/instruments/uves/inst.html and https://www.eso.org/sci/facilities/paranal/instruments/uves/overview.html. Dekker H. et al. 2000, Proc. SPIE 4008, 534, doi:10.1117/12.395512.
* Keck HIRES: https://www2.keck.hawaii.edu/inst/hires/slitres.html and https://www2.keck.hawaii.edu/inst/hires/binning.html. Vogt S. S. et al. 1994, Proc. SPIE 2198, 362–375, doi:10.1117/12.176725.
* Subaru HDS: https://subarutelescope.org/Instruments/HDS/ and https://subarutelescope.org/en/about/instrument/hds/index.html. Noguchi K. et al. 2002, PASJ 54, 855–864, doi:10.1093/pasj/54.6.855.
* Dorn R. J. et al. 2023, "CRIRES⁺ on sky at the ESO Very Large Telescope", A&A 671, A24. doi:10.1051/0004-6361/202245217, arXiv:2301.08048. ESO CRIRES overview: https://www.eso.org/sci/facilities/paranal/instruments/crires/overview.html
* IGRINS: Gemini page https://www.gemini.edu/instrumentation/igrins and McDonald page https://mcdonald.utexas.edu/research-facilities/HJST/IGRINS. Park C. et al. 2014, Proc. SPIE 9147, 91471D, doi:10.1117/12.2056431. Le H. A. N. et al. 2015, Adv. Space Res. 55, 2509–2518, doi:10.1016/j.asr.2015.03.007, arXiv:1501.03249.
* Wilson J. C. et al. 2019, "The APOGEE Spectrographs", PASP 131, 055001. doi:10.1088/1538-3873/ab0075, arXiv:1902.00928.
* Cropper M. et al. 2018, "Gaia Data Release 2: Gaia Radial Velocity Spectrometer", A&A 616, A5. doi:10.1051/0004-6361/201832763, arXiv:1804.09369.
* DESI Collaboration (Abareshi B. et al.) 2022, "Overview of the Instrumentation for the Dark Energy Spectroscopic Instrument", AJ 164, 207. doi:10.3847/1538-3881/ac882b, arXiv:2205.10939.
* Smee S. A. et al. 2013, "The Multi-object, Fiber-fed Spectrographs for SDSS and BOSS", AJ 146, 32. doi:10.1088/0004-6256/146/2/32, arXiv:1208.2233.
* LAMOST DR11 v1.1 LRS data-production description, https://www.lamost.org/dr11/v1.1/doc/lr-data-production-description. The page was read through tavily_extract because the site timed out to curl. Cui X.-Q. et al. 2012, RAA 12, 1197–1242, doi:10.1088/1674-4527/12/9/003. Liu C. et al. 2020, "LAMOST-MRS: scientific goals and survey plan", arXiv:2005.07210.
* de Jong R. S. et al. 2019, "4MOST: Project overview and information for the First Call for Proposals", The Messenger 175, 3. doi:10.18727/0722-6691/5117 (from the arXiv metadata; Crossref returns 404 because it is not a Crossref DOI), arXiv:1903.02464.
* Subaru PFS instrument parameters, https://pfs.ipmu.jp/research/parameters.html (values as of July 2024; page updated 2026-10-21). Tamura N. et al. 2016, Proc. SPIE 9908, 99081M, doi:10.1117/12.2232103, arXiv:1608.01075.
* Jin S. et al. 2024, "The wide-field, multiplexed, spectroscopic facility WEAVE", MNRAS 530, 2688–2730. doi:10.1093/mnras/stad557, arXiv:2212.03981.

**Calibration, laboratory and solar data**

* Lovis C., Pepe F. 2007, "A new list of thorium and argon spectral lines in the visible", A&A 468, 1115–1121. doi:10.1051/0004-6361:20077249, arXiv:astro-ph/0703412.
* Murphy M. T., Berengut J. C. 2014, "Laboratory atomic transition data for precise optical quasar absorption spectroscopy", MNRAS 438, 388–411. doi:10.1093/mnras/stt2204, arXiv:1311.2949. Text at `data/raw/relativity_cosmology/arxiv_pdf/1311.2949.txt`; q values in `data/notes/murphy_berengut_2014_mm_transitions.json`.
* Brown R. C. et al. 2013, PRA 87, 032504. doi:10.1103/PhysRevA.87.032504. Local copy `papers/brown2013_pra87_032504.pdf`, from https://www.agpohl.physik.uni-mainz.de/files/2022/01/PRA_Brown.pdf.
* Sansonetti C. J. et al. 2011, PRL 107, 023001, doi:10.1103/PhysRevLett.107.023001. Erratum: PRL 109, 259901 (2012), doi:10.1103/PhysRevLett.109.259901. The erratum table was read via exa.
* Li R. et al. 2020, PRL 124, 063002 (⁶Li D lines), doi:10.1103/PhysRevLett.124.063002. Local copy `papers/li6_dlines_ecnu.pdf`.
* NIST Atomic Spectra Database (Kramida A., Ralchenko Yu., Reader J. and NIST ASD Team), https://physics.nist.gov/asd. Windows queried 2026-11-23 for H I, D I, Li I, Hg I, Ba II, Ca II and Fe I; the exact URLs are in `nist_asd/*.csv.url`.
* Sansonetti J. E., Martin W. C. 2005, "Handbook of Basic Atomic Spectroscopic Data", J. Phys. Chem. Ref. Data 34, 1559. doi:10.1063/1.1800011.
* Reiners A. et al. 2016, "The IAG solar flux atlas: Accurate wavelengths and absolute convective blueshift in standard solar spectra", A&A 587, A65. doi:10.1051/0004-6361/201527530, arXiv:1511.03014. Text at `data/raw/relativity_cosmology/arxiv_pdf/`.
* Löhner-Böttcher J. et al. 2019, "Convective blueshifts in the solar atmosphere III", A&A 624, A57. doi:10.1051/0004-6361/201834925, arXiv:1901.07606.
* González Hernández J. I. et al. 2020, "The solar gravitational redshift from HARPS-LFC Moon spectra", A&A 643, A146. doi:10.1051/0004-6361/202038937, arXiv:2009.10558.
* Fischer D. A. et al. 2016, "State of the Field: Extreme Precision Radial Velocities", PASP 128, 066001. doi:10.1088/1538-3873/128/964/066001, arXiv:1602.07939.
* Probst R. A. et al. 2020, "A crucial test for astronomical spectrograph calibration with frequency combs", Nature Astronomy 4, 603–608. doi:10.1038/s41550-020-1010-x.
* Li C.-H. et al. 2008, Nature 452, 610–612, doi:10.1038/nature06854.
* Liske J. et al. 2008, "Cosmic dynamics in the era of Extremely Large Telescopes", MNRAS 386, 1192–1218. doi:10.1111/j.1365-2966.2008.13090.x, arXiv:0802.1532.
* Prochaska J. X., McWilliam A. 2000, ApJ 537, L57–L60, doi:10.1086/312749. Cited for context only; its text was not read.
* Planck Collaboration 2020, "Planck 2018 results. VI", A&A 641, A6, doi:10.1051/0004-6361/201833910. Used through `astropy.cosmology.Planck18`.

**Refractive index of air (formulas in `scripts/research/airvac_formulas.py`)**

* Edlén B. 1953, JOSA 43, 339, doi:10.1364/JOSA.43.000339.
* Edlén B. 1966, Metrologia 2, 71, doi:10.1088/0026-1394/2/2/002.
* Peck E. R., Reeder K. 1972, JOSA 62, 958, doi:10.1364/JOSA.62.000958.
* Birch K. P., Downs M. J. 1994, Metrologia 31, 315, doi:10.1088/0026-1394/31/4/006.
* Ciddor P. E. 1996, Appl. Opt. 35, 1566, doi:10.1364/AO.35.001566.

---

## 8. Caveats

* **Centroid precisions are photon-noise limits** for a single, isolated, unresolved Gaussian line at SNR 100 per pixel. They exclude:
  * read noise (important at SNR ≲ 10);
  * blends, intrinsic line width and asymmetry, and LSF mismatch;
  * flat-field and pixel-size errors (HARPS floor ≈ 0.01 px);
  * wavelength-solution errors.

  Real per-line accuracies are therefore worse. RV precision of ensembles improves as √N_lines only against photon noise, not against systematics.
* **Assumed and derived instrument values.**
  * DESI R values are mid-points of the quoted ranges.
  * LAMOST-MRS and WEAVE samplings (3 px) are assumptions; their papers give no sampling.
  * UVES and CRIRES+ samplings are slit-width/pixel-scale estimates.
  * HIRES sampling comes from the binning page.
  * LAMOST-LRS sampling is on the 10⁻⁴-dex output grid, not on the CCD.
* **Reference wavelengths are conventions.** Each instrument row uses one wavelength. FWHM scales with λ across each instrument's range.
* **R conventions differ.** "R" is design or median as quoted. HARPS has 115 000 (Mayor 2003) against 110 000 measured by Lovis & Pepe (2007). ESPRESSO HR has 140 000 (ESO) against 138 000 (Pepe 2021 Table 3). IGRINS has 45 000 (Gemini) against 40 000 used in the ETC paper.
* **The Handbook–ASD statistics are a snapshot.** They use `data/handbook_asd.json` as of 2026-11-23 17:48, which another process may still update. Rerun the script to refresh them.
* **Air–vacuum comparisons use each formula's own standard air.** Edlén 1953 and 1966 specify 300 ppm CO₂; Peck & Reeder about 300–330 ppm; Birch & Downs and Ciddor 450 ppm. Ciddor and Peck & Reeder eq. 2 are extrapolated at 2000 Å and 22 000 Å.
* **The Li vacuum wavelengths are derived** as λ = c/ν from the published frequencies; the air values use Peck & Reeder eq. 3, the formula NIST ASD uses.
* **Bibcodes.** Only the Bouchy et al. (2001) and Mayor et al. (2003) bibcodes are given. The others were not checked against ADS, so DOIs are the primary identifiers.
