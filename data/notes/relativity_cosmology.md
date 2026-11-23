# Relativity and cosmology in wavelength / redshift calibration

Research notes for the Spectral Line Wiki, wavelength-calibration section.

**Provenance rules used in these notes**

* **[lit]** marks a number quoted from a paper. The paper is cited, and the text was read from the arXiv PDF or a Crossref/publisher abstract stored under `data/raw/relativity_cosmology/`.
* **[derived]** marks a number computed by the scripts in `scripts/research/` from the stated formula. Constants come from `astropy.constants` (IAU 2015 Res. B3 nominal values, IAU 2012 au), ERFA (`erfa.ELG`, `erfa.ELB`), astropy `Planck18`, or the DE440s kernel at `~/JPLDE/de440s.bsp`.
* Wavelength strings copied from catalogues are verbatim and are **vacuum** unless stated.

**Scripts and outputs**

| file | content |
|---|---|
| `scripts/research/relcos_numbers.py` → `data/notes/relativity_cosmology_numbers.json` | every [derived] number in Sections 1–8 |
| `scripts/research/relcos_ephemeris.py` → `data/notes/relativity_cosmology_ephemeris.json` | DE440s magnitudes: Sun about the SSB, lunar reflex, observer clock term |
| `scripts/research/parse_murphy_berengut_2014.py` → `data/notes/murphy_berengut_2014_mm_transitions.json` | the 43 many-multiplet transitions of Murphy & Berengut (2014): ω0, λ0 (vacuum), δv_lab, q, and Δv for Δα/α = 1 ppm |
| `scripts/research/fetch_arxiv_meta*.py`, `parse_arxiv_meta.py`, `fetch_arxiv_pdfs.sh`, `crossref_lookup.py` | bibliographic metadata and full texts (raw data under `data/raw/relativity_cosmology/`) |

---

## 1. Why 1e-3 Å? The velocity scale of the wavelength grid

The velocity equivalent of a wavelength step is

$$\delta v = c\,\frac{\delta\lambda}{\lambda}, \qquad \delta\lambda = \lambda\,\frac{\delta v}{c}.$$

[derived] **1e-3 Å corresponds to:**

| λ (Å) | δv for 1e-3 Å | δv for 1e-4 Å | δλ for 20 m/s | δλ for 1 m/s |
|---|---|---|---|---|
| 1215.67 (Lyα) | 246.6 m/s | 24.7 m/s | 8.1e-5 Å | 4.1e-6 Å |
| 1608.451 (Fe II) | 186.4 m/s | 18.6 m/s | 1.07e-4 Å | 5.4e-6 Å |
| 2796.354 (Mg II) | 107.2 m/s | 10.7 m/s | 1.87e-4 Å | 9.3e-6 Å |
| 5000 | 60.0 m/s | 6.0 m/s | 3.34e-4 Å | 1.67e-5 Å |
| 6562.8 (Hα) | 45.7 m/s | 4.6 m/s | 4.38e-4 Å | 2.19e-5 Å |
| 10830 (He I) | 27.7 m/s | 2.8 m/s | 7.2e-4 Å | 3.6e-5 Å |
| 16000 (H band) | 18.7 m/s | 1.9 m/s | 1.07e-3 Å | 5.3e-5 Å |
| 21655 (Brγ) | 13.8 m/s | 1.4 m/s | 1.45e-3 Å | 7.2e-5 Å |

**Ladder of physical shifts.** Each shift is expressed as a velocity, and as a Δλ at 5000 Å, where Δλ = v/c × 5000 Å.

| effect | size | Δλ at 5000 Å | source |
|---|---|---|---|
| air vs vacuum (n − 1 ≈ 3e-4) | ≈ 90 km/s | ≈ 1.5 Å | [lit] Davis et al. 2011 footnote 21 (uses n = 1.0003) |
| Earth orbital velocity | 29.26–30.31 km/s (mean 29.78) | ≈ 0.5 Å | [derived] DE440s |
| Earth rotation at the equator | 465.1 m/s | 7.8 mÅ | [derived] ω_E R_E |
| solar gravitational redshift (surface → ∞) | 636.31 m/s | 10.6 mÅ | [derived] GM☉/(R☉c) |
| solar convective blueshift | about −550 to 0 m/s, depending on line depth | up to −9 mÅ | [lit] Reiners et al. 2016 |
| lunar reflex of the Earth (Earth about the Earth–Moon barycentre) | 11.7–13.4 m/s | 0.2 mÅ | [derived] DE440s |
| Sun about the solar-system barycentre (heliocentric vs barycentric frame) | 8.5–15.9 m/s (2000–2040) | 0.2 mÅ | [derived] DE440s |
| transverse Doppler of a star with v_t = 100 km/s | 16.7 m/s | 0.28 mÅ | [derived] v²/2c |
| observer clock (solar + terrestrial potential + Earth's speed) | 4.65 m/s mean, 0.20 m/s peak-to-peak | 0.08 mÅ | [derived]; [lit] Lindegren & Dravins 2003 |
| cross term of a multiplicative barycentric correction | ≈ 3 m/s | 0.05 mÅ | [lit] Wright & Eastman 2014 |
| Δα/α = 1 ppm in the many-multiplet method | −30 to +19 m/s per line | ≈ 2e-4 Å (rest frame, λ ≈ 2000–2600 Å) | [derived] from Murphy & Berengut 2014 q values |
| requirement on lab rest frequencies for α work | δv ≲ 20 m/s | ≈ 1.3e-4 Å at 2000 Å | [lit] Murphy & Berengut 2014 |
| Earth-twin exoplanet signal | ≈ 9 cm/s | 1.5e-6 Å | [lit] Wright & Eastman 2014 |
| redshift drift (Sandage–Loeb), Planck18 | 0.1–5.4 cm/s per decade for 1 < z < 4 | ≤ 1e-6 Å over 20 yr | [derived] |

**Conclusion.** A wavelength scale quoted to **1e-3 Å** is precise to about 15–250 m/s across the UV to near-IR, and to about 60 m/s at 5000 Å. That precision has three consequences:

1. It is the **natural precision of good laboratory and Ritz wavelengths**. Lovis & Pepe (2007) give the typical internal uncertainty of HARPS Th-Ar lines as ≈ 10 m/s, or 0.18 mÅ. Redman et al. (2014) find a standard deviation of 8.2e-5 nm (0.82 mÅ) between Palmer & Engleman (1983) and their Ritz wavelengths. The best Ritz wavelengths reach about 2e-5 nm, or 0.2 mÅ.
2. It **resolves every astrophysical and relativistic shift above about 50 m/s**. These include the solar gravitational redshift (636 m/s, about 10.6 mÅ), convective blueshifts (hundreds of m/s), Earth rotation (465 m/s) and the radio/optical/relativistic velocity-convention differences. A coarser grid (0.01 Å = 600 m/s at 5000 Å) blurs the gravitational redshift and the convective shift, so 1e-3 Å is about the minimum precision at which line positions are physically meaningful at the 100 m/s level.
3. It is **not sufficient** for:
   * Many-multiplet α/μ work, which needs ≲ 20 m/s, i.e. ≈ 1–2e-4 Å in the rest-frame UV. Murphy & Berengut (2014) therefore tabulate 6–7 decimals.
   * Exoplanet radial velocities at 1 m/s (≈ 1.7e-5 Å).
   * Redshift drift at cm/s (≈ 1e-7 Å).

   These measurements are differential and rely on an instrumental calibrator (Section 10), not on the absolute rest wavelength of each line.

When a source gives fewer than three decimals in Å, the wiki keeps the original string and flags it. Examples are OH airglow lines computed from energy levels, "secure to about one hundredth of Å" (Rousselot et al. 2000), i.e. ≈ 190 m/s at 1.6 µm, and many Fraunhofer-line compilations.

---

## 2. Redshifts compose multiplicatively, not additively

Each process maps a frequency to a frequency, so the factors (1+z) multiply:

$$1+z_{\rm obs} = (1+\bar z_{\rm cos})\,(1+z_{\rm pec})\,(1+z_{\rm grav})\,(1+z_{\odot})\,(1+z_{\oplus})\cdots$$

* Davis et al. (2011, ApJ 741, 67) give this as eq. 1: (1+z) = (1+z̄)(1+z_pec)(1+z_grav).
* Carr et al. (2022, PASA 39, e046) give eqs. 2, 3, 6, 7:
  * $(1+z_{\rm hel})=(1+z^{\times}_{\rm CMB})(1+z_{\rm Sun})$, so $z^{\times}_{\rm CMB} = (1+z_{\rm hel})/(1+z_{\rm Sun})-1$.
  * $(1+z_{\rm hel})=(1+z_{\rm HD})(1+z_{\rm Sun})(1+z_p)$ and $z_{\rm HD}=(1+z_{\rm CMB})/(1+z_p)-1$.
  * $z_{\rm Sun}=\sqrt{(1-v_{\rm Sun}/c)/(1+v_{\rm Sun}/c)}-1$, where $v_{\rm Sun}=\mathbf v_{\rm Sun}\cdot\hat n_{\rm obj}=v^{\max}_{\rm Sun}\cos\alpha$. The Sun's motion toward an object makes it appear blueshifted.
* Wright & Eastman (2014, PASP 126, 838) give eq. 10, bary(z_meas) = (1+z_meas)(1+z_B) − 1, and eq. 11, z_true − z_meas = z_B + z_B z_meas. They stress that "the correction to z_meas is multiplicative, not additive". For stellar RVs the cross term z_B z_meas is ≈ 3 m/s.

**The additive approximation** is $z_{\rm hel}= z^{+}_{\rm CMB}+z_{\rm Sun}$ (Carr et al. eq. 1). Its error is exactly $z_{\rm CMB}\,z_{\rm Sun}$, and $|z_{\rm Sun}|\le v_{\rm Sun}/c = 1.234\times10^{-3}$ for v_Sun = 369.82 km/s.

[derived] Error of the additive inversion $z_{\rm CMB}\approx z_{\rm hel}+\beta\cos\theta$:

| z_CMB | object at the apex (cos θ = +1) | object at the antapex (cos θ = −1) | at 90° (cos θ = 0, pure γ term) |
|---|---|---|---|
| 0.01 | −1.16e-5 | +1.31e-5 | −7.7e-7 |
| 0.1 | −1.23e-4 | +1.24e-4 | −8.4e-7 |
| 0.5 | −6.16e-4 | +6.18e-4 | −1.1e-6 |
| 1.0 | −1.23e-3 | +1.24e-3 | −1.5e-6 |

[derived] The same error for a peculiar velocity of 300 km/s, $(1+z_c)(1+z_p)-1-(z_c+v/c)$: 1.05e-5 at z = 0.01, 1.0e-4 at z = 0.1, 1.0e-3 at z = 1.

**How large are such errors in cosmology?** [lit]

* Carr et al. 2022: "by z ∼ 1, the error is on the order of 10⁻³, which is an order of magnitude larger than most reported statistical uncertainties in redshifts". Their full set of redshift updates changed $H_0$ by −0.11 km s⁻¹ Mpc⁻¹ and w by −0.001. In their sample, the replaced redshifts changed by between 1e-8 and 0.038, with RMS ≈ 3e-3.
* Calcino & Davis 2017 (JCAP 01, 038): "a small systematic redshift shift (Δz ∼ 10⁻⁵) … can cause a significant bias (∼1%) in the recovery of cosmological parameters".
* Davis et al. 2011: neglecting the CMB dipole gives Δw ≈ 0.04 at fixed Ω_m. Neglecting coherent peculiar velocities gives a ≈ 2% systematic shift in w.
* Davis et al. 2011, Appendix A: the peculiar-velocity error must be propagated as σ_z = (1+z̄)σ_zpec (eq. A1). A redshift error turns into a magnitude error σ_µ = σ_z (5/ln10)[1/(1+z̄) + c/(χ̃H(z̄))] (eq. A3). The empty-universe form is σ_µ ≈ σ_z (5/ln10)(1+z̄)/[z̄(1+z̄/2)] (eq. A4).

[derived] Distance-modulus change for a redshift error δz, computed as Planck18 μ(z+δz) − μ(z):

| z | δz = 1e-4 | δz = 1e-5 |
|---|---|---|
| 0.01 | 0.0218 mag | 0.0022 mag |
| 0.023 | 0.0096 mag | 0.00096 mag |
| 0.1 | 0.0023 mag | 0.00023 mag |
| 0.5 | 0.00052 mag | 5.2e-5 mag |
| 1.0 | 0.00027 mag | 2.7e-5 mag |

**Do not use special-relativistic formulae for cosmological redshift.** Davis et al. 2011 (footnote 20) state that the SR velocity–redshift relation "is only appropriate for peculiar velocities … special relativistic corrections should never be applied" to recession velocities. Davis & Lineweaver (2004) rule out the SR Doppler interpretation of cosmological redshifts at 23σ using supernova magnitudes. Davis & Scrimgeour (2014) show that the common approximation v_p = cz − H₀D overestimates peculiar velocities by ≈ 100 km/s at z ≈ 0.04. The correct form (their eq. 16) is $v_p = c\,(z-\bar z)/(1+\bar z)$.

**Recipe for the wiki.**

* Store redshifts, not velocities, and combine them as products of (1+z).
* Convert a frame velocity v to $(1+z)=\sqrt{(1+\beta)/(1-\beta)}$ along the line of sight. For an oblique direction use γ(1+β cos θ).

---

## 3. Velocity conventions and reference frames (optical, radio, relativistic, FITS)

Greisen et al. (2006, A&A 446, 747; "FITS Paper III") define the following, with a receding object positive:

$$V_{\rm radio}=c\,\frac{\nu_0-\nu}{\nu_0}=c\,\frac{\lambda-\lambda_0}{\lambda},\qquad Z_{\rm opt}=c\,\frac{\lambda-\lambda_0}{\lambda_0}=cz,$$
$$v_{\rm app}=c\,\frac{\nu_0^2-\nu^2}{\nu_0^2+\nu^2}=c\,\frac{\lambda^2-\lambda_0^2}{\lambda^2+\lambda_0^2}\quad\text{(eqs. 14, 16: "apparent radial velocity", VELO)}.$$

The FITS CTYPE codes are FREQ, ENER, WAVN, VRAD (radio), WAVE (vacuum wavelength), VOPT, ZOPT, AWAV (air wavelength), VELO and BETA (v/c). Non-linear algorithm codes include F2W, W2F, V2W, W2A and so on, and GRI/GRA for grisms. Greisen et al. call the apparent radial velocity "itself a conventional velocity" and refer readers to Lindegren & Dravins (2003).

[derived] Three conventions for the same redshift:

| z | cz (optical) | radio cz/(1+z) | relativistic v_app |
|---|---|---|---|
| 0.001 | 299.792 km/s | 299.493 km/s | 299.643 km/s |
| 0.01 | 2997.925 | 2968.242 | 2982.936 |
| 0.023 | 6895.227 | 6740.202 | 6815.952 |
| 0.1 | 29979.246 | 27253.860 | 28487.066 |
| 1.0 | 299792.458 | 149896.229 | 179875.475 |

At z = 0.001 the radio and optical conventions already differ by 0.3 km/s. Quoting a "velocity" without its convention is therefore an error of order cz²/c.

**Transverse motion** (Greisen et al. eq. 2, after Lang 1974): $\lambda=\lambda_0 (c+v_r)/\sqrt{c^2-v_r^2-v_t^2}$. Greisen et al. state that "a transverse velocity of 100 km/s produces an apparent velocity of 15 m/s". [derived] v_t²/2c = 16.7 m/s.

**Spectral reference frames.** [lit] Greisen et al. 2006, Table 12, gives the SPECSYS values with conventional magnitudes:

| SPECSYS | definition | magnitude | apex (l, b) | reference given |
|---|---|---|---|---|
| TOPOCENT | topocentric | 0 | – | |
| GEOCENTR | geocentric | 0.5 km/s | – | |
| BARYCENT | barycentric | 30 km/s | – | Stumpff (1980) |
| HELIOCEN | heliocentric | 30 km/s | – | Stumpff (1980) |
| LSRK | kinematic LSR | 20 km/s | (56, 23); normally quoted as RA 18h, Dec +30° (1900) | Delhaye (1965); Gordon (1976) |
| LSRD | dynamical LSR | 16.6 km/s | (53, 25) | Delhaye (1965) |
| GALACTOC | Galactocentric | 220 km/s | (90, 0) | Kerr & Lynden-Bell (1986) |
| LOCALGRP | Local Group | 300 km/s | (90, 0) | de Vaucouleurs (1976) |
| CMBDIPOL | CMB dipole | 368 km/s | (263.85, 48.25) | Bennett et al. (2003) |
| SOURCE | source rest frame | any | – | |

* The radio LSR convention is Gordon (1976), "Computer programs for radio astronomy", in *Methods of Experimental Physics* 12C, ed. M. L. Meeks (Academic Press, New York), p. 277, as cited by Greisen et al.
* Modern solar motion relative to the LSR: (U, V, W)☉ = (11.1 ± 0.74, 12.24 ± 0.47, 7.25 ± 0.37) km/s, with systematic uncertainties of about (1, 2, 0.5) km/s (Schönrich, Binney & Dehnen 2010).
* Modern CMB dipole: 369.82 ± 0.11 km/s toward (l, b) = (264.021° ± 0.011°, 48.253° ± 0.005°) (Planck 2018, as adopted by Carr et al. 2022). The earlier COBE value was 371 ± 1 km/s toward (264.14°, 48.26°) (Fixsen et al. 1996).
* The conventional FITS magnitudes are fixed definitions. A quoted LSRK velocity means "20 km/s toward 18h, +30° (B1900)", not the modern solar motion, so the definition used must be recorded.

**Heliocentric versus barycentric.** [derived] The Sun moves 8.5–15.9 m/s about the solar-system barycentre (DE440s, 2000–2040), so the two frames differ by up to ≈ 16 m/s, or ≈ 0.27 mÅ at 5000 Å. IAU 2000 Resolution C1 specifies the barycentre, not the Sun.

**Other points from Greisen et al. (2006):**

* A frame correction "is a true velocity and thus the relativistic formula … should be used".
* Differential Doppler errors across a field are 8.7 m/s for a 1° offset (diurnal, worst case), 17.5 m/s over 12 h, and 1.0 km/s for exposures months apart. "They may amount to 5 km/s for the Local Group correction over a field of one degree."

---

## 4. Gravitational redshift and convective shifts

### 4.1 Weak-field formulae

The weak-field shift is $v_g \simeq GM/(Rc)$. Between emitter potential Φ_e and observer potential Φ_o:

$$1+z_{\rm grav}=\frac{1-\Phi_o/c^2}{1-\Phi_e/c^2}\;\;\Rightarrow\;\; c z_{\rm grav}\simeq\frac{\Phi_o-\Phi_e}{c}.$$

For a static Schwarzschild emitter observed from infinity, $1+z=(1-2GM/Rc^2)^{-1/2}$.

### 4.2 The Sun

* [derived] With IAU 2015 nominal GM☉ = 1.3271244e20 m³ s⁻² and R☉ = 6.957e8 m:
  * GM☉/(R☉c) = **636.31 m/s**. The Schwarzschild form gives 636.312 m/s.
  * Observer at 1 au: −GM☉/(1 au · c) = −2.959 m/s, giving **633.35 m/s**.
  * Earth's potential: L_G c = W₀/c = 0.209 m/s, giving **633.14 m/s** at the geoid. GM_E/(R_E c) = 0.208 m/s also gives 633.14 m/s. González Hernández et al. quote 633.10 m/s for the same three terms; the 0.04 m/s difference presumably comes from rounding or their constants.
  * An Earth observer also has transverse-Doppler blueshift v²/2c = 1.48 m/s relative to the Sun, since the Sun lies roughly perpendicular to the orbital velocity. The barycentric correction normally removes that term.
* [lit] Lindegren & Dravins 2003, using R☉ = 6.95508e8 m: 636.5 m/s to infinity and 633.5 m/s at the Earth's mean distance (their footnote 3 gives 636.486 ± 0.024 m/s). A chromospheric line formed at 30 Mm has about 20 m/s less; a coronal line perhaps 100 m/s less.
* [lit] González Hernández et al. 2020, HARPS-LFC spectra of the Moon:
  * Theory: 636.31 → 633.35 → 633.10 m/s.
  * 15 Fe I lines with 150 < EW < 550 mÅ give a measured solar gravitational redshift of **639 ± 14 m/s**. In such strong lines the cores are "not affected" by convection.
  * About 97 weak lines (EW < 180 mÅ), fitted with 3D CO5BOLD profiles, give 638 ± 6 m/s.
  * Earlier full-disk values quoted there: 612 ± 58 m/s (Allende Prieto & García López 1998), 601 ± 43 m/s (Molaro & Monai 2012), 698 ± 133 m/s (Takeda & Ueno 2012), and 600.4 ± 0.8 m/s from the K I 7699 line, 1976–2013 (Roca Cortés & Pallé 2014).

### 4.3 Stars and compact objects

* [lit] Lindegren & Dravins 2003: the shift scales as (M/M☉)(R/R☉)⁻¹. It ranges from white dwarfs ("some 30 km/s") to supergiants ("some 30 m/s"). For single stars it is normally not possible to compute it to better than 50 m/s.
* [derived] Illustrative values of 636.3 m/s × M/R (inputs are illustrative, not catalogue values):

  | object | M | R | v_g |
  |---|---|---|---|
  | M dwarf | 0.2 | 0.22 | 578 m/s |
  | subgiant | 1.2 | 2 | 382 m/s |
  | giant | 1.5 | 10 | 95 m/s |
  | supergiant | 15 | 500 | 19 m/s |
  | white dwarf | 0.6 | 0.0125 | 30.5 km/s |
  | white dwarf | 1.0 | 0.008 | 79.5 km/s |

  A 1.4 M☉, 12 km neutron star has z = 0.235 (strong field).
* [lit] White dwarfs:
  * Mean for 449 DA white dwarfs: ⟨v_g⟩ = **32.57 ± 1.17 km/s**, corresponding to ⟨M⟩ = 0.647 M☉ (Falcon et al. 2010).
  * Sirius B: **80.65 ± 0.77 km/s**, measured relative to Hα of Sirius A, giving M = 1.017 ± 0.025 M☉ (Joyce et al. 2018).
  * The classic survey is Greenstein & Trimble (1967).
* [lit] Solar-type dwarfs are expected to be redshifted about 500 m/s more than giants. Pasquini et al. (2011) found that M67 dwarfs and giants share the same apparent RV to within ≈ 100 m/s. Their explanation is that stronger lines are more blueshifted in dwarfs, so convection compensates the gravitational redshift.
* [lit] S2 near Sgr A*: the combined gravitational redshift and transverse Doppler effect is z ≈ 200 km/s / c at pericentre, where the orbital speed is ≈ 7650 km/s. The fitted parameter is f = 0.90 ± 0.09 (stat) ± 0.15 (sys) (GRAVITY Collaboration 2018).

### 4.4 Galaxy clusters and large-scale structure

* [lit] Wojtak, Hansen & Hjorth (2011) estimate cz_gr ≈ 10 km/s for ~10¹⁴ M☉ clusters, after refs. 8–10 of their paper. Within a 6 Mpc aperture they measure Δ = **−7.7 ± 3.0 km/s**, i.e. galaxies are blueshifted relative to the brightest cluster galaxy. Simulations predict −(5–10) km/s.
* [derived] GM/(Rc) for 10¹⁵ M☉ within 1 Mpc is 14.3 km/s (order of magnitude only).
* [lit] Calcino & Davis 2017 note that gravitational redshifts of order 4e-5 have been measured in galaxy clusters, and that any monopole Δz ~ 1e-5 biases supernova cosmology at the 1% level.

### 4.5 The observer's potential

* [derived] The Earth's surface potential is absorbed in TT through L_G = 6.969290134e-10 (IAU 2000 B1.9; erfa.ELG), equivalent to 0.209 m/s.
* [derived] Observatory altitude adds g h/c ≈ 3.3e-5 m/s per km: 7.9e-5 m/s at 2400 m and 1.4e-4 m/s at 4200 m. This is negligible.
* [lit] Pound & Rebka (1960) is the laboratory demonstration of gravitational redshift.

### 4.6 Convective blueshift (Sun and stars)

Photospheric absorption lines are blueshifted because bright rising granules contribute more photons than the dark sinking lanes (Dravins 1982; Dravins, Lindegren & Nordlund 1981).

* [lit] Lindegren & Dravins 2003:
  * The solar shift is ≈ 400 m/s after removing the gravitational redshift.
  * 3D models predict about 1000 m/s for F stars and about 200 m/s for K dwarfs.
  * Differential shifts between strong and weak solar lines in the visible are ≈ 200 m/s, reaching ≈ 1000 m/s for Procyon.
  * Shifts can change sign to convective redshift in the cores of very strong lines and in the far UV.
* [lit] Reiners et al. 2016 (IAG FTS atlas, accurate to ±10 m/s over 405–1065 nm), Table 5, absolute convective blueshift of 1252 Fe I lines after subtracting 636 m/s:

  | relative line depth | mean shift (m/s) | σ (m/s) |
  |---|---|---|
  | 0.1 | −547 | 303 |
  | 0.2 | −513 | 199 |
  | 0.3 | −515 | 232 |
  | 0.4 | −491 | 155 |
  | 0.5 | −450 | 284 |
  | 0.6 | −419 | 224 |
  | 0.7 | −334 | 215 |
  | 0.8 | −207 | 193 |
  | 0.9 | +2 | 128 |

  The fitted polynomial (eq. 2) is Δv = −504.891 − 43.7963x − 145.560x² + 884.308x³ m/s, for depth x between 0.05 and 0.95. The PDF text layer drops the "x" on the linear term.
* [lit] Kitt Peak FTS atlases show wavelength-scale deviations of several tens to 100 m/s (Reiners et al. 2016). Kurucz et al. (1984) are accurate to about 10 m/s only in 530–570 nm.
* [lit] Löhner-Böttcher et al. 2018 (LARS, LFC-calibrated to 1 m/s, Fe lines at 6302 Å): at disk centre the bisector convective blueshift reaches **−450 m/s** in the line wings. Toward the limb the line cores saturate at a redshift of +100 m/s. The 2019 paper on 26 lines confirms the systematic dependence on line depth.
* Consequence for the wiki: the solar (Fraunhofer) wavelengths measured in the solar spectrum differ from laboratory wavelengths by +636 m/s (gravitational) plus −0 to −550 m/s (convective). At 5000 Å this is roughly **+1 to +10 mÅ**, well above the 1e-3 Å precision of the page. Solar-atlas wavelengths therefore must not be mixed with laboratory rest wavelengths, and each Fraunhofer entry must state which one it gives.

---

## 5. Kinematic relativistic terms: time dilation, Shapiro delay, aberration, perspective

### 5.1 Full observer–source relation

Lindegren & Dravins (2003), eqs. 40–43, to O(c⁻²):

$$1+z_{\rm obs}=\Big(1-\frac{\Phi_{\rm obs}}{c^2}-\frac{|\mathbf v_{\rm obs}|^2}{2c^2}\Big)\Big(1+\frac{\mathbf k'\mathbf v_{\rm obs}}{c}\Big)^{-1}\Big(1+\frac{\mathbf k'\mathbf v_*}{c}\Big)\Big(1-\frac{\Phi_*}{c^2}-\frac{|\mathbf v_*|^2}{2c^2}\Big)^{-1},$$

$$1+z_B=(1+z_{\rm obs})\Big(1-\frac{\Phi_{\rm obs}}{c^2}-\frac{|\mathbf v_{\rm obs}|^2}{2c^2}\Big)^{-1}\Big(1+\frac{\mathbf k'\mathbf v_{\rm obs}}{c}\Big)=\Big(1+\frac{\mathbf k'\mathbf v_*}{c}\Big)\Big(1-\frac{\Phi_*}{c^2}-\frac{|\mathbf v_*|^2}{2c^2}\Big)^{-1}(1+X).$$

Here **k** is the *coordinate* (BCRS) direction to the star, "not the observed (aberrated and refracted) direction". X collects the effects that cannot be separated: convection, pressure shifts, cosmological redshift and so on. Wright & Eastman (2014), eqs. 3 and 25, give the same structure with z_GR and the γ factors.

### 5.2 Observer clock term

* [lit] Lindegren & Dravins 2003 eq. 42: for an Earth observer, ⟨Φ_obs⟩ ≈ GM☉/A + GM⊕/R⊕ ≈ 8.934e8 m² s⁻² and ⟨v_obs⟩ ≈ 29785 m/s. The factor $(1-\Phi/c^2-v^2/2c^2)^{-1}\simeq 1+1.550\times10^{-8}$.
* The three contributions are **3.0 m/s** (solar potential), **0.2 m/s** (Earth's potential) and **1.5 m/s** (Earth's velocity). The orbital eccentricity modulates the factor by 2eGM☉/c² ≈ 3.3e-10, or about **0.1 m/s**.
* [derived] DE440s 2000–2040: mean **4.648 m/s**, peak-to-peak 0.203 m/s. The mean equals L_B c = 4.648 m/s (IAU 2006 B3). It splits into L_C c = 4.439 m/s (TCG–TCB) and L_G c = 0.209 m/s.
* [lit] Wright & Eastman 2014: the GR term at the Earth is ≈ 3 m/s with ≈ 3 cm/s annual modulation. The time-dilation term varies by a similar amount, ≈ 3 cm/s. Both are "constant at the 1 m/s level".

### 5.3 Transverse Doppler

[derived] v²/2c for representative speeds:

| motion | v | v²/2c |
|---|---|---|
| Earth rotation | 465 m/s | 0.36 mm/s |
| Earth orbit | 29.78 km/s | 1.48 m/s |
| star | 50 km/s | 4.2 m/s |
| star | 100 km/s | 16.7 m/s |
| halo star | 300 km/s | 150 m/s |
| Sun relative to the CMB | 369.82 km/s | 228 m/s |

### 5.4 Shapiro delay

* [derived] For a source at infinity and solar elongation ε, Δt = −(2GM☉/c³) ln(1 − cos ε). The induced RV is

$$\delta v=c\,\frac{d\Delta t}{dt}=-\frac{2GM_\odot}{c^2}\cot\frac{\varepsilon}{2}\,\dot\varepsilon,$$

  with 2GM☉/c² = 2953 m and |ε̇| ≤ 2π per sidereal year.
* [derived] Maxima: 0.59 mm/s at ε = 90°, 2.2 mm/s at 30°, 6.7 mm/s at 10°, 6.7 cm/s at 1°.
* [lit] Wright & Eastman 2014: "at most 1.5 mm/s" for τ Ceti, reached when the Sun is 23° away. Targets closer than about 32° to the Sun are not observable from the ground. For a star passing within 1′ of Jupiter the effect is ≪ 1 mm/s.

### 5.5 Aberration and direction errors

* Aberration has no separate effect if the relativistic Doppler formula is used with the BCRS (catalogue) direction.
* [derived] If the apparent (aberrated) direction n′ is inserted into the first-order formula, then n′·v − n·v = (v²/c) sin²θ. The maximum error is **v²/c = 2.96 m/s**.
* [derived] A position error δθ projects the orbital velocity: v_orb δθ = 0.144 m/s per arcsec.
* [lit] Wright & Eastman 2014: stellar positions must be better than 1″ for work at or below 10 cm/s. Dividing 1 cm/s by 30 km/s gives a required accuracy of 69 mas. Earth-orientation errors up to 4.4″ are acceptable at 1 cm/s, because 1 cm/s divided by 460 m/s is 4.4″. The flux-weighted exposure midpoint must be known to 250 ms.

### 5.6 Perspective, secular acceleration and other terms

* [lit] Wright & Eastman 2014:
  * Secular acceleration is ∼ r₀µ²Δt/c. For a star at 10 ly with µ = 2″/yr it is 0.3 m/s per year.
  * The parallax–proper-motion cross term is up to 10 cm/s for the nearest stars.
  * Relativistic addition of the systemic and planetary velocities is a ~1e-14 effect.
* [lit] Kürster et al. 2003 made the first detection of RV secular acceleration, for Barnard's star: 5.15 ± 0.89 m/s/yr measured against 4.50 m/s/yr predicted.

### 5.7 Time scales

* IAU 2000 C1 refers z_B to an observer "at zero gravitational potential … at rest with respect to the solar-system barycenter". This is the TCB frame.
* [derived] Expressing velocities in TDB-compatible rather than TCB units rescales them by 1 − L_B. That is 4.7e-4 m/s at 30 km/s, which is negligible. The *frequency* offset L_B c = 4.65 m/s is not negligible and must be part of the correction.

---

## 6. Redshift drift (Sandage–Loeb test)

For a comoving source (Sandage 1962; McVittie 1962; Loeb 1998; Liske et al. 2008, eq. 4):

$$\dot z\equiv\frac{dz}{dt_0}=(1+z)H_0-H(z),\qquad \dot v=\frac{c\,\dot z}{1+z}=cH_0\Big[1-\frac{E(z)}{1+z}\Big].$$

[derived] Planck18 in astropy (H₀ = 67.66 km s⁻¹ Mpc⁻¹, Ω_m = 0.30966, with radiation and neutrinos). The velocity-drift scale is cH₀ = 2.074 cm/s per year.

| z | E(z) | ż (yr⁻¹) | v̇ (cm/s/yr) | Δv in 10 yr (cm/s) | Δv in 20 yr (cm/s) | Δλ in 20 yr at Lyα(1+z) (Å) |
|---|---|---|---|---|---|---|
| 0.5 | 1.3188 | +1.254e-11 | +0.251 | +2.51 | +5.01 | +3.0e-7 |
| 1.0 | 1.7829 | +1.502e-11 | +0.225 | **+2.25** | +4.50 | +3.7e-7 |
| 1.5 | 2.3564 | +9.94e-12 | +0.119 | +1.19 | +2.38 | +2.4e-7 |
| 2.0 | 3.0157 | −1.08e-12 | −0.011 | **−0.11** | −0.22 | −2.6e-8 |
| 2.5 | 3.7467 | −1.707e-11 | −0.146 | −1.46 | −2.92 | −4.2e-7 |
| 3.0 | 4.5407 | −3.742e-11 | −0.280 | **−2.80** | −5.61 | −9.1e-7 |
| 3.5 | 5.3915 | −6.169e-11 | −0.411 | −4.11 | −8.22 | −1.5e-6 |
| 4.0 | 6.2947 | −8.959e-11 | −0.537 | **−5.37** | −10.74 | −2.2e-6 |
| 5.0 | 8.2452 | −1.554e-10 | −0.776 | −7.76 | −15.53 | −3.8e-6 |

* [derived] The drift is zero at **z = 1.960** and reaches its maximum positive value, **+2.58 cm/s per decade**, at z ≈ 0.64.
* [derived] For comparison, a flat Ω_m = 0.3, H₀ = 70 model gives +2.57, +0.24, −2.47 and −5.07 cm/s per decade at z = 1, 2, 3, 4. An Einstein–de Sitter model gives −8.6, −15.2, −20.7 and −25.6 cm/s per decade. Liske et al. (2008) call a redshift region with ż > 0 "the hallmark of Ω_Λ ≠ 0".
* [lit] Liske et al. 2008: "At z = 4 the redshift drift is of order 10⁻⁹ or 6 cm s⁻¹ per decade" (Ω_m, Ω_Λ = 0.3, 0.7). Martins et al. 2024 (ANDES): "Δz = żΔt ≈ H₀Δt ≈ 10⁻⁹ or ∼6 cm/s over a decade".
* [lit] Liske et al. 2008, eq. 16, photon-noise accuracy using all lines down to Lyγ:

  $$\sigma_v = 1.35\Big(\frac{S/N}{2370}\Big)^{-1}\Big(\frac{N_{\rm QSO}}{30}\Big)^{-1/2}\Big(\frac{1+z_{\rm QSO}}{5}\Big)^{-1.7}\ {\rm cm\,s^{-1}}.$$

  The exponent becomes −0.9 at z > 4. With Lyα only, the normalisation is 2 instead of 1.35 (eq. 15).
* [derived] Example: S/N = 3000, N = 30, z = 3 gives σ_v = 1.56 cm/s.
* [lit] Liske et al. 2008: a 42 m telescope can detect the drift "over a period of ~20 yr using 4000 h". Peculiar accelerations add only ~1e-3 cm/s per decade per line.
* [lit] Cooke 2020: Lyα-forest peculiar accelerations are more than an order of magnitude below the cosmological signal.
* [lit] Steinmetz et al. 2008: the experiment needs Doppler drifts of about 1 cm/s per year.
* [lit] Martins et al. 2024, ANDES requirements: total precision at the cm/s level. The residual bias of the wavelength solution "may not drift by more than ∼1 cm/s" over about two decades. The design goals are precision σ < 0.5 m/s, accuracy < 1 m/s and stability < 1 cm/s. The redshift-drift case is a *goal*, not a strict requirement.
* Wavelength scale: the drift is ≈ 1e-7 to 1e-6 Å over 20 years, i.e. 10⁻⁴ to 10⁻³ of the 1e-3 Å grid. It is only accessible differentially against a laser frequency comb.

---

## 7. Fine-structure constant α and proton-to-electron mass ratio µ

### 7.1 The many-multiplet method

The many-multiplet method was introduced by Webb et al. (1999) and Dzuba, Flambaum & Webb (1999a, PRL 82, 888; 1999b, PRA 59, 230):

$$\omega_z=\omega_0+q\,x,\qquad x=\Big(\frac{\alpha_z}{\alpha_0}\Big)^2-1\simeq 2\frac{\Delta\alpha}{\alpha}.$$

Murphy & Berengut (2014), eq. 1:

$$\frac{\Delta\alpha}{\alpha}\approx-\frac{\Delta v_i\,\omega_i}{2c\,q_i}\quad\Longleftrightarrow\quad \Delta v_i\approx-2c\,\frac{\Delta\alpha}{\alpha}\,\frac{q_i}{\omega_i}.$$

Errors in q cannot create a spurious non-zero Δα/α. Only velocity shifts *between* transitions can, which come from lab-wavelength errors or calibration distortions.

### 7.2 Selected transitions

[lit] values are verbatim from Murphy & Berengut 2014, Tables 1–11. The last two columns are [derived] for Δα/α = +1 ppm (rest frame).

| ion | λ0 (vac, Å) | ω0 (cm⁻¹) | δv_lab (m/s) | q (cm⁻¹) | Δv (m/s) | Δλ (Å) |
|---|---|---|---|---|---|---|
| Mg II 2796 | 2796.353790(16) | 35760.85414(20) | 1.7 | 212(2) | −3.56 | −3.3e-5 |
| Mg II 2803 | 2803.530982(16) | 35669.30440(20) | 1.7 | 121(2) | −2.03 | −1.9e-5 |
| Mg I 2852 | 2852.962797(15) | 35051.28076(19) | 1.6 | 90(10) | −1.54 | −1.5e-5 |
| Al II 1670 | 1670.78861(11) | 59851.976(4) | 20.0 | 270(10) | −2.71 | −1.5e-5 |
| Si II 1526 | 1526.706980(16) | 65500.4538(7) | 3.2 | 47(4) | −0.43 | −2.2e-6 |
| Fe II 1608 | 1608.450852(78) | 62171.623(3) | 14.5 | −1165(300) | +11.24 | +6.0e-5 |
| Fe II 2344 | 2344.212747(76) | 42658.2443(14) | 9.7 | 1375(300) | −19.33 | −1.5e-4 |
| Fe II 2382 | 2382.763995(80) | 41968.0674(14) | 10.1 | 1505(100) | −21.50 | −1.7e-4 |
| Fe II 2600 | 2600.172114(88) | 38458.9926(13) | 10.1 | 1370(100) | −21.36 | −1.9e-4 |
| Cr II 2066 | 2066.163819(60) | 48398.8729(14) | 8.8 | −1421(70) | +17.60 | +1.2e-4 |
| Zn II 2026 | 2026.136946(31) | 49355.00544(76) | 4.6 | 2470(25) | −30.01 | −2.0e-4 |
| Ni II 1741 | 1741.55289(12) | 57420.0191(40) | 20.9 | −1400(250) | +14.62 | +8.5e-5 |

* The full set of 43 transitions is in `data/notes/murphy_berengut_2014_mm_transitions.json`.
* [derived] Zn II 2026 against Cr II 2066 shifts by 47.6 m/s per ppm, the most sensitive pair in this table.
* At observed redshift z, the wavelength shifts scale by (1+z).

### 7.3 Laboratory precision requirement

* [lit] Murphy & Berengut 2014: a transition is "useful for measuring Δα/α" only if its rest-frame frequency has "been measured with velocity precision δv ≲ 20 m s⁻¹". The tabulated uncertainties range "from just 0.04 m s⁻¹ for the frequency comb-calibrated Mg II isotopes … up to ≈20 m s⁻¹ for the Ni II transitions".
* [derived] 20 m/s corresponds to **1.1e-4 Å at 1608 Å, 1.35e-4 Å at 2026 Å and 1.9e-4 Å at 2796 Å**. The required lab precision is therefore ≈ 1e-4 Å, about 10 times finer than the wiki's 1e-3 Å grid. Isotopic and hyperfine structure, up to ~0.85 km/s for Mg II (Murphy & Berengut 2014), must also be included.

### 7.4 Measurements and calibration systematics

* [lit] Historical results:
  * Webb et al. 1999: Δα/α = (−1.1 ± 0.4)e-5 for 30 Keck absorbers.
  * Murphy, Webb & Flambaum 2003: (−0.543 ± 0.116)e-5 for 128 absorbers.
  * Webb et al. 2011 and King et al. 2012: a spatial dipole toward RA 17.3 h, Dec −61°, with amplitude 0.97 (+0.22/−0.20) e-5 at 4.1σ.
  * Murphy et al. 2022 (ESPRESSO with LFC calibration, HE 0515−4414, z = 1.15): **Δα/α = 1.3 ± 1.3 (stat) ± 0.4 (sys) ppm**.
* [lit] Calibration systematics:
  * Murphy et al. 2007 (MNRAS 378, 221): UVES spectra calibrated with the older Th-Ar list contain ±30–75 m/s distortions.
  * Schmidt et al. 2021: ESPRESSO Th-Ar+FP against LFC calibration differs by up to 24 m/s, depending on wavelength.
  * Martins et al. 2024: ANDES aims at about 0.3 ppm, corresponding to ∼6 m/s relative line shifts. Relative line positions must then be accurate to ∼1 m/s, whereas "ESPRESSO's accuracy in this context is ∼20 m s⁻¹".

### 7.5 Proton-to-electron mass ratio µ

* The relation is λ_i/λ_i⁰ = (1+z_Q)(1 + K_i Δµ/µ), so Δv_i ≈ c K_i Δµ/µ (Ubachs et al. 2007, eq. 3).
* [lit] Ubachs et al. 2007, Table 8, lists H₂ Lyman/Werner K_i values including K = −0.01872 for L0 P(4) and +0.04417 for L11 R(0). This excerpt shows K_i of order −0.02 to +0.05; the full extremes of the table were not extracted.
* [derived] For Δµ/µ = 1 ppm these two lines shift by −5.6 and +13.2 m/s.
* [lit] Radio methanol in PKS 1830−211 (z = 0.89): Δµ/µ = (−1.0 ± 0.8 stat ± 1.0 sys)e-7 (Bagdonaite et al. 2013, PRL 111, 231101; see also Science 339, 46). Kanekar et al. (2015) give a VLA methanol constraint for the same lens.
* [lit] ANDES: ∼0.5 ppm on µ needs ∼5 m/s relative-velocity systematics across the U band (Martins et al. 2024). The review is Ubachs et al. (2016, RMP 88, 021003).

---

## 8. Cosmological redshift versus peculiar velocity; H₀ error budget

### 8.1 Recession velocity is not a Doppler velocity

* The cosmological redshift is 1+z = a(t₀)/a(t_e). It is not a Doppler shift.
* The proper recession velocity today is v_rec = H₀ D_C(z), which exceeds c beyond the Hubble sphere (Davis & Lineweaver 2004).
* [derived] Planck18:

  | z | v_rec/c | SR "Doppler" v/c | cz/c |
  |---|---|---|---|
  | 0.1 | 0.098 | 0.095 | 0.1 |
  | 0.5 | 0.439 | 0.385 | 0.5 |
  | 1 | 0.766 | 0.600 | 1 |
  | 2 | 1.198 | 0.800 | 2 |
  | 3 | 1.468 | 0.882 | 3 |
  | 1089.8 | 3.13 | 1.000 | 1089.8 |

* [derived] v_rec = c at **z = 1.476**.

### 8.2 Low-z Hubble diagram

* [lit] Riess et al. 2022, eqs. 4 and 17, use a cosmographic expansion:

  $$a_B=\log\Big\{cz\big[1+\tfrac12(1-q_0)z-\tfrac16(1-q_0-3q_0^2+j_0)z^2\big]\Big\}-0.2m_B^0$$

  with q₀ = −0.55, j₀ = 1 and Hubble-flow SNe at **0.0233 < z < 0.15**.
* [lit] Riess et al. 2022: dropping the flow corrections lowers H₀ by 0.5 km s⁻¹ Mpc⁻¹. An H₀ from Cepheid-host redshifts alone (mean cz ≈ 2000 km/s) would be limited to 3–4% by local flows.
* [derived] H₀D_L/(cz) − 1 = +1.8% at z = 0.0233 and +11% at z = 0.15 (Planck18). The (1+z) and q₀ terms are therefore mandatory at these redshifts.

### 8.3 Peculiar-velocity noise

* [derived] Using Davis et al. 2011 eqs. A1 and A4:

  | σ_v | z | σ_µ | distance error per object |
  |---|---|---|---|
  | 250 km/s | 0.0233 | 0.080 mag | 3.7% |
  | 250 km/s | 0.1 | 0.021 mag | 0.96% |
  | 150 km/s | 0.0233 | 0.048 mag | 2.2% |

  This is why the H₀ ladder starts at z ≈ 0.023.
* Systematic redshift errors (frame, dipole, Δz monopole) matter at the 1e-5 level (Calcino & Davis 2017). That is about 3 km/s, or ≈ 0.05 Å at 5000 Å, far coarser than 1e-3 Å. Redshift-measurement precision is therefore not limited by the wiki's rest-wavelength precision, but the conversion formulae must be exact.

---

## 9. "Radial velocity" is not a kinematic velocity: IAU 2000 Resolutions C1 and C2

The following text is verbatim from the Appendix of Lindegren & Dravins (2003), which reproduces the resolutions (Rickman 2001, Trans. IAU XXIVB).

**Resolution C1**, on the Definition of a Spectroscopic "Barycentric Radial-Velocity Measure":

> Recognising … 2. That, due to relativistic effects, measurements being made inside gravitational fields, and alternative choices of coordinate frames, the naive concept of radial velocity being equal to the time derivative of distance, becomes ambiguous at accuracy levels around 100 m/s;
>
> Considering 1. That, although many effects may influence the precise shifts of spectroscopic wavelengths and frequencies, only local ones (i.e. arising within the solar system, and depending on the gravitational potential of the observer, and the observer's position and motion relative to the solar-system barycenter) can in general be reliably calculated; 2. That, although the wavelength displacement (or frequency shift) corrected for such local effects can thus be derived from spectroscopic measurements, the resulting quantity cannot unambiguously be interpreted as a radial motion of the object;
>
> Therefore recommend That, whenever radial velocities are considered to a high accuracy, the spectroscopic result from a measurement of shifts in wavelength or frequency be given as the "barycentric radial-velocity measure" cz_B, after correcting for gravitational effects caused by solar-system objects, and effects by the observer's displacement and motion relative to the solar-system barycenter. Here, c equals the conventional speed of light = 299,792,458 m/s, and z_B = (λ − λ0)/λ0, where λ0 is the rest-frame wavelength and λ the wavelength observed by a hypothetical observer at zero gravitational potential, located at, and being at rest with respect to, the solar-system barycenter. The epoch of the observation equals the barycentric time of light arrival.

**Resolution C2** defines the geometric (astrometric) radial velocity as dr/dt_B, with r the barycentric coordinate distance and t_B the TCB of light arrival at the barycentre.

**Interpretation** (Lindegren & Dravins 2003):

* cz_B is "just a quantification of the spectroscopic shift, not of physical velocity". It contains the source's gravitational redshift, the transverse-Doppler term |v*|²/2c and the factor (1+X): convection, line asymmetries, pressure shifts, rotation and activity, plus the cosmological redshift.
* For a solar-type star the sum is of order +636 − 400 m/s ≈ +0.2 km/s, and it differs line by line by ~200 m/s.
* RV standard stars are therefore unlikely to be "true" to better than 0.5 (perhaps 0.3) km/s without a physical model of the star. Kinematic v_r follows from cz_B only if Φ*, v_t and X are known.
* Note 6 of the paper: the "B" refers to the solar-system barycentre, "not the rest-frame defined by the cosmological microwave background".

**Consequences for the wiki:**

* Laboratory rest wavelengths are λ0.
* Solar or stellar line positions are *apparent* wavelengths that include +633 m/s (solar gravitational redshift as seen from Earth) and the line-dependent convective shift.
* Telluric lines are in the topocentric frame, with pressure shifts of their own.
* Each list must carry its frame (SPECSYS), medium and convention.

---

## 10. Instrumental wavelength calibrators and their precision

| calibrator | range | precision / accuracy | source |
|---|---|---|---|
| Th-Ar hollow-cathode lamp, HARPS list | 3785–6915 Å | per-line internal uncertainty ≈ 10 m/s (3.3e-8, 0.18 mÅ); global HARPS solution precise to 20 cm/s (0.0037 mÅ); absolute scale tied to Palmer & Engleman 1983 | [lit] Lovis & Pepe 2007 |
| Th Ritz wavelengths (NIST FTS) | 250–5500 nm, 19874 lines | from globally optimised Th I/II/III levels; σ(PE83 − Ritz) = 8.2e-5 nm, σ(LP07 − Ritz) = 13.2e-5 nm; "large scatter" in LP07; best individual Ritz values ~2e-5 nm (e.g. 873.482334 ± 0.000017 nm) | [lit] Redman, Nave & Sansonetti 2014 |
| Th-Ar line selection for UVES | | old list caused ±30–75 m/s distortions; selection algorithm improves residuals by more than 3× | [lit] Murphy et al. 2007 |
| Th-Ar in the near-IR | 691–5804 nm | Kerber et al. 2008 atlas "is sufficient for precisions of ≈ 50 m/s on CRIRES"; telluric CO₂ lines have reached 6–10 m/s (Figueira et al. 2010) | [lit] as quoted by Redman et al. 2011 |
| U lamps and U-Ne | 850–4000 nm; H band 1454–1638 nm | 10,100 U I/II lines, "six times as many calibration lines as thorium in the NIR"; U/Ne residuals 50% larger than a modal-noise-limited LFC | [lit] Redman et al. 2011, 2012 |
| I₂ absorption cell | visible (the cell's absorption band; range not extracted here) | 3 m/s Doppler precision | [lit] Marcy & Butler 1992; Butler et al. 1996 |
| Fabry–Pérot etalon combined with a hollow-cathode lamp | | combined solution removes 50 m/s distortions; differs from LFC by ≤ 10 m/s; chromatic drifts must be modelled (HPF) | [lit] Wildi et al. 2010; Bauer et al. 2015; Cersullo et al. 2019; Terrien et al. 2021 |
| laser frequency comb | visible and near-IR | mode frequencies f_n = n f_rep + f_0, known to better than 1e-12 relative; simulated photon limit ≈ 1 cm/s over 4000 Å; see list below | [lit] Murphy et al. 2007 (MNRAS 380, 839); Steinmetz et al. 2008; Li et al. 2008; Wilken et al. 2012; Probst et al. 2020; Milaković et al. 2020 |
| telluric O₂ lines | | stable to 10 m/s rms over 6 years (HARPS); 5 m/s "out of the box"; 1–2 m/s with simple atmospheric modelling | [lit] Figueira et al. 2010 |
| OH (and O₂) airglow lines | 0.6–2.6 µm | wavelengths computed from Abrams et al. 1994 OH energy levels, "secure to about one hundredth of Å" (≈ 190 m/s at 1.6 µm); O₂ lines "a few hundredths of Å"; vacuum, with an air table given | [lit] Oliva & Origlia 1992; Osterbrock et al. 1996, 1997; Rousselot et al. 2000 |
| solar atlases as references | 405–1065 nm | IAG FTS atlas accurate to ±10 m/s; Kitt Peak atlases deviate by tens to 100 m/s; LARS (LFC-calibrated) solar spectra accurate to 1 m/s | [lit] Reiners et al. 2016; Löhner-Böttcher et al. 2018, 2019 |

Laser frequency comb results in detail:

* Steinmetz et al. 2008: first astronomical use, ≈ 9 m/s at 1.5 µm.
* Li et al. 2008: an "astro-comb" should allow 1 cm/s.
* Wilken et al. 2012: short-term repeatability of 2.5 cm/s on HARPS.
* Probst et al. 2020: two independent LFCs confirm 1 cm/s stability.
* Milaković et al. 2020: ≲ 4 cm/s per exposure and 1 cm/s averaged over a few hours, but ≈ 60 cm/s zero-point offsets between combs in the raw data.
* Martins et al. 2024 (ANDES): ≈ 50 cm/s inconsistency between two HARPS LFCs, and LFC flux-dependent shifts up to ±2 m/s on ESPRESSO.

**Hierarchy (precision per line, then per exposure):**

| calibrator | precision |
|---|---|
| OH airglow | ~100 m/s |
| Th-Ar or U-Ne line lists | ~10 m/s per line, 0.1–1 m/s per exposure |
| telluric O₂ | ~10 m/s |
| I₂ cell | ~1–3 m/s |
| FP etalon | relative precision at the cm/s level, but needs an absolute anchor |
| LFC | ~1–3 cm/s, absolute |

Only the LFC is traceable to an atomic clock, which the redshift-drift and α experiments need. For lamp-based work, accuracy is limited by the laboratory wavelengths, i.e. by the same 1e-4 to 1e-3 Å precision this wiki tabulates. Martins et al. (2024) note that HCL laboratory data "currently seem not accurate enough to achieve the 1 m/s accuracy goal".

---

## 11. Recommended reduction chain for the wiki's calibration section

1. **Measure** λ_obs in **vacuum**, or convert from air with a stated formula.
2. **Remove the observer terms** with the BCRS direction **k** of the source (not the apparent direction):
   $$1+z_B=(1+z_{\rm obs})\,\gamma_{\rm obs}\big(1+\mathbf k\cdot\boldsymbol\beta_{\rm obs}\big)\big/\big(1+z_{\rm GR,obs}\big),\qquad z_{\rm GR,obs}\simeq-\textstyle\sum_i GM_i/(c^2 d_i).$$
   The observer terms are Earth rotation, orbital motion, lunar reflex, the Sun about the SSB, and the 4.65 m/s mean observer-clock offset. Tools are barycorrpy (Kanodia & Wright 2018) and the Wright & Eastman (2014) algorithm, which reach 1 cm/s. The required inputs are timing to 250 ms, position to 69 mas per cm/s, and proper motion and parallax.
3. **Change frame** to LSR, GSR, LG or CMB with a multiplicative factor. Here β_f is the velocity of the barycentre (or Sun) relative to frame f, and n̂ is the unit vector toward the source:
   $$1+z_{\rm frame}=(1+z_B)\,\gamma_f\big(1+\boldsymbol\beta_f\cdot\hat{\mathbf n}\big).$$
   Motion toward the source blueshifts it, so z_frame > z_B. In line-of-sight form this is $z_{\rm frame}=(1+z_B)/(1+z_{\odot,f})-1$ with $z_{\odot,f}=\sqrt{(1-\beta_f\cos\alpha)/(1+\beta_f\cos\alpha)}-1$ (Carr et al. 2022, eqs. 3 and 5, for the CMB). The γ_f factor differs from that form only at order β². Record which apex is used: FITS conventional or modern.
4. **Separate the source terms** with $(1+z_{\rm frame})=(1+\bar z)(1+z_{\rm pec})(1+z_{\rm grav,*})(1+z_{\rm conv})$. Do not interpret c z̄ as a velocity in the SR sense.
5. **Publish z or cz_B**, with SPECSYS, SSYSOBS and the velocity convention (VOPT, VRAD or VELO), following Greisen et al. (2006).

---

## 12. References

All DOIs were checked against arXiv metadata or Crossref. Journal, volume and page for 47 DOIs were re-checked with `scripts/research/crossref_doi_check.py`; the results are in `data/raw/relativity_cosmology/crossref_doi_check.json`.

* Abrams, M. C., Davis, S. P., Rao, M. L. P., Engleman, R., Jr., Brault, J. W. 1994, ApJS 93, 351. doi:10.1086/192058. High-resolution FTS of the OH Meinel system, used for the OH line positions.
* Allende Prieto, C. & García López, R. J. 1998, A&AS 129, 41. doi:10.1051/aas:1998173. Fe I line shifts in the solar spectrum.
* Bagdonaite, J., Daprà, M., Jansen, P., Bethlem, H. L., Ubachs, W., et al. 2013, PRL 111, 231101. doi:10.1103/PhysRevLett.111.231101. arXiv:1311.3438.
* Bagdonaite, J., Jansen, P., Henkel, C., Bethlem, H. L., Menten, K. M., Ubachs, W. 2013, Science 339, 46. doi:10.1126/science.1224898.
* Barstow, M. A., Bond, H. E., Holberg, J. B., et al. 2005, MNRAS 362, 1134. doi:10.1111/j.1365-2966.2005.09359.x.
* Bauer, F. F., Zechmeister, M., Reiners, A. 2015, A&A 581, A117. doi:10.1051/0004-6361/201526462. arXiv:1506.07887.
* Bennett, C. L., et al. 2003, ApJS 148, 1, as cited by Greisen et al. 2006 for CMBDIPOL. Not independently checked.
* Butler, R. P., Marcy, G. W., Williams, E., McCarthy, C., Dosanjh, P., Vogt, S. S. 1996, PASP 108, 500. doi:10.1086/133755.
* Calcino, J. & Davis, T. 2017, JCAP 01, 038. doi:10.1088/1475-7516/2017/01/038. arXiv:1610.07695.
* Carr, A., Davis, T. M., Scolnic, D., Said, K., et al. 2022, PASA 39, e046. doi:10.1017/pasa.2022.41. arXiv:2112.01471.
* Carr, A. & Davis, T. M. 2021, ApJ 914, 97. doi:10.3847/1538-4357/abfb6e. Corrections to the NED redshift calculator.
* Cersullo, F., Coffinet, A., Chazelas, B., Lovis, C., Pepe, F. 2019, A&A 624, A122. doi:10.1051/0004-6361/201833852. arXiv:1903.01942.
* Cooke, R. 2020, MNRAS 492, 2044. doi:10.1093/mnras/stz3465. arXiv:1912.04983.
* Cristiani, S., et al. 2007, Nuovo Cim. B 122, 1159. arXiv:0712.4152. CODEX-ESPRESSO. The DOI given by arXiv, 10.1393/ncb/i2008-10450-9, returns 404 on Crossref.
* Davis, T. M. & Lineweaver, C. H. 2004, PASA 21, 97. doi:10.1071/AS03040. arXiv:astro-ph/0310808.
* Davis, T. M., Hui, L., Frieman, J. A., Haugbølle, T., et al. 2011, ApJ 741, 67. doi:10.1088/0004-637X/741/1/67. arXiv:1012.2912.
* Davis, T. M. & Scrimgeour, M. I. 2014, MNRAS 442, 1117. doi:10.1093/mnras/stu920. arXiv:1405.0105.
* Delhaye, J. 1965, in Stars and Stellar Systems V, ed. Blaauw & Schmidt (Univ. Chicago Press), 61, as cited by Greisen et al. 2006.
* Dravins, D. 1982, ARA&A 20, 61. doi:10.1146/annurev.aa.20.090182.000425.
* Dravins, D., Lindegren, L., Nordlund, Å. 1981, A&A 96, 345. Cited from the literature; no DOI was found via Crossref.
* Dravins, D., Lindegren, L., Madsen, S. 1999, A&A 348, 1040. arXiv:astro-ph/9907145. Astrometric radial velocities I.
* Dzuba, V. A., Flambaum, V. V., Webb, J. K. 1999, PRL 82, 888. doi:10.1103/PhysRevLett.82.888. Also 1999, PRA 59, 230. doi:10.1103/PhysRevA.59.230.
* Falcon, R. E., Winget, D. E., Montgomery, M. H., Williams, K. A. 2010, ApJ 712, 585. doi:10.1088/0004-637X/712/1/585. arXiv:1002.2009.
* Figueira, P., Pepe, F., Lovis, C., Mayor, M. 2010, A&A 515, A106. doi:10.1051/0004-6361/201014005. arXiv:1003.0541.
* Fixsen, D. J., et al. 1996, ApJ 473, 576. doi:10.1086/178173. arXiv:astro-ph/9605054.
* González Hernández, J. I., Rebolo, R., Pasquini, L., Lo Curto, G., et al. 2020, A&A 643, A146. doi:10.1051/0004-6361/202038937. arXiv:2009.10558.
* Gordon, M. A. 1976, in Methods of Experimental Physics 12C (Astrophysics: Radio Observations), ed. M. L. Meeks (Academic Press, New York), p. 277, as cited by Greisen et al. 2006.
* GRAVITY Collaboration (Abuter, R., et al.) 2018, A&A 615, L15. doi:10.1051/0004-6361/201833718. arXiv:1807.09409.
* Greenstein, J. L. & Trimble, V. L. 1967, ApJ 149, 283. doi:10.1086/149254.
* Greisen, E. W., Calabretta, M. R., Valdes, F. G., Allen, S. L. 2006, A&A 446, 747. doi:10.1051/0004-6361:20053818. arXiv:astro-ph/0507293.
* Joyce, S. R. G., Barstow, M. A., Holberg, J. B., Bond, H. E., et al. 2018, MNRAS 481, 2361. doi:10.1093/mnras/sty2404. arXiv:1809.01240.
* Kanekar, N., Ubachs, W., Menten, K. M., Bagdonaite, J., Brunthaler, A., et al. 2015, MNRAS 448, L104. doi:10.1093/mnrasl/slu206. arXiv:1412.7757.
* Kanodia, S. & Wright, J. 2018, RNAAS 2, 4. doi:10.3847/2515-5172/aaa4b7. arXiv:1801.01634. barycorrpy.
* Kerber, F., Nave, G., Sansonetti, C. J. 2008, ApJS 178, 374. doi:10.1086/590111.
* Kerr, F. J. & Lynden-Bell, D. 1986, MNRAS 221, 1023. doi:10.1093/mnras/221.4.1023.
* King, J. A., Webb, J. K., Murphy, M. T., Flambaum, V. V., Carswell, R. F., et al. 2012, MNRAS 422, 3370. doi:10.1111/j.1365-2966.2012.20852.x. arXiv:1202.4758.
* Kürster, M., Endl, M., Rouesnel, F., Els, S., et al. 2003, A&A 403, 1077. doi:10.1051/0004-6361:20030396. arXiv:astro-ph/0303528.
* Li, C.-H., Benedick, A. J., Fendel, P., Glenday, A. G., Kärtner, F. X., Phillips, D. F., Sasselov, D., Szentgyorgyi, A., Walsworth, R. L. 2008, Nature 452, 610. doi:10.1038/nature06854. arXiv:0804.0955.
* Lindegren, L. & Dravins, D. 2003, A&A 401, 1185. doi:10.1051/0004-6361:20030181. arXiv:astro-ph/0302522. Includes the full text of IAU 2000 Resolutions C1 and C2 (Rickman 2001, Trans. IAU XXIVB).
* Liske, J., Grazian, A., Vanzella, E., Dessauges, M., et al. 2008, MNRAS 386, 1192. doi:10.1111/j.1365-2966.2008.13090.x. arXiv:0802.1532.
* Loeb, A. 1998, ApJ 499, L111. doi:10.1086/311375. arXiv:astro-ph/9802122.
* Löhner-Böttcher, J., Schmidt, W., Stief, F., Steinmetz, T., Holzwarth, R. 2018, A&A 611, A4. doi:10.1051/0004-6361/201732107. Also 2019, A&A 624, A57. doi:10.1051/0004-6361/201834925.
* Lovis, C. & Pepe, F. 2007, A&A 468, 1115. doi:10.1051/0004-6361:20077249. arXiv:astro-ph/0703412.
* Marcy, G. W. & Butler, R. P. 1992, PASP 104, 270. doi:10.1086/132989.
* Martins, C. J. A. P., Cooke, R., Liske, J., Murphy, M. T., et al. 2024, "Cosmology and fundamental physics with the ELT-ANDES spectrograph", arXiv:2311.16274v2. The journal reference was not checked.
* Milaković, D., Pasquini, L., Webb, J. K., Lo Curto, G. 2020, MNRAS 493, 3997. doi:10.1093/mnras/staa356. arXiv:2002.05182.
* Murphy, M. T., Webb, J. K., Flambaum, V. V. 2003, MNRAS 345, 609. doi:10.1046/j.1365-8711.2003.06970.x.
* Murphy, M. T., Tzanavaris, P., Webb, J. K., Lovis, C. 2007, MNRAS 378, 221. doi:10.1111/j.1365-2966.2007.11768.x.
* Murphy, M. T., Udem, Th., Holzwarth, R., Sizmann, A., Pasquini, L., Araujo-Hauck, C., et al. 2007, MNRAS 380, 839. doi:10.1111/j.1365-2966.2007.12147.x.
* Murphy, M. T. & Berengut, J. C. 2014, MNRAS 438, 388. doi:10.1093/mnras/stt2204. arXiv:1311.2949. Not in VizieR; checked on the CfA mirror. The table was parsed from the arXiv v3 PDF.
* Murphy, M. T., Molaro, P., Leite, A. C. O., Cupani, G., Cristiani, S., et al. 2022, A&A 658, A123. doi:10.1051/0004-6361/202142257. arXiv:2112.05819.
* Oliva, E. & Origlia, L. 1992, A&A 254, 466. ADS 1992A&A...254..466O.
* Osterbrock, D. E., Fulbright, J. P., Martel, A. R., Keane, M. J., Trager, S. C., Basri, G. 1996, PASP 108, 277. doi:10.1086/133722. Also Osterbrock, Fulbright & Bida 1997, PASP 109, 614. doi:10.1086/133920.
* Pasquini, L., Melo, C., Chavero, C., Dravins, D., et al. 2011, A&A 526, A127. doi:10.1051/0004-6361/201015337. arXiv:1011.4635.
* Planck Collaboration 2020, A&A 641, A1 (Planck 2018 I). doi:10.1051/0004-6361/201833880. Also A&A 641, A6 (Planck 2018 VI). doi:10.1051/0004-6361/201833910.
* Pound, R. V. & Rebka, G. A. 1960, PRL 4, 337. doi:10.1103/PhysRevLett.4.337.
* Prša, A., Harmanec, P., Torres, G., Mamajek, E., Asplund, M., et al. 2016, AJ 152, 41. doi:10.3847/0004-6256/152/2/41. IAU 2015 Res. B3 nominal values.
* Probst, R. A., Milaković, D., Toledo-Padrón, B., Lo Curto, G., et al. 2020, Nature Astronomy 4, 603. doi:10.1038/s41550-020-1010-x. arXiv:2002.08868.
* Redman, S. L., Lawler, J. E., Nave, G., Ramsey, L. W., Mahadevan, S. 2011, ApJS 195, 24. doi:10.1088/0067-0049/195/2/24. arXiv:1107.4091.
* Redman, S. L., Ycas, G. G., Terrien, R., Mahadevan, S., Ramsey, L. W., Bender, C. F., et al. 2012, ApJS 199, 2. doi:10.1088/0067-0049/199/1/2. arXiv:1112.1062.
* Redman, S. L., Nave, G., Sansonetti, C. J. 2014, ApJS 211, 4. doi:10.1088/0067-0049/211/1/4. arXiv:1308.5229.
* Reiners, A., Mrotzek, N., Lemke, U., Hinrichs, J., Reinsch, K. 2016, A&A 587, A65. doi:10.1051/0004-6361/201527530. arXiv:1511.03014.
* Riess, A. G., Yuan, W., Macri, L. M., Scolnic, D., et al. 2022, ApJL 934, L7. doi:10.3847/2041-8213/ac5c5b. arXiv:2112.04510.
* Roca Cortés, T. & Pallé, P. L. 2014, MNRAS 443, 1837. doi:10.1093/mnras/stu1238.
* Rousselot, P., Lidman, C., Cuby, J.-G., Moreels, G., Monnet, G. 2000, A&A 354, 1134–1150. ADS 2000A&A...354.1134R.
* Sandage, A. 1962, ApJ 136, 319. doi:10.1086/147385.
* Schmidt, T. M., Molaro, P., Murphy, M. T., Lovis, C., Cupani, G., Cristiani, S., et al. 2021, A&A 646, A144. doi:10.1051/0004-6361/202039345. arXiv:2011.13963.
* Schönrich, R., Binney, J., Dehnen, W. 2010, MNRAS 403, 1829. doi:10.1111/j.1365-2966.2010.16253.x. arXiv:0912.3693.
* Steinmetz, T., Wilken, T., Araujo-Hauck, C., Holzwarth, R., Hänsch, T. W., Pasquini, L., et al. 2008, Science 321, 1335. doi:10.1126/science.1161030.
* Stumpff, P. 1980, A&AS 41, 1, as cited by Greisen et al. 2006.
* Terrien, R. C., Ninan, J. P., Diddams, S. A., Mahadevan, S., et al. 2021, AJ 161, 252. doi:10.3847/1538-3881/abef68.
* Ubachs, W., Buning, R., Eikema, K. S. E., Reinhold, E. 2007, J. Mol. Spectrosc. 241, 155. doi:10.1016/j.jms.2006.12.004.
* Ubachs, W., Bagdonaite, J., Salumbides, E. J., Murphy, M. T., Kaper, L. 2016, Rev. Mod. Phys. 88, 021003. doi:10.1103/RevModPhys.88.021003.
* Webb, J. K., Flambaum, V. V., Churchill, C. W., Drinkwater, M. J., Barrow, J. D. 1999, PRL 82, 884. doi:10.1103/PhysRevLett.82.884. arXiv:astro-ph/9803165.
* Webb, J. K., King, J. A., Murphy, M. T., Flambaum, V. V., Carswell, R. F., Bainbridge, M. B. 2011, PRL 107, 191101. doi:10.1103/PhysRevLett.107.191101.
* Wildi, F., Pepe, F., Chazelas, B., Lo Curto, G., Lovis, C. 2010, Proc. SPIE 7735, 77354X. doi:10.1117/12.857951.
* Wilken, T., Lo Curto, G., Probst, R. A., Steinmetz, T., Manescau, A., Pasquini, L., González Hernández, J. I., Rebolo, R., Hänsch, T. W., Udem, T., Holzwarth, R. 2012, Nature 485, 611. doi:10.1038/nature11092.
* Wojtak, R., Hansen, S. H., Hjorth, J. 2011, Nature 477, 567. doi:10.1038/nature10445. arXiv:1109.6571.
* Wright, J. T. & Eastman, J. D. 2014, PASP 126, 838. doi:10.1086/678541. arXiv:1409.4774.
* Time-scale constants: IAU 2000 Res. B1.9 (L_G = 6.969290134e-10) and IAU 2006 Res. B3 (L_B = 1.550519768e-8), both from ERFA (`erfa.ELG`, `erfa.ELB`). The IAU web site (www.iau.org) could not be reached from this machine, so the resolution PDFs were not downloaded.
