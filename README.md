# Spectral Line Wiki (1 mÅ edition)

A self-contained, offline HTML reference for precise spectroscopy:

* **Element emission lines** — 11 968 strong lines of the neutral and singly-ionised atoms of H–Es (Z = 1–99) from the
  NIST *Handbook of Basic Atomic Spectroscopic Data* (Sansonetti & Martin 2005; local e-book copy in
  `~/Desktop/ebookbasd`; 12 008 printed, 12 duplicated Pt lines dropped, 28 lines printed without an ionisation stage
  kept aside), cross-matched line-by-line with the NIST Atomic Spectra Database (ASD v5.12, queried 2026-11-23).
* **Earth-atmosphere emission (airglow / aurora / geocorona / meteors)** — 214 atomic lines (NIST ASD), 171 molecular
  band heads/origins, the PALACE v1.0 airglow model line list (Noll et al. 2025; 3000–10 500 Å plus all 183 O I lines to
  2.13 µm; laboratory positions), the UVES/VLT observed night-sky atlas (Hanuschik 2003; identifications re-derived with
  the Cosby et al. 2006 criterion) and 15 802 near-IR OH Meinel / O₂ lines (0.74–4.2 µm; HITRAN2024 positions, PALACE
  intensities to 2.5 µm, including 2336 OH Δv = 1 fundamental and Δv = 2 lines at 2.5–4.2 µm without intensities).
* **Fraunhofer absorption lines** — 62 classic letter designations (lab values from NIST ASD / HITRAN, solar values from
  Moore, Minnaert & Houtgast 1966, Allende Prieto & García López 1998, the IAG atlas), 8455 solar photospheric lines
  (3374–12 818 Å; none in the H and K bands) and
  23 976 telluric lines (HITRAN2024, depth ≥ 1 % at R = 10⁵ from sea level, with pressure shifts).
* **Why 1 mÅ** — velocity equivalents, spectrograph centroid precision (31 instrument modes), a 45-row table of physical
  effects at 5000 Å, and a verdict on where 1 mÅ is necessary, insufficient or overkill.
* **Wavelength calibration** — air ↔ vacuum (seven standard formulas + the full Ciddor 1996 model), Doppler
  conventions and FITS frames, Earth rotation / orbit / barycentric correction (in-browser JPL DE440s calculator,
  Wright & Eastman 2014), time scales and BJD, LSR and Galactic rotation, Local Group, CMB dipole and large-scale flows,
  cosmological redshift composition, redshift drift, α-variation, the Galactic acceleration of the barycentre,
  instrumental calibrators and a line-shift checklist; a forward model gives the expected observed wavelength of a line
  for a given target, site, time, source velocity/frame and spectrograph medium, and the explorer can Doppler-shift the
  laboratory lines to a source velocity.

**Live page:** <https://exoplanet5.github.io/spectral-line-wiki/> — or open `index.html` locally in any modern browser
(Chrome, Safari, Firefox); it needs no network access.

Raw downloads (`data/raw/`: NIST ASD query dumps, HITRAN line files, PALACE archives, paper PDFs and scans) are not
included in this repository; the fetch scripts under `scripts/` and `scripts/research/` recreate them. Everything the
page is built from — curated catalogues (`data/curated/`), normalised data blocks (`data/build/`), parsed Handbook and
ASD match tables, research notes and the ephemeris fit — is included.

## Conventions

* Wavelengths in Å. **Standard air** for λ ≥ 2000 Å, **vacuum** below (IAU/NIST convention); the top-bar toggle switches
  the explorer, finder and sorting to vacuum. Tables always show both.
* Values a source publishes are shown verbatim. Converted values show the digits the source supports; digits beyond
  that are faded. Air ↔ vacuum conversions use Peck & Reeder (1972) eq. 3, the formula NIST ASD itself uses (identified
  empirically: 99.99 % of 35 988 ASD air/vacuum pairs reproduced within rounding).
* σ<sub>λ</sub> is the uncertainty stated by the source (or the bound implied by a source code, e.g. HITRAN position codes),
  otherwise an implied 2.5 units of the last published digit (NIST ASD: implied uncertainties are "generally between 2.5 and
  25 units of the least significant figure"); photographic solar atlases at least 4 mÅ. Grades: A ≤ 1 mÅ, B ≤ 5 mÅ,
  C ≤ 10 mÅ, D > 10 mÅ. *Quoted* precision and *demonstrated* uncertainty are therefore distinguished: 78 % of the element
  lines are quoted to ≥ 3 decimals after the ASD upgrade, but 35 % have σ ≤ 1 mÅ.

## Directory layout

```
index.html                 the page (generated — do not edit by hand)
src/template.html          page skeleton and hard-coded text
src/content/*.html         prose fragments substituted into {{PLACEHOLDERS}}
src/style.css              styles (light/dark themes)
src/spectro_core.js        computations: time scales, DE440s Chebyshev evaluation, Earth orientation, barycentric
                           correction, Galactic frames, air/vacuum formulas, Ciddor full model, Doppler conventions
src/ui.js                  DataTable, SpectrumView (zoomable stick spectra), SVG line charts, formatting
src/app_core.js            boot, element section, periodic table, explorer, finder
src/app_calc.js            calculators (precision, air/vacuum, barycentric, frames, cosmology)
src/app_sets.js            generic catalogue sets and static tables
scripts/                   data pipeline (below);  scripts/research/  research-agent scripts and table generators
data/handbook.json         parsed Handbook;  data/raw/nist_asd/  raw ASD CSV per spectrum
data/curated/*.csv         curated catalogues ('#' comment headers document every column)
data/notes/*.md|json       research notes with cited numbers
data/build/*.json          normalised data blocks embedded in the page
data/refs_*.json           reference lists;  data/ephem_cheb.json  DE440s Chebyshev fit
tests/                     headless-Chrome smoke test and numerical cross-checks
```

## Rebuilding

Python: `~/.venvs/astro313/bin/python` (astropy, numpy, pandas, jplephem). Node is used only for tests.

```sh
P=~/.venvs/astro313/bin/python
$P scripts/parse_handbook.py      # e-book HTML  -> data/handbook.json
$P scripts/fetch_asd.py           # NIST ASD per spectrum -> data/raw/nist_asd/ (network; ~193 queries, show_av=4,
                                  #   plus 33 vacuum queries (show_av=3) for 1850-2000 A, where show_av=4 already gives air)
$P scripts/match_asd.py           # Handbook <-> ASD matching (1850-2000 A vacuum lines matched to the ASD vacuum strings)
$P scripts/build_elements.py      # adopted wavelengths, uncertainties, flags -> data/build/elements.json
$P scripts/build_sets.py          # curated catalogues -> data/build/{ag_*,fh_classic,solar,telluric}.json
$P scripts/build_frames.py        # LSR/GSR/LG/CMB frame definitions -> data/build/frames.json
$P scripts/build_ephem.py         # DE440s (~/JPLDE/de440s.bsp) Chebyshev fit 1990-2066 -> data/ephem_cheb.json
$P scripts/build_refs.py          # reference list -> src/content/REFS.html (reports cited keys without entries)
$P scripts/build_page.py          # -> index.html   (--out <path> writes elsewhere)
```

## Validation

| check | result |
|---|---|
| barycentric correction vs astropy 7.2 `radial_velocity_correction` (DE440s), 300 random cases 1990–2030 | rms 0.26 mm/s, max 0.68 mm/s (`tests/test_bary.js`) |
| DE440s Chebyshev fit | Earth velocity ≤ 0.28 mm/s, position ≤ 0.27 km |
| TDB / BJD vs astropy | ≤ 40 µs (float64 JD resolution) |
| air/vacuum formulas vs validated reference implementation (284 published-value checks) | identical to machine precision (`tests/check_airvac.js`) |
| LSRK / LSRD / LSR / Galactocentric projections vs astropy transforms | ≤ 0.06 m/s (`tests/check_frames.py`) |
| redshift drift vs astropy Planck18 | ≤ 0.07 cm/s per decade (the page omits radiation/neutrinos) |
| page boot, all tables and calculators | no JS errors (`tests/headless.py`) |

Content was researched, written and adversarially fact-checked by independent agents against primary sources; the
research notes (`data/notes/`) record the source of every number.

## Known limitations

* 2527 element lines remain at 0.01 Å or coarser with no stated uncertainty: neither the Handbook nor ASD has a better
  laboratory value (mostly lanthanides, actinides and heavy transition metals). 302 lines where Handbook and ASD differ by
  more than twice their combined uncertainty are flagged ⚠; 118 unresolved multiplets keep a blend value (flag "mult").
  Seven ASD lines are still claimed by two Handbook lines (matching is not one-to-one).
* OH Meinel positions are formally good to 2.5–250 mÅ by HITRAN codes (internal consistency of independent lists is
  1–6 mÅ at low J); O₂, CO and most CO₂ lines are better than 1 mÅ.
* Observed UVES sky-line positions are accurate to ≈ 17 mÅ with setting-dependent offsets of ≈ +8 to +15 mÅ; use the
  laboratory/model positions for calibration.
* Telluric positions are zero-pressure HITRAN values; the sea-level pressure shift (+2 mÅ for O₂ up to +22 mÅ for CH₄)
  is given per line but line mixing and speed dependence are not modelled.
* Air wavelengths above 1.69 µm extrapolate the dispersion formulas.
* The published Cosby et al. (2006) identification table could not be downloaded (publisher blocks automated access);
  the UVES identifications are a reproduction of their procedure (97.5 % identified, as in the original).
* The barycentric calculator neglects polar motion (< 1 mm/s) and needs UT1−UTC from the user (default 0); the
  leap-second table is current through 2026-12-31.
