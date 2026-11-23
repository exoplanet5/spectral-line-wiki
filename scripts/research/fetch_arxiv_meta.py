#!/usr/bin/env python
"""Fetch arXiv Atom metadata (title, authors, journal_ref, DOI, abstract) for the
relativity / cosmology / calibration references used in
data/notes/relativity_cosmology.md.

Usage:
    python fetch_arxiv_meta.py OUTDIR

Writes one raw Atom XML file per query into OUTDIR (untrusted download) and
nothing else.  Parsing is done by parse_arxiv_meta.py (run with python -I).
arXiv API etiquette: >= 3 s between requests.
"""
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

# key -> arXiv search_query (title phrase, optionally author)
QUERIES = {
    "davis_lineweaver_2004": 'ti:"Expanding Confusion" AND au:Davis',
    "davis_2011_pecvel_sn": 'ti:"effect of peculiar velocities on supernova cosmology"',
    "carr_2022_pantheonplus_z": 'ti:"Improving the Redshifts and Peculiar Velocities"',
    "calcino_davis_2017": 'ti:"need for accurate redshifts in supernova cosmology"',
    "davis_scrimgeour_2014": 'ti:"Deriving accurate peculiar velocities"',
    "greisen_2006_fits3": 'ti:"Representations of spectral coordinates in FITS"',
    "lindegren_dravins_2003": 'ti:"fundamental definition of radial velocity"',
    "dravins_1999_astrometric_rv": 'ti:"Astrometric radial velocities" AND au:Dravins',
    "liske_2008_codex": 'ti:"Cosmic dynamics in the era of Extremely Large Telescopes"',
    "loeb_1998": 'ti:"Direct Measurement of Cosmological Parameters from the Cosmic Deceleration"',
    "webb_1999_alpha": 'ti:"Search for time variation of the fine structure constant" AND au:Webb',
    "dzuba_1999_q": 'ti:"Space-time variation of physical constants and relativistic corrections in atoms"',
    "murphy_berengut_2014": 'ti:"Laboratory atomic transition data for precise optical quasar absorption spectroscopy"',
    "murphy_2003_alpha": 'ti:"Further evidence for a variable fine-structure constant from Keck"',
    "webb_2011_dipole": 'ti:"Indications of a spatial variation of the fine structure constant"',
    "murphy_2022_espresso": 'ti:"Fundamental physics with ESPRESSO" AND ti:"HE 0515"',
    "murphy_2007_thar": 'ti:"Selection of ThAr lines for wavelength calibration of echelle spectra"',
    "ubachs_2016_mu": 'ti:"drifting proton-electron mass ratio from H2"',
    "bagdonaite_2013_methanol": 'ti:"stringent constraint on a drifting proton-to-electron mass ratio from alcohol"',
    "redman_2014_thorium": 'ti:"spectrum of thorium from 250 nm to 5500 nm"',
    "lovis_pepe_2007": 'ti:"new list of thorium and argon spectral lines in the visible"',
    "redman_2011_uranium": 'ti:"infrared spectrum of uranium hollow cathode lamps"',
    "redman_2012_une": 'ti:"atlas of uranium-neon in the H band"',
    "steinmetz_2008_lfc": 'ti:"Laser frequency combs for astronomical observations"',
    "wilken_2012_lfc": 'ti:"spectrograph for exoplanet observations calibrated at the centimetre-per-second level"',
    "murphy_2007_lfc": 'ti:"High-precision wavelength calibration of astronomical spectrographs with laser frequency combs"',
    "li_2008_lfc": 'ti:"laser frequency comb that enables radial velocity measurements"',
    "probst_2020_lfc": 'ti:"crucial test for astronomical spectrograph calibration with frequency combs"',
    "milakovic_2020_lfc": 'ti:"Precision and consistency of astrocombs"',
    "figueira_2010_o2": 'ti:"Evaluating the stability of atmospheric lines with HARPS"',
    "bauer_2015_fp": 'ti:"Calibrating echelle spectrographs with Fabry-Perot etalons"',
    "wright_eastman_2014": 'ti:"Barycentric Corrections at 1 cm/s for precise Doppler velocities"',
    "kanodia_wright_2018": 'ti:"barycorrpy"',
    "wojtak_2011_cluster_gz": 'ti:"Gravitational redshift of galaxies in clusters as predicted by general relativity"',
    "falcon_2010_wd": 'ti:"gravitational redshift determination of the mean mass of white dwarfs"',
    "joyce_2018_siriusb": 'ti:"gravitational redshift of Sirius B"',
    "gonzalez_hernandez_2020_sun_gz": 'ti:"solar gravitational redshift from HARPS-LFC Moon spectra"',
    "reiners_2016_iag": 'ti:"IAG solar flux atlas"',
    "lohner_bottcher_2019_lars": 'ti:"Convective blueshifts in the solar atmosphere"',
    "pasquini_2011_gz_m67": 'ti:"Gravitational redshifts in main-sequence and giant stars"',
    "planck2018_vi": 'ti:"Planck 2018 results. VI. Cosmological parameters"',
    "planck2018_i": 'ti:"Planck 2018 results. I. Overview"',
    "riess_2022_shoes": 'ti:"Comprehensive Measurement of the Local Value of the Hubble Constant with 1 km"',
    "schoenrich_2010_lsr": 'ti:"Local kinematics and the local standard of rest"',
    "kuerster_2003_barnard": 'ti:"low-level radial velocity variability in Barnard"',
    "prsa_2016_nominal": 'ti:"Nominal values for selected solar and planetary quantities"',
    "martins_2024_andes_cosmo": 'ti:"Cosmology and fundamental physics with the ELT-ANDES spectrograph"',
    "marconi_2022_andes": 'ti:"ANDES, the high resolution spectrograph for the ELT"',
    "pepe_2021_espresso": 'ti:"ESPRESSO at VLT. On-sky performance and first results"',
    "rousselot_2000_oh": 'ti:"Night-sky spectral atlas of OH emission lines in the near-infrared"',
    "cristiani_2007_codex": 'ti:"CODEX-ESPRESSO experiment"',
    "cersullo_2019_fp": 'ti:"infrared Fabry-Perot-based radial-velocity-reference module"',
    "schmidt_2021_espresso_lfc_alpha": 'ti:"Fundamental physics with ESPRESSO: Towards an accurate wavelength calibration"',
    "chaplin_2019_granulation_rv": 'ti:"Filtering solar-like oscillations for exoplanet detection in radial velocity observations"',
}


def main():
    outdir = Path(sys.argv[1])
    outdir.mkdir(parents=True, exist_ok=True)
    only = set(sys.argv[2:])
    for key, q in QUERIES.items():
        if only and key not in only:
            continue
        dest = outdir / f"{key}.xml"
        if dest.exists() and dest.stat().st_size > 1500:
            continue
        url = ("https://export.arxiv.org/api/query?" +
               urllib.parse.urlencode({"search_query": q, "max_results": 5}))
        req = urllib.request.Request(url, headers={"User-Agent": "spectral-line-wiki-research/0.1"})
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                data = r.read()
            dest.write_bytes(data)
            print(f"ok   {key}  {len(data)} bytes")
        except Exception as exc:  # noqa: BLE001
            print(f"FAIL {key}: {exc}")
        time.sleep(3.2)


if __name__ == "__main__":
    main()
