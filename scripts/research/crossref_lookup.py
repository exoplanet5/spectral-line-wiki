#!/usr/bin/env python
"""Look up bibliographic metadata on Crossref for references that are not
(well) covered by arXiv.  Two modes:

  python crossref_lookup.py fetch OUTDIR      # download raw JSON (untrusted)
  python -I crossref_lookup.py parse OUTDIR   # print top matches

Generic User-Agent, no personal information.
"""
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

Q = {
    "wilken_2012": "Wilken A spectrograph for exoplanet observations calibrated at the centimetre-per-second level Nature 2012",
    "steinmetz_2008": "Steinmetz Laser frequency combs for astronomical observations Science 2008",
    "ubachs_2016": "Ubachs Colloquium: Search for a drifting proton-electron mass ratio from H2 Reviews of Modern Physics 2016",
    "bagdonaite_2013_sci": "Bagdonaite A stringent limit on a drifting proton-to-electron mass ratio from alcohol in the early universe Science 2013",
    "rousselot_2000": "Rousselot Night-sky spectral atlas of OH emission lines in the near-infrared Astronomy and Astrophysics 2000",
    "oliva_origlia_1992": "Oliva Origlia The OH airglow spectrum - A calibration source for infrared spectrometers 1992",
    "butler_1996": "Butler Marcy Attaining Doppler Precision of 3 m/s PASP 1996",
    "marcy_butler_1992": "Marcy Butler Precision radial velocities with an iodine absorption cell PASP 1992",
    "sandage_1962": "Sandage The change of redshift and apparent luminosity of galaxies due to the deceleration of selected expanding universes 1962",
    "dravins_1982": "Dravins Photospheric spectrum line asymmetries and wavelength shifts Annual Review 1982",
    "dravins_1981": "Dravins Lindegren Nordlund Solar granulation influence of convection on spectral line asymmetries and wavelength shifts 1981",
    "pepe_2021": "Pepe ESPRESSO at VLT On-sky performance and first results Astronomy Astrophysics 2021",
    "kerr_lyndenbell_1986": "Kerr Lynden-Bell Review of galactic constants MNRAS 1986",
    "carr_davis_2021": "Carr Davis heliocentric redshift correction CMB frame NED 2021",
    "mcvittie_1962": "McVittie Appendix to the change of redshift and apparent luminosity of galaxies 1962",
    "allende_prieto_1998": "Allende Prieto Garcia Lopez Fe I line shifts in the optical spectrum of the Sun 1998",
    "kerber_2008": "Kerber Nave Sansonetti The spectrum of Th-Ar hollow cathode lamps in the 691-5804 nm region 2008",
    "palmer_engleman_1983": "Palmer Engleman Atlas of the thorium spectrum Los Alamos 1983",
    "stumpff_1980": "Stumpff Two self-consistent FORTRAN subroutines for the computation of the Earth's motion 1980",
    "irwin_fukushima_1999": "Irwin Fukushima A numerical time ephemeris of the Earth 1999",
    "wildi_2010": "Wildi Pepe Chazelas A Fabry-Perot calibrator of the HARPS radial velocity spectrograph SPIE 2010",
    "terrien_2021_hpf_lfc": "Terrien Habitable-zone Planet Finder laser frequency comb wavelength calibration 2021",
    "molaro_2013_alpha_uves": "Molaro UVES large programme for testing fundamental physics 2013",
    "murphy_cooksey_2017": "Murphy Cooksey Subaru Telescope limits on cosmological variations in the fine-structure constant 2017",
    "roca_cortes_palle_2014": "Roca Cortes Palle The Mark-I helioseismic spectrometer solar gravitational redshift 2014",
    "dzuba_1999_pra": "Dzuba Flambaum Webb Calculations of the relativistic effects in many-electron atoms and space-time variation of fundamental constants Physical Review A 1999",
    "savedoff_1956": "Savedoff Physical constants in extra-galactic nebulae Nature 1956",
    "greenstein_trimble_1967": "Greenstein Trimble The Einstein redshift in white dwarfs 1967",
    "barstow_2005_siriusb": "Barstow Bond Holberg Hubble Space Telescope spectroscopy of the Balmer lines in Sirius B 2005",
    "pound_rebka_1960": "Pound Rebka Apparent weight of photons 1960",
    "kurucz_1984": "Kurucz Furenlid Brault Testerman Solar flux atlas from 296 to 1300 nm 1984",
    "karachentsev_makarov_1996": "Karachentsev Makarov The galaxy motion relative to nearby galaxies and the local velocity field 1996",
    "ubachs_2007_jms": "Ubachs Buning Eikema Reinhold On a possible variation of the proton-to-electron mass ratio H2 spectra in the line of sight of high-redshift quasars and in the laboratory 2007",
    "osterbrock_1996": "Osterbrock Fulbright Martel Night-sky high-resolution spectral atlas of OH and O2 emission lines for echelle spectrograph wavelength calibration 1996",
    "abrams_1994": "Abrams Davis Rao Engleman Brault High-resolution Fourier transform spectroscopy of the Meinel system of OH 1994",
    "bennett_2003_wmap": "Bennett First-Year Wilkinson Microwave Anisotropy Probe observations preliminary maps and basic results 2003",
    "reinhold_2006": "Reinhold Indication of a cosmological variation of the proton-electron mass ratio based on laboratory measurement and reanalysis of H2 spectra 2006",
    "yahil_1977": "Yahil Tammann Sandage The local group: the solar motion relative to its centroid 1977",
}


def fetch(outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    for k, q in Q.items():
        d = outdir / f"{k}.json"
        if d.exists() and d.stat().st_size > 200:
            continue
        url = ("https://api.crossref.org/works?" + urllib.parse.urlencode({
            "query.bibliographic": q, "rows": 3,
            "select": "DOI,title,author,container-title,volume,page,issued,article-number"}))
        try:
            with urllib.request.urlopen(urllib.request.Request(
                    url, headers={"User-Agent": "spectral-line-wiki-research/0.1"}), timeout=60) as r:
                d.write_bytes(r.read())
            print("ok", k)
        except Exception as e:  # noqa: BLE001
            print("FAIL", k, e)
        time.sleep(1.0)


def parse(outdir):
    res = {}
    for f in sorted(outdir.glob("*.json")):
        items = json.loads(f.read_text())["message"]["items"]
        print("===", f.stem)
        res[f.stem] = []
        for it in items:
            au = it.get("author", [])
            names = ", ".join(f"{a.get('family','')} {a.get('given','')[:1]}." for a in au[:4])
            yr = it.get("issued", {}).get("date-parts", [[None]])[0][0]
            rec = dict(doi=it.get("DOI"), title=(it.get("title") or [""])[0],
                       journal=(it.get("container-title") or [""])[0], volume=it.get("volume"),
                       page=it.get("page") or it.get("article-number"), year=yr, authors=names)
            res[f.stem].append(rec)
            print(f"  {rec['authors']} ({yr}) {rec['title'][:110]} | {rec['journal']} {rec['volume']} {rec['page']} | doi:{rec['doi']}")
    (outdir.parent / "crossref_parsed.json").write_text(json.dumps(res, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    mode, od = sys.argv[1], Path(sys.argv[2])
    fetch(od) if mode == "fetch" else parse(od)
