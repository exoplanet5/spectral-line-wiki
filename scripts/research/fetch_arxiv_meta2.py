#!/usr/bin/env python
"""Second-round arXiv metadata queries (same conventions as fetch_arxiv_meta.py).
Usage: python fetch_arxiv_meta2.py OUTDIR"""
import sys, time, urllib.parse, urllib.request
from pathlib import Path
Q = {
 "wilken_2012_lfc_b": 'au:Wilken AND ti:"centimetre"',
 "ubachs_2016_mu_b": 'au:Ubachs AND ti:"proton-electron mass ratio" AND ti:Colloquium',
 "bagdonaite_2013_b": 'au:Bagdonaite AND ti:methanol',
 "bagdonaite_2013_c": 'au:Bagdonaite AND ti:"mass ratio" AND ti:alcohol',
 "pepe_2021_espresso_b": 'au:Pepe AND ti:"ESPRESSO at VLT"',
 "rousselot_2000_oh_b": 'au:Rousselot AND ti:"OH emission lines"',
 "cersullo_2019_fp_b": 'au:Cersullo AND ti:"Fabry"',
 "marconi_2022_andes_b": 'au:Marconi AND ti:ANDES',
 "carr_davis_2021": 'au:Carr AND au:Davis AND ti:"heliocentric"',
 "gravity_2018_s2": 'ti:"Detection of the gravitational redshift in the orbit of the star S2"',
 "kerber_2008_thar": 'au:Kerber AND ti:"thorium"',
 "nave_2012_feii": 'au:Nave AND ti:"Fe II"',
 "king_2012_alpha": 'au:King AND ti:"Spatial variation in the fine-structure constant"',
 "dumusque_2021_hcl": 'au:Dumusque AND ti:"Th-Ar"',
 "fixsen_1996_cobe": 'au:Fixsen AND ti:"Cosmic Microwave Background Spectrum from the Full COBE"',
 "cooke_2020_drift": 'au:Cooke AND ti:"redshift drift"',
 "kim_2015_drift_design": 'ti:"Sandage-Loeb" AND ti:"test"',
 "molaro_2022_espresso_sun": 'au:Molaro AND ti:"ESPRESSO" AND ti:"solar"',
 "dravins_2008_lab_wavelengths": 'au:Dravins AND ti:"spectral line"',
}
out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
for k, q in Q.items():
    d = out / f"{k}.xml"
    if d.exists() and d.stat().st_size > 1500: continue
    url = "https://export.arxiv.org/api/query?" + urllib.parse.urlencode({"search_query": q, "max_results": 6})
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "spectral-line-wiki-research/0.1"}), timeout=60) as r:
            b = r.read()
        d.write_bytes(b); print("ok", k, len(b))
    except Exception as e:
        print("FAIL", k, e)
    time.sleep(3.2)
