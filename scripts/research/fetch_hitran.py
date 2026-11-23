#!/usr/bin/env python
"""Download HITRAN line-by-line data (160-char .par records) for the main telluric
absorbers over 3000-25000 A (vacuum), i.e. 4000-33334 cm^-1.

Uses the HITRANonline line-by-line API endpoint that HAPI v1 used
(https://hitran.org/lbl/api).  Data = HITRAN2024 edition (as served Nov 2026).

usage: python fetch_hitran.py OUTDIR
"""
import json
import os
import sys
import time
import urllib.request
import http.client
import datetime

API = "https://hitran.org/lbl/api"

# HITRAN global isotopologue IDs (from https://hitran.org/docs/iso-meta/ , molparam.txt)
MOLS = {
    "H2O": dict(mol=1, iso=[1, 2, 3, 4, 5, 6, 129]),
    "CO2": dict(mol=2, iso=[7, 8, 9, 10, 11, 12, 13, 14, 121, 15, 120, 122]),
    "O3": dict(mol=3, iso=[16, 17, 18, 19, 20]),
    "N2O": dict(mol=4, iso=[21, 22, 23, 24, 25]),
    "CO": dict(mol=5, iso=[26, 27, 28, 29, 30, 31]),
    "CH4": dict(mol=6, iso=[32, 33, 34, 35]),
    "O2": dict(mol=7, iso=[36, 37, 38]),
}

# wavenumber chunk edges (cm^-1); 25000 A -> 4000 cm^-1, 3000 A -> 33333.3 cm^-1
EDGES = [3990, 4500, 5000, 5500, 6000, 6500, 7000, 8000, 9000, 10000, 11000,
         12000, 13000, 14000, 16000, 18000, 20000, 25000, 33400]


def fetch(isos, numin, numax):
    url = (f"{API}?iso_ids_list={','.join(str(i) for i in isos)}"
           f"&numin={numin}&numax={numax}&head=0&fixwidth=0&request_params=par_line")
    req = urllib.request.Request(url, headers={"User-Agent": "python-urllib/research-script"})
    with urllib.request.urlopen(req, timeout=600) as r:
        try:
            data = r.read()
        except http.client.IncompleteRead as e:
            # hitran.org sets Content-Length as if lines ended in CRLF (162 B/record)
            # but sends LF (161 B/record).  Accept only if the shortfall is exactly
            # one byte per complete 160-char record.
            data = e.partial
            nrec = data.count(b"\n")
            if not (len(data) == 161 * nrec and e.expected == nrec):
                raise
    return url, data


def main(outdir):
    os.makedirs(outdir, exist_ok=True)
    manifest = []
    for name, d in MOLS.items():
        for a, b in zip(EDGES[:-1], EDGES[1:]):
            fn = os.path.join(outdir, f"{name}_{a}_{b}.par")
            if os.path.exists(fn) and os.path.getsize(fn) >= 0 and os.path.exists(fn + ".ok"):
                continue
            t0 = time.time()
            try:
                url, data = fetch(d["iso"], a, b)
            except Exception as e:  # noqa
                print(name, a, b, "ERROR", e, flush=True)
                manifest.append(dict(molecule=name, numin=a, numax=b, error=str(e)))
                continue
            with open(fn, "wb") as f:
                f.write(data)
            open(fn + ".ok", "w").write("ok\n")
            n = data.count(b"\n")
            print(f"{name} {a}-{b}: {n} lines, {len(data)/1e6:.1f} MB, {time.time()-t0:.1f}s", flush=True)
            manifest.append(dict(molecule=name, numin=a, numax=b, url=url, nlines=n,
                                 bytes=len(data),
                                 utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
            time.sleep(1.0)
    mf = os.path.join(outdir, "manifest.json")
    old = json.load(open(mf)) if os.path.exists(mf) else []
    json.dump(old + manifest, open(mf, "w"), indent=1)


if __name__ == "__main__":
    main(sys.argv[1])
