#!/usr/bin/env python
"""Download HITRAN2024 total internal partition sums Q(T) (TIPS-2025, 1 K steps)
for the isotopologues used in the telluric list.  usage: fetch_hitran_q.py OUTDIR"""
import os
import sys
import time
import urllib.request

GIDS = [1, 2, 3, 4, 5, 6, 129,
        7, 8, 9, 10, 11, 12, 13, 14, 121, 15, 120, 122,
        16, 17, 18, 19, 20,
        21, 22, 23, 24, 25,
        26, 27, 28, 29, 30, 31,
        32, 33, 34, 35,
        36, 37, 38]

out = sys.argv[1]
os.makedirs(out, exist_ok=True)
for g in GIDS:
    fn = os.path.join(out, f"q{g}.txt")
    if os.path.exists(fn) and os.path.getsize(fn) > 1000:
        continue
    url = f"https://hitran.org/data/Q/q{g}.txt"
    req = urllib.request.Request(url, headers={"User-Agent": "python-urllib/research-script"})
    with urllib.request.urlopen(req, timeout=120) as r:
        open(fn, "wb").write(r.read())
    print(g, os.path.getsize(fn), flush=True)
    time.sleep(0.3)
