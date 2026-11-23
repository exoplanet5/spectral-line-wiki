#!/usr/bin/env python
"""Crossref metadata for the references of precision_rationale.md.
  python    prec_crossref.py fetch OUTDIR   (downloads raw JSON = untrusted data)
  python -I prec_crossref.py parse OUTDIR   (writes OUTDIR/parsed.json, prints a table)
Bibliographic queries (no DOI known) are listed in QUERIES; their first hit is only
accepted by the human reader after checking the printed title.
"""
import json, sys, time, urllib.request, urllib.parse
from pathlib import Path

DOIS = """10.1051/0004-6361:20010730 10.1051/0004-6361/202038306 10.3847/1538-3881/ab811d 10.1117/12.2233002
10.1002/asna.201512172 10.1051/0004-6361/202245217 10.1088/1538-3873/ab0075 10.1051/0004-6361/201832763
10.3847/1538-3881/ac882b 10.1088/0004-6256/146/2/32 10.18727/0722-6691/5117 10.1093/mnras/stad557
10.1088/1538-3873/128/964/066001 10.1117/12.2232103 10.1051/0004-6361:20077249 10.1051/0004-6361/201834925
10.1016/j.asr.2015.03.007 10.1093/mnras/stt2204 10.1051/0004-6361/202038937 10.1038/nature11092
10.1093/pasj/54.6.855 10.1117/12.176725 10.1117/12.395512 10.1117/12.2056431 10.1088/1674-4527/12/9/003
10.1103/PhysRevLett.107.023001 10.1364/AO.35.001566 10.1088/0026-1394/2/2/002 10.1364/JOSA.62.000958
10.1051/0004-6361/201527530 10.1111/j.1365-2966.2008.13090.x 10.1038/s41550-020-1010-x 10.1038/nature06854
10.1364/JOSA.43.000339 10.1088/0026-1394/31/4/006 10.1063/1.1800011 10.1103/PhysRevA.87.032504
10.1103/PhysRevLett.109.259901 10.1103/PhysRevLett.124.063002 10.1051/0004-6361/201833910""".split()

QUERIES = {
    "schwab2016_neid": "Design of NEID, an extreme precision Doppler spectrograph for WIYN Schwab 2016",
    "lenz_ayres1992": "Errors associated with fitting Gaussian profiles to noisy emission-line spectra Lenz Ayres",
    "milakovic2020": "Precision and consistency of astrocombs Milakovic",
    "kurucz_hfs_prochaska2000": "On the hyperfine splitting corrections of Sc V Mn Co Prochaska McWilliam",
}


def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "spectral-line-wiki-research/0.1"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read()


def fetch(od):
    od.mkdir(parents=True, exist_ok=True)
    for d in DOIS:
        f = od / (d.replace("/", "_").replace(":", "_") + ".json")
        if f.exists() and f.stat().st_size > 200:
            continue
        try:
            f.write_bytes(get("https://api.crossref.org/works/" + urllib.parse.quote(d)))
        except Exception as e:
            print("FAIL", d, e)
        time.sleep(1.0)
    for k, q in QUERIES.items():
        f = od / ("query_" + k + ".json")
        if f.exists() and f.stat().st_size > 200:
            continue
        try:
            f.write_bytes(get("https://api.crossref.org/works?rows=3&query.bibliographic=" + urllib.parse.quote(q)))
        except Exception as e:
            print("FAIL", k, e)
        time.sleep(1.0)


def summarise(m):
    au = m.get("author", [])
    names = [a.get("family", "") + (", " + a.get("given", "")[:1] + "." if a.get("given") else "") for a in au]
    yr = (m.get("published-print") or m.get("published-online") or m.get("issued") or {}).get("date-parts", [[None]])[0][0]
    return {"doi": m.get("DOI"), "title": " ".join((m.get("title") or [""])[0].split()),
            "authors": names[:6], "n_authors": len(names), "year": yr,
            "journal": (m.get("container-title") or [""])[0], "volume": m.get("volume"),
            "issue": m.get("issue"), "page": m.get("page") or m.get("article-number")}


def parse(od):
    out = {}
    for f in sorted(od.glob("*.json")):
        if f.name == "parsed.json":
            continue
        j = json.loads(f.read_text())
        msg = j.get("message", {})
        if "items" in msg:
            out[f.stem] = [summarise(it) for it in msg["items"]]
        else:
            out[f.stem] = summarise(msg)
    (od / "parsed.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    for k, v in out.items():
        for s in (v if isinstance(v, list) else [v]):
            print(f"{k[:28]:28s} | {s['year']} | {'; '.join(s['authors'][:3])}{' et al.' if s['n_authors'] > 3 else ''} | "
                  f"{s['journal'][:40]} {s['volume']}, {s['page']} | {s['title'][:90]} | {s['doi']}")


if __name__ == "__main__":
    (fetch if sys.argv[1] == "fetch" else parse)(Path(sys.argv[2]))
