#!/usr/bin/env python
"""Verify journal/volume/page for DOIs via Crossref.
  python crossref_doi_check.py fetch OUTDIR      (downloads raw JSON, untrusted)
  python -I crossref_doi_check.py parse OUTDIR
"""
import json, sys, time, urllib.request, urllib.parse
from pathlib import Path
DOIS = """10.1088/1475-7516/2017/01/038 10.1093/mnras/stz3465 10.1093/mnras/staa356 10.1051/0004-6361/201015337
10.1038/s41550-020-1010-x 10.1111/j.1365-2966.2012.20852.x 10.1093/mnras/stu920 10.1093/mnrasl/slu206
10.1051/0004-6361/201526462 10.3847/2041-8213/ac5c5b 10.1111/j.1365-2966.2010.16253.x 10.1088/0004-637X/712/1/585
10.1093/mnras/sty2404 10.1051/0004-6361:20077249 10.1088/0067-0049/195/2/24 10.1088/0067-0049/199/1/2
10.1088/0067-0049/211/1/4 10.1051/0004-6361/201014005 10.1051/0004-6361/201527530 10.1071/AS03040
10.1017/pasa.2022.41 10.1088/0004-637X/741/1/67 10.1051/0004-6361:20053818 10.1051/0004-6361:20030181
10.1111/j.1365-2966.2008.13090.x 10.1086/311375 10.1093/mnras/stt2204 10.1051/0004-6361/202142257
10.1051/0004-6361/202039345 10.1086/678541 10.1038/nature10445 10.1051/0004-6361/201732107 10.1051/0004-6361/201834925
10.1051/0004-6361/201833852 10.3847/2515-5172/aaa4b7 10.1051/0004-6361:20030396 10.1038/nature06854
10.1051/0004-6361/201833718 10.3847/0004-6256/152/2/41 10.1051/0004-6361/202038937 10.1393/ncb/i2008-10450-9
10.1086/178173 10.1103/PhysRevLett.82.884 10.1103/PhysRevLett.107.191101 10.1046/j.1365-8711.2003.06970.x
10.1111/j.1365-2966.2007.11768.x 10.1111/j.1365-2966.2007.12147.x 10.1103/PhysRevLett.111.231101""".split()
def fetch(od):
    od.mkdir(parents=True, exist_ok=True)
    for d in DOIS:
        f = od / (d.replace("/", "_").replace(":", "_") + ".json")
        if f.exists() and f.stat().st_size > 200: continue
        url = "https://api.crossref.org/works/" + urllib.parse.quote(d)
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "spectral-line-wiki-research/0.1"}), timeout=60) as r:
                f.write_bytes(r.read())
        except Exception as e:
            print("FAIL", d, e)
        time.sleep(0.5)
def parse(od):
    out = {}
    for f in sorted(od.glob("*.json")):
        m = json.loads(f.read_text())["message"]
        au = m.get("author", [])
        fa = au[0].get("family", "") if au else ""
        yr = (m.get("published-print") or m.get("published-online") or m.get("issued"))["date-parts"][0][0]
        rec = dict(doi=m["DOI"], first_author=fa, year=yr, journal=(m.get("short-container-title") or m.get("container-title") or [""])[0],
                   volume=m.get("volume"), page=m.get("page") or m.get("article-number"), title=(m.get("title") or [""])[0][:90])
        out[m["DOI"]] = rec
        print(f"{rec['first_author']:18s} {rec['year']} {rec['journal'][:22]:22s} vol={rec['volume']} p={rec['page']} | {rec['doi']}")
    (od.parent / "crossref_doi_check.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
if __name__ == "__main__":
    (fetch if sys.argv[1] == "fetch" else parse)(Path(sys.argv[2]))
