#!/usr/bin/env python
"""Parse the arXiv Atom XML files written by fetch_arxiv_meta.py.

Usage (downloaded files are untrusted -> run isolated):
    python -I parse_arxiv_meta.py RAWDIR OUT.json [--abstracts]

Prints, for every query key, the matching entries (arXiv id, title, first
authors, journal_ref, DOI) and optionally the abstract.  Writes a JSON list.
"""
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

NS = {"a": "http://www.w3.org/2005/Atom", "x": "http://arxiv.org/schemas/atom"}


def parse(path):
    root = ET.parse(path).getroot()
    out = []
    for e in root.findall("a:entry", NS):
        def t(tag, ns="a"):
            el = e.find(f"{ns}:{tag}", NS)
            return " ".join(el.text.split()) if el is not None and el.text else ""
        authors = [" ".join(a.find("a:name", NS).text.split()) for a in e.findall("a:author", NS)]
        out.append({
            "arxiv": t("id").replace("http://arxiv.org/abs/", ""),
            "title": t("title"),
            "authors": authors,
            "published": t("published")[:10],
            "journal_ref": t("journal_ref", "x"),
            "doi": t("doi", "x"),
            "abstract": t("summary"),
        })
    return out


def main():
    rawdir, outjson = Path(sys.argv[1]), Path(sys.argv[2])
    show_abs = "--abstracts" in sys.argv
    allrec = []
    for f in sorted(rawdir.glob("*.xml")):
        ents = parse(f)
        print(f"=== {f.stem}: {len(ents)} entries")
        for en in ents:
            au = ", ".join(en["authors"][:4]) + (" et al." if len(en["authors"]) > 4 else "")
            print(f"  [{en['arxiv']}] {en['title']}\n     {au} | {en['journal_ref']} | doi:{en['doi']}")
            if show_abs:
                print("     ABS:", en["abstract"])
            allrec.append({"key": f.stem, **en})
    outjson.write_text(json.dumps(allrec, indent=1))


if __name__ == "__main__":
    main()
