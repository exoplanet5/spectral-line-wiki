"""Validate ADS bibcodes via the public link_gateway (PUB_HTML / EPRINT_HTML redirect) and compare to expected DOI.
Usage: python -I check_bibcodes.py   (prints bibcode, status, redirect)"""
import sys, urllib.parse, urllib.request, http.client

PAIRS = [
 ("2014PASP..126..838W", "10.1086/678541"), ("2003A&A...401.1185L", "10.1051/0004-6361:20030181"),
 ("2010MNRAS.403.1829S", "10.1111/j.1365-2966.2010.16253.x"), ("1998MNRAS.298..387D", "10.1046/j.1365-8711.1998.01600.x"),
 ("2018A&A...615L..15G", "10.1051/0004-6361/201833718"), ("2019A&A...625L..10G", "10.1051/0004-6361/201935656"),
 ("2021A&A...647A..59G", "10.1051/0004-6361/202040208"), ("2022A&A...657L..12G", "10.1051/0004-6361/202142465"),
 ("2019ApJ...885..131R", "10.3847/1538-4357/ab4a11"), ("2020ApJ...892...39R", "10.3847/1538-4357/ab76cd"),
 ("2004ApJ...616..872R", "10.1086/424960"), ("2020A&A...641A...1P", "10.1051/0004-6361/201833880"),
 ("2020A&A...641A...6P", "10.1051/0004-6361/201833910"), ("1996ApJ...473..576F", "10.1086/178173"),
 ("1993ApJ...419....1K", "10.1086/173453"), ("1999AJ....118..337C", "10.1086/300942"),
 ("2000ApJ...529..786M", "10.1086/308304"), ("2008ApJ...676..184T", "10.1086/527428"),
 ("2012ApJ...753....8V", "10.1088/0004-637X/753/1/8"), ("2014MNRAS.443.1688D", "10.1093/mnras/stu1210"),
 ("2011ApJ...741...67D", "10.1088/0004-637X/741/1/67"), ("2022PASA...39...46C", "10.1017/pasa.2022.41"),
 ("2023ApJ...944...94T", "10.3847/1538-4357/ac94d8"), ("2019MNRAS.482.1417B", "10.1093/mnras/sty2813"),
 ("2021A&A...649A...9G", "10.1051/0004-6361/202039734"), ("2016ARA&A..54..529B", "10.1146/annurev-astro-081915-023441"),
 ("2018RNAAS...2....4K", "10.3847/2515-5172/aaa4b7"), ("2019MNRAS.489.2395T", "10.1093/mnras/stz2181"),
 ("2017ApJ...837...18B", "10.3847/1538-4357/aa5ead"), ("2019AJ....158...40B", "10.3847/1538-3881/ab24c3"),
 ("1986MNRAS.221.1023K", "10.1093/mnras/221.4.1023"), ("1977ApJ...217..903Y", "10.1086/155636"),
 ("1996AJ....111..794K", "10.1086/117825"), ("2023A&A...674A..37G", "10.1051/0004-6361/202243797"),
 ("2018RNAAS...2..210D", "10.3847/2515-5172/aaef8b"), ("2019RAA....19...68D", "10.1088/1674-4527/19/5/68"),
 ("2009ApJ...700..137R", "10.1088/0004-637X/700/1/137"), ("2010PASP..122..935E", "10.1086/655938"),
 ("1962ApJ...136..319S", "10.1086/147385"), ("1998ApJ...499L.111L", "10.1086/311375"),
 ("1965gast.book...61D", None),
]

def head(bib, kind):
    path = "/link_gateway/" + urllib.parse.quote(bib) + "/" + kind
    c = http.client.HTTPSConnection("ui.adsabs.harvard.edu", timeout=30)
    c.request("GET", path, headers={"User-Agent": "research-script/0.1"})
    r = c.getresponse(); loc = r.getheader("Location"); st = r.status; c.close()
    return st, loc

for bib, doi in PAIRS:
    st, loc = head(bib, "PUB_HTML")
    if st != 302:
        st2, loc2 = head(bib, "ADS_PDF")
        print(f"{bib:22s} PUB_HTML {st}  ADS_PDF {st2} {loc2}")
        continue
    ok = (doi is not None and loc and urllib.parse.unquote(loc).lower().endswith(doi.lower()))
    print(f"{bib:22s} {st} {'DOI-MATCH' if ok else 'CHECK'} {loc}")
