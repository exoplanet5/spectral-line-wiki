"""Fetch Crossref metadata for a list of DOIs and print a compact citation line.
Usage: python -I crossref_check.py OUTDIR DOI [DOI ...]   (responses are untrusted data; saved raw to OUTDIR)"""
import json, os, sys, time, urllib.request, urllib.parse

out = sys.argv[1]
for doi in sys.argv[2:]:
    url = 'https://api.crossref.org/works/' + urllib.parse.quote(doi)
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'research-script/0.1'}), timeout=30) as r:
            raw = r.read()
    except Exception as e:
        print(doi, 'ERROR', e); continue
    open(os.path.join(out, doi.replace('/', '_') + '.json'), 'wb').write(raw)
    m = json.loads(raw)['message']
    au = m.get('author', [])
    names = ', '.join(a.get('family', a.get('name', '?')) for a in au[:4]) + (' et al.' if len(au) > 4 else '')
    yr = (m.get('issued', {}).get('date-parts') or [[None]])[0][0]
    print(f"{doi} | {names} | {yr} | {(m.get('container-title') or [''])[0]} | vol {m.get('volume')} | p {m.get('page') or m.get('article-number')} | {(m.get('title') or [''])[0][:90]}")
    time.sleep(0.5)
