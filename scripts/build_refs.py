#!/usr/bin/env python
"""Compile the reference list (src/content/REFS.html) from
   * data/refs_base.json            -- references used by the pipeline / calculators / methods text
   * data/refs_writers.json         -- refs returned by the content-writing agents
and check that every <a class="cite" href="#ref-KEY"> in src/content/*.html and in the static tables
has an entry (and report unused entries)."""
import glob, html, json, os, re, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))


def norm_doi(d):
    d = (d or '').strip()
    d = re.sub(r'^(https?://(dx\.)?doi\.org/|doi:)', '', d, flags=re.I)
    return d


def main():
    refs = {}
    for fn in ('refs_base.json', 'refs_writers.json'):
        p = os.path.join(ROOT, 'data', fn)
        if not os.path.exists(p):
            continue
        for r in json.load(open(p)):
            k = r['key'].strip().lower()
            if k in refs:
                # keep the longer citation, fill missing ids
                old = refs[k]
                if fn == 'refs_base.json' or len(r.get('citation', '')) > len(old.get('citation', '')):
                    old['citation'] = r['citation']
                for f in ('doi', 'bibcode', 'url'):
                    if not old.get(f) and r.get(f):
                        old[f] = r[f]
            else:
                refs[k] = dict(r, key=k)
    # cited keys
    cited = set()
    texts = [open(f, encoding='utf-8').read() for f in glob.glob(os.path.join(ROOT, 'src', 'content', '*.html'))]
    texts += [open(f, encoding='utf-8').read() for f in glob.glob(os.path.join(ROOT, 'data', 'build', 'st_*.json'))]
    texts += [open(os.path.join(ROOT, 'src', 'template.html'), encoding='utf-8').read()]
    for t in texts:
        cited |= {k.lower() for k in re.findall(r'#ref-([A-Za-z0-9_\-]+)', t)}
    missing = sorted(cited - set(refs))
    unused = sorted(set(refs) - cited)
    if missing:
        print('MISSING reference entries for cited keys:', missing, file=sys.stderr)
    print(f'{len(refs)} references, {len(cited)} cited, {len(unused)} uncited (listed anyway)')

    def sort_key(k):
        m = re.match(r'([a-z\-_]+)(\d{4})', k)
        return (m.group(1), m.group(2)) if m else (k, '')

    items = []
    for k in sorted(refs, key=sort_key):
        r = refs[k]
        links = []
        doi = norm_doi(r.get('doi'))
        if doi:
            links.append(f'<a href="https://doi.org/{html.escape(doi)}">doi:{html.escape(doi)}</a>')
        if r.get('bibcode'):
            links.append(f'<a href="https://ui.adsabs.harvard.edu/abs/{html.escape(r["bibcode"])}">{html.escape(r["bibcode"])}</a>')
        if r.get('url') and not doi and not (r.get('bibcode') and 'adsabs' in r['url']):
            links.append(f'<a href="{html.escape(r["url"])}">link</a>')
        items.append(f'<li id="ref-{html.escape(k)}"><span class="muted small">[{html.escape(k)}]</span> {html.escape(r["citation"])}'
                     + (f' — {" · ".join(links)}' if links else '') + '</li>')
    out = ('<p class="small muted">Inline citations link here; every catalogue table also names its source per row. '
           'Data files, scripts and research notes used to build this page are in the project directory '
           '(<code>data/</code>, <code>scripts/</code>, <code>data/notes/</code>).</p>\n<ol class="refs">\n' + '\n'.join(items) + '\n</ol>\n')
    open(os.path.join(ROOT, 'src', 'content', 'REFS.html'), 'w', encoding='utf-8').write(out)
    return missing


if __name__ == '__main__':
    main()
