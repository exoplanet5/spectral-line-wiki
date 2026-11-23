#!/usr/bin/env python
"""Assemble the self-contained page  index.html  from src/ (template, CSS, JS) and data/build/*.json.

Every dataset is embedded as <script type="application/json" id="data-...">; JSON text is
escaped so that "</script" cannot terminate the block.
"""
import json, os, re, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
SRC = os.path.join(ROOT, 'src')
BUILD = os.path.join(ROOT, 'data', 'build')


def columnar(rows, cols=None):
    if cols is None:
        cols = []
        for r in rows:
            for k in r:
                if k not in cols:
                    cols.append(k)
    return {'cols': cols, 'rows': [[r.get(c) for c in cols] for r in rows]}


def jdump(obj):
    s = json.dumps(obj, separators=(',', ':'), ensure_ascii=False, allow_nan=False)
    return s.replace('</', '<\\/')


def main():
    data = {}
    E = json.load(open(os.path.join(BUILD, 'elements.json')))
    keep = ['Z', 'el', 'ion', 'hb', 'med', 'w', 'u', 'us', 'src', 'air', 'vac', 'wn', 'int', 'fl', 'ref', 'm', 'dl',
            'ai', 'lo', 'up', 'A', 'acc', 'Ei', 'Ek', 'aw', 'ak', 'au', 'ty', 'lr', 'P', 'lv', 'hA', 'dx', 'bl']
    for r in E['lines']:
        # the page re-derives: air (= native string when med 'a'), vac when med 'v', wn = 1e8/vac
        r.pop('air', None); r.pop('wn', None)
        if r['med'] == 'v':
            r.pop('vac', None)
        else:
            r['vac'] = round(r['vac'], 5)
        if r.get('aw') and r['aw'] == r['w']:
            r.pop('aw')
        if r.get('lo'):
            r.pop('lv', None)           # ASD classification supersedes the Handbook level text
        for k in list(r):
            if r[k] in ('', None):
                r.pop(k)
    keep = [k for k in keep if k not in ('air', 'wn')]
    el = columnar(E['lines'], keep)
    el['dicts'] = {}
    for j, c in enumerate(el['cols']):
        vals = [t[j] for t in el['rows'] if t[j] is not None]
        if vals and all(isinstance(v, str) for v in vals) and len(set(vals)) < 0.35 * len(vals):
            uniq = sorted(set(vals)); ix = {v: i for i, v in enumerate(uniq)}
            for t in el['rows']:
                if t[j] is not None:
                    t[j] = ix[t[j]]
            el['dicts'][c] = uniq
    el['elements'] = E['elements']
    el['stats'] = E['stats']
    data['elements'] = el
    for name in sorted(os.listdir(BUILD)):
        if name.endswith('.json') and name != 'elements.json':
            data[name[:-5]] = json.load(open(os.path.join(BUILD, name)))
    eph = json.load(open(os.path.join(ROOT, 'data', 'ephem_cheb.json')))
    eph.pop('test_vectors', None)
    data['ephem'] = eph

    tpl = open(os.path.join(SRC, 'template.html'), encoding='utf-8').read()
    css = open(os.path.join(SRC, 'style.css'), encoding='utf-8').read()
    js = '\n'.join(open(os.path.join(SRC, f), encoding='utf-8').read()
                   for f in ('spectro_core.js', 'ui.js', 'app_core.js', 'app_calc.js', 'app_sets.js')
                   if os.path.exists(os.path.join(SRC, f)))
    blocks = '\n'.join(f'<script type="application/json" id="data-{k}">{jdump(v)}</script>' for k, v in data.items())
    page = tpl.replace('/*{{STYLE}}*/', css).replace('<!--{{DATA}}-->', blocks).replace('/*{{SCRIPT}}*/', js.replace('</script', '<\\/script'))
    meta = json.load(open(os.path.join(ROOT, 'data', 'build_meta.json'))) if os.path.exists(os.path.join(ROOT, 'data', 'build_meta.json')) else {}
    for k, v in meta.items():
        page = page.replace('{{' + k + '}}', str(v))
    cdir = os.path.join(SRC, 'content')
    for k in set(re.findall(r'\{\{([A-Z_]+)\}\}', page)):
        fn = os.path.join(cdir, k + '.html')
        if os.path.exists(fn):
            page = page.replace('{{' + k + '}}', open(fn, encoding='utf-8').read())
    left = re.findall(r'\{\{[A-Z_]+\}\}', page)
    if left:
        print('WARNING unreplaced placeholders:', sorted(set(left)), file=sys.stderr)
    out = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else os.path.join(ROOT, 'index.html')
    open(out, 'w', encoding='utf-8').write(page)
    print('wrote', out, f'{os.path.getsize(out) / 1e6:.2f} MB;', ', '.join(f'{k}:{len(jdump(v)) / 1e3:.0f}kB' for k, v in data.items()))


if __name__ == '__main__':
    main()
