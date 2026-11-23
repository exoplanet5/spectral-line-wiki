#!/usr/bin/env python
"""Parse the NIST Handbook of Basic Atomic Spectroscopic Data (Sansonetti & Martin 2005,
J. Phys. Chem. Ref. Data 34, 1559) HTML e-book into machine-readable JSON.

Input : ~/Desktop/ebookbasd/*.htm
Output: data/handbook.json  {elements:[...], lines:[...], persistent:[...], unstaged:[...], parse_notes:{...}}

The wavelength strings are kept verbatim (wl_str) so that the published number of
decimals -- the only uncertainty indicator in the Handbook -- is preserved.

Strong lines (table2):
  * lines printed without an ionization stage (Os 5 lines MCS75, Po 1 CHPT55, Es 22 WLGC74) cannot be
    matched to an ASD spectrum; they are NOT put in `lines` but listed in `unstaged` (ion '').
  * lines printed twice in the e-book (identical element, stage, wavelength string, intensity, flags,
    medium and reference; 12 Pt lines) are kept once; the duplicates are listed in parse_notes.

Ground states (table1), per spectrum in elements[].spectra[ion]:
  ground  original text with the markup stripped, as before ("1s22s22p63s23p63d64s2 5D4"; for Pa II,
          which has no ionization-energy line, taken from the ground-state row itself)
  conf    subshells separated by '.', occupancy taken from the e-book's explicit <Sup> markup
          (none printed = 1, so it is omitted, as in ASD: "3d6.4s"); a parent term in parentheses is
          attached to the subshell it belongs to ("4f13(2F°7/2).6s"); the e-book's grouping
          parentheses around the closed core are dropped.  Token grammar:
          (\\d+)([spdfghik])(\\d*)(\\(parent term\\))?
          The e-book's j subscripts of jj-coupled subshells (Pb I 6p2 with subscript 1/2) are dropped
          here; the term carries the jj designation.
  term    the level designation as printed ("5D4", "2P°3/2", "10D°5/2", "(7/2,1/2)°4"; Th II "3/2"),
          corrected only where GROUND_FIX says so (missing parity marks; Hf II)
  core    noble-gas core abbreviation ("[Ar]", "" if none) and
  val     the remaining subshells in printed order ("3d6.4s2"), so the page can show [Ar]3d⁶4s²
  note    only where the e-book is internally inconsistent or disagrees with NIST ASD (GROUND_FIX)
"""
import glob, html, json, os, re, sys
from collections import Counter

SRC = os.path.expanduser('~/Desktop/ebookbasd')
OUT = os.path.join(os.path.dirname(__file__), '..', 'data', 'handbook.json')


def txt(s):
    return html.unescape(re.sub(r'<[^>]+>', '', s))


# ------------------------------------------------------------------------------------------ ground states
CORES = [('[Rn]', '1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6 4f14 5d10 6s2 6p6'),
         ('[Xe]', '1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6'),
         ('[Kr]', '1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6'), ('[Ar]', '1s2 2s2 2p6 3s2 3p6'),
         ('[Ne]', '1s2 2s2 2p6'), ('[He]', '1s2')]
GTOK = re.compile(r'\s*(?:(?P<sub>(?P<n>\d+)<i>(?P<l>[a-z])</i>(?P<occ>(?:<sup>\d+</sup>)*)(?:<sub>[^<]*</sub>)*)'
                  r'|(?P<par>\((?=<sup>|\d+[A-Z])(?P<pt>.*?)\))|\(|\))', re.I)
# Printed ground designations that are internally inconsistent (electron count, parity of the configuration)
# or that NIST ASD designates differently (checked against the ASD lines whose lower level has E = 0).
# value: ({printed subshell token: corrected token} or None, corrected term or None, note)
GROUND_FIX = {
    ('Dy', 'II'): ({'4f9(5I8)': '4f10(5I8)'}, None, 'e-book prints 4f9(5I8)6s (64 electrons); corrected to 4f10(5I8)6s as in NIST ASD'),
    ('Hf', 'II'): (None, '2D3/2', 'e-book prints 2P°3/2, impossible for the even 5d6s2 configuration; NIST ASD: 5d6s2 a 2D3/2'),
    ('Al', 'I'): (None, '2P°1/2', 'e-book omits the odd-parity mark (3p); NIST ASD: 2P°1/2'),
    ('Tl', 'I'): (None, '2P°1/2', 'e-book omits the odd-parity mark (6p); NIST ASD: 2P°1/2'),
    ('I', 'I'): (None, '2P°3/2', 'e-book omits the odd-parity mark (5p5); NIST ASD: 2P°3/2'),
    ('Bi', 'I'): (None, '4S°3/2', 'e-book omits the odd-parity mark (6p3); NIST ASD: 4S°3/2'),
    ('Th', 'II'): (None, None, 'e-book gives 6d2 7s with J = 3/2 but no term; NIST ASD designates the ground level 6d7s2 2D3/2'),
}


def parse_ground_conf(raw):
    """e-book HTML configuration -> (list of subshell tokens, electron count)."""
    s = re.sub(r'<(/?)(sup|sub|i)>', lambda m: f'<{m.group(1)}{m.group(2).lower()}>', raw.strip(), flags=re.I)
    pos, toks, nel = 0, [], 0
    while pos < len(s):
        m = GTOK.match(s, pos)
        if not m or m.end() == pos:
            raise ValueError(f'unparsed ground configuration at {s[pos:]!r}')
        pos = m.end()
        if m.group('sub'):
            occ = ''.join(re.findall(r'<sup>(\d+)</sup>', m.group('occ')))
            toks.append(f"{m.group('n')}{m.group('l')}{occ}")
            nel += int(occ) if occ else 1
        elif m.group('par'):
            toks[-1] += '(' + txt(m.group('pt')).replace(' ', '') + ')'
    return toks, nel


def split_core(toks):
    have = set(toks)
    for k, v in CORES:
        cs = v.split()
        if all(c in have for c in cs) and len(toks) > len(cs):
            return k, [t for t in toks if t not in cs]
    return '', toks


def parse_general(name):
    s = open(f'{SRC}/{name}table1.htm', encoding='latin-1').read()
    title = re.search(r'Atomic Data for\s+(.*?)\s*\((.*?)\)', txt(s))
    Z = int(re.search(r'Atomic Number\s*=\s*(\d+)', s).group(1))
    w = re.search(r'Atomic Weight\s*=\s*([^<]+)', s)
    sym = title.group(2).strip()
    ies = {}
    body = txt(s)
    for m in re.finditer(r'\b' + re.escape(sym) + r'\s+(I|II)\s+Ground State\s+(.*?)\n\s*Ionization energy\s+([\d.()\[\]]+)\s*cm-1\s*\(([\d.()\[\]]+)\s*eV\)', body, re.S):
        ies[m.group(1)] = {'ground': ' '.join(m.group(2).split()), 'ie_cm': m.group(3), 'ie_eV': m.group(4)}
    # explicit superscript markup of the same rows
    for m in re.finditer(r'^' + re.escape(sym) + r'\s+(I|II)\s+Ground State\s+(.*)$', s, re.M):
        ion = m.group(1)
        parts = re.split(r'\s{2,}', m.group(2).strip())
        toks, nel = parse_ground_conf(' '.join(parts[:-1]))
        term = txt(parts[-1]).strip()
        fix = GROUND_FIX.get((sym, ion))
        if fix and fix[0]:
            for a, b in fix[0].items():
                assert a in toks, (sym, ion, a)
                toks = [b if t == a else t for t in toks]
                nel += int(re.match(r'\d+[a-z](\d*)', b).group(1) or 1) - int(re.match(r'\d+[a-z](\d*)', a).group(1) or 1)
        if fix and fix[1]:
            term = fix[1]
        core, val = split_core(toks)
        if nel != Z - (0 if ion == 'I' else 1):
            print('WARN ground electron count', sym, ion, nel, file=sys.stderr)
        rec = ies.setdefault(ion, {})
        if not rec.get('ground'):          # row without an ionization-energy line (Pa II): original text from the row
            rec['ground'] = ' '.join(txt(m.group(2)).split())
        rec.update({'conf': '.'.join(toks), 'term': term, 'core': core, 'val': '.'.join(val)})
        if fix:
            rec['note'] = fix[2]
    return {'name': title.group(1).strip(), 'symbol': sym, 'Z': Z,
            'weight': w.group(1).strip() if w else None, 'spectra': ies, 'file': name}


# ------------------------------------------------------------------------------------------ strong lines
LINE_RE = re.compile(r'^\s*(\d+)?(\*)?\s*([A-Za-z*]+(?:,[A-Za-z*]+)*)?\s+(\d+\.\d*|\d+)\s+([A-Z][a-z]?)(?:\s+(I|II))?\s+(\S+)?\s*$')


def parse_strong(name, sym):
    s = open(f'{SRC}/{name}table2.htm', encoding='latin-1').read()
    body = s.split('</h3>', 1)[1]
    medium, out = None, []
    for raw in body.split('\n'):
        t = txt(raw).rstrip()
        st = t.strip()
        if st in ('Air', 'Vacuum'):
            medium = st.lower(); continue
        # intensity, optional flags, wavelength, symbol, [ion], ref
        m = LINE_RE.match(t + ' ')
        if not m:
            continue
        inten, star, flags, wl, el, ion, ref = m.groups()
        if el != sym:
            print('WARN symbol mismatch', name, t, file=sys.stderr)
        fl = (flags or '')
        if star: fl = ('*,' + fl) if fl else '*'
        dec = len(wl.split('.')[1]) if '.' in wl else 0
        out.append({'el': el, 'ion': ion or '', 'int': int(inten) if inten else None,
                    'flags': fl, 'wl_str': wl, 'wl': float(wl), 'dec': dec,
                    'medium': medium, 'ref': ref})
    return out


def parse_persistent(name, sym):
    res = []
    for tab, ion, more in ((3, 'I', 'moredata'), (4, 'II', 'moredata2')):
        f = f'{SRC}/{name}table{tab}.htm'
        if not os.path.exists(f):
            continue
        s = open(f, encoding='latin-1').read()
        hdr = re.search(r'Persistent Lines of (.*?)</h3>', s)
        if not hdr:
            continue
        ion = 'II' if 'Singly Ionized' in hdr.group(1) else 'I'
        more = 'moredata2' if ion == 'II' else 'moredata'
        # configurations
        conf = {}
        mf = f'{SRC}/{name}{more}.htm'
        if os.path.exists(mf):
            ms = open(mf, encoding='latin-1').read()
            for blk in re.split(r'<a name="', ms)[1:]:
                key = blk.split('"', 1)[0]
                b = txt(blk.split('>', 1)[1]).split('return')[0]
                lines = [l for l in b.split('\n') if l.strip()]
                conf[key] = [' '.join(l.split()) for l in lines]
        for raw in s.split('\n'):
            if not raw.startswith('<a name='):
                continue
            key = re.match(r'<a name="([^"]+)"', raw).group(1)
            if not txt(raw).strip():
                continue
            # split the raw HTML row into: head text | level anchor | reference anchors
            head = raw.split('</a>', 1)[1]
            lev = re.search(r'<a href="[^"]*more[^"]*">(.*?)</a>', head)
            refs = re.findall(r'<a href="[^"]*table7\.htm#[^"]*">(.*?)</a>', head)
            pre = txt(head[:lev.start()] if lev else re.split(r'<a ', head)[0])
            m = re.match(r'^\s*(\d+\*?)?\s*([A-Za-z*]+(?:,[A-Za-z*]+)*)?\s*(\d+\.\d*)\s*([\d.]+(?:[eE][+-]?\d+)?\*?)?\s*$', pre)
            if not m:
                print('WARN persistent unparsed', f, repr(pre), file=sys.stderr); continue
            inten, flags, wl, aki = m.groups()
            ei = ek = None
            if lev:
                lv = txt(lev.group(1))
                mm = re.match(r'^\s*(\S+)\s*-\s*(\S+)\s*$', lv)
                if mm:
                    ei, ek = mm.groups()
            if inten and inten.endswith('*'):
                inten = inten[:-1]; flags = ('*,' + flags) if flags else '*'
            res.append({'el': sym, 'ion': ion, 'int': int(inten) if inten else None, 'flags': flags or '',
                        'wl_str': wl, 'wl': float(wl), 'Aki_1e8': aki, 'Ei': ei, 'Ek': ek,
                        'line_ref': refs[0] if refs else None, 'A_ref': refs[1] if len(refs) > 1 else None,
                        'levels': conf.get(key)})
    return res


def main():
    names = sorted({os.path.basename(f)[:-len('table1.htm')] for f in glob.glob(f'{SRC}/*table1.htm')})
    elements, lines, pers, unstaged, dups = [], [], [], [], []
    n_printed = 0
    for n in names:
        g = parse_general(n)
        elements.append(g)
        seen = set()
        for l in parse_strong(n, g['symbol']):
            n_printed += 1
            l['Z'] = g['Z']
            k = (l['el'], l['ion'], l['wl_str'], l['int'], l['flags'], l['medium'], l['ref'])
            if k in seen:
                dups.append(f"{l['el']} {l['ion']} {l['wl_str']}")
                continue
            seen.add(k)
            (lines if l['ion'] else unstaged).append(l)
        pers += parse_persistent(n, g['symbol'])
    elements.sort(key=lambda e: e['Z'])
    lines.sort(key=lambda l: l['wl'])
    unstaged.sort(key=lambda l: (l['Z'], l['wl']))
    notes = {'strong_lines_printed': n_printed, 'duplicates_dropped': dups,
             'unstaged': f'{len(unstaged)} strong lines printed without an ionization stage '
                         f'({dict(Counter(l["el"] for l in unstaged))}); kept in `unstaged`, not matched to ASD',
             'ground_notes': {f'{e["symbol"]} {ion}': s['note'] for e in elements for ion, s in e['spectra'].items() if s.get('note')}}
    print('elements', len(elements), 'printed', n_printed, 'lines', len(lines), 'unstaged', len(unstaged),
          'duplicates', len(dups), 'persistent', len(pers))
    print('medium', Counter(l['medium'] for l in lines))
    print('decimals', sorted(Counter(l['dec'] for l in lines).items()))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump({'source': 'J.E. Sansonetti & W.C. Martin, Handbook of Basic Atomic Spectroscopic Data, '
                         'J. Phys. Chem. Ref. Data 34, 1559 (2005); NIST e-book edition',
               'elements': elements, 'lines': lines, 'persistent': pers, 'unstaged': unstaged, 'parse_notes': notes},
              open(OUT, 'w'), indent=0)


if __name__ == '__main__':
    main()
