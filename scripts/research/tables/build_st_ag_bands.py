#!/usr/bin/env python
"""Build data/build/st_ag_bands.json (static table, host #agBands) from data/curated/airglow_bands.csv.

usage: python -I build_st_ag_bands.py <airglow_bands.csv> <out.json>
"""
import csv, html, json, math, re, sys

REFS = {  # key in the curated CSV -> (page cite key, link text)
    'HH79-WebBook': ('huber1979', 'Huber &amp; Herzberg 1979'),
    'Barthelemy25': ('barthelemy2025', 'Barthélemy et al. 2025'),
    'Gattinger10': ('gattinger2010', 'Gattinger et al. 2010'),
    'Krupenie72': ('krupenie1972', 'Krupenie 1972'),
    'HITRAN2020': ('gordon2022', 'Gordon et al. 2022'),
    'Whiter23': ('whiter2023', 'Whiter et al. 2023'),
    'Hunten67': ('hunten1967', 'Hunten 1967'),
    'Noxon68': ('noxon1968', 'Noxon 1968'),
    'Meier91': ('meier1991', 'Meier 1991'),
    'Evans10': ('evans2010', 'Evans et al. 2010'),
    'PALACE': ('noll2025', 'Noll et al. 2025'),
    'LK77': ('lofthus1977', 'Lofthus &amp; Krupenie 1977'),
}
KEYRE = '|'.join(re.escape(k) for k in sorted(REFS, key=len, reverse=True))


def link(k):
    key, txt = REFS[k]
    return f'<a class="cite" href="#ref-{key}">{txt}</a>'


def uni(s):
    """plain-text tidy-up (before HTML escaping)"""
    s = s.replace('->', '→').replace('<=', '≤').replace('>=', '≥').replace('>~', '≳ ')
    s = re.sub(r'\bH-alpha\b', 'Hα', s)
    s = re.sub(r'cm-1\b', 'cm⁻¹', s)
    s = re.sub(r'cm-2\b', 'cm⁻²', s)
    s = re.sub(r'\bs-1\b', 's⁻¹', s)
    s = re.sub(r'(?<=[\d)]) A\b', ' Å', s)
    s = s.replace('several A', 'several Å')
    s = re.sub(r'(?<=\d) ?um\b', ' µm', s)
    s = re.sub(r'(?<=\d)-(?=\d)', '–', s)
    s = s.replace('lambda_H', 'λ_H').replace('lambda_0', 'λ₀').replace('sigma_0', 'σ₀')
    s = s.replace('16OH', '¹⁶OH').replace('O(1D)', 'O(¹D)').replace('N(4S)', 'N(⁴S)').replace('O(3P)', 'O(³P)')
    s = s.replace(' + hv', ' + hν').replace('v≤', 'v ≤ ').replace('-several', '–several').replace('@1 keV', 'at 1 keV')
    s = s.replace('(PALACE R1)', '[PALACE reaction R1]').replace('NO-gamma', 'NO γ').replace('C 3Delta_u', 'C3Delta_u')
    s = re.sub(r'N2\+', 'N₂⁺', s)
    s = re.sub(r'O2\+', 'O₂⁺', s)
    s = re.sub(r'\bHO2\b', 'HO₂', s)
    s = re.sub(r'\bNO2', 'NO₂', s)
    s = re.sub(r'(?<![A-Za-z0-9])O3\b', 'O₃', s)
    s = re.sub(r'(?<![A-Za-z0-9])O2', 'O₂', s)
    s = re.sub(r'(?<![A-Za-z0-9])N2(?![0-9])', 'N₂', s)
    s = re.sub(r"\b([abcAXd]'?)-([XabA])\b", r'\1–\2', s)
    s = s.replace('T.III', 'Table III').replace('T.II', 'Table II').replace('Fig.11', 'Fig. 11').replace('p.', 'p. ')
    return s


TERM = re.compile(r"(?<![A-Za-z])([A-Za-z]'?)(\d)(Sigma|Pi|Delta)(?:_([gui]))?([+-])?")
GREEK = {'Sigma': 'Σ', 'Pi': 'Π', 'Delta': 'Δ'}


def terms(s):
    def rep(m):
        out = f'{m.group(1).replace(chr(39), "′")}<sup>{m.group(2)}</sup>{GREEK[m.group(3)]}'
        if m.group(4):
            out += f'<sub>{m.group(4)}</sub>'
        if m.group(5):
            out += f'<sup>{"+" if m.group(5) == "+" else "−"}</sup>'
        return out
    return TERM.sub(rep, s)


def linkify(s):
    s = re.sub(r'\[(' + KEYRE + r')(?: ([^\]]*))?\]', lambda m: '(' + link(m.group(1)) + (', ' + m.group(2) if m.group(2) else '') + ')', s)
    s = re.sub(r'(?<![#\w-])(' + KEYRE + r')(?![\w-])', lambda m: link(m.group(1)), s)
    return s


def rich(s):
    return linkify(terms(html.escape(uni(s), quote=False))).replace('λ_H', 'λ<sub>H</sub>')


def transition(s):
    s = s.replace("X2Pi (v') - X2Pi (v'')", 'X2Pi (v′) - X2Pi (v″)')
    s = s.replace("X2A'' (200-000), A2A'(001)-X2A''(000)", 'X²A″ (200–000); A²A′ (001) – X²A″ (000)')
    s = uni(s).replace(' - ', ' – ').replace('NO + O → NO₂* (chemiluminescence)', 'NO + O → NO₂* (chemiluminescence)')
    return terms(html.escape(s, quote=False))


MOL = {'N2+': 'N₂⁺', 'N2': 'N₂', 'O2': 'O₂', 'NO': 'NO', 'OH': 'OH', 'FeO': 'FeO', 'NO2': 'NO₂', 'HO2': 'HO₂'}


def num(s):
    s2 = s.strip().lstrip('~').strip('()')
    return float(s2) if re.fullmatch(r'\d+(\.\d+)?', s2) else None


def origin(r):
    """band-origin wavelength from the tabulated sigma_0, rounded to the precision of sigma (derived)"""
    if r['sigma_kind'] != 'origin' or not r['sigma_cm']:
        return None
    w = r['origin_wl_air_A'] if (r['medium'] == 'air' and r['origin_wl_air_A']) else r['origin_wl_vac_A']
    if not w:
        return None
    w = float(w)
    sig = r['sigma_cm']
    k = len(sig.split('.')[1]) if '.' in sig else 0
    dl = w * w * 10.0 ** (-k) / 1e8          # one unit of the last digit of sigma, in Angstrom
    dec = min(3, max(0, math.ceil(-math.log10(dl))))
    return round(w, dec)


def box(h, em):
    """give long-text cells a minimum width so rows stay short in a horizontally scrolling table"""
    return f'<span style="display:inline-block;min-width:{em}em">{h}</span>' if h else ''


def prec(p):
    p = p.replace('(LK77: measurement uncertainty about 0.5 A)', '(measurement uncertainty)')
    p = p.replace('(McEwen calc. vs obs. within 0.1 A)', '(McEwen calculated vs observed)')
    p = p.replace('(from fit; definition-dependent)', '(fit; definition-dependent)')
    p = p.replace('(15 A resolution spectra)', '(15 Å-resolution spectra)')
    p = p.replace('from its own sigma', 'from its own σ')
    p = re.sub(r'^~', '≈ ', p)
    return rich(p)


def main(src, out):
    rows = []
    with open(src, newline='', encoding='utf-8') as f:
        rd = csv.DictReader(l for l in f if not l.startswith('#'))
        for r in rd:
            wl = num(r['wl_A'])
            kind = uni(r['wl_kind'])
            if not re.fullmatch(r'\d+(\.\d+)?', r['wl_A'].strip()):
                kind += f'; printed: {uni(r["wl_A"])}' + (' Å' if any(c.isdigit() for c in r['wl_A']) else '')
            note = r['note'].replace(' A band origin is a model quantity: no line lies at it; individual lines are listed in HITRAN / sky atlases.', '').strip()
            srcs = r['wl_source']
            if srcs.startswith('derived:'):
                srcs_h = ('derived here: ²Π Hamiltonian fit (T, A, B, D, γ) per v to ' + link('HITRAN2020') +
                          ' ¹⁶OH X-state level energies (J ≤ 12.5); ν₀ = T(v′) − T(v″); cross-check ' + link('HH79-WebBook'))
            else:
                head, _, tail = srcs.partition(': ')
                extra = ''
                m = re.search(r'band origins \(b\) from Wilkinson / Miller \(unc\. 0\.003 cm-1 for Miller\)', tail)
                if m:
                    extra = '; origins: Wilkinson / Miller (unc. 0.003 cm-1 for Miller)'
                elif 'Babcock & Herzberg' in tail:
                    extra = '; origins: Babcock & Herzberg, refit Albritton et al.'
                elif 'Herzberg & Herzberg' in tail:
                    extra = '; origins: Herzberg & Herzberg; 0-1 is a twilight head'
                elif 'Broida & Gaydon' in tail:
                    extra = '; emission heads (Broida & Gaydon; Barth & Kaplan)'
                elif 'Benesch' in tail:
                    extra = '; origins: Benesch et al.'
                srcs_h = rich(head + extra)
            rows.append([
                MOL.get(r['molecule'], r['molecule']),
                {'gamma': 'γ', 'delta': 'δ'}.get(r['system'], uni(r['system']).replace('dv=', 'Δv = ')),
                transition(r['transition']),
                uni(r['band']),
                wl,
                kind,
                r['medium'],
                box(prec(r['precision_A']), 9),
                origin(r),
                r['category'].replace('|', ', '),
                box(rich(r["mechanism"]), 20),
                box(rich(r["altitude_km"]), 14),
                box(rich(r["typical_R"]), 22),
                box(srcs_h, 16),
                box(rich(note), 18),
            ])
    cols = [
        dict(k='mol', label='Molecule', num=False, html=False),
        dict(k='sys', label='System', num=False, html=False),
        dict(k='tr', label='Transition', num=False, html=True, title='upper – lower electronic state'),
        dict(k='band', label='v′–v″', num=False, html=False, title='vibrational band (upper–lower)'),
        dict(k='wl', label='λ (Å)', num=True, html=False,
             title='head or origin wavelength as tabulated by the source (air above 2000 Å, vacuum below; see Medium). Parenthesised/approximate values are marked under Feature'),
        dict(k='kind', label='Feature', num=False, html=False, title='which feature the wavelength refers to: band head (branch and degradation), band origin, nominal position or range'),
        dict(k='med', label='Medium', num=False, html=False),
        dict(k='prec', label='Precision (Å)', num=False, html=True, title='stated or estimated precision of the tabulated wavelength'),
        dict(k='org', label='Origin λ₀ (Å)', num=True, html=False,
             title='derived: band origin from the tabulated σ₀ (λ = 10⁸/σ₀, then Peck & Reeder 1972 to air when the row is in air), rounded to the precision of σ₀ (max. 3 decimals)'),
        dict(k='cat', label='Category', num=False, html=False),
        dict(k='mech', label='Excitation', num=False, html=True),
        dict(k='alt', label='Altitude (km)', num=False, html=True),
        dict(k='I', label='Typical intensity', num=False, html=True, title='with viewing condition: zenith, nadir from 600 km, Paranal climatology, or auroral model per erg cm⁻² s⁻¹'),
        dict(k='src', label='Position source', num=False, html=True),
        dict(k='note', label='Note', num=False, html=True),
    ]
    caption = ('Molecular airglow and auroral band systems (171 rows). Heads and origins are <i>different</i> features: a head is the '
               'turning point of one branch, the origin ν₀ = T(v′) − T(v″) lies between branches and coincides with no line. '
               'N₂ and N₂⁺: <a class="cite" href="#ref-lofthus1977">Lofthus &amp; Krupenie 1977</a>; O₂: '
               '<a class="cite" href="#ref-krupenie1972">Krupenie 1972</a> (both read from the table images of the JPCRD reprints); '
               'NO: nominal positions from <a class="cite" href="#ref-meier1991">Meier 1991</a>; OH Meinel origins derived here from '
               '<a class="cite" href="#ref-gordon2022">HITRAN2020</a> level energies; continua after '
               '<a class="cite" href="#ref-noll2025">Noll et al. 2025</a>. “Origin λ₀” is derived from the tabulated wavenumber. '
               'Line-resolved OH and O₂ positions are in §5.3–5.4.')
    T = {'host': '#agBands', 'caption': caption, 'tall': True, 'search': True, 'cols': cols, 'rows': rows}
    json.dump(T, open(out, 'w', encoding='utf-8'), ensure_ascii=False, indent=0)
    print('wrote', out, len(rows), 'rows')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
