"""Build data/build/st_frame_apex.json (static table for #frameApex) from data/notes/velocity_frames.json.
All apex RA/Dec and (U,V,W) are derived with astropy from the published (V,l,b) or (U,V,W)."""
import os
import json, sys, re, numpy as np, astropy.units as u, warnings
from astropy.coordinates import SkyCoord, Galactocentric, galactocentric_frame_defaults
from astropy.coordinates.builtin_frames.galactocentric import get_matrix_vectors
warnings.filterwarnings('ignore')
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
N = json.load(open(ROOT + '/data/notes/velocity_frames.json'))
FR = N['frames']

def C(key, text):
    return f'<a class="cite" href="#ref-{key}">{text}</a>'

def uvw_of(V, l, b):
    L, B = np.radians(l), np.radians(b)
    return np.array([V*np.cos(B)*np.cos(L), V*np.cos(B)*np.sin(L), V*np.sin(B)])

def lb_of(v):
    V = float(np.linalg.norm(v)); return V, float(np.degrees(np.arctan2(v[1], v[0])) % 360), float(np.degrees(np.arcsin(v[2]/V)))

def sexa(l, b):
    s = SkyCoord(l=l*u.deg, b=b*u.deg, frame='galactic').icrs
    ra = s.ra.to_string(unit=u.hour, sep=('ʰ', 'ᵐ', 'ˢ'), precision=1, pad=True)
    de = s.dec.to_string(sep=('°', '′', '″'), precision=0, alwayssign=True, pad=True)
    return s.ra.deg, s.dec.deg, f'{ra} {de}'.replace('-', '−')

def f(x, n=2):
    return f'{x:.{n}f}'.replace('-', '−')

def uvwstr(v, n=2):
    return '(' + ', '.join(f(x, n) for x in v) + ')'

rows = []
def add(typ, name, V, l, b, uvw, unc, notes, refs, ndec=2):
    ra, dec, sx = sexa(l, b)
    rows.append([name, round(float(V), 2), round(float(l), 3), round(float(b), 3), round(ra, 3), round(dec, 3),
                 uvwstr(uvw, ndec), unc, ', '.join(refs), re.sub(r'</?code>', '', notes)])

def fr(k):
    d = FR[k]; return d['V'], d['l'], d['b'], np.array(d['uvw'])

# ---------------- LSR
V, l, b, w = fr('LSRK_Gordon1975')
add('LSR', 'LSRK — standard solar motion', V, l, b, w, 'exact (definition)',
    '20 km/s toward 18ʰ, +30° (B1900); FITS SPECSYS <code>LSRK</code>; astropy <code>LSRK</code>; SLALIB <code>sla_RVLSRK</code>. Reid et al. 2009 quote (10.3, 15.3, 7.7) km/s.',
    [C('gordon1976', 'Gordon 1976'), C('greisen2006', 'Greisen et al. 2006'), C('reid2009', 'Reid et al. 2009')])
V, l, b, w = fr('LSRD_Delhaye1965')
add('LSR', 'LSRD — Delhaye 1965', V, l, b, w, 'exact (definition)',
    '(U, V, W) = (9, 12, 7) km/s; FITS <code>LSRD</code> (16.6 km/s); astropy <code>LSRD</code>; SLALIB <code>sla_RVLSRD</code>.',
    [C('delhaye1965', 'Delhaye 1965'), C('greisen2006', 'Greisen et al. 2006')])
V, l, b, w = fr('LSR_Schoenrich2010')
add('LSR', 'LSR — Schönrich et al. 2010', V, l, b, w,
    '(+0.69/−0.75, ±0.47, +0.37/−0.36) stat; (1, 2, 0.5) sys',
    'astropy <code>LSR</code> / <code>GalacticLSR</code> default; Planck 2018 I Table 3 quotes 17.9 ± 2.0 km/s toward (48° ± 7°, 23° ± 4°).',
    [C('schoenrich2010', 'Schönrich et al. 2010')])
V, l, b, w = fr('LSR_DehnenBinney1998')
add('LSR', 'LSR — Dehnen & Binney 1998', V, l, b, w, '±(0.36, 0.62, 0.38)',
    'Published V⊙ = 5.25 (arXiv v1 abstract: 5.23).', [C('dehnen1998', 'Dehnen &amp; Binney 1998')])
V, l, b, w = fr('LSR_BHG2016')
add('LSR', 'LSR — Bland-Hawthorn & Gerhard 2016', V, l, b, w, '±(1, 2, 0.5)',
    'Review consensus values.', [C('blandhawthorn2016', 'Bland-Hawthorn &amp; Gerhard 2016')], ndec=2)
V, l, b, w = fr('LSR_Ding2019_GaiaDR2')
add('LSR', 'LSR — Ding et al. 2019 (Gaia DR2)', V, l, b, w, '±(0.64, 0.49, 0.36)',
    'Gaia DR2 kinematics.', [C('ding2019', 'Ding et al. 2019')])
L = N['literature']['Sun_LSR_Planck2018']
add('LSR', 'Sun–LSR, Planck 2018 Table 3', L['V'], L['l'], L['b'], uvw_of(L['V'], L['l'], L['b']), '±2.0 km/s; ±7°, ±4°',
    'Rounded from Schönrich et al. 2010 (Planck 2018 I Table 3).', [C('planck2020', 'Planck Collab. 2020'), C('schoenrich2010', 'Schönrich et al. 2010')])

# ---------------- GSR
V, l, b, w = fr('GSR_IAU1985_Delhaye')
add('GSR', 'GSR — IAU 1985 Θ₀ + LSRD (derived)', V, l, b, w, 'exact (definition)',
    'Θ₀ = 220 km/s, R₀ = 8.5 kpc (IAU 1985) added to Delhaye (9, 12, 7); reproduces the NED/RC3 vector.',
    [C('kerr1986', 'Kerr &amp; Lynden-Bell 1986'), C('delhaye1965', 'Delhaye 1965')])
V, l, b, w = fr('NED_GSR_RC3')
add('GSR', 'GSR — NED (RC3)', V, l, b, w, 'NED: 4 % of correction',
    '232.3 km/s toward (87.8°, +1.7°).', [C('devaucouleurs1991', 'de Vaucouleurs et al. 1991'), C('ned2026', 'NED')])
add('GSR', 'GALACTOC — LSR rotation (component)', 220.0, 90.0, 0.0, uvw_of(220, 90, 0), 'exact (definition)',
    'Conventional magnitude of FITS SPECSYS <code>GALACTOC</code>: Θ₀ = 220 km/s toward (90°, 0°); add a solar motion to obtain the Sun’s velocity.',
    [C('greisen2006', 'Greisen et al. 2006'), C('kerr1986', 'Kerr &amp; Lynden-Bell 1986')])
L = N['literature']['LSR_GC_Planck2018']
add('GSR', 'LSR–GC, Planck 2018 (component)', L['V'], L['l'], L['b'], uvw_of(L['V'], L['l'], L['b']), '±5 km/s',
    'Circular speed at the Sun from McMillan 2011 (Planck 2018 I Table 3).', [C('planck2020', 'Planck Collab. 2020'), C('mcmillan2011', 'McMillan 2011')])
# astropy Galactocentric: exact Galactic-axis vectors (galcen_v_sun is given in Galactocentric axes)
gb = SkyCoord(l=[0, 90, 0]*u.deg, b=[0, 0, 90]*u.deg, frame='galactic').icrs.cartesian.xyz.value
for par, label, refs, extra in [
        ('v4.0', 'GSR — astropy “v4.0” (default)',
         [C('astropy2022', 'Astropy Collab. 2022'), C('drimmel2018', 'Drimmel &amp; Poggio 2018'), C('gravity2018', 'GRAVITY Collab. 2018'), C('reid2004', 'Reid &amp; Brunthaler 2004'), C('bennett2019', 'Bennett &amp; Bovy 2019')],
         'galcen_v_sun = (12.9, 245.6, 7.78) km/s in Galactocentric axes; R₀ = 8.122 kpc; z⊙ = 20.8 pc; centre (266.4051°, −28.936175°). (U, V, W) shown in Galactic axes (derived).'),
        ('pre-v4.0', 'GSR — astropy “pre-v4.0”',
         [C('astropy2022', 'Astropy Collab. 2022'), C('schoenrich2010', 'Schönrich et al. 2010')],
         'galcen_v_sun = (11.1, 232.24, 7.25) km/s in Galactocentric axes; R₀ = 8.3 kpc; z⊙ = 27 pc. (U, V, W) shown in Galactic axes (derived).')]:
    with galactocentric_frame_defaults.set(par):
        gc = Galactocentric()
        vs = gc.galcen_v_sun.d_xyz.to_value(u.km/u.s)
        M, _ = get_matrix_vectors(gc)
        vg = gb.T @ (M.T @ vs)
    V, l, b = lb_of(vg)
    add('GSR', label, V, l, b, vg, 'convention (software default)', extra, refs)
V, l, b, w = fr('GSR_Reid2019_composite')
add('GSR', 'GSR — Reid et al. 2019 composite (derived)', V, l, b, w, '±4 km/s (V)',
    'Θ₀ + V⊙ = 247 ± 4 km/s (Θ₀ = 236 ± 7 km/s, R₀ = 8.15 ± 0.15 kpc) with U, W from Schönrich et al. 2010; a composite, not a published vector.',
    [C('reid2019', 'Reid et al. 2019'), C('schoenrich2010', 'Schönrich et al. 2010')])
for k, r0key, r0txt in [('GSR_SgrApm_GRAVITY2019', 'gravity2019', 'R₀ = 8.178 kpc'), ('GSR_SgrApm_GRAVITY2021', 'gravity2021', 'R₀ = 8.275 kpc')]:
    V, l, b, w = fr(k); d = FR[k]
    add('GSR', f'GSR — Sgr A* reflex, {r0txt} (derived)', V, l, b, w, f'±{d["sVy"]:.2f} km/s (V)',
        f'V = 4.74047 |μ<sub>l</sub>| R₀, W = 4.74047 |μ<sub>b</sub>| R₀ with μ = (−6.411, −0.219) mas/yr; U = 11.1 km/s from Schönrich et al. 2010; {r0txt}.',
        [C('reid2020', 'Reid &amp; Brunthaler 2020'), C(r0key, 'GRAVITY Collab. ' + r0key[-4:]), C('schoenrich2010', 'Schönrich et al. 2010')])

# ---------------- Local Group
V, l, b, w = fr('LG_RC2_300sinl')
add('LG', 'LG — RC2 / “IAU” 300 sin l cos b', V, l, b, np.array([0.0, 300.0, 0.0]), 'exact (definition)',
    'Δv = 300 sin l cos b; FITS SPECSYS <code>LOCALGRP</code> uses the same 300 km/s toward (90°, 0°).',
    [C('devaucouleurs1976', 'de Vaucouleurs et al. 1976'), C('mould2000', 'Mould et al. 2000'), C('greisen2006', 'Greisen et al. 2006')])
V, l, b, w = fr('LG_YTS1977')
add('LG', 'LG — Yahil, Tammann & Sandage 1977', V, l, b, w, 'exact (formula)',
    'Δv = −79 cos l cos b + 296 sin l cos b − 36 sin b (Mould et al. 2000 eq. A1); apex derived.',
    [C('yahil1977', 'Yahil et al. 1977'), C('mould2000', 'Mould et al. 2000')], ndec=0)
for k, unc, note, refs in [
        ('LG_KM1996', '±5 km/s; ±2°, ±2°', 'NED’s Local Group frame (NED error model 6 %).', [C('karachentsev1996', 'Karachentsev &amp; Makarov 1996'), C('ned2026', 'NED')]),
        ('LG_CvdB1999', '±18 km/s; ±5°, ±4°', '', [C('courteau1999', 'Courteau &amp; van den Bergh 1999')]),
        ('LG_Diaz2014', '±15 km/s; ±3.6°, ±3.0°', 'Adopted in Planck 2018 I Table 3.', [C('diaz2014', 'Diaz et al. 2014'), C('planck2020', 'Planck Collab. 2020')])]:
    V, l, b, w = fr(k)
    name = {'LG_KM1996': 'LG — Karachentsev & Makarov 1996 (NED)', 'LG_CvdB1999': 'LG — Courteau & van den Bergh 1999', 'LG_Diaz2014': 'LG — Diaz et al. 2014'}[k]
    add('LG', name, V, l, b, w, unc, note, refs)

# ---------------- CMB
V, l, b, w = fr('CMB_Planck2018')
add('CMB', 'CMB — Planck 2018', V, l, b, w, '±0.11 km/s; ±0.011°, ±0.005°',
    'β = (1.23357 ± 0.00036)×10⁻³; dipole amplitude 3362.08 ± 0.99 μK; published (α, δ) = (167.942°, −6.944°).',
    [C('planck2020', 'Planck Collab. 2020')])
V, l, b, w = fr('CMB_Fixsen1996_NED')
add('CMB', 'CMB — COBE/FIRAS 1996, NED “3K”', V, l, b, w, '±0.30°, ±0.30° (95 % CL)',
    '3.372 ± 0.007 mK, T₀ = 2.728 ± 0.004 K ⇒ cΔT/T₀ = 370.56 km/s (derived); NED adopts 371.0 km/s (error model 7 %).',
    [C('fixsen1996', 'Fixsen et al. 1996'), C('ned2026', 'NED')])
L = N['literature']['CMB_Kogut1993_dipole']
Vk = N['kogut1993_v_from_dT_with_FIRAS_T0_kms']
add('CMB', 'CMB — COBE/DMR 1993', Vk, L['l'], L['b'], uvw_of(Vk, L['l'], L['b']), '3.365 ± 0.027 mK; ±0.3°, ±0.5°',
    '|V| from ΔT with T₀ = 2.728 K (derived).', [C('kogut1993', 'Kogut et al. 1993')])
add('CMB', 'CMB — FITS CMBDIPOL', 368.0, 263.85, 48.25, uvw_of(368, 263.85, 48.25), 'exact (definition)',
    'Conventional magnitude of FITS SPECSYS <code>CMBDIPOL</code>, from WMAP (Bennett et al. 2003).',
    [C('greisen2006', 'Greisen et al. 2006'), C('bennett2003', 'Bennett et al. 2003')])
for k, name, unc, note, refs in [
        ('LG_CMB_Kogut1993', 'LG–CMB, Kogut et al. 1993 (component)', '±22 km/s; ±3°, ±3°', 'Velocity of the LG, not of the Sun.', [C('kogut1993', 'Kogut et al. 1993')]),
        ('LG_CMB_Planck2018', 'LG–CMB, Planck 2018 (component)', '±15 km/s; ±2.0°, ±1.4°', 'Non-relativistic vector difference of the Planck dipole and Diaz et al. 2014 Sun–LG.', [C('planck2020', 'Planck Collab. 2020'), C('diaz2014', 'Diaz et al. 2014')]),
        ('GC_CMB_Planck2018', 'GC–CMB, Planck 2018 (component)', '±5 km/s; ±0.20°, ±0.28°', 'From Schönrich et al. 2010 + McMillan 2011 (Planck 2018 I Table 3).', [C('planck2020', 'Planck Collab. 2020')])]:
    V, l, b, w = fr(k)
    add('CMB', name, V, l, b, w, unc, note, refs)

cols = [
    {"k": "frame", "label": "Frame / definition", "num": False, "html": False, "title": "Sun's (barycentre's) velocity relative to the frame unless marked 'component'"},
    {"k": "v", "label": "|V| (km/s)", "num": True, "html": False, "title": "speed of the Sun relative to the frame"},
    {"k": "l", "label": "l<sub>apex</sub> (°)", "num": True, "html": False, "title": "Galactic longitude of the apex"},
    {"k": "b", "label": "b<sub>apex</sub> (°)", "num": True, "html": False, "title": "Galactic latitude of the apex"},
    {"k": "ra", "label": "α<sub>apex</sub> (°)", "num": True, "html": False, "title": "ICRS (≈ J2000) right ascension of the apex; derived with astropy unless published"},
    {"k": "dec", "label": "δ<sub>apex</sub> (°)", "num": True, "html": False, "title": "ICRS (≈ J2000) declination of the apex; derived with astropy unless published"},
    {"k": "uvw", "label": "(U, V, W) (km/s)", "num": False, "html": False, "title": "Galactic Cartesian components: U toward l = 0°, V toward l = 90°, W toward b = +90° (derived where the source gives (V, l, b))"},
    {"k": "unc", "label": "Uncertainty", "num": False, "html": False, "title": "as published (1σ unless stated)"},
    {"k": "ref", "label": "Reference", "num": False, "html": True, "title": ""},
    {"k": "notes", "label": "Convention / notes", "num": False, "html": True, "title": ""},
]
caption = ('Apex of the Sun’s (solar-system barycentre’s) velocity relative to each frame; rows marked “component” give the motion of one frame relative to the next. '
           'The correction added to a barycentric radial velocity is Δv = |V| cos θ, θ being the angle between target and apex (eq. 7.6.2). '
           'Apex RA/Dec and (U, V, W) are derived with astropy 7.2 from the published (|V|, l, b) or (U, V, W) except where the source gives them; “derived” rows are combinations made for this page, not published vectors. '
           'Values from <code>data/notes/velocity_frames.json</code>.')
out = {"host": "#frameApex", "caption": caption, "tall": False, "search": True, "cols": cols, "rows": rows}
json.dump(out, open(ROOT + '/data/build/st_frame_apex.json', 'w'), ensure_ascii=False, indent=1)
for r in rows: print(r[0].ljust(42), r[1:7], r[7])

# ---- concise notes (full detail lives in the §7.6/§7.7 text); rewrite the notes column and trim long reference lists
SHORT = {
 'LSRK — standard solar motion': '20 km/s toward 18ʰ, +30° (B1900); FITS, astropy and SLALIB LSRK',
 'LSRD — Delhaye 1965': '(9, 12, 7) km/s; FITS LSRD (16.6 km/s); astropy, SLALIB LSRD',
 'LSR — Schönrich et al. 2010': 'astropy LSR default',
 'LSR — Dehnen & Binney 1998': 'Hipparcos; arXiv v1 abstract gives V = 5.23',
 'LSR — Bland-Hawthorn & Gerhard 2016': 'review-adopted values',
 'LSR — Ding et al. 2019 (Gaia DR2)': 'Gaia DR2',
 'Sun–LSR, Planck 2018 Table 3': 'rounded from Schönrich et al. 2010',
 'GSR — IAU 1985 Θ₀ + LSRD (derived)': 'Θ₀ = 220 km/s, R₀ = 8.5 kpc, plus (9, 12, 7); equals the NED vector',
 'GSR — NED (RC3)': 'NED “GSR”; error model 4 % of the correction',
 'GALACTOC — LSR rotation (component)': 'FITS GALACTOC magnitude; add a solar motion',
 'LSR–GC, Planck 2018 (component)': 'circular speed from McMillan 2011',
 'GSR — astropy “v4.0” (default)': 'R₀ = 8.122 kpc, z⊙ = 20.8 pc; galcen_v_sun = (12.9, 245.6, 7.78)',
 'GSR — astropy “pre-v4.0”': 'R₀ = 8.3 kpc, z⊙ = 27 pc; galcen_v_sun = (11.1, 232.24, 7.25)',
 'GSR — Reid et al. 2019 composite (derived)': 'Θ₀ + V⊙ = 247 ± 4 km/s; U, W from Schönrich et al. 2010',
 'GSR — Sgr A* reflex, R₀ = 8.178 kpc (derived)': 'V, W = 4.74047 |μ| R₀; μ = (−6.411, −0.219) mas/yr',
 'GSR — Sgr A* reflex, R₀ = 8.275 kpc (derived)': 'V, W = 4.74047 |μ| R₀; μ = (−6.411, −0.219) mas/yr',
 'LG — RC2 / “IAU” 300 sin l cos b': 'Δv = 300 sin l cos b; FITS LOCALGRP',
 'LG — Yahil, Tammann & Sandage 1977': 'Δv = −79 cos l cos b + 296 sin l cos b − 36 sin b',
 'LG — Karachentsev & Makarov 1996 (NED)': 'NED Local Group frame; error model 6 %',
 'LG — Courteau & van den Bergh 1999': 'measured apex',
 'LG — Diaz et al. 2014': 'used in Planck 2018 Table 3',
 'CMB — Planck 2018': 'β = (1.23357 ± 0.00036)×10⁻³; 3362.08 ± 0.99 μK',
 'CMB — COBE/FIRAS 1996, NED “3K”': '3.372 mK, cΔT/T₀ = 370.56 km/s; NED uses 371.0',
 'CMB — COBE/DMR 1993': '|V| from 3.365 mK with T₀ = 2.728 K (derived)',
 'CMB — FITS CMBDIPOL': 'FITS conventional value, after WMAP',
 'LG–CMB, Kogut et al. 1993 (component)': 'motion of the LG, not of the Sun',
 'LG–CMB, Planck 2018 (component)': 'vector difference: Planck dipole − Diaz et al. Sun–LG',
 'GC–CMB, Planck 2018 (component)': 'from Schönrich et al. 2010 + McMillan 2011',
}
for r in rows:
    r[-1] = SHORT[r[0]]
    if r[0].startswith('GSR — astropy “v4.0”'):
        r[-2] = ', '.join([C('astropy2022', 'Astropy Collab. 2022'), C('drimmel2018', 'Drimmel &amp; Poggio 2018'), C('gravity2018', 'GRAVITY Collab. 2018')])
out['rows'] = rows
out['cols'][-1]['html'] = False
json.dump(out, open(ROOT + '/data/build/st_frame_apex.json', 'w'), ensure_ascii=False, indent=1)
print('rewrote with short notes;', len(rows), 'rows')

# ---- merge (l, b) and (α, δ) into single columns to keep the table narrow enough to read
def pm(x, n=3):
    return (f'{x:+.{n}f}').replace('-', '−')
new_rows = []
for r in rows:
    name, V, l, b, ra, dec, uvw, unc, ref, notes = r
    new_rows.append([name, V, f'{l:.3f}, {pm(b)}', f'{ra:.3f}, {pm(dec)}', uvw, unc, ref, notes])
out['rows'] = new_rows
out['cols'] = [
    {"k": "frame", "label": "Frame / definition", "num": False, "html": False, "title": "Sun's (barycentre's) velocity relative to the frame unless marked 'component'"},
    {"k": "v", "label": "|V| (km/s)", "num": True, "html": False, "title": "speed of the Sun relative to the frame"},
    {"k": "lb", "label": "Apex l, b (°)", "num": False, "html": False, "title": "Galactic longitude and latitude of the apex"},
    {"k": "radec", "label": "Apex α, δ (°, J2000)", "num": False, "html": False, "title": "ICRS (≈ J2000) apex; derived with astropy unless published"},
    {"k": "uvw", "label": "(U, V, W) (km/s)", "num": False, "html": False, "title": "Galactic Cartesian components: U toward l = 0°, V toward l = 90°, W toward b = +90° (derived where the source gives |V|, l, b)"},
    {"k": "unc", "label": "Uncertainty", "num": False, "html": False, "title": "as published (1σ unless stated)"},
    {"k": "ref", "label": "Reference", "num": False, "html": True, "title": ""},
    {"k": "notes", "label": "Convention / notes", "num": False, "html": False, "title": ""},
]
json.dump(out, open(ROOT + '/data/build/st_frame_apex.json', 'w'), ensure_ascii=False, indent=1)
print('merged coordinate columns;', len(new_rows), 'rows,', len(out['cols']), 'cols')
