#!/usr/bin/env python
"""Kinematic frame definitions (solar motion relative to each frame) for the frame calculator.
Apex directions are converted between Galactic and equatorial (ICRS) with astropy."""
import json, os, warnings
import numpy as np, astropy.units as u
from astropy.coordinates import SkyCoord, FK4
warnings.filterwarnings('ignore')
ROOT = os.path.join(os.path.dirname(__file__), '..')

def lb_of_uvw(U, V, W):
    v = np.sqrt(U*U + V*V + W*W)
    return v, np.degrees(np.arctan2(V, U)) % 360, np.degrees(np.arcsin(W / v))

def icrs(l, b):
    s = SkyCoord(l=l*u.deg, b=b*u.deg, frame='galactic').icrs
    return round(s.ra.deg, 4), round(s.dec.deg, 4)

k0 = SkyCoord(ra=270*u.deg, dec=30*u.deg, frame=FK4(equinox='B1900'))
k = k0.galactic
lk, bk = k.l.deg, k.b.deg
k_icrs = k0.icrs
F = []
def add(key, name, desc, ref, **kw):
    d = dict(key=key, name=name, desc=desc, ref=ref, **kw)
    if 'uvw' in kw:
        v, l, b = lb_of_uvw(*kw['uvw']); d.update(v=round(v, 3), apex=[round(l, 3), round(b, 3)], apex_icrs=icrs(l, b))
    elif 'apex' in kw:
        d['apex_icrs'] = icrs(*kw['apex'])
    F.append(d)

add('lsrk', 'LSRK (kinematic LSR)', 'Standard solar motion: 20 km/s toward α = 18ʰ, δ = +30° (B1900) ⇒ (l, b) = (%.3f°, %.3f°); radio-astronomy / FITS LSRK convention' % (lk, bk),
    'Gordon 1976; Greisen et al. 2006 (FITS LSRK); astropy LSRK', v=20.0, apex=[round(lk, 6), round(bk, 6)])
F[-1]['apex_icrs'] = (round(k_icrs.ra.deg, 6), round(k_icrs.dec.deg, 6))
add('lsrd', 'LSRD (dynamical, Delhaye 1965)', '(U, V, W) = (9, 12, 7) km/s (FITS LSRD; astropy LSRD)', 'Delhaye 1965; Greisen et al. 2006', uvw=[9.0, 12.0, 7.0])
add('lsr', 'LSR (Schönrich et al. 2010)', '(U, V, W) = (11.1 +0.69/−0.75, 12.24 +0.47/−0.47, 7.25 +0.37/−0.36) km/s (± 1, 2, 0.5 km/s syst.; astropy LSR default)', 'Schönrich, Binney & Dehnen 2010, MNRAS 403, 1829', uvw=[11.1, 12.24, 7.25])
# astropy Galactocentric v4.0 galcen_v_sun = (12.9, 245.6, 7.78) km/s is expressed in Galactocentric axes (centred on
# Sgr A*, tilted by z_sun = 20.8 pc); on Galactic axes it is (12.8798, 245.6000, 7.8126) km/s (astropy 7.2 matrix).
add('gsr', 'Galactocentric (astropy v4.0 defaults)', 'v⊙ = (12.9, 245.6, 7.78) km/s in Galactocentric axes = (12.880, 245.600, 7.813) km/s on Galactic axes (Drimmel & Poggio 2018; R₀ = 8.122 kpc, z⊙ = 20.8 pc)',
    'Drimmel & Poggio 2018; GRAVITY Collaboration 2018; Reid & Brunthaler 2004; Bennett & Bovy 2019', uvw=[12.8798, 245.6000, 7.8126])
add('gsr_mod', 'Galactocentric (Sgr A* proper motion × GRAVITY 2021 R₀)', 'v⊙ = (11.1, 251.49, 8.59) km/s: U from Schönrich et al. 2010, V and W from the Sgr A* proper motion (Reid & Brunthaler 2020) at R₀ = 8.275 kpc (derived combination)',
    'Reid & Brunthaler 2020; GRAVITY Collaboration 2021; Schönrich et al. 2010', uvw=[11.1, 251.4868, 8.5908])
add('gsr_iau', 'GSR (IAU 1985 Θ₀ = 220 km/s + LSRD) = NED/RC3', 'Θ₀ = 220 km/s toward l = 90° plus Delhaye (9, 12, 7) ⇒ 232.3 km/s toward (87.8°, 1.7°)',
    'Kerr & Lynden-Bell 1986; Delhaye 1965; de Vaucouleurs et al. 1991 (RC3); NED', uvw=[9.0, 232.0, 7.0])
add('lg_km96', 'Local Group (Karachentsev & Makarov 1996)', '316 ± 5 km/s toward (l, b) = (93° ± 2°, −4° ± 2°)', 'Karachentsev & Makarov 1996, AJ 111, 794', v=316.0, apex=[93.0, -4.0])
add('lg_cvdb99', 'Local Group (Courteau & van den Bergh 1999)', '306 ± 18 km/s toward (l, b) = (99° ± 5°, −4° ± 4°)', 'Courteau & van den Bergh 1999, AJ 118, 337', v=306.0, apex=[99.0, -4.0])
add('lg_diaz14', 'Local Group (Diaz et al. 2014)', '299 km/s toward (l, b) = (98.4°, −5.9°)', 'Diaz et al. 2014, MNRAS 443, 1688', v=299.0, apex=[98.4, -5.9])
add('lg_yts77', 'Local Group (Yahil, Tammann & Sandage 1977)', 'Δv = −79 cos l cos b + 296 sin l cos b − 36 sin b km/s', 'Yahil, Tammann & Sandage 1977, ApJ 217, 903', formula='yts77', v=round(float(np.sqrt(79**2+296**2+36**2)), 2))
add('lg_iau76', 'Local Group (RC2 / “IAU” 300 sin l cos b)', 'Δv = 300 sin l cos b km/s (FITS LOCALGRP)', 'de Vaucouleurs, de Vaucouleurs & Corwin 1976 (RC2); de Vaucouleurs 1976, Trans. IAU 16B, 201; Mould et al. 2000; Greisen et al. 2006', formula='iau76', v=300.0)
add('cmb', 'CMB rest frame (Planck 2018)', '369.82 ± 0.11 km/s toward (l, b) = (264.021° ± 0.011°, 48.253° ± 0.005°)', 'Planck Collaboration 2020, A&A 641, A1', v=369.82, apex=[264.021, 48.253])
add('cmb_cobe', 'CMB rest frame (COBE/FIRAS, Fixsen et al. 1996; NED “3K”)', '371.0 km/s toward (l, b) = (264.14°, 48.26°)', 'Fixsen et al. 1996, ApJ 473, 576; NED', v=371.0, apex=[264.14, 48.26])
json.dump({'frames': F}, open(os.path.join(ROOT, 'data', 'build', 'frames.json'), 'w'), indent=1)
for f in F: print(f['key'], f.get('v'), f.get('apex'), f.get('apex_icrs'))
