"""Compare the page's frame projections with astropy transforms (radial velocity of a source at rest in the
barycentric frame, seen in LSRK / LSRD / LSR / Galactocentric)."""
import json, subprocess, warnings, numpy as np, astropy.units as u
from astropy.coordinates import SkyCoord, ICRS, LSRK, LSRD, LSR, Galactocentric, galactocentric_frame_defaults
warnings.filterwarnings('ignore')
F = {f['key']: f for f in json.load(open('data/build/frames.json'))['frames']}
rng = np.random.default_rng(3)
ra = rng.uniform(0, 360, 40); dec = np.degrees(np.arcsin(rng.uniform(-1, 1, 40)))
js = subprocess.run(['node', 'tests/frames_js.js'], input=json.dumps([F, list(ra), list(dec)]), capture_output=True, text=True).stdout
J = json.loads(js)
res = {k: [] for k in ['lsrk', 'lsrd', 'lsr', 'gsr']}
with galactocentric_frame_defaults.set('v4.0'):
    for i in range(len(ra)):
        c = ICRS(ra=ra[i]*u.deg, dec=dec[i]*u.deg, distance=1*u.kpc, pm_ra_cosdec=0*u.mas/u.yr, pm_dec=0*u.mas/u.yr, radial_velocity=0*u.km/u.s)
        for k, fr in [('lsrk', LSRK()), ('lsrd', LSRD()), ('lsr', LSR())]:
            res[k].append(c.transform_to(fr).radial_velocity.to(u.km/u.s).value - J[i][k])
        g = c.transform_to(Galactocentric())
        # radial velocity in the Galactocentric frame as seen from the Sun's position: project the
        # source's galactocentric velocity on the Sun->source unit vector
        sun = Galactocentric().realize_frame(ICRS(ra=0*u.deg, dec=0*u.deg, distance=0*u.kpc).transform_to(Galactocentric()).cartesian)
        sv = ICRS(ra=0*u.deg, dec=0*u.deg, distance=0*u.kpc).transform_to(Galactocentric()).cartesian.xyz.to(u.kpc).value
        p = g.cartesian.xyz.to(u.kpc).value; n = (p - sv) / np.linalg.norm(p - sv)
        v = g.velocity.d_xyz.to(u.km/u.s).value
        res['gsr'].append(np.dot(v, n) - J[i]['gsr'])
print('ICRS->Gal->ICRS round trip max (deg):', max(x['rt'] for x in J))
for k, v in res.items():
    print(f'{k:5s} page - astropy: max |dv| = {np.max(np.abs(v)) * 1000:.2f} m/s')
