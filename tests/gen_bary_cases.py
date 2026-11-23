"""Generate barycentric-correction reference values with astropy + DE440s for the JS core."""
import os
import json, sys, warnings
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, solar_system_ephemeris
warnings.filterwarnings('ignore')
solar_system_ephemeris.set('file://' + os.path.expanduser('~/JPLDE/de440s.bsp'))
rng = np.random.default_rng(42)
sites = [(-155.4681, 19.8283, 4160), (-70.4045, -24.6272, 2635), (-70.7345, -29.2563, 2347),
         (117.5750, 40.3958, 960), (0.0, 0.0, 0), (17.88, 28.76, 2396), (-17.8792, 28.7569, 2326), (100.0302, 26.6951, 3200)]
N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
mjd = rng.uniform(48000, 62500, N)   # 1990 .. 2030 (IERS-A predicted beyond)
out = []
for i in range(N):
    lon, lat, h = sites[i % len(sites)]
    ra = rng.uniform(0, 360); dec = np.degrees(np.arcsin(rng.uniform(-1, 1)))
    t = Time(Time(mjd[i], format='mjd', scale='utc').isot, scale='utc')
    loc = EarthLocation.from_geodetic(lon * u.deg, lat * u.deg, h * u.m)
    sc = SkyCoord(ra * u.deg, dec * u.deg)
    v = sc.radial_velocity_correction(kind='barycentric', obstime=t, location=loc).to(u.km / u.s).value
    vh = sc.radial_velocity_correction(kind='heliocentric', obstime=t, location=loc).to(u.km / u.s).value
    tdb = Time(t, location=loc).tdb
    try:
        dut1 = float(t.delta_ut1_utc)
    except Exception:
        dut1 = 0.0
    lt = t.light_travel_time(sc, kind='barycentric', location=loc, ephemeris='jpl').to(u.s).value \
        if False else (t.tdb + Time(t).light_travel_time(sc, 'barycentric', location=loc)).jd
    out.append(dict(ra=ra, dec=dec, lon=lon, lat=lat, h=h, iso=t.isot, dut1=dut1, vcorr=v, vhelio=vh,
                    jd_tdb=t.tdb.jd, bjd=lt))
json.dump(out, open('tests/bary_cases.json', 'w'))
print('wrote', len(out))
