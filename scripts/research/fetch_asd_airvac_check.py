#!/usr/bin/env python
"""Fetch NIST ASD line lists for a few well-measured spectra, in three wavelength-medium
modes, to reverse-engineer / validate the air<->vacuum conversion that ASD applies.

  show_av=3 : vacuum, all wavelengths
  show_av=2 : vacuum <2000 A, air 2000-20000 A, vacuum >20000 A   (ASD default)
  show_av=4 : vacuum <1850 A, air >1850 A  (per the ASD lines form label)

Usage: python fetch_asd_airvac_check.py <outdir>
Output: <outdir>/<El>_<ion>_av<N>.csv  (raw ASD CSV, untrusted data)
"""
import os, sys, time, urllib.parse, urllib.request

OUTD = sys.argv[1]
os.makedirs(OUTD, exist_ok=True)
BASE = 'https://physics.nist.gov/cgi-bin/ASD/lines1.pl?'
COMMON = dict(unit=0, format=2, line_out=0, output=0, show_obs_wl=1, show_calc_wl=1, unc_out=1,
              show_wn=1, intens_out='on', allowed_out=1, forbid_out=1, no_spaces=1,
              enrg_out='on', J_out='on', term_out='on', conf_out='on', A_out=0, bibrefs=1, order_out=0, remove_js='on',
              page_size=15)
SPECTRA = [('Fe', 'I', 1800, 26000), ('Ne', 'I', 1800, 26000), ('Ar', 'I', 1800, 26000),
           ('Kr', 'I', 1800, 26000), ('Hg', 'I', 1800, 26000), ('Na', 'I', 1800, 26000),
           ('H', 'I', 1800, 26000), ('Th', 'I', 1800, 26000)]
for el, ion, lo, hi in SPECTRA:
    for av in (3, 2, 4):
        fn = os.path.join(OUTD, f'{el}_{ion}_av{av}.csv')
        if os.path.exists(fn) and os.path.getsize(fn) > 200:
            continue
        q = dict(COMMON, spectra=f'{el} {ion}', low_w=str(lo), upp_w=str(hi), show_av=av)
        url = BASE + urllib.parse.urlencode(q)
        for attempt in range(4):
            try:
                req = urllib.request.Request(url, headers={'User-Agent': 'research-script/1.0'})
                txt = urllib.request.urlopen(req, timeout=240).read().decode('utf-8', 'replace')
                break
            except Exception as e:  # noqa
                print('retry', el, ion, av, e, file=sys.stderr)
                time.sleep(10 * (attempt + 1))
        else:
            print('FAILED', el, ion, av)
            continue
        open(fn, 'w').write(txt)
        print(f'{el} {ion} av={av}: {txt.count(chr(10))} rows', flush=True)
        time.sleep(2)
print('DONE')
