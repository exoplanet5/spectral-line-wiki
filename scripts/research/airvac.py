"""Air <-> vacuum wavelength conversions (standard dry air, 15 C, 101325 Pa, 0.03-0.045% CO2).

All functions take wavelengths in Angstrom.

* Peck & Reeder (1972), J. Opt. Soc. Am. 62, 958  -- used by NIST ASD for 200 nm - 2 um
      n - 1 = 1e-8 [ 8060.51 + 2480990/(132.274 - s2) + 17455.7/(39.32957 - s2) ],  s = 1/lambda_vac[um]
* Edlen (1966), Metrologia 2, 71 (as quoted by Morton 2000, ApJS 130, 403 -- the IAU standard)
      n - 1 = 1e-8 [ 8342.13 + 2406030/(130 - s2) + 15997/(38.9 - s2) ]
* Ciddor (1996), Appl. Opt. 35, 1566 (standard air, 450 ppm CO2) -- as written in the SDSS / VALD notes
      n - 1 = 0.05792105/(238.0185 - s2) + 0.00167917/(57.362 - s2)
"""
import numpy as np


def n_peck_reeder(lam_vac_A):
    s2 = (1e4 / np.asarray(lam_vac_A, float)) ** 2
    return 1 + 1e-8 * (8060.51 + 2480990.0 / (132.274 - s2) + 17455.7 / (39.32957 - s2))


def n_edlen66(lam_vac_A):
    s2 = (1e4 / np.asarray(lam_vac_A, float)) ** 2
    return 1 + 1e-8 * (8342.13 + 2406030.0 / (130.0 - s2) + 15997.0 / (38.9 - s2))


def n_ciddor96(lam_vac_A):
    s2 = (1e4 / np.asarray(lam_vac_A, float)) ** 2
    return 1 + 0.05792105 / (238.0185 - s2) + 0.00167917 / (57.362 - s2)


_N = {'peck_reeder': n_peck_reeder, 'edlen66': n_edlen66, 'ciddor96': n_ciddor96}


def vac_to_air(lam_vac_A, formula='peck_reeder'):
    lam = np.asarray(lam_vac_A, float)
    return lam / _N[formula](lam)


def air_to_vac(lam_air_A, formula='peck_reeder', iters=6):
    lam_air = np.asarray(lam_air_A, float)
    lam = lam_air.copy()
    for _ in range(iters):
        lam = lam_air * _N[formula](lam)
    return lam


if __name__ == '__main__':
    # check against NIST ASD Na I D2: vac 5891.583264, air 5889.95095
    for f in _N:
        print(f, '%.5f' % vac_to_air(5891.583264, f), '%.6f' % air_to_vac(5889.95095, f))
