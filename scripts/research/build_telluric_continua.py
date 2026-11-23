#!/usr/bin/env python
"""Telluric continuum / collision-induced absorption bands 3000-25000 A.

Derived from:
  * HITRAN2024 CIA files (Karman et al. 2019 format): O2-air, O2-O2, N2-air
  * Serdyuchenko-Gorshelev O3 cross sections (2014; Zenodo 10.5281/zenodo.5793207)
  * AFGL-1986 U.S. Standard atmosphere (number density, T, O2, O3), N2 vmr from MIPAS-2007 RFM file
Vertical (airmass 1) optical depth from sea level:
  CIA:  tau(nu) = sum_layers k(nu, T_layer) * n_X * n_Y * dz        [k in cm^5 molecule^-2]
  O3 :  tau(lambda) = sigma(lambda, T) * N_O3                     [sigma in cm^2 molecule^-1]
Peak positions are grid maxima of the tabulated data (precision = table grid, ~1 cm^-1 or 0.01 nm),
NOT line positions: these are broad continua (FWHM tens to hundreds of A).

usage: build_telluric_continua.py CIADIR O3FILE AFGLDIR MIPASFILE OUTCSV OUTJSON
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.signal import find_peaks

CIADIR, O3FILE, AFGLDIR, MIPAS, OUTCSV, OUTJSON = sys.argv[1:7]


def read_cia(fn):
    sets = []
    with open(fn) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        h = lines[i]
        if not h.strip():
            i += 1
            continue
        pair = h[:20].strip()
        wmin, wmax = float(h[20:30]), float(h[30:40])
        n = int(h[40:47]); T = float(h[47:54])
        # HITRAN CIA header nominally A20,F10.3,F10.3,I7,F7.1,E10.3,F6.3,A21,I4 but comments
        # overflow in some files -> parse tokens: 7 numeric fields, free comment, integer ref last
        tok = h.split()
        ref = tok[-1] if tok[-1].isdigit() and len(tok) > 7 else ""
        comment = " ".join(tok[7:-1] if ref else tok[7:])
        a = np.array([[float(x) for x in l.split()[:2]] for l in lines[i + 1:i + 1 + n]])
        sets.append(dict(pair=pair, wmin=wmin, wmax=wmax, n=n, T=T, comment=comment, ref=ref,
                         nu=a[:, 0], k=a[:, 1]))
        i += 1 + n
    return sets


# atmosphere (sea level, 0.1 km layers to 60 km)
t1 = pd.read_csv(os.path.join(AFGLDIR, "table_1f.csv"))
t2 = pd.read_csv(os.path.join(AFGLDIR, "table_2a.csv"))
zf = np.arange(0.05, 60, 0.1)
n_air = np.exp(np.interp(zf, t1["z"], np.log(t1["n"])))
T = np.interp(zf, t1["z"], t1["t"])
o2 = np.interp(zf, t2["z"], t2["O2"]) * 1e-6
o3 = np.exp(np.interp(zf, t1["z"], np.log(t1["O3"]))) * 1e-6
# N2 vmr from MIPAS RFM file
txt = open(MIPAS).read().split("*")
n2v = None
for blk in txt:
    if blk.startswith("N2 "):
        n2v = np.array([float(x) for x in blk.split("\n", 1)[1].split()]) * 1e-6
        break
hgt = None
for blk in txt:
    if blk.startswith("HGT"):
        hgt = np.array([float(x) for x in blk.split("\n", 1)[1].split()])
        break
n2 = np.interp(zf, hgt, n2v)
dz = 0.1e5
nO2, nN2 = o2 * n_air, n2 * n_air
NO3 = (o3 * n_air).sum() * dz
DU = 2.6867e16  # molecules cm^-2 per Dobson unit (definition: 0.01 mm at STP)

rows = []
meta = dict(N_O3_cm2=NO3, O3_DU=NO3 / DU, N2_vmr_surface=float(n2[0]), O2_vmr_surface=float(o2[0]))


def tau_cia(sets, nx, ny, wmin_key):
    # sets of one band at several T: choose per layer the set with nearest T, interpolate to common grid
    ref = max(sets, key=lambda s: s["n"])
    grid = ref["nu"]
    Ts = np.array([s["T"] for s in sets])
    kk = np.array([np.interp(grid, s["nu"], s["k"], left=0, right=0) for s in sets])
    tau = np.zeros_like(grid)
    for j in range(len(zf)):
        it = np.argmin(abs(Ts - T[j]))
        tau += kk[it] * nx[j] * ny[j] * dz
    return grid, tau


def add_band(pair, comment, grid, tau, refcode, note, prominence_frac=0.2):
    m = (grid > 0)
    pk, prop = find_peaks(tau, prominence=prominence_frac * tau.max())
    if len(pk) == 0:
        pk = [int(np.argmax(tau))]
    for p in pk:
        half = tau[p] / 2
        lo = p
        while lo > 0 and tau[lo] > half:
            lo -= 1
        hi = p
        while hi < len(tau) - 1 and tau[hi] > half:
            hi += 1
        lam = 1e8 / grid[p]
        rows.append(dict(species=pair, band=comment, source=refcode,
                         peak_wavenumber_cm1=round(float(grid[p]), 3),
                         peak_lambda_vac_A=round(float(lam), 1),
                         fwhm_lambda_vac_A_lo=round(float(1e8 / grid[hi]), 1),
                         fwhm_lambda_vac_A_hi=round(float(1e8 / grid[lo]), 1),
                         tau_peak_airmass1_sealevel=float("%.3g" % tau[p]),
                         absorption_peak_percent=float("%.3g" % (100 * (1 - np.exp(-tau[p])))),
                         note=note))


oair = read_cia(os.path.join(CIADIR, "O2-air_2024.cia"))
ooo = read_cia(os.path.join(CIADIR, "O2-O2_2024.cia"))
nair_ = read_cia(os.path.join(CIADIR, "N2-air_2018.cia"))
# O2-air bands within 4000-33333 cm^-1
groups = {}
for s in oair:
    if s["wmax"] < 4000:
        continue
    key = round(s["wmin"], -2) if s["wmin"] < 8400 else round(s["wmin"])
    groups.setdefault(("O2-Air", s["comment"] or "a(0)-X(0) 1.27um", round(s["wmax"], -2)), []).append(s)
for (pair, comment, _), sets in groups.items():
    grid, tau = tau_cia(sets, nO2, n_air, None)
    add_band("O2-air CIA", comment, grid, tau, f"HITRAN2024 CIA O2-air_2024.cia ref {sets[0]['ref']}",
             "tau = k(O2-air) n_O2 n_air dz")
# O2-O2 bands not covered by O2-air (B band, double transitions)
groups = {}
for s in ooo:
    if s["wmax"] < 4000 or s["comment"] in ("a(1)-X(0)", "a(2)-X(0)") or s["comment"].startswith("A Band"):
        continue
    if s["wmin"] < 8400:
        continue  # 1.27 um band covered by O2-air
    groups.setdefault((s["comment"], round(s["wmin"], -3), s["ref"]), []).append(s)
for (comment, _, ref), sets in groups.items():
    grid, tau = tau_cia(sets, nO2, nO2, None)
    add_band("O2-O2 CIA", comment, grid, tau, f"HITRAN2024 CIA O2-O2_2024.cia ref {ref}",
             "tau = k(O2-O2) n_O2^2 dz", prominence_frac=0.05)
# N2-air 2.15 um
sets = [s for s in nair_ if s["wmin"] >= 4000]
if sets:
    grid, tau = tau_cia(sets, nN2, n_air, None)
    add_band("N2-air CIA", "N2 CIA near 2.16 um (" + sets[0]["comment"] + ")", grid, tau,
             f"HITRAN2024 CIA N2-air_2018.cia ref {sets[0]['ref']}", "tau = k(N2-air) n_N2 n_air dz")

# O3 continuum (Serdyuchenko-Gorshelev), 223 K column (stratospheric) -> column 9
a = np.loadtxt(O3FILE, skiprows=45)
lam_nm, sig223, sig293 = a[:, 0], a[:, 8], a[:, 1]
tau3 = sig223 * NO3
lamA = lam_nm * 10
sel = (lamA >= 3000) & (lamA <= 11000)
# Chappuis peak: smooth with 5 nm boxcar on 0.01 nm grid
w = 500
sm = np.convolve(tau3, np.ones(w) / w, mode="same")
m = (lamA > 4500) & (lamA < 7500)
ip = np.where(m)[0][np.argmax(sm[m])]
half = sm[ip] / 2
lo = ip
while sm[lo] > half:
    lo -= 1
hi = ip
while sm[hi] > half:
    hi += 1
rows.append(dict(species="O3", band="Chappuis band (continuum)",
                 source="Serdyuchenko et al. 2014 (223 K), 5 nm boxcar",
                 peak_wavenumber_cm1=round(1e8 / lamA[ip], 1), peak_lambda_vac_A=round(lamA[ip], 1),
                 fwhm_lambda_vac_A_lo=round(lamA[lo], 1), fwhm_lambda_vac_A_hi=round(lamA[hi], 1),
                 tau_peak_airmass1_sealevel=float("%.3g" % sm[ip]),
                 absorption_peak_percent=float("%.3g" % (100 * (1 - np.exp(-sm[ip])))),
                 note=f"tau = sigma(223K) N_O3; N_O3 = {NO3:.3e} cm^-2 ({NO3 / DU:.0f} DU, AFGL US Std)"))
# Huggins band: wavelengths where tau crosses 3, 1, 0.1, 0.01 (going redward from 3000 A)
mm = (lamA >= 3000) & (lamA <= 4000)
for thr in [3.0, 1.0, 0.1, 0.01]:
    idx = np.where(mm & (tau3 < thr))[0]
    # first wavelength beyond which tau stays below thr (in 1 nm smoothed curve)
    s1 = np.convolve(tau3, np.ones(100) / 100, mode="same")
    ok = np.where(mm & (s1 >= thr))[0]
    lam_cross = lamA[ok.max()] if len(ok) else np.nan
    rows.append(dict(species="O3", band=f"Hartley-Huggins edge: tau(1-nm smoothed) last >= {thr}",
                     source="Serdyuchenko et al. 2014 (223 K)", peak_wavenumber_cm1=np.nan,
                     peak_lambda_vac_A=round(float(lam_cross), 1), fwhm_lambda_vac_A_lo=np.nan,
                     fwhm_lambda_vac_A_hi=np.nan, tau_peak_airmass1_sealevel=thr,
                     absorption_peak_percent=float("%.3g" % (100 * (1 - np.exp(-thr)))),
                     note="UV cutoff of ground-based spectra; Huggins vibronic structure 3100-3600 A"))
# Huggins band local maxima 3100-3500 A
s1 = np.convolve(tau3, np.ones(50) / 50, mode="same")
mh = np.where((lamA >= 3100) & (lamA <= 3600))[0]
pk, _ = find_peaks(s1[mh], prominence=0.02 * s1[mh].max(), distance=100)
meta["huggins_peaks_vac_A_223K_0p5nm_smoothed"] = [round(float(lamA[mh][p]), 1) for p in pk]
meta["o3_tau_at"] = {f"{l}": float("%.3g" % np.interp(l, lamA, tau3)) for l in
                     [3000, 3100, 3200, 3300, 3400, 3500, 4500, 5000, 5500, 6000, 6500, 7000, 8000, 9000, 10000]}
meta["o3_tau_293K_vs_223K_6000A"] = float(np.interp(6000, lamA, sig293) / np.interp(6000, lamA, sig223))

REFS = {
    "4": "Tran, Boulet & Hartmann 2006, JGR 111, D15210, doi:10.1029/2005JD006869",
    "24": "Hartmann, Boulet & Toon 2017, JGR Atmos 122, 2419, doi:10.1002/2016JD025677",
    "27": "Karman et al. 2018, Nature Chem 10, 549, doi:10.1038/s41557-018-0015-x",
    "29": "Spiering & van der Zande 2012, PCCP 14, 9923, doi:10.1039/c2cp40961e",
    "30": "Spiering et al. 2011, PCCP 13, 9616, doi:10.1039/C1CP20403C",
    "32": "Thalman & Volkamer 2013, PCCP 15, 15371, doi:10.1039/c3cp50968k",
    "53": "Mondelain et al. 2019, JGR Atmos 124, 414, doi:10.1029/2018JD029317; Kassi et al. 2021, JGR Atmos 126, e2021JD034860, doi:10.1029/2021JD034860",
    "54": "Finkenzeller & Volkamer 2022, JQSRT 279, 108063, doi:10.1016/j.jqsrt.2021.108063",
    "58": "Adkins et al. 2023, JQSRT 310, 108732, doi:10.1016/j.jqsrt.2023.108732",
}
for r in rows:
    code = r["source"].split(" ref ")[-1].strip() if " ref " in r["source"] else ""
    r["original_reference"] = REFS.get(code, "Serdyuchenko et al. 2014, AMT 7, 625, doi:10.5194/amt-7-625-2014; data doi:10.5281/zenodo.5793207" if r["species"] == "O3" else "")
# O2-O2 double transitions: two overlapping data sets (ref 32 to 630 nm; ref 54 297-500 nm, newer).
# keep ref 54 where both give the same peak (within 10 A)
keep = []
for r in rows:
    if r["species"] == "O2-O2 CIA" and r["source"].endswith(" ref 32"):
        if any(q["species"] == "O2-O2 CIA" and q["source"].endswith(" ref 54")
               and abs(q["peak_lambda_vac_A"] - r["peak_lambda_vac_A"]) < 10 for q in rows):
            continue
    keep.append(r)
rows = keep
df = pd.DataFrame(rows).sort_values("peak_lambda_vac_A")
df.to_csv(OUTCSV, index=False)
json.dump(meta, open(OUTJSON, "w"), indent=1, default=float)
pd.set_option("display.width", 250)
print(df.to_string())
print(json.dumps(meta, indent=1, default=float))
