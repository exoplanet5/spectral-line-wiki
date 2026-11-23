#!/usr/bin/env python
"""Build the curated telluric absorption-line list (3000-25000 A) from HITRAN2024.

Inputs (all downloaded, treated as data only):
  PARDIR    HITRAN 160-char .par chunks (fetch_hitran.py)
  QDIR      HITRAN Q(T) files q<gid>.txt (fetch_hitran_q.py)
  MOLPARAM  HITRAN molparam.txt
  AFGLDIR   AFGL-1986 tables (joseki/data/afgl_1986: table_1f.csv, table_2a.csv)
  NOAADIR   NOAA GML global annual means (co2/ch4/n2o_annmean_gl.txt)
Outputs:
  OUTCSV    curated line list
  OUTMETA   JSON metadata (columns, model, counts, precision statistics)
  OUTBANDS  per-band summary CSV

Method (everything derived is documented in OUTMETA):
  * lambda_vac = 1e8 / nu  (nu = HITRAN vacuum wavenumber, cm^-1)
  * lambda_air = lambda_vac / n, n from Peck & Reeder (1972, JOSA 62, 958) eq. for
    standard air (15 C, 101325 Pa, dry) -- the formula used by NIST ASD:
        (n-1)*1e8 = 8060.51 + 2480990/(132.274 - s^2) + 17455.7/(39.32957 - s^2),
        s = 1/lambda_vac[um]
  * line-centre optical depth / depth / equivalent width for an isolated line through
    a layered AFGL-1986 U.S. Standard atmosphere (CO2, CH4, N2O rescaled to NOAA 2025
    global means), airmass 1, Voigt profile, HITRAN S(T) with TIPS Q(T), gamma(p,T).
  * pressure shift: nu* = nu + delta_air * p  (HITRAN definition), evaluated at 1 atm
    and at the column-intensity-weighted mean pressure of the absorbing column.

usage: build_telluric.py PARDIR QDIR MOLPARAM AFGLDIR NOAADIR OUTCSV OUTMETA OUTBANDS
"""
import csv
import glob
import json
import math
import os
import sys

import numpy as np
import pandas as pd
from scipy.special import voigt_profile

PARDIR, QDIR, MOLPARAM, AFGLDIR, NOAADIR, OUTCSV, OUTMETA, OUTBANDS = sys.argv[1:9]

# ---------------------------------------------------------------- constants (CODATA 2018 exact / HITRAN defs)
C2 = 1.4387769          # cm K   (HITRAN definitions page)
KB = 1.380649e-23       # J/K    (exact SI)
AMU = 1.66053906660e-27  # kg    (CODATA 2018)
CLIGHT = 299792458.0    # m/s    (exact)
NA = 6.02214076e23      # 1/mol  (exact SI)
TREF = 296.0
ATM_HPA = 1013.25

MOLS = {1: "H2O", 2: "CO2", 3: "O3", 4: "N2O", 5: "CO", 6: "CH4", 7: "O2"}
GID = {  # HITRAN global isotopologue IDs in local-ID order (hitran.org/docs/iso-meta)
    1: [1, 2, 3, 4, 5, 6, 129],
    2: [7, 8, 9, 10, 11, 12, 13, 14, 121, 15, 120, 122],
    3: [16, 17, 18, 19, 20],
    4: [21, 22, 23, 24, 25],
    5: [26, 27, 28, 29, 30, 31],
    6: [32, 33, 34, 35],
    7: [36, 37, 38],
}
ISOCHAR = {**{str(i): i for i in range(1, 10)}, "0": 10, "A": 11, "B": 12}

# ---------------------------------------------------------------- molparam
molparam = {}  # (mol, localiso) -> dict
cur = None
for line in open(MOLPARAM):
    s = line.split()
    if not s or s[0] == "Molecule":
        continue
    if len(s) == 2 and s[1].startswith("("):
        cur = int(s[1].strip("()"))
        k = 0
        continue
    if cur in MOLS:
        k += 1
        molparam[(cur, k)] = dict(afgl=s[0], abund=float(s[1]), q296=float(s[2]),
                                  gj=int(s[3]), mass=float(s[4]))

# ---------------------------------------------------------------- partition sums
QT = {}


def qfunc(gid):
    if gid not in QT:
        a = np.loadtxt(os.path.join(QDIR, f"q{gid}.txt"))
        QT[gid] = a
    a = QT[gid]
    return lambda T: np.interp(T, a[:, 0], a[:, 1])


# ---------------------------------------------------------------- parse HITRAN
recs = set()
rows = []
for fn in sorted(glob.glob(os.path.join(PARDIR, "*.par"))):
    with open(fn, "rb") as f:
        for raw in f:
            raw = raw.rstrip(b"\r\n")
            if len(raw) != 160:
                continue
            s = raw.decode("ascii")
            if s in recs:
                continue
            recs.add(s)
            rows.append(s)
print("parsed unique records:", len(rows), flush=True)

cols = dict(mol=(0, 2), iso=(2, 3), nu=(3, 15), S=(15, 25), A=(25, 35), gair=(35, 40),
            gself=(40, 45), elow=(45, 55), nair=(55, 59), dair=(59, 67), Vu=(67, 82),
            Vl=(82, 97), Qu=(97, 112), Ql=(112, 127), ierr=(127, 133), iref=(133, 145),
            lm=(145, 146), gu=(146, 153), gl=(153, 160))
df = pd.DataFrame({k: [r[a:b] for r in rows] for k, (a, b) in cols.items()})
del rows, recs
df["mol"] = df["mol"].astype(int)
df["liso"] = df["iso"].map(ISOCHAR)
for k in ["nu", "S", "A", "gair", "gself", "elow", "nair", "dair"]:
    df[k + "_f"] = df[k].astype(float)
df = df[df["mol"].isin(MOLS)].reset_index(drop=True)
df["gid"] = [GID[m][i - 1] for m, i in zip(df["mol"], df["liso"])]
df = df[(df["nu_f"] >= 4000.0) & (df["nu_f"] <= 1e8 / 3000.0)].reset_index(drop=True)
print("records in 4000-33333.3 cm-1:", len(df), flush=True)
nrec_by_mol = df.groupby("mol").size().rename(index=MOLS).to_dict()

# ---------------------------------------------------------------- atmosphere
t1 = pd.read_csv(os.path.join(AFGLDIR, "table_1f.csv"))   # US Standard: z,p(mb),t,n(cm-3),H2O,O3,N2O,CO,CH4 ppmv
t2 = pd.read_csv(os.path.join(AFGLDIR, "table_2a.csv"))   # z,...,CO2,O2 ppmv


def noaa_last(fn):
    yrs, vals = [], []
    for line in open(fn):
        if line.startswith("#") or not line.strip():
            continue
        s = line.split()
        yrs.append(int(s[0])); vals.append(float(s[1]))
    return yrs[-1], vals[-1]


yco2, co2 = noaa_last(os.path.join(NOAADIR, "co2_annmean_gl.txt"))   # ppm
ych4, ch4 = noaa_last(os.path.join(NOAADIR, "ch4_annmean_gl.txt"))   # ppb
yn2o, n2o = noaa_last(os.path.join(NOAADIR, "n2o_annmean_gl.txt"))   # ppb
scale = {"CO2": co2 / t2["CO2"].iloc[0], "CH4": (ch4 / 1000) / t1["CH4"].iloc[0],
         "N2O": (n2o / 1000) / t1["N2O"].iloc[0]}

zt = t1["z"].values
prof = {"p": t1["p"].values, "t": t1["t"].values, "n": t1["n"].values,
        "H2O": t1["H2O"].values, "O3": t1["O3"].values, "N2O": t1["N2O"].values * scale["N2O"],
        "CO": t1["CO"].values, "CH4": t1["CH4"].values * scale["CH4"],
        "CO2": np.interp(zt, t2["z"].values, t2["CO2"].values) * scale["CO2"],
        "O2": np.interp(zt, t2["z"].values, t2["O2"].values)}


def make_layers(z0, ztop=80.0):
    edges = [z0]
    for zmax, dz in [(3, 0.25), (12, 0.5), (30, 1.0), (50, 2.5), (ztop, 5.0)]:
        while edges[-1] + dz <= zmax + 1e-9:
            edges.append(round(edges[-1] + dz, 6))
        if edges[-1] < zmax - 1e-9 and zmax <= ztop:
            edges.append(zmax)
    edges = np.array(sorted(set(edges)))
    edges = edges[edges >= z0]
    zf = np.arange(z0, ztop + 1e-9, 0.005)
    lp = np.interp(zf, zt, np.log(prof["p"]))
    ln = np.interp(zf, zt, np.log(prof["n"]))
    tt = np.interp(zf, zt, prof["t"])
    n = np.exp(ln)
    p = np.exp(lp)
    L = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (zf >= a - 1e-9) & (zf <= b + 1e-9)
        z_, n_, p_, t_ = zf[m], n[m], p[m], tt[m]
        u_air = np.trapezoid(n_, z_ * 1e5)  # molecules/cm2 of air
        lay = dict(z0=a, z1=b, u_air=u_air,
                   p_atm=np.trapezoid(n_ * p_, z_ * 1e5) / u_air / ATM_HPA,
                   T=np.trapezoid(n_ * t_, z_ * 1e5) / u_air)
        for g in ["H2O", "CO2", "O3", "N2O", "CO", "CH4", "O2"]:
            v = np.exp(np.interp(z_, zt, np.log(prof[g]))) * 1e-6
            lay["u_" + g] = np.trapezoid(n_ * v, z_ * 1e5)
            lay["vmr_" + g] = lay["u_" + g] / u_air
        L.append(lay)
    return pd.DataFrame(L)


SCEN = {"sl": 0.0, "2p5km": 2.5}
LAY = {k: make_layers(z) for k, z in SCEN.items()}
atm_summary = {}
for k, L in LAY.items():
    pwv_mm = L["u_H2O"].sum() * 18.015 / NA * 10.0  # g/cm2 -> mm (rho_w = 1 g/cm3)
    atm_summary[k] = dict(z_obs_km=SCEN[k], n_layers=len(L), PWV_mm=round(pwv_mm, 2),
                          p_lowest_layer_mean_hPa=round(L["p_atm"].iloc[0] * ATM_HPA, 1),
                          columns_cm2={g: float("%.4e" % L["u_" + g].sum())
                                       for g in ["H2O", "CO2", "O3", "N2O", "CO", "CH4", "O2"]})
print(json.dumps(atm_summary, indent=1), flush=True)

# ---------------------------------------------------------------- prefilter (linear EW at 296 K)
umol_sl = {MOLS[m]: LAY["sl"]["u_" + MOLS[m]].sum() for m in MOLS}
df["Wlin296"] = df["S_f"] * df["mol"].map(lambda m: umol_sl[MOLS[m]])
cand = df[df["Wlin296"] > 1e-4].reset_index(drop=True)
print("prefilter candidates:", len(cand), flush=True)


def line_layer_params(sub, L):
    """return c (nl,nlay) = S(T)*u, sigma_G, gamma_L, p (nlay) for lines in sub."""
    mol = sub["mol"].values
    nu = sub["nu_f"].values[:, None]
    T = L["T"].values[None, :]
    p = L["p_atm"].values[None, :]
    S296 = sub["S_f"].values[:, None]
    El = sub["elow_f"].values[:, None]
    Q296 = np.array([molparam[(m, i)]["q296"] for m, i in zip(mol, sub["liso"])])[:, None]
    QTl = np.vstack([qfunc(g)(L["T"].values) for g in sub["gid"].values])
    with np.errstate(over="ignore", invalid="ignore"):
        ST = (S296 * (Q296 / QTl) * np.exp(-C2 * El / T) / np.exp(-C2 * El / TREF)
              * (1 - np.exp(-C2 * nu / T)) / (1 - np.exp(-C2 * nu / TREF)))
    u = np.vstack([L["u_" + MOLS[m]].values for m in mol])
    vmr = np.vstack([L["vmr_" + MOLS[m]].values for m in mol])
    psf = vmr * p
    gam = (TREF / T) ** sub["nair_f"].values[:, None] * (
        sub["gair_f"].values[:, None] * (p - psf) + sub["gself_f"].values[:, None] * psf)
    mass = np.array([molparam[(m, i)]["mass"] for m, i in zip(mol, sub["liso"])])[:, None] * AMU
    sig = nu / CLIGHT * np.sqrt(KB * T / mass)  # Gaussian sigma in cm^-1
    return ST * u, sig, gam, np.broadcast_to(p, ST.shape)


def tau0(sub, L):
    c, sig, gam, p = line_layer_params(sub, L)
    v0 = voigt_profile(0.0, sig, gam)
    return (c * v0).sum(1), (c * p).sum(1) / c.sum(1)


CH = 20000
res = {k: [] for k in SCEN}
pbar = []
for i0 in range(0, len(cand), CH):
    sub = cand.iloc[i0:i0 + CH]
    for k, L in LAY.items():
        t0, pb = tau0(sub, L)
        res[k].append(t0)
        if k == "sl":
            pbar.append(pb)
for k in SCEN:
    cand["tau0_" + k] = np.concatenate(res[k])
cand["pbar_sl"] = np.concatenate(pbar)
cand["depth_sl"] = 1 - np.exp(-cand["tau0_sl"])
cand["depth_2p5km"] = 1 - np.exp(-cand["tau0_2p5km"])

DEPTH_MIN = 0.01   # selection threshold on the R=100000 central depth (sea level, airmass 1)
RES = 100000.0     # reference resolving power (Gaussian LSF, FWHM = lambda/R)
pre = cand[cand["depth_sl"] >= DEPTH_MIN].sort_values("nu_f").reset_index(drop=True)
print("infinite-resolution depth_sl >= %.2f:" % DEPTH_MIN, len(pre), flush=True)


# ---------------------------------------------------------------- equivalent widths and R=1e5 depths
def ew(sub, L, npts=601):
    """isolated-line EW (cm^-1) and central depth after Gaussian LSF of FWHM nu/RES."""
    c, sig, gam, p = line_layer_params(sub, L)
    t0 = (c * voigt_profile(0.0, sig, gam)).sum(1)
    gmax = gam.max(1)
    smin = np.minimum((sig * 1.1774).min(1), gam.min(1) + 1e-12)
    a = np.maximum(0.3 * smin, 1e-5)
    sR = sub["nu_f"].values / RES / (2 * np.sqrt(2 * np.log(2)))
    H = np.maximum(40 * gmax * np.sqrt(np.maximum(t0, 1.0)), 12 * sig.max(1) * 1.1774)
    H = np.maximum(np.maximum(H, 20 * gmax), 8 * sR)
    W = np.empty(len(sub)); DR = np.empty(len(sub)); WL = (c.sum(1))
    for j in range(len(sub)):
        tmax = np.arcsinh(H[j] / a[j])
        x = a[j] * np.sinh(np.linspace(-tmax, tmax, npts))
        tau = (c[j][:, None] * voigt_profile(x[None, :], sig[j][:, None], gam[j][:, None])).sum(0)
        ab = -np.expm1(-tau)
        # analytic Lorentz far-wing tail beyond +-H (optically thin there): sum_l c_l (2/pi) gamma_l / H
        W[j] = np.trapezoid(ab, x) + (c[j] * gam[j]).sum() * 2 / np.pi / H[j]
        g = np.exp(-0.5 * (x / sR[j]) ** 2) / (sR[j] * np.sqrt(2 * np.pi))
        DR[j] = np.trapezoid(ab * g, x) / np.trapezoid(g, x)
    return W, DR, WL


def run_ew(df_, L):
    out = [ew(df_.iloc[i0:i0 + 2000], L) for i0 in range(0, len(df_), 2000)]
    return [np.concatenate([o[k] for o in out]) for k in range(3)]


W, DR, WL = run_ew(pre, LAY["sl"])
pre["EW_cm_sl"], pre["depthR_sl"], pre["Wlin_sl"] = W, DR, WL
weak = pre["tau0_sl"] < 0.02
ew_check = float(np.median(pre.loc[weak, "EW_cm_sl"] / pre.loc[weak, "Wlin_sl"])) if weak.any() else None
print("check: median EW/(sum S u) for tau0<0.02 lines =", ew_check, flush=True)
sel = pre[pre["depthR_sl"] >= DEPTH_MIN].reset_index(drop=True)
print("selected (R=1e5 depth_sl >= %.2f):" % DEPTH_MIN, len(sel), flush=True)
W, DR, WL = run_ew(sel, LAY["2p5km"])
sel["EW_cm_2p5km"], sel["depthR_2p5km"] = W, DR
print("EW done", flush=True)


# ---------------------------------------------------------------- wavelengths, shifts, labels
def n_peck_reeder(lam_vac_A):
    s2 = (1e4 / lam_vac_A) ** 2
    return 1 + 1e-8 * (8060.51 + 2480990.0 / (132.274 - s2) + 17455.7 / (39.32957 - s2))


def n_ciddor450(lam_vac_A):
    s2 = (1e4 / lam_vac_A) ** 2
    return 1 + 1e-8 * (5792105.0 / (238.0185 - s2) + 167917.0 / (57.362 - s2))


def n_morton2000(lam_vac_A):
    s2 = (1e4 / lam_vac_A) ** 2
    return 1 + 8.34254e-5 + 2.406147e-2 / (130 - s2) + 1.5998e-4 / (38.9 - s2)


nu = sel["nu_f"].values
lv = 1e8 / nu
la = lv / n_peck_reeder(lv)
sel["lambda_vac_A"] = lv
sel["lambda_air_A"] = la
dla = lv / n_ciddor450(lv) - la
dlm = lv / n_morton2000(lv) - la

code = sel["ierr"].str[0].astype(int)
nu_unc = np.where(code > 0, 10.0 ** (-(code - 1)), np.nan)
sel["nu_unc_code"] = code
sel["nu_unc_max_cm"] = nu_unc
sel["lambda_unc_max_mA"] = 1e8 * nu_unc / nu ** 2 * 1e3
sel["S_unc_code"] = sel["ierr"].str[1].astype(int)
sel["delta_unc_code"] = sel["ierr"].str[5].astype(int)
sel["dlam_pshift_1atm_mA"] = (1e8 / (nu + sel["dair_f"].values) - lv) * 1e3
sel["dlam_pshift_eff_mA_sl"] = (1e8 / (nu + sel["dair_f"].values * sel["pbar_sl"].values) - lv) * 1e3


def clean(s):
    return " ".join(s.split())


def band_label(row):
    m = MOLS[row.mol]
    Vu, Vl = clean(row.Vu), clean(row.Vl)
    if m == "O2":
        # HITRAN O2 global quanta: electronic state letter + v
        eu, vu = Vu.split()[0], Vu.split()[-1]
        el, vl = Vl.split()[0], Vl.split()[-1]
        lab = f"{eu}({vu})-{el}({vl})"
        names = {"b(0)-X(0)": "O2 A band b-X(0-0)", "b(1)-X(0)": "O2 B band b-X(1-0)",
                 "b(2)-X(0)": "O2 gamma band b-X(2-0)", "b(3)-X(0)": "O2 delta band b-X(3-0)",
                 "a(0)-X(0)": "O2 1.27um a-X(0-0)", "a(1)-X(0)": "O2 1.06um a-X(1-0)",
                 "a(0)-X(1)": "O2 1.58um a-X(0-1)", "b(0)-X(1)": "O2 b-X(0-1)", "b(1)-X(1)": "O2 b-X(1-1)"}
        return lab, names.get(lab, "O2 " + lab)
    if m == "H2O":
        try:
            u = [int(x) for x in Vu.replace("-", " -").split()]
            l_ = [int(x) for x in Vl.replace("-", " -").split()]
            if len(u) == 3 and len(l_) == 3 and min(u) >= 0 and min(l_) >= 0:
                P = lambda v: 2 * v[0] + v[1] + 2 * v[2]
                dP = P(u) - P(l_)
                return f"({u[0]}{u[1]}{u[2]})-({l_[0]}{l_[1]}{l_[2]})", f"H2O polyad dP={dP}"
        except ValueError:
            pass
        return f"{Vu} - {Vl}", "H2O (unassigned/other)"
    if m == "CO2":
        return f"{Vu.replace(' ', '')}-{Vl.replace(' ', '')}", None
    return f"{Vu} - {Vl}", None


labs = [band_label(r) for r in sel.itertuples()]
sel["band"] = [a for a, b in labs]
sel["band_group"] = [b for a, b in labs]
# spectral regions: per molecule, cluster selected lines in wavenumber (new region when the
# gap to the previous selected line exceeds 40 cm^-1); name = molecule + median lambda_vac (um)
sel = sel.sort_values(["mol", "nu_f"]).reset_index(drop=True)
reg = np.zeros(len(sel), int)
k = 0
for i in range(1, len(sel)):
    if sel["mol"].iat[i] != sel["mol"].iat[i - 1] or sel["nu_f"].iat[i] - sel["nu_f"].iat[i - 1] > 40.0:
        k += 1
    reg[i] = k
sel["reg_id"] = reg
rname = {}
for r, g in sel.groupby("reg_id"):
    rname[r] = f"{MOLS[g['mol'].iat[0]]} {np.median(g['lambda_vac_A']) / 1e4:.3f}um"
sel["region"] = sel["reg_id"].map(rname)
# for molecules without a conventional group name, group by spectral region
m = sel["band_group"].isna()
sel.loc[m, "band_group"] = sel.loc[m, "region"]
sel = sel.sort_values("nu_f", ascending=False).reset_index(drop=True)   # increasing wavelength

sel["molecule"] = sel["mol"].map(MOLS)
sel["isotopologue"] = [f"{MOLS[m_]}({molparam[(m_, i)]['afgl']})" for m_, i in zip(sel["mol"], sel["liso"])]

out = pd.DataFrame({
    "molecule": sel["molecule"],
    "isotopologue": sel["isotopologue"],
    "hitran_mol": sel["mol"],
    "hitran_iso": sel["liso"],
    "wavenumber_vac_cm-1": sel["nu"].str.strip(),
    "lambda_vac_A": sel["lambda_vac_A"].map(lambda v: f"{v:.4f}"),
    "lambda_air_A": sel["lambda_air_A"].map(lambda v: f"{v:.4f}"),
    "nu_unc_code": sel["nu_unc_code"],
    "lambda_unc_max_mA": sel["lambda_unc_max_mA"].map(lambda v: "" if np.isnan(v) else f"{v:.3g}"),
    "S296_cm-1/(molec.cm-2)": sel["S"].str.strip(),
    "S_unc_code": sel["S_unc_code"],
    "A_s-1": sel["A"].str.strip(),
    "E_lower_cm-1": sel["elow"].str.strip(),
    "gamma_air_cm-1/atm": sel["gair"].str.strip(),
    "gamma_self_cm-1/atm": sel["gself"].str.strip(),
    "n_air": sel["nair"].str.strip(),
    "delta_air_cm-1/atm": sel["dair"].str.strip(),
    "delta_unc_code": sel["delta_unc_code"],
    "dlam_pshift_1atm_mA": sel["dlam_pshift_1atm_mA"].map(lambda v: f"{v:.3f}"),
    "p_eff_atm_sl": sel["pbar_sl"].map(lambda v: f"{v:.3f}"),
    "dlam_pshift_eff_mA_sl": sel["dlam_pshift_eff_mA_sl"].map(lambda v: f"{v:.3f}"),
    "lambda_vac_obs_sl_A": (sel["lambda_vac_A"] + sel["dlam_pshift_eff_mA_sl"] / 1e3).map(lambda v: f"{v:.4f}"),
    "lambda_air_obs_sl_A": (sel["lambda_air_A"] + sel["dlam_pshift_eff_mA_sl"] / 1e3
                            / n_peck_reeder(sel["lambda_vac_A"])).map(lambda v: f"{v:.4f}"),
    "depth_R1e5_sl": sel["depthR_sl"].map(lambda v: f"{v:.4f}"),
    "depth_inf_sl": sel["depth_sl"].map(lambda v: f"{v:.4f}"),
    "EW_mA_sl": (sel["EW_cm_sl"] * sel["lambda_vac_A"] ** 2 / 1e8 * 1e3).map(lambda v: f"{v:.3f}"),
    "depth_R1e5_2p5km": sel["depthR_2p5km"].map(lambda v: f"{v:.4f}"),
    "depth_inf_2p5km": sel["depth_2p5km"].map(lambda v: f"{v:.4f}"),
    "EW_mA_2p5km": (sel["EW_cm_2p5km"] * sel["lambda_vac_A"] ** 2 / 1e8 * 1e3).map(lambda v: f"{v:.3f}"),
    "band": sel["band"],
    "band_group": sel["band_group"],
    "region": sel["region"],
    "upper_global_quanta": sel["Vu"].map(clean),
    "lower_global_quanta": sel["Vl"].map(clean),
    "upper_local_quanta": sel["Qu"].map(clean),
    "lower_local_quanta": sel["Ql"].map(clean),
    "g_upper": sel["gu"].str.strip(),
    "g_lower": sel["gl"].str.strip(),
    "source": "HITRAN2024",
})
os.makedirs(os.path.dirname(OUTCSV), exist_ok=True)
out.to_csv(OUTCSV, index=False, quoting=csv.QUOTE_MINIMAL)
print("wrote", OUTCSV, len(out), flush=True)

# ---------------------------------------------------------------- band summary
sel["EW_mA_sl"] = sel["EW_cm_sl"] * sel["lambda_vac_A"] ** 2 / 1e8 * 1e3
bands = []
for (mol, grp), g in sel.groupby(["molecule", "band_group"]):
    j = g["depthR_sl"].idxmax()
    bands.append(dict(molecule=mol, band_group=grp, n_lines=len(g),
                      n_lines_depthR_ge_0p1=int((g["depthR_sl"] >= 0.1).sum()),
                      n_lines_depthR_ge_0p5=int((g["depthR_sl"] >= 0.5).sum()),
                      n_lines_depthR_ge_0p01_2p5km=int((g["depthR_2p5km"] >= 0.01).sum()),
                      lambda_vac_min_A=f"{g['lambda_vac_A'].min():.3f}",
                      lambda_vac_max_A=f"{g['lambda_vac_A'].max():.3f}",
                      lambda_vac_median_A=f"{g['lambda_vac_A'].median():.1f}",
                      strongest_lambda_vac_A=f"{g.loc[j, 'lambda_vac_A']:.4f}",
                      strongest_lambda_air_A=f"{g.loc[j, 'lambda_air_A']:.4f}",
                      strongest_depthR_sl=f"{g.loc[j, 'depthR_sl']:.3f}",
                      sum_EW_mA_sl=f"{g['EW_mA_sl'].sum():.1f}",
                      median_dlam_pshift_1atm_mA=f"{g['dlam_pshift_1atm_mA'].median():.2f}",
                      median_dlam_pshift_eff_mA_sl=f"{g['dlam_pshift_eff_mA_sl'].median():.2f}",
                      max_abs_dlam_pshift_1atm_mA=f"{g['dlam_pshift_1atm_mA'].abs().max():.2f}",
                      median_lambda_unc_max_mA=f"{g['lambda_unc_max_mA'].median():.3g}",
                      frac_unc_le_1mA=f"{(g['lambda_unc_max_mA'] <= 1.0).mean():.3f}"))
bdf = pd.DataFrame(bands).sort_values(["molecule", "lambda_vac_median_A"], key=lambda s: s if s.name == "molecule" else s.astype(float))
bdf.to_csv(OUTBANDS, index=False)
print(bdf.to_string(), flush=True)

# ---------------------------------------------------------------- metadata / stats
stats = dict(
    n_records_downloaded_4000_33333=nrec_by_mol,
    n_prefilter_candidates=int(len(cand)),
    n_selected=int(len(sel)),
    n_selected_by_molecule=sel.groupby("molecule").size().to_dict(),
    n_selected_by_isotopologue=sel.groupby("isotopologue").size().to_dict(),
    n_infinite_res_depth_ge_0p01=int(len(pre)),
    ew_weak_line_check_median_ratio=ew_check,
    n_depthR_ge_0p1_sl=int((sel["depthR_sl"] >= 0.1).sum()),
    n_depthR_ge_0p5_sl=int((sel["depthR_sl"] >= 0.5).sum()),
    n_depthR_ge_0p01_2p5km=int((sel["depthR_2p5km"] >= 0.01).sum()),
    n_selected_by_band_group=sel.groupby("band_group").size().to_dict(),
    lambda_unc_max_mA_median_by_molecule=sel.groupby("molecule")["lambda_unc_max_mA"].median().to_dict(),
    frac_unc_le_1mA_by_molecule=sel.groupby("molecule")["lambda_unc_max_mA"].apply(lambda v: float((v <= 1.0).mean())).to_dict(),
    frac_unc_le_1mA_depthR_ge_0p1=float((sel.loc[sel["depthR_sl"] >= 0.1, "lambda_unc_max_mA"] <= 1.0).mean()),
    nu_unc_code_counts=sel["nu_unc_code"].value_counts().sort_index().to_dict(),
    frac_lambda_unc_le_1mA=float((sel["lambda_unc_max_mA"] <= 1.0).mean()),
    frac_lambda_unc_le_0p1mA=float((sel["lambda_unc_max_mA"] <= 0.1).mean()),
    frac_delta_air_zero=float((sel["dair_f"] == 0).mean()),
    pshift_1atm_mA=dict(median=float(sel["dlam_pshift_1atm_mA"].median()),
                        p5=float(sel["dlam_pshift_1atm_mA"].quantile(0.05)),
                        p95=float(sel["dlam_pshift_1atm_mA"].quantile(0.95)),
                        min=float(sel["dlam_pshift_1atm_mA"].min()),
                        max=float(sel["dlam_pshift_1atm_mA"].max())),
    pshift_eff_by_molecule_median_mA=sel.groupby("molecule")["dlam_pshift_eff_mA_sl"].median().to_dict(),
    p_eff_by_molecule_median_atm=sel.groupby("molecule")["pbar_sl"].median().to_dict(),
    air_formula_diff_mA=dict(
        ciddor450_minus_peckreeder=dict(min=float(dla.min() * 1e3), max=float(dla.max() * 1e3)),
        morton2000_minus_peckreeder=dict(min=float(dlm.min() * 1e3), max=float(dlm.max() * 1e3))),
    lambda_vac_range_A=[float(lv.min()), float(lv.max())],
)
# formula differences on a reference grid (independent of selection)
grid = np.array([3000, 4000, 5000, 6000, 7000, 7600, 8000, 9000, 10000, 12700, 15000, 16900, 20000, 25000.0])
stats["air_formula_diff_grid_mA"] = {
    f"{g:.0f}": dict(vac_minus_air_PR_A=round(float(g - g / n_peck_reeder(g)), 4),
                     ciddor450_minus_PR_mA=round(float((g / n_ciddor450(g) - g / n_peck_reeder(g)) * 1e3), 3),
                     morton2000_minus_PR_mA=round(float((g / n_morton2000(g) - g / n_peck_reeder(g)) * 1e3), 3))
    for g in grid}

meta = dict(
    title="Telluric absorption lines 3000-25000 A (HITRAN2024)",
    source="HITRAN2024 via https://hitran.org/lbl/api (par_line 160-char records)",
    selection=f"isolated-line central depth >= {DEPTH_MIN} after Gaussian LSF R={RES:.0f}, airmass 1, observer at sea level, "
              "AFGL-1986 U.S. Standard atmosphere (CO2, CH4, N2O profiles scaled to NOAA GML "
              f"{yco2} global annual means: CO2 {co2} ppm, CH4 {ch4} ppb, N2O {n2o} ppb), "
              "Voigt profile, line mixing and continua ignored",
    atmosphere=atm_summary,
    scale_factors=scale,
    columns={
        "wavenumber_vac_cm-1": "HITRAN transition wavenumber in vacuum at zero pressure, verbatim F12.6 string",
        "lambda_vac_A": "derived: 1e8 / wavenumber (vacuum, zero pressure), 4 decimals",
        "lambda_air_A": "derived: lambda_vac / n_air, n from Peck & Reeder 1972 standard air (15 C, 101325 Pa, dry); "
                        "formula fitted 0.185-1.69 um, extrapolated beyond 1.69 um",
        "nu_unc_code": "HITRAN position uncertainty code (0 unreported; k: 10^-k <= err < 10^-(k-1) cm^-1)",
        "lambda_unc_max_mA": "derived upper bound of the position uncertainty from nu_unc_code, in mA",
        "S296_cm-1/(molec.cm-2)": "HITRAN line intensity at 296 K, weighted by terrestrial isotopic abundance",
        "S_unc_code": "HITRAN relative uncertainty code for S (3: >=20%, 4: 10-20%, 5: 5-10%, 6: 2-5%, 7: 1-2%, 8: <1%...)",
        "A_s-1": "Einstein A coefficient",
        "E_lower_cm-1": "lower-state energy",
        "gamma_air_cm-1/atm": "air-broadened HWHM at 296 K, 1 atm",
        "gamma_self_cm-1/atm": "self-broadened HWHM at 296 K, 1 atm",
        "n_air": "temperature exponent of gamma_air",
        "delta_air_cm-1/atm": "air pressure shift at 296 K (nu* = nu + delta_air*p); 0.000000 often means unknown",
        "delta_unc_code": "HITRAN absolute uncertainty code for delta_air (0 = unreported)",
        "dlam_pshift_1atm_mA": "derived: vacuum-wavelength shift (mA) produced by delta_air at p = 1 atm",
        "p_eff_atm_sl": "derived: S(T)*column-weighted mean pressure of the absorbing column (sea-level model)",
        "dlam_pshift_eff_mA_sl": "derived: wavelength shift (mA) of an unsaturated line seen through the sea-level model atmosphere",
        "lambda_vac_obs_sl_A": "derived: lambda_vac_A + dlam_pshift_eff_mA_sl/1000 = expected (pressure-shifted) vacuum position seen from sea level",
        "lambda_air_obs_sl_A": "derived: the same converted to standard air (Peck & Reeder 1972)",
        "depth_R1e5_sl": "derived: central depth of the isolated line after convolution with a Gaussian LSF of FWHM lambda/100000, airmass 1, sea level (selection quantity)",
        "depth_inf_sl": "derived: central depth 1-exp(-tau0) of the isolated line, airmass 1, sea level, infinite resolution",
        "EW_mA_sl": "derived: equivalent width (mA) of the isolated line, airmass 1, sea level",
        "depth_R1e5_2p5km": "derived: as depth_R1e5_sl but observer at 2.5 km altitude (PWV ~4.4 mm)",
        "depth_inf_2p5km": "derived: as depth_inf_sl but observer at 2.5 km altitude",
        "region": "derived: contiguous spectral region of selected lines of the molecule (gap > 40 cm^-1 starts a new one), named by median lambda_vac in um",
        "EW_mA_2p5km": "derived: same as EW_mA_sl but observer at 2.5 km altitude",
        "band": "vibrational/electronic band from HITRAN global quanta (upper-lower)",
        "band_group": "band grouping (O2 electronic bands; H2O polyad change dP=2v1+v2+2v3; others by band)",
        "upper/lower_global/local_quanta": "HITRAN quantum-number strings (whitespace collapsed)",
    },
    stats=stats,
    citations=[
        "Gordon, I.E., Rothman, L.S., Hargreaves, R.J., Gomez, F.M., Bertin, T., Hill, C., et al. 2026, "
        "The HITRAN2024 molecular spectroscopic database, JQSRT 353, 109807, doi:10.1016/j.jqsrt.2026.109807 "
        "(data: https://hitran.org, line-by-line API https://hitran.org/lbl/api, downloaded 2026-11-23)",
        "Kochanov, R.V., Gordon, I.E., Rothman, L.S., Wcislo, P., Hill, C., Wilzewski, J.S. 2016, HAPI, "
        "JQSRT 177, 15, doi:10.1016/j.jqsrt.2016.03.005",
        "Gamache, R.R., Vispoel, B., Tennyson, J., Yurchenko, S.N., Polyansky, O.L., Gordon, I.E., et al. 2025, "
        "Partition sums for molecules and their isotopologues for HITRAN2024 (TIPS-2025), JQSRT 345, 109568, doi:10.1016/j.jqsrt.2025.109568 "
        "(Q(T) files https://hitran.org/data/Q/q<id>.txt)",
        "Anderson, G.P., Clough, S.A., Kneizys, F.X., Chetwynd, J.H., Shettle, E.P. 1986, AFGL Atmospheric "
        "Constituent Profiles (0-120 km), AFGL-TR-86-0110, ADS 1986afgl.rept.....A (tables as packaged in "
        "joseki 2.7.0, https://pypi.org/project/joseki/, data/afgl_1986/table_1f.csv, table_2a.csv)",
        "NOAA Global Monitoring Laboratory, global annual mean CO2/CH4/N2O (Lan, X., Tans, P., Thoning, K.W., "
        "et al.), https://gml.noaa.gov/webdata/ccgg/trends/{co2,ch4,n2o}/*_annmean_gl.txt (downloaded 2026-11-23)",
        "Peck, E.R., Reeder, K. 1972, Dispersion of air, JOSA 62, 958, doi:10.1364/JOSA.62.000958",
        "Ciddor, P.E. 1996, Refractive index of air: new equations for the visible and near infrared, "
        "Appl. Opt. 35, 1566, doi:10.1364/AO.35.001566 (comparison only)",
        "Morton, D.C. 2000, ApJS 130, 403, doi:10.1086/317349 (IAU standard formula after Birch & Downs 1994; comparison only)",
    ],
)
json.dump(meta, open(OUTMETA, "w"), indent=1, default=float)
print(json.dumps(stats, indent=1, default=float))
