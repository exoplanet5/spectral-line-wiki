"""Validate decoding of PALACE OH IDs against Keck/HIRES (Osterbrock+1996/97, VizieR III/211) OH labels.
Converts PALACE vacuum wavelengths to air with Edlen (1966) as printed in Hanuschik (2003) eq. 5,
then compares with the Osterbrock air wavelengths (Abrams+1994 based, Edlen 1953 air).
Usage: python -I validate_palace_vs_hires.py <pmd_reflineint_mol.dat> <III_211_oh.tsv>
"""
import sys, re
import numpy as np

def n_edlen66(lam_vac_A):
    s2 = (1e4 / np.asarray(lam_vac_A, float)) ** 2   # sigma in um^-1 (vacuum wavenumber)
    return 1 + 1e-8 * (8342.13 + 2406030.0 / (130.0 - s2) + 15997.0 / (38.9 - s2))

def vac2air(lv):
    return np.asarray(lv, float) / n_edlen66(lv)

DJ = {'O': -2, 'P': -1, 'Q': 0, 'R': 1, 'S': 2}

def oh_label(idc, parity_mode='upper'):
    vu = 10 if idc[4] == 'X' else int(idc[4]); vl = int(idc[5])
    br = idc[6:8]; Fu = int(idc[8]); Fl = int(idc[9]); Nu = int(idc[10:12]); par = idc[12]
    Ju = Nu + 0.5 if Fu == 1 else Nu - 0.5
    Jl = Ju - DJ[br[1]]
    if parity_mode == 'lower' and br[1] == 'Q':
        par = {'e': 'f', 'f': 'e'}[par]
    if Fu == Fl and br[0] == br[1]:
        lab = f"{br[1]}{Fu}{par}({Jl:.1f})"
    else:
        lab = f"{br[0]}{br[1]}{Fu}{Fl}{par}({Jl:.1f})"
    return vu, vl, lab

pal = {}
for line in open(sys.argv[1]):
    if line.startswith('#') or line.strip().startswith('lam'):
        continue
    p = line.split()
    if len(p) < 11 or not p[4].startswith('OHXX'):
        continue
    lv = float(p[0]) * 1e4
    for mode in ('upper', 'lower'):
        vu, vl, lab = oh_label(p[4], mode)
        pal.setdefault(mode, {})[(vu, vl, lab)] = (lv, p[0], p[4])

hires = []
for line in open(sys.argv[2]):
    if line.startswith('#') or not line.strip():
        continue
    q = line.rstrip('\n').split('\t')
    if len(q) >= 5 and q[0].strip().isdigit():
        hires.append((q[1].strip(), int(q[2]), int(q[3]), q[4].strip()))

for mode in ('upper', 'lower'):
    d = []; miss = []
    for wl, vu, vl, tr in hires:
        key = (vu, vl, tr)
        if key in pal[mode]:
            lv = pal[mode][key][0]
            d.append(float(wl) - float(vac2air(lv)))
        else:
            miss.append((wl, vu, vl, tr))
    d = np.array(d)
    print(f"parity={mode}: matched {len(d)}/{len(hires)}; HIRES(air) - PALACE(vac->air Edlen66): "
          f"median {np.median(d)*1e3:.2f} mA, mean {d.mean()*1e3:.2f} mA, std {d.std()*1e3:.2f} mA, "
          f"|d|>5mA: {(abs(d)>0.005).sum()}, |d|>20mA: {(abs(d)>0.02).sum()}")
    print('   first misses:', miss[:8])
