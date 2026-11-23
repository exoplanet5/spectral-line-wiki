"""Build the curated optical airglow line tables for the Spectral Line Wiki.

Inputs (all downloaded, treated as untrusted data; run with `python -I`):
  --viz    data/raw/airglow_optical/vizier_J_AA_407_1157   (Hanuschik 2003, CDS tables 4-9, asu-tsv)
  --eso    data/raw/airglow_optical/eso_uves/ident_extracted (ESO QC web gident_*.dat line tables)
  --hires  data/raw/airglow_optical/vizier_III_211          (Osterbrock+1996/97 Keck/HIRES, VizieR III/211)
  --palace data/raw/airglow_optical/palace_zenodo/PMD_extracted/PMD (Noll+2025 PALACE v1.0 data release)
  --out    data/curated

Outputs:
  airglow_optical.csv                  UVES/VLT observed night-sky emission lines (Hanuschik 2003),
                                       all settings merged, with derived identifications
  airglow_optical_hires_oh.csv         Keck/HIRES OH/O2 lines (Osterbrock+1996/97) with lab-based air wavelengths
  airglow_model_palace_3000_10500.csv  PALACE v1.0 model airglow line list 3000-10500 A (vacuum)
  airglow_model_palace_atomic_ir.csv   PALACE v1.0 atomic lines above 10500 A (all O I recombination lines, to 2.14 um)
  airglow_optical_qa.json              cross-match statistics

Formulas used (documented, derived quantities are flagged as such in the CSV headers):
  * Air <-> vacuum: Edlen (1966) standard-air dispersion as printed in Hanuschik (2003) eq. 5:
      n = 1 + 1e-8 * (8342.13 + 2406030/(130 - s^2) + 15997/(38.9 - s^2)),  s = 1/lambda_vac [um^-1]
      lambda_air = lambda_vac / n ; air->vac by fixed-point iteration lambda_vac = lambda_air * n(lambda_vac).
  * Flux -> aeronomic units: F [1e-16 erg s^-1 cm^-2 A^-1 arcsec^-2] -> R/nm with
      R/nm = F*1e-16 * lambda/(h c) * (206264.806...)^2 * 4 pi * 1e-6 * 10   (CODATA via astropy.constants)
    (checked: 0.5313 at 3550.3833 A -> 50.77 R/nm = Cosby+2006 Table 3 value).
"""
import argparse, csv, json, math, os, re
from collections import defaultdict
import numpy as np
from astropy import constants as const

# ----------------------------------------------------------------------------- helpers
def n_edlen66(lam_vac_A):
    s2 = (1e4 / np.asarray(lam_vac_A, float)) ** 2
    return 1.0 + 1e-8 * (8342.13 + 2406030.0 / (130.0 - s2) + 15997.0 / (38.9 - s2))

def vac2air(lv):
    lv = np.asarray(lv, float)
    return lv / n_edlen66(lv)

def air2vac(la):
    la = np.asarray(la, float)
    lv = la * n_edlen66(la)
    for _ in range(6):
        lv = la * n_edlen66(lv)
    return lv

ARCSEC2_PER_SR = (180.0 / math.pi * 3600.0) ** 2
HC_ERG_CM = (const.h * const.c).cgs.value

def flux_to_R_per_nm(F, lam_A):
    """F in 1e-16 erg/s/cm2/A/arcsec2 -> Rayleigh per nm."""
    e_ph = HC_ERG_CM / (lam_A * 1e-8)
    ph_per_A_sr = F * 1e-16 / e_ph * ARCSEC2_PER_SR
    return 4 * math.pi * 1e-6 * ph_per_A_sr * 10.0

def fmt_sig(x, n=4):
    return '' if x is None or not np.isfinite(x) else f"{x:.{n}g}"

# ----------------------------------------------------------------------------- UVES inputs
SETTINGS = {  # setting: (CDS table, nominal range from Hanuschik 2003 Table 3 / ESO QC page, typical continuum 1e-16 cgs)
    '346': ('table4', 3140.0, 3760.0, 0.17),
    '437': ('table5', 3740.0, 4860.0, 0.14),
    '580L': ('table6', 4800.0, 5770.0, 0.09),
    '580U': ('table7', 5830.0, 6790.0, 0.10),
    '860L': ('table8', 6700.0, 8560.0, 0.08),
    '860U': ('table9', 8600.0, 10430.0, 0.07),
    '800U': (None, 8530.0, 8630.0, 0.07),   # ESO web only (added 2003); continuum value of 860U assumed for the rough ratio only
}

def read_vizier_tsv(path, mincols=5):
    rows = []
    for line in open(path, encoding='utf-8', errors='replace'):
        if line.startswith('#') or not line.strip():
            continue
        p = line.rstrip('\n').split('\t')
        if len(p) >= mincols and p[0].strip().isdigit():
            rows.append([x.strip() for x in p])
    return rows

def read_nist_hei(path):
    """NIST ASD CSV (format=2). Keep He I lines whose lower level is the metastable 1s.2s 3S (Cosby+2006 Table 1)."""
    out = []
    with open(path, encoding='utf-8', errors='replace') as f:
        rd = csv.reader(f)
        hdr = next(rd)
        def clean(x):
            return x.strip().lstrip('=').strip('"')
        for row in rd:
            if len(row) < len(hdr) - 1:
                continue
            d = {h: clean(v) for h, v in zip(hdr, row)}
            if d.get('conf_i') == '1s.2s' and d.get('term_i') == '3S':
                wl = d.get('ritz_wl_air(A)') or d.get('obs_wl_air(A)')
                if wl:
                    out.append(dict(wl=wl, unc=d.get('unc_ritz_wl', ''), lab=f"He I 1s2s 3S{d['J_i']} - {d['conf_k'].replace('.', '')} {d['term_k']}{d['J_k']}"))
    return out

ESO_RE = re.compile(r'^\s+(\d+)\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$')

def read_eso_dat(path):
    rows = []
    for line in open(path, encoding='ascii', errors='replace'):
        m = ESO_RE.match(line)
        if m:
            rows.append(list(m.groups()))
    return rows

# ----------------------------------------------------------------------------- PALACE decoding
DN = {'N': -3, 'O': -2, 'P': -1, 'Q': 0, 'R': 1, 'S': 2, 'T': 3}
ISO = {'OH': 'OH', 'QH': '18OH', 'O2': 'O2', 'OG': '16O17O', 'OQ': '16O18O'}
ESTATE = {'X': 'X', 'A': 'A', 'a': 'a', 'b': 'b', 'c': 'c'}

def vib(c):
    return 10 + 'ABCDEF'.index(c) if c in 'ABCDEF' else (10 if c == 'X' else int(c))

def decode_mol(idc):
    """Return (species, system, band, label) from a 13-character PALACE molecular ID."""
    iso = ISO.get(idc[0:2], idc[0:2])
    up, lo = idc[2], idc[3]
    vu, vl = vib(idc[4]), vib(idc[5])
    br = idc[6:8]; Fu = int(idc[8]); Fl = int(idc[9]); Nu = int(idc[10:12]); tag = idc[12]
    if iso in ('OH', '18OH'):
        Ju = Nu + 0.5 if Fu == 1 else Nu - 0.5
        Jl = Ju - DN[br[1]]
        if Fu == Fl and br[0] == br[1]:
            lab = f"{br[1]}{Fu}{tag}({Jl:.1f})"
        else:
            lab = f"{br[0]}{br[1]}{Fu}{Fl}{tag}({Jl:.1f})"
        return iso, 'OH Meinel X-X', f"{vu}-{vl}", f"{iso} {vu}-{vl} {lab}", (iso, vu, vl, br, Fu, Fl, Nu)
    upn = "A'" if (up == 'A' and lo == 'a') else up
    system = {('b', 'X'): 'O2 b-X Atmospheric', ('A', 'X'): 'O2 A-X Herzberg I', ('A', 'a'): "O2 A'-a Chamberlain",
              ('c', 'X'): 'O2 c-X Herzberg II', ('c', 'b'): 'O2 c-b', ('a', 'X'): 'O2 a-X IR Atmospheric'}.get((up, lo), f'O2 {up}-{lo}')
    Nl = Nu - DN[br[0]]
    # J'' : X 3Sigma_g-  F1: J=N+1, F2: J=N, F3: J=N-1 ; singlet lower states (a, b): J=N
    Jl = Nl + {1: 1, 2: 0, 3: -1}[Fl] if lo == 'X' else Nl
    # Cosby+2006 style: DeltaJ-branch F'F''(J''), followed by [DeltaN DeltaJ, N''] as in PALACE
    lab = f"{br[1]}{Fu}{Fl}({Jl:.1f}) [{br} N''={Nl}]" + (' quadrupole' if tag == 'q' else '')
    return iso, system, f"{vu}-{vl}", f"{iso} {upn}-{lo}({vu}-{vl}) {lab}", None

def decode_atom(idc):
    sp = idc[0:2].replace('_', '')
    cu, tu, cl, tl, ju, jl = idc[2:4], idc[4:6], idc[6:8], idc[8:10], idc[10], idc[11]
    def jstr(c):
        v = int(c)
        if sp == 'O':
            return str(v)
        return f"{v}/2" if v % 2 else str(v // 2)
    return sp, f"{sp} I", f"{sp} I {cu} {tu}{jstr(ju)} - {cl} {tl}{jstr(jl)}"

def read_palace(pmd):
    lines = []
    for fname, kind in (('pmd_reflineint_mol.dat', 'mol'), ('pmd_reflineint_atom.dat', 'atom')):
        for line in open(os.path.join(pmd, fname), encoding='utf-8', errors='replace'):
            if line.startswith('#') or line.strip().startswith('lam') or not line.strip():
                continue
            p = line.split()
            if len(p) < 11:
                continue
            lam_um = p[0]
            lv = float(lam_um) * 1e4
            if kind == 'mol':
                sp, system, band, label, ohkey = decode_mol(p[4])
            else:
                sp, system, label = decode_atom(p[4]); band = ''; ohkey = None
            lines.append(dict(lam_um=lam_um, lv=lv, trans=float(p[1]), fH2O=p[2], I=float(p[3]), I_str=p[3], id=p[4],
                              chem=p[5], Tkin=p[7], hlay=p[8], src=p[9], vclass=p[10], kind=kind,
                              species=sp, system=system, band=band, label=label, ohkey=ohkey))
    la = vac2air([d['lv'] for d in lines])
    for d, a in zip(lines, la):
        d['la'] = float(a)
    lines.sort(key=lambda d: d['la'])
    return lines

def read_hires(hdir):
    oh = []
    for r in read_vizier_tsv(os.path.join(hdir, 'oh.tsv')):
        oh.append(dict(wl=r[1], vu=int(r[2]), vl=int(r[3]), trans=r[4]))
    o2 = []
    for r in read_vizier_tsv(os.path.join(hdir, 'o2.tsv'), mincols=3):
        o2.append(dict(wl=r[1], trans=r[2]))
    return oh, o2

# ----------------------------------------------------------------------------- identification (Cosby+2006 criterion)
def identify(lam, fwhm, pal_la, pal, Imin=1e-3):
    w = 0.5 * max(lam / 40000.0, fwhm)
    i0 = np.searchsorted(pal_la, lam - w); i1 = np.searchsorted(pal_la, lam + w)
    cands = [pal[i] for i in range(i0, i1) if pal[i]['I'] * pal[i]['trans'] >= Imin]
    if not cands:
        return w, 0, 'U', []
    # group OH Lambda doublets (identical ID except the e/f character)
    units = defaultdict(list)
    for c in cands:
        key = c['id'][:12] if c['chem'] == 'OH' else c['id']
        units[key].append(c)
    ulist = sorted(units.values(), key=lambda u: -sum(c['I'] * c['trans'] for c in u))
    ui = [sum(c['I'] * c['trans'] for c in u) for u in ulist]
    if len(ulist) == 1 or ui[0] >= 5 * ui[1]:
        sel = [ulist[0]]
        flag = 'L2' if len(ulist[0]) == 2 else 'S1'
    else:
        tot = sum(ui)
        sel = [u for u, x in zip(ulist, ui) if x >= 0.1 * tot][:4]
        ntr = sum(len(u) for u in sel)
        if len(sel) == 1:
            flag = 'L2' if len(sel[0]) == 2 else 'S1'
        elif ntr == 2 and len(sel) == 2 and ui[0] < 3 * ui[1]:
            flag = 'B2'
        else:
            flag = f"M{min(ntr, 9)}"
    return w, len(cands), flag, sel

def unit_label(u):
    if len(u) == 2 and u[0]['chem'] == 'OH':
        a = u[0]['label']
        # merge e/f labels: "OH 6-2 P1e(8.5)" + "OH 6-2 P1f(8.5)" -> "OH 6-2 P1e/f(8.5)"
        m = re.match(r'^(.*\d)([ef])(\(.*\))$', a)
        return f"{m.group(1)}e/f{m.group(3)}" if m else a + '+' + u[1]['label']
    return u[0]['label']

# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--viz', required=True); ap.add_argument('--eso', required=True)
    ap.add_argument('--hires', required=True); ap.add_argument('--palace', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--nist-hei', required=True, help='NIST ASD He I CSV (format=2, air)')
    a = ap.parse_args()

    pal = read_palace(a.palace)
    pal_la = np.array([d['la'] for d in pal])
    hires_oh, hires_o2 = read_hires(a.hires)
    hei = read_nist_hei(a.nist_hei)

    # ---- UVES rows: CDS (A&A 2003) version first, then ESO-web-only lines
    rows = []
    for s, (tab, lo, hi, cont) in SETTINGS.items():
        eso = read_eso_dat(os.path.join(a.eso, f'gident_{s}.dat'))
        eso_w = np.array([float(r[1]) for r in eso]) if eso else np.array([])
        used = set()
        if tab:
            for r in read_vizier_tsv(os.path.join(a.viz, f'{tab}.tsv')):
                wl = float(r[2])
                # match the ESO-web row: same running number first, else nearest unused row within 6e-4 A
                j = next((k for k, x in enumerate(eso) if x[0] == r[1] and abs(eso_w[k] - wl) < 6e-4 and k not in used), -1)
                if j < 0 and len(eso_w):
                    order = np.argsort(abs(eso_w - wl))
                    j = next((int(k) for k in order if int(k) not in used and abs(eso_w[k] - wl) < 6e-4), -1)
                e = eso[j] if j >= 0 else None
                if e is not None:
                    used.add(j)
                rows.append(dict(setting=s, src='Hanuschik2003_CDS_' + tab, idx=r[1], wl=r[2], fwhm=r[3], flux=r[4],
                                 eso_wl=e[1] if e else '', eso_int=e[2] if e else '', cds=True))
        for j, e in enumerate(eso):
            if j in used:
                continue
            rows.append(dict(setting=s, src=f'Hanuschik2003_ESOweb_gident_{s}', idx=e[0], wl=e[1], fwhm=e[3], flux=e[4],
                             eso_wl=e[1], eso_int=e[2], cds=False))
    for r in rows:
        r['lam'] = float(r['wl']); r['fw'] = float(r['fwhm']); r['F'] = float(r['flux'])
    rows.sort(key=lambda r: (r['lam'], r['setting']))

    # ---- duplicates across overlapping settings (mutual nearest neighbours within the Cosby window)
    by_set = defaultdict(list)
    for k, r in enumerate(rows):
        by_set[r['setting']].append(k)
    gid = 0
    for r in rows:
        r['dup'] = ''; r['primary'] = 1; r['dupnote'] = ''
    sets = list(SETTINGS)
    for ia in range(len(sets)):
        for ib in range(ia + 1, len(sets)):
            A, B = sets[ia], sets[ib]
            loA, hiA = SETTINGS[A][1:3]; loB, hiB = SETTINGS[B][1:3]
            ka = [k for k in by_set[A]]; kb = [k for k in by_set[B]]
            if not ka or not kb:
                continue
            wa = np.array([rows[k]['lam'] for k in ka]); wb = np.array([rows[k]['lam'] for k in kb])
            if wa.max() < wb.min() - 1 or wb.max() < wa.min() - 1:
                continue
            for i, k in enumerate(ka):
                j = int(np.argmin(abs(wb - wa[i])))
                if int(np.argmin(abs(wa - wb[j]))) != i:
                    continue
                r1, r2 = rows[k], rows[kb[j]]
                tol = 0.5 * max(r1['lam'] / 40000.0, r1['fw'], r2['fw'])
                if abs(wa[i] - wb[j]) >= tol:
                    continue
                g = r1['dup'] or r2['dup']
                if not g:
                    gid += 1; g = f"D{gid:03d}"
                r1['dup'] = r2['dup'] = g
    # entries listed twice in the same source table (identical wavelength string, slightly different FWHM/flux)
    srcdup = defaultdict(list)
    for k, r in enumerate(rows):
        srcdup[(r['src'], r['wl'])].append(k)
    for key, ks in srcdup.items():
        if len(ks) > 1:
            gid += 1; g = f"S{gid:03d}"
            for n, k in enumerate(ks):
                rows[k]['srcdup'] = g
                rows[k]['primary_srcdup'] = 1 if n == 0 else 0
    groups = defaultdict(list)
    for k, r in enumerate(rows):
        if r['dup']:
            groups[r['dup']].append(k)
    for g, ks in groups.items():
        def score(k):
            r = rows[k]; lo, hi = SETTINGS[r['setting']][1:3]
            return (1 if r['cds'] else 0, min(r['lam'] - lo, hi - r['lam']))
        best = max(ks, key=score)
        desc = '; '.join(f"{rows[k]['setting']}:{rows[k]['wl']}" for k in ks)
        for k in ks:
            rows[k]['primary'] = 1 if k == best else 0
            rows[k]['dupnote'] = (f"same feature measured in overlapping settings ({desc}); primary = "
                                  f"{rows[best]['setting']} (CDS version preferred, then larger distance from setting edge)")

    for r in rows:
        if r.get('srcdup'):
            note = (f"entry listed twice in the source table at the same wavelength (group {r['srcdup']}); "
                    + ('first entry kept as primary' if r['primary_srcdup'] else 'second entry marked non-primary'))
            r['dupnote'] = (note + '; ' + r['dupnote']).strip('; ')
            if not r['primary_srcdup']:
                r['primary'] = 0
            r['dup'] = r['dup'] or r['srcdup']
    # ---- identifications
    qa = defaultdict(list)
    for r in rows:
        lam, fw = r['lam'], r['fw']
        w, ncand, flag, sel = identify(lam, fw, pal_la, pal)
        r['w'] = w; r['ncand'] = ncand
        # rough observed integrated intensity (R): continuum-subtracted peak (typical continuum per setting,
        # Hanuschik 2003 Table 3) x FWHM x 1.0645 ; floor at 20% of the peak; only used for the ratio flag
        cont = SETTINGS[r['setting']][3]
        peak = max(r['F'] - cont, 0.2 * r['F'])
        Iobs = flux_to_R_per_nm(peak, lam) * 0.1 * fw * 1.0645
        r['R_per_nm'] = flux_to_R_per_nm(r['F'], lam)
        if sel:
            comps = [c for u in sel for c in u]
            Imod = sum(c['I'] * c['trans'] for c in comps)
            ratio = Imod / Iobs if Iobs > 0 else float('nan')
            if ratio < 0.1:
                flag += 'w'
            top = sel[0]
            top_la = sum(c['la'] for c in top) / len(top)
            r.update(flag=flag, ident=' + '.join(unit_label(u) for u in sel),
                     species=top[0]['species'], system=top[0]['system'],
                     mod_la=f"{top_la:.3f}", dlam=f"{lam - top_la:+.4f}", Imod=fmt_sig(Imod, 3), ratio=fmt_sig(ratio, 3),
                     comps='; '.join(f"{c['id']}|vac={float(c['lam_um'])*1e4:.3f}|air={c['la']:.3f}|I={c['I_str']}|T={c['trans']:.3f}" for c in comps))
            r['method'] = 'xmatch:PALACE_v1.0 lines, Cosby+2006 window/labels (derived)'
            r['top_label'] = unit_label(top)
            if flag in ('S1', 'L2'):
                qa[(top[0]['species'], r['setting'])].append(lam - top_la)
        else:
            r.update(flag='U', ident='', species='', system='', mod_la='', dlam='', Imod='', ratio='', comps='', top_label='',
                     method='xmatch:PALACE_v1.0 lines, Cosby+2006 window/labels (derived)')
            hh = [h for h in hei if abs(lam - float(h['wl'])) < w]
            if hh:
                mla = sum(float(h['wl']) for h in hh) / len(hh)
                r.update(flag='S1', ident=hh[0]['lab'].rsplit(' ', 1)[0] + ' ' + hh[0]['lab'].rsplit(' ', 1)[1][:2] + ' multiplet (NIST)',
                         species='He', system='He I 2s 3S - np 3P (metastable)', mod_la=f"{mla:.3f}", dlam=f"{lam - mla:+.4f}",
                         comps='; '.join(f"{h['lab']}|air={h['wl']}|unc={h['unc']}" for h in hh),
                         method='xmatch:NIST ASD He I lines with lower level 1s2s 3S (Ritz, air), Cosby+2006 window (derived; no model intensity)')
        # HIRES matches
        hm = []
        for h in hires_oh:
            d = lam - float(h['wl'])
            if abs(d) < w:
                hm.append((abs(d), h['wl'], f"OH {h['vu']}-{h['vl']} {h['trans']}", d))
        for h in hires_o2:
            if h['wl']:
                d = lam - float(h['wl'])
                if abs(d) < w:
                    hm.append((abs(d), h['wl'], f"O2 b-X(0-1) {h['trans'].strip()}", d))
        hm.sort()
        r['h_wl'] = ';'.join(x[1] for x in hm); r['h_id'] = ';'.join(x[2] for x in hm)
        r['h_d'] = ';'.join(f"{x[3]:+.4f}" for x in hm)
        if hm and r['flag'] in ('S1', 'L2') and r['species'] == 'OH':
            norm = lambda t: re.sub(r'(\d)(e/f|e|f)\(', r'\1(', t)
            same = [float(x[1]) for x in hm if norm(x[2]) == norm(r['top_label'])]
            if same:
                qa[('UVES-HIRES_sameOHtransition', r['setting'])].append(lam - sum(same) / len(same))

    # ---- write airglow_optical.csv
    os.makedirs(a.out, exist_ok=True)
    out1 = os.path.join(a.out, 'airglow_optical.csv')
    hdr = [
        '# Spectral Line Wiki -- optical night-sky (airglow) emission lines, 3140-10430 A, UVES/VLT Paranal atlas.',
        '# Observed lines: Hanuschik R.W. 2003, A&A 407, 1157 (2003A&A...407.1157H, doi:10.1051/0004-6361:20030885);',
        '#   CDS/VizieR J/A+A/407/1157 tables 4-9 (doi:10.26093/cds/vizier.34071157, fetched via vizier.cfa.harvard.edu)',
        '#   plus lines present only in the later ESO QC web line tables (gident_*.dat, https://www.eso.org/qc/UVES/txt/sky/UVES_ident.tar.Z).',
        '# wl_air_A: verbatim source string, air wavelength (UVES ThAr-calibrated pipeline scale). CDS strings have 4 decimals;',
        '#   ESO-web-only rows keep the ESO .dat string (6 decimals printed from FITS float32 -> digits beyond ~1e-3 A are not significant).',
        '# uncert_A: 0.017 = author-stated likely error of a single line position (Hanuschik 2003 sect. 5.2: UVES-Keck <d> = +8 +- 17 mA,',
        '#   N=529 lines >5200 A). Internal dispersion-solution rms 0.0047 A (346 nm setting) to 0.010 A (860 nm setting).',
        '#   Known systematic: UVES - NIST = +0.0013 +- 0.0002 nm for 13 atomic lines (Cosby et al. 2006 Table 5); UVES - HIRES +0.00116 nm.',
        '# flux: verbatim FLUX column = Gaussian-fit central (peak) flux density INCLUDING the local continuum,',
        '#   unit 1e-16 erg s-1 cm-2 A-1 arcsec-2 (ESO QC page; VizieR lists 1e-16 erg/s/cm2/A), extinction-corrected.',
        '# flux_peak_R_per_nm_derived: flux converted to Rayleigh/nm (derived; R/nm = F*1e-16*lambda/(hc)*206264.806^2*4pi*1e-6*10).',
        '# fwhm_A: verbatim FWHM of the fitted Gaussian (A). R ~ 43000-45000 (1 arcsec slit).',
        '# wl_vac_A: blank (not given by the source). wl_vac_A_derived: DERIVED from wl_air_A with Edlen (1966) as printed in',
        '#   Hanuschik 2003 eq. 5, n = 1 + 1e-8*(8342.13 + 2406030/(130-s^2) + 15997/(38.9-s^2)), s = 1/lambda_vac [um^-1], iterated.',
        '# ident: DERIVED identification (NOT the Cosby et al. 2006 published list, whose machine-readable supplement could not be',
        '#   retrieved). Each line was cross-matched to the PALACE v1.0 model line list (Noll et al. 2025, GMD 18, 4353;',
        '#   doi:10.5281/zenodo.14064023; positions from HITRAN2020 for OH and O2 b-X/a-X, Cosby+2006 models for O2 A-X/A\'-a/c-X/c-b,',
        '#   NIST for atoms; vacuum -> air with Edlen 1966) using the Cosby+2006 criterion |lam_obs-lam_model| < 0.5*max(lam/40000, FWHM),',
        '#   model intensity weighted by PALACE Paranal transmission, I*T >= 1e-3 R. Labels (after Cosby+2006): S1 single (or >=5x dominant);',
        '#   L2 OH Lambda doublet pair; B2 two transitions within a factor 3; Mn n blended transitions (each >=10% of the total, <=4 units);',
        '#   U no model line in window; suffix w = assigned model intensity <0.1 x rough observed intensity (tentative).',
        '#   OH labels: band v\'-v\'\', branch (DeltaN DeltaJ for satellites) F\'F\'\', e/f = parity of upper level (PALACE/HITRAN), (J\'\').',
        '#   O2 labels: system(v\'-v\'\') branch DeltaN DeltaJ F\'F\'\' (N\'\'). ident_model_wl_air_A = mean model air wavelength of the strongest unit.',
        '# hires_*: Keck/HIRES lines (Osterbrock et al. 1996 PASP 108, 277; 1997 PASP 109, 614; VizieR III/211) within the same window;',
        '#   their air wavelengths are laboratory/calculated (Abrams et al. 1994 ApJS 93, 351; Edlen 1953 air), not sky-measured.',
        '# dup_group/is_primary: the same feature measured in two overlapping settings; both rows kept, one marked primary.',
        '# Lines in 5770-5830 A and 8560-8600 A are missing or sparse (CCD mosaic gaps; 800U fills part of the latter).',
    ]
    cols = ['wl_air_A', 'wl_vac_A', 'medium_of_source', 'uncert_A', 'flux', 'fwhm_A', 'ident', 'source_key',
            'setting', 'src_index', 'wl_air_A_eso_dat', 'int_counts_eso', 'flux_peak_R_per_nm_derived', 'wl_vac_A_derived',
            'species_main', 'system_main', 'ident_flag', 'ident_window_A', 'ident_n_candidates', 'ident_model_wl_air_A',
            'ident_dlam_A', 'ident_model_I_R', 'ident_I_ratio_model_over_obs', 'ident_components', 'ident_method',
            'hires_wl_air_A', 'hires_ident', 'hires_dlam_A', 'dup_group', 'is_primary', 'notes']
    with open(out1, 'w', newline='') as f:
        for h in hdr:
            f.write(h + '\n')
        wr = csv.writer(f, lineterminator='\n')
        wr.writerow(cols)
        for r in rows:
            notes = r['dupnote']
            if not r['cds']:
                notes = ('line only in ESO QC web table (not in A&A/CDS 2003 version); ' + notes).strip('; ')
            if r['F'] == 0.0:
                notes = ('zero flux in source table; ' + notes).strip('; ')
            if r['fw'] < 0.0:
                notes = ('negative FWHM as published; fit artefact; ' + notes).strip('; ')
            nd = 4
            wr.writerow([r['wl'], '', 'air', '0.017', r['flux'], r['fwhm'], r['ident'], r['src'], r['setting'], r['idx'],
                         r['eso_wl'], r['eso_int'], fmt_sig(r['R_per_nm'], 4), f"{float(air2vac(r['lam'])):.{nd}f}",
                         r['species'], r['system'], r['flag'], f"{r['w']:.4f}", r['ncand'], r['mod_la'], r['dlam'],
                         r['Imod'], r['ratio'], r['comps'], r['method'],
                         r['h_wl'], r['h_id'], r['h_d'], r['dup'], r['primary'], notes])

    # ---- HIRES curated table
    pal_oh = {}
    for d in pal:
        if d['chem'] == 'OH' and d['species'] == 'OH':
            m = re.match(r'^OH (\d+)-(\d+) (\S+)$', d['label'])
            pal_oh[(int(m.group(1)), int(m.group(2)), m.group(3))] = d
    out2 = os.path.join(a.out, 'airglow_optical_hires_oh.csv')
    with open(out2, 'w', newline='') as f:
        for h in ['# Keck/HIRES night-sky OH (and O2 b-X 0-1) line list, Osterbrock et al. 1996 PASP 108, 277 (1996PASP..108..277O)',
                  '#   and 1997 PASP 109, 614 (1997PASP..109..614O); CDS/VizieR III/211 (fetched via vizier.cfa.harvard.edu).',
                  '# wl_air_A: verbatim. OH: laboratory/calculated wavelengths (Abrams et al. 1994 ApJS 93, 351, Table 29, or computed by',
                  '#   Osterbrock+ from Abrams term values for higher J), converted to air with Edlen (1953); given to 0.001 A.',
                  '#   Lines are those identified in the HIRES sky spectra (Mauna Kea). O2 rows: computed wavelengths, 0.01 A (Paper 1 Table 3).',
                  '# uncert_A: OH tabulated Abrams lines are accurate to ~1e-5 A (Osterbrock 1996 sect. 3); blends of e/f components',
                  '#   separated by <0.2 A have errors "well below 0.01 A" (Osterbrock 1996). O2: ~0.01 A.',
                  '# palace_*: DERIVED cross-check - the same transition in PALACE v1.0 (HITRAN2020), vacuum->air with Edlen 1966;',
                  '#   check_flag set when |HIRES - PALACE| > 0.02 A (possible typo in III/211, e/f label convention, or high-J model difference).',
                  '#   VizieR column headers: nu\' = upper, nu" = lower vibrational level (VizieR text says the opposite but values are v\'-v\'\').']:
            f.write(h + '\n')
        wr = csv.writer(f, lineterminator='\n')
        wr.writerow(['wl_air_A', 'wl_vac_A', 'medium_of_source', 'uncert_A', 'ident', 'v_upper', 'v_lower', 'transition',
                     'source_key', 'wl_vac_A_derived', 'palace_wl_vac_A', 'palace_wl_air_A_derived', 'hires_minus_palace_A', 'check_flag'])
        dd = []
        for h in hires_oh:
            d = pal_oh.get((h['vu'], h['vl'], h['trans']))
            diff = float(h['wl']) - d['la'] if d else float('nan')
            if d and abs(diff) < 0.005:
                dd.append(diff)
            flag = '' if (d and abs(diff) <= 0.02) else ('no PALACE match' if not d else f"|diff|>0.02 A")
            wr.writerow([h['wl'], '', 'air (Edlen 1953 standard air)', '~1e-5 (lab); see header', f"OH {h['vu']}-{h['vl']} {h['trans']}",
                         h['vu'], h['vl'], h['trans'], 'Osterbrock1996_1997_HIRES_VizieR_III_211_oh',
                         f"{float(air2vac(float(h['wl']))):.3f}", f"{float(d['lam_um'])*1e4:.3f}" if d else '',
                         f"{d['la']:.3f}" if d else '', f"{diff:+.4f}" if d else '', flag])
        for h in hires_o2:
            wr.writerow([h['wl'], '', 'air', '~0.01 (computed)', f"O2 b-X(0-1) {h['trans'].strip()}", 0, 1, h['trans'].strip(),
                         'Osterbrock1996_HIRES_VizieR_III_211_o2 (computed)',
                         f"{float(air2vac(float(h['wl']))):.2f}" if h['wl'] else '', '', '', '', '' if h['wl'] else 'no wavelength in source'])
    qa_hires = dict(n=len(dd), median_mA=float(np.median(dd) * 1e3), std_mA=float(np.std(dd) * 1e3))

    # ---- PALACE model lists: 3000-10500 A (all species) and atomic lines above 10500 A
    PAL_COLS = ['wl_vac_um_verbatim', 'wl_vac_A', 'wl_air_A_derived', 'medium_of_source', 'species', 'system', 'band', 'ident',
                'palace_id', 'I_R', 'trans', 'I_times_trans_R', 'vclass', 'src', 'hlay_km', 'source_key']
    PAL_HDR0 = ['# PALACE v1.0 model airglow line list (Noll, Schmidt, Hannawald, Kausch & Kimeswenger 2025, Geosci. Model Dev. 18, 4353,',
                '#   doi:10.5194/gmd-18-4353-2025; data doi:10.5281/zenodo.14064023, files PMD/pmd_reflineint_mol.dat + pmd_reflineint_atom.dat).']
    PAL_HDR1 = ['# wl_vac_A = wl_vac_um_verbatim x 1e4 (decimal shift). wl_air_A_derived: DERIVED with Edlen (1966) (Hanuschik 2003 eq. 5).',
                '# I_R: PALACE reference intensity in Rayleigh (Paranal climatological mean, zenith, before extinction).',
                '# trans: Paranal atmospheric transmission at zenith, PWV 2.5 mm (Noll+2012 sky model); I_R*trans ~ intensity seen at ground.',
                '# src: 4-bit code (1: X-shooter-based, 2: UVES Noll+2020, 3: other data Cosby+06/Oliva+15, 4: model=1/measurement=0).',
                '# ident: decoded from the 13/12-character PALACE ID (OH e/f = upper-level parity; J\'\' in parentheses; O2 N\'\').',
                '# uncert_A: not stated per line by the source; rounding step 0.001 A; real accuracy follows the parent database.']

    def write_palace(path, hdr_mid, keep):
        n = 0
        with open(path, 'w', newline='') as f:
            for h in PAL_HDR0 + hdr_mid + PAL_HDR1:
                f.write(h + '\n')
            wr = csv.writer(f, lineterminator='\n')
            wr.writerow(PAL_COLS)
            for d in sorted(pal, key=lambda d: d['lv']):
                if not keep(d):
                    continue
                s = d['lam_um']; ip, fp = s.split('.')
                wl_vac_A = f"{int(ip)*10000 + int(fp[:4])}.{fp[4:]}" if len(fp) == 7 else f"{float(s)*1e4:.3f}"
                wr.writerow([s, wl_vac_A, f"{d['la']:.3f}", 'vacuum', d['species'], d['system'], d['band'], d['label'], d['id'],
                             d['I_str'], f"{d['trans']:.4f}", fmt_sig(d['I'] * d['trans'], 4), d['vclass'], d['src'], d['hlay'],
                             'Noll2025_PALACE_v1.0'])
                n += 1
        return n

    nmod = write_palace(os.path.join(a.out, 'airglow_model_palace_3000_10500.csv'), [
        '# Subset 3000 <= lambda_vac <= 10500 A. wl_vac_um_verbatim: source string (vacuum, um, 7 decimals = 0.001 A).',
        '#   Positions: HITRAN2020 (OH; O2 a-X; O2 b-X bands 0-0, 0-1, 1-0, 1-1, 2-0, 2-1 = 469 16O2 lines, and all 16O18O/16O17O lines),',
        '#   Cosby et al. 2006 (the other 465 16O2 b-X lines without a HITRAN2024 counterpart: bands with v\' >= 3 or v\'\' >= 2 plus five',
        '#   N\' = 24-34 lines of 0-1/1-1; O2 A-X, A\'-a, c-X, c-b models), NIST ASD (atoms). Atoms > 10500 A: airglow_model_palace_atomic_ir.csv.'],
        lambda d: 3000.0 <= d['lv'] <= 10500.0)
    atom_ir = [d for d in pal if d['kind'] == 'atom' and d['lv'] > 10500.0]
    sp_ir = sorted(set(d['species'] for d in atom_ir))
    n_ir = write_palace(os.path.join(a.out, 'airglow_model_palace_atomic_ir.csv'), [
        f'# Subset: ALL atomic lines (pmd_reflineint_atom.dat) with lambda_vac > 10500 A: {len(atom_ir)} lines, species {"/".join(sp_ir)}, '
        f'{min(d["lv"] for d in atom_ir):.1f}-{max(d["lv"] for d in atom_ir):.1f} A, summed I_R = {sum(d["I"] for d in atom_ir):.2f} R.',
        '#   PALACE has no atomic lines above 2.14 um and none of other species above 10500 A; companion file airglow_model_palace_3000_10500.csv.',
        '#   These are O I radiative-recombination lines (PALACE class Orc, O+ + e -> O*, F region, hlay 300 km); PALACE scales them from its',
        '#   measured 777/845 nm features with Escalante & Victor (1992) recombination coefficients and NIST A values (no IR line was measured).',
        '# wl_vac_um_verbatim: source string (vacuum, um, 7 decimals = 0.001 A); position = NIST ASD (Kramida et al.) as rounded by PALACE.'],
        lambda d: d['kind'] == 'atom' and d['lv'] > 10500.0)
    print(f'PALACE atomic lines > 10500 A: {n_ir}')

    # ---- QA summary
    def stats(v):
        v = np.array(v)
        if len(v) == 0:
            return None
        med = float(np.median(v)); mad = float(np.median(abs(v - med)) * 1.4826)
        return dict(n=int(len(v)), median_mA=round(med * 1e3, 2), robust_sigma_mA=round(mad * 1e3, 2), mean_mA=round(float(v.mean()) * 1e3, 2))
    flags = defaultdict(int)
    for r in rows:
        flags[r['flag']] += 1
    q = dict(n_rows=len(rows), n_cds=sum(r['cds'] for r in rows), n_eso_only=sum(not r['cds'] for r in rows),
             n_primary=sum(r['primary'] for r in rows), n_dup_groups=len(groups),
             wl_range=[min(r['lam'] for r in rows), max(r['lam'] for r in rows)],
             ident_flags=dict(sorted(flags.items())),
             n_identified_primary=sum(1 for r in rows if r['primary'] and not r['flag'].startswith('U')),
             species_counts_primary={k: v for k, v in sorted(
                 ((sp, sum(1 for r in rows if r['primary'] and r['species'] == sp)) for sp in set(r['species'] for r in rows)),
                 key=lambda kv: (-kv[1], kv[0]))},
             offsets_obs_minus_model={f"{k[0]}|{k[1]}": stats(v) for k, v in sorted(qa.items())},
             hires_minus_palace_OH=qa_hires, n_model_lines_3000_10500=nmod, n_hires_oh=len(hires_oh), n_hires_o2=len(hires_o2))
    with open(os.path.join(a.out, 'airglow_optical_qa.json'), 'w') as f:
        json.dump(q, f, indent=1)
    print(json.dumps(q, indent=1))

if __name__ == '__main__':
    main()
