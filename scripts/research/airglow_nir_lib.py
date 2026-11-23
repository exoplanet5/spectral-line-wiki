"""Parsers for NIR airglow line sources (HITRAN .par, PALACE FITS, Oliva+2015 VizieR TSV,
Rousselot+2000 ESO ISAAC list). Pure functions; paths are passed in by the caller."""
import re, glob, os
import numpy as np

# ---------------- HITRAN 160-character .par records ----------------
def parse_par_line(s):
    return dict(
        molec=int(s[0:2]), iso=s[2:3], nu_str=s[3:15].strip(), nu=float(s[3:15]),
        S=float(s[15:25]), A=float(s[25:35]), gamma_air=s[35:40], gamma_self=s[40:45],
        Elow=float(s[45:55]), n_air=s[55:59], delta_air=s[59:67],
        Vu=s[67:82], Vl=s[82:97], Qu=s[97:112], Ql=s[112:127],
        ierr=s[127:133], iref=s[133:145], flag=s[145:146],
        gu=float(s[146:153]), gl=float(s[153:160]))

def read_par_files(pattern):
    recs = []
    for f in sorted(glob.glob(pattern)):
        with open(f) as fh:
            for line in fh:
                line = line.rstrip('\n').rstrip('\r')
                if len(line) != 160:
                    continue
                r = parse_par_line(line); r['file'] = os.path.basename(f); recs.append(r)
    return recs

# HITRAN position uncertainty code -> [lo, hi) in cm-1 (HITRAN docs / Rothman+2005 Table 5)
IERR_NU = {0: (1.0, np.inf), 1: (0.1, 1.0), 2: (0.01, 0.1), 3: (1e-3, 1e-2), 4: (1e-4, 1e-3),
           5: (1e-5, 1e-4), 6: (1e-6, 1e-5), 7: (1e-7, 1e-6), 8: (1e-8, 1e-7), 9: (1e-9, 1e-8)}

DJ = {'M': -4, 'N': -3, 'O': -2, 'P': -1, 'Q': 0, 'R': 1, 'S': 2, 'T': 3, 'U': 4}

def oh_quanta(r):
    """OH X2Pi - X2Pi (Meinel). V field: 'X3/2   v' (Omega, v). Ql: ' Br J'' sym sym'."""
    om_u, v_u = r['Vu'].split(); om_l, v_l = r['Vl'].split()
    Fu = 1 if om_u == 'X3/2' else 2   # OH X state inverted: F1 = 2Pi3/2, F2 = 2Pi1/2
    Fl = 1 if om_l == 'X3/2' else 2
    q = r['Ql']
    br = q[1:3]; Jl = float(q[3:8]); sym = q[8:10]
    dN, dJ = br[0], br[1]
    Ju = Jl + DJ[dJ]
    return dict(v_u=int(v_u), v_l=int(v_l), F_u=Fu, F_l=Fl, br=br, dN=dN, dJ=dJ,
                J_l=Jl, J_u=Ju, sym=sym)

def o2_quanta(r):
    """O2: V field 'b      0' (state, v). Ql: ' dN N'' dJ J''   d|q'."""
    st_u, v_u = r['Vu'].split(); st_l, v_l = r['Vl'].split()
    q = r['Ql']
    m = re.match(r'^\s*([MNOPQRSTU])\s*(\d+)\s*([MNOPQRSTU])\s*(\d+)\s*([dq]?)\s*$', q)
    if not m:
        return None
    dN, Nl, dJ, Jl, typ = m.groups()
    Nl = int(Nl); Jl = int(Jl)
    return dict(st_u=st_u, st_l=st_l, v_u=int(v_u), v_l=int(v_l), dN=dN, N_l=Nl,
                dJ=dJ, J_l=Jl, N_u=Nl + DJ[dN], J_u=Jl + DJ[dJ], mtype=typ)

# ---------------- VizieR ASU-TSV ----------------
def read_vizier_tsv(path):
    """Return list of dict rows (strings verbatim) from a VizieR asu-tsv file."""
    with open(path) as fh:
        lines = [l.rstrip('\n') for l in fh]
    i = 0
    while i < len(lines) and (lines[i].startswith('#') or lines[i].strip() == ''):
        i += 1
    header = lines[i].split('\t'); units = lines[i+1]; dashes = lines[i+2]
    rows = []
    for l in lines[i+3:]:
        if not l.strip() or l.startswith('#'):
            continue
        parts = l.split('\t')
        rows.append({h: p for h, p in zip(header, parts)})
    return rows

OLIVA_ID = re.compile(r'^\[(\d+)-(\d+)\]([PQR])([12])([ef])\((\d+\.5)\)\s*$')

def parse_oliva_id(s):
    m = OLIVA_ID.match(s.strip())
    if not m:
        return None
    vu, vl, br, F, par, Jl = m.groups()
    return dict(v_u=int(vu), v_l=int(vl), dJ=br, F=int(F), parity=par, J_l=float(Jl))

def read_rousselot(path):
    out = []
    with open(path) as fh:
        for l in fh:
            if l.startswith('#') or not l.strip():
                continue
            a, b = l.split()[:2]
            out.append((a, float(a), b, float(b)))
    return out

# ---------------- air wavelength (derived) ----------------
def vac_to_air_morton2000(lam_vac_A):
    """IAU standard (Morton 2000, ApJS 130, 403, eq. 8; Edlen 1966 / Birch & Downs 1994 based):
    n - 1 = 8.34254e-5 + 2.406147e-2/(130 - s^2) + 1.5998e-4/(38.9 - s^2), s = 1e4/lam_vac[A] (um^-1).
    Standard dry air at 15 C and 101325 Pa (Edlen 1966 standard air, revised by Birch & Downs 1994)."""
    s2 = (1e4 / np.asarray(lam_vac_A, float))**2
    n = 1 + 8.34254e-5 + 2.406147e-2 / (130 - s2) + 1.5998e-4 / (38.9 - s2)
    return np.asarray(lam_vac_A, float) / n


def palace_id_from_hitran(r):
    """Build the 13-character PALACE v1.0 line ID (pmd_reflineint_mol.dat header) from HITRAN quanta."""
    vc = lambda v: ('X' if v == 10 else str(v)) if r['mol'] == 'OH' else ('ABCDEF'[v - 10] if v >= 10 else str(v))
    if r['mol'] == 'OH':
        Nu = r['J_u'] - 0.5 if r['F_u'] == 1 else r['J_u'] + 0.5
        return 'OHXX%s%s%s%d%d%02d%s' % (vc(r['v_u']), vc(r['v_l']), r['br'], r['F_u'], r['F_l'], int(round(Nu)), r['sym'][0])
    pre = {'16O2': 'O2', '16O18O': 'OQ', '16O17O': 'OG'}[r['isotop']]
    Fl = {1: 1, 0: 2, -1: 3}.get(r['J_l'] - r['N_l'], 9)
    return '%s%s%s%s%s%s%s1%d%02d%s' % (pre, r['st_u'], r['st_l'], vc(r['v_u']), vc(r['v_l']), r['dN'], r['dJ'], Fl, r['N_u'], r['mtype'])
