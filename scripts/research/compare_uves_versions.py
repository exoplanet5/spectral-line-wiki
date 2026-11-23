"""Compare ESO-web UVES gident_*.dat line tables with the CDS/VizieR (A&A 2003) tables.
Usage: python -I compare_uves_versions.py <eso_dat_dir> <vizier_dir>
"""
import sys, os, re
eso_dir, viz_dir = sys.argv[1], sys.argv[2]
pairs = {'346': 'table4', '437': 'table5', '580L': 'table6', '580U': 'table7', '860L': 'table8', '860U': 'table9'}

def read_eso(p):
    rows = []
    for line in open(p):
        m = re.match(r'^\s+(\d+)\s+([\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*$', line)
        if m:
            rows.append(m.groups())
    return rows

def read_viz(p):
    rows = []
    for line in open(p):
        if line.startswith('#') or not line.strip():
            continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) >= 5 and parts[0].strip().isdigit():
            rows.append([x.strip() for x in parts])
    return rows

for s, t in pairs.items():
    e = read_eso(os.path.join(eso_dir, f'gident_{s}.dat'))
    v = read_viz(os.path.join(viz_dir, f'{t}.tsv'))
    ew = [float(r[1]) for r in e]
    vw = [float(r[2]) for r in v]
    # match each viz line to nearest eso line
    unmatched_v = [w for w in vw if min(abs(w - x) for x in ew) > 0.0006]
    unmatched_e = [w for w in ew if min(abs(w - x) for x in vw) > 0.0006]
    # check flux/FWHM agreement for matched
    nd = 0
    for r in v:
        w = float(r[2])
        j = min(range(len(ew)), key=lambda k: abs(ew[k] - w))
        if abs(ew[j] - w) <= 0.0006:
            if abs(float(e[j][3]) - float(r[3])) > 0.00006 or abs(float(e[j][4]) - float(r[4])) > 0.00006:
                nd += 1
    print(f'{s} vs {t}: eso={len(e)} viz={len(v)}  viz-not-in-eso={len(unmatched_v)} eso-not-in-viz={len(unmatched_e)} matched-with-diff-FWHM/flux={nd}')
    if unmatched_v:
        print('   viz only:', unmatched_v[:40])
    if unmatched_e:
        print('   eso only:', unmatched_e[:40])
