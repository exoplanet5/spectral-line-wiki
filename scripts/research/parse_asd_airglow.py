#!/usr/bin/env python
"""Parse the NIST ASD CSV downloads (av2 = default air/vac display + wavenumber, av3 = vacuum)
and select the transitions listed in airglow_targets.py.

Usage: python -I parse_asd_airglow.py <asd_dir> <scripts_dir> <out_json>
(run with -I because it reads downloaded files; scripts_dir is added to sys.path explicitly)

For every selected transition, writes ASD strings verbatim (observed / Ritz, air / vacuum,
uncertainties, wavenumber), A_ki, configurations, terms, J, transition type, ASD bibliographic
codes.  M1 and E2 rows of the same forbidden transition are merged (A values summed, parts kept).
"""
import csv, glob, json, os, re, sys

asd_dir, scripts_dir, out_json = sys.argv[1:4]
sys.path.insert(0, scripts_dir)
from airglow_targets import T  # noqa: E402

def clean(v):
    v = (v or '').strip()
    m = re.fullmatch(r'="(.*)"', v)
    return m.group(1) if m else v

def load(fn):
    with open(fn, newline='') as f:
        rows = list(csv.reader(f))
    if not rows:
        return []
    hdr = [h.strip() for h in rows[0]]
    out = []
    for r in rows[1:]:
        if not any(x.strip() for x in r):
            continue
        d = {h: clean(r[i]) if i < len(r) else '' for i, h in enumerate(hdr) if h}
        out.append(d)
    return out

def num(s):
    s = re.sub(r'[^0-9.eE+-]', '', s or '')
    try:
        return float(s)
    except ValueError:
        return None

def wlkeys(d):
    """return (obs, unc_obs, ritz, unc_ritz, medium) strings."""
    med = None
    obs = ritz = ''
    for k in d:
        if k.startswith('obs_wl_'):
            obs = d[k]; med = 'air' if '_air' in k else 'vac'
        if k.startswith('ritz_wl_'):
            ritz = d[k]; med = 'air' if '_air' in k else 'vac'
    return obs, d.get('unc_obs_wl', ''), ritz, d.get('unc_ritz_wl', ''), med

def key(d):
    return (d.get('Ei(cm-1)', ''), d.get('Ek(cm-1)', ''), d.get('J_i', ''), d.get('J_k', ''),
            d.get('term_i', ''), d.get('term_k', ''), d.get('conf_i', ''), d.get('conf_k', ''))

# load all files
tab = {}
for fn in sorted(glob.glob(os.path.join(asd_dir, '*_av2.csv'))):
    base = os.path.basename(fn)[:-8]
    spec = ' '.join(base.split('_')[:2])
    rows2 = load(fn)
    rows3 = load(fn.replace('_av2.csv', '_av3.csv')) if os.path.exists(fn.replace('_av2.csv', '_av3.csv')) else []
    tab.setdefault(spec, []).append((base, rows2, rows3))

def term_ok(want, have):
    if want is None:
        return True
    return have.replace(' ', '') == want or have.split()[-1:] == [want]

results = []
for t in T:
    spec = t['spec']
    cands = []; seen = set()
    for base, rows2, rows3 in tab.get(spec, []):
        for d in rows2:
            obs, uo, ritz, ur, med = wlkeys(d)
            w = num(ritz) or num(obs)
            if w is None or abs(w - t['nominal']) > t['tol']:
                continue
            if spec != 'H I' and not d.get('term_i') and not d.get('conf_i'):
                continue
            if spec == 'H I':
                # keep the n-n' multiplet row (no term) and E1 fine-structure components only
                if not d.get('conf_i'):
                    continue
                if d.get('Type') not in ('', None) or not d.get('Aki(s^-1)'):
                    continue
            if not (term_ok(t['ti'], d.get('term_i', '')) and term_ok(t['tk'], d.get('term_k', ''))):
                continue
            sig = (key(d), d.get('Type', ''), obs, ritz)
            if sig in seen:          # same ASD row present in two overlapping download windows
                continue
            seen.add(sig)
            cands.append((base, d, rows3))
    # merge rows sharing the same transition (M1+E2 parts)
    groups = {}
    for base, d, rows3 in cands:
        groups.setdefault(key(d), []).append((base, d, rows3))
    if not groups:
        results.append(dict(target=t, status='NOT FOUND'))
        continue
    # for non-H targets choose the transition closest to nominal; H I keeps all components
    def dist(g):
        d = g[0][1]; obs, uo, ritz, ur, med = wlkeys(d)
        return abs((num(ritz) or num(obs)) - t['nominal'])
    glist = sorted(groups.values(), key=dist)
    if not t.get('all') and (spec not in ('H I', 'He I') or t['group'] in ('HeI_10830', 'HeI_584', 'HeI_537')):
        glist = glist[:1]
    for g in glist:
        base, d, rows3 = g[0]
        obs, uo, ritz, ur, med = wlkeys(d)
        # vacuum partner(s)
        vac = [r for r in rows3 if key(r) == key(d)]
        vobs = vritz = vuo = vur = ''
        if vac:
            vobs, vuo, vritz, vur, _ = wlkeys(vac[0])
        parts = []
        A_tot = 0.0; A_any = False
        for _, dd, _ in g:
            a = num(dd.get('Aki(s^-1)', ''))
            parts.append(dict(type=dd.get('Type', '') or 'E1', Aki=dd.get('Aki(s^-1)', ''), Acc=dd.get('Acc', ''),
                              fik=dd.get('fik', ''), tp_ref=dd.get('tp_ref', ''), line_ref=dd.get('line_ref', ''),
                              intens=dd.get('intens', '')))
            if a is not None:
                A_tot += a; A_any = True
        results.append(dict(
            target=t, status='ok', asd_file=base,
            medium_default=med,
            obs_default=obs, unc_obs_default=uo, ritz_default=ritz, unc_ritz_default=ur,
            obs_vac=vobs, unc_obs_vac=vuo, ritz_vac=vritz, unc_ritz_vac=vur,
            wn_cm=d.get('wn(cm-1)', ''),
            Ei=d.get('Ei(cm-1)', ''), Ek=d.get('Ek(cm-1)', ''),
            conf_i=d.get('conf_i', ''), term_i=d.get('term_i', ''), J_i=d.get('J_i', ''),
            conf_k=d.get('conf_k', ''), term_k=d.get('term_k', ''), J_k=d.get('J_k', ''),
            parts=parts, Aki_total=(f'{A_tot:.3e}' if A_any else ''), n_parts=len(g)))

json.dump(results, open(out_json, 'w'), indent=1)
bad = [r for r in results if r['status'] != 'ok']
print(len(results), 'records;', len(bad), 'not found')
for r in bad:
    print('  NOT FOUND', r['target'])
for r in results:
    if r['status'] == 'ok':
        print(f"{r['target']['spec']:6s} {r['target']['group']:18s} air/def={r['ritz_default'] or r['obs_default']:>14s} "
              f"({r['medium_default']}) vac={r['ritz_vac'] or r['obs_vac']:>14s} unc={r['unc_ritz_default'] or r['unc_obs_default'] or '-':>8s} "
              f"{r['term_i']}{r['J_i']}-{r['term_k']}{r['J_k']} A={r['Aki_total']} n={r['n_parts']}")
