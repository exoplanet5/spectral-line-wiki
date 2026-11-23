#!/usr/bin/env python
"""Match every NIST-Handbook strong line to the current NIST ASD line list of the same spectrum.

Inputs: data/handbook.json, data/raw/nist_asd/<El>_<ion>.csv (show_av=4: vacuum < 1850 A, air >= 1850 A)
and, for Handbook vacuum lines in 1850-2000 A, <El>_<ion>_vac1850.csv (show_av=3, vacuum everywhere; its
vacuum strings are used directly).  Only if that file is missing are ASD air values in 1850-2000 A converted
to vacuum with Peck & Reeder (1972) (fallback, flagged 'asd_conv').

Notation: d = decimals of the Handbook string, qh = 10^-d; da = decimals of an ASD string after removing
non-digit characters ('+' and '*' mark ASD's online-computed Ritz values, '?', brackets ...), qa = 10^-da;
u = uncertainty stated by ASD (None if not given).

1. Candidates: ASD lines of the same spectrum with an observed or Ritz value within
   lim = max(5*qh, 0.02 A) of the Handbook value.
2. Best candidate: the nearest (distance d0 = closest of its observed/Ritz values).  Candidates within
   d0 + min(0.5*qh, 0.5 mA) are ties, resolved by: allowed E1 line (empty ASD Type) first, then the largest
   A_ki, then the higher ASD intensity, then distance.
3. Match class from the closest value of the best candidate, with q = max(qh, qa):
       exact  |dl| <= 0.5*q + 0.1 mA      (Handbook value = ASD value rounded)
       close  |dl| <= max(3*q, 3 mA)      (same line, re-measured or re-derived)
       loose  |dl| <= lim                  (tentative cross-reference)
       none   no ASD line within lim
4. Adoption (exact/close).  Values consistent with the Handbook line: those of the best candidate within
   max(3*qh, 3 mA) of it (window on the Handbook scale).  A value with a stated u is eligible only if
   da >= d or u <= 3*qh.  Adopt the eligible value with the smallest u; else an ASD value without stated u
   only if it has more decimals than the Handbook (da > d); else keep the Handbook value.
5. loose-unique: a loose match whose window (lim) holds exactly one ASD line, and that line has a stated
   u <= qh, adopts that ASD value (class 'loose-unique').  Other loose matches keep the Handbook value.
6. Unresolved multiplets (tested when the best candidate is an allowed, designated line, and always for
   H I / He II).  Components S = allowed (E1), designated lines of the best candidate's multiplet
   (same lower and upper configuration and term, J ignored; for H I and He II the same n - n') within lim
   whose position (Ritz, else observed) is nearer to this Handbook line than to any other Handbook line of
   the spectrum (components the Handbook lists separately are not part of the blend).  With >= 2
   components, centroid c = sum(w*p)/sum(w), w = g_k*A_ki (g_k*intensity if an A_ki is missing),
   the Handbook line is a blend ('blend' = number of components) if
     a. 'shared-obs': ASD gives one observed value for >= 2 components, the Handbook value lies within the
        close window of it and is at least as close to it as to any component position; or
     b. 'centroid': c differs from the strongest component by more than that component's uncertainty
        (stated, else 2.5 units of its last digit), so the blend matters, and either
          - both c and the nearest component lie within one unit of the Handbook's last digit (<= qh) and
            ASD has not observed every component separately: the Handbook value cannot single out the
            component, so it is taken as the blend; or
          - the Handbook value is closer to c than to any component position and is the rounding (within
            0.5*qh + 0.1 mA) of c, or of an ASD observed value that is itself closer to c than to its own
            component's Ritz value (an observation of the blend); or
     c. "n-n'": H I / He II values with <= 2 decimals (the whole n - n' group).
   A blend never takes a single component's Ritz value.  Its value pool: ASD observed values shared by
   >= 2 components, observed values closer to c than to their own component's Ritz value (measured on
   the blend), and ASD composite entries (rows listing several J; for <= 2-decimal H I / He II values the
   bare "n - n'" rows); step 4 is applied
   to that pool, else the Handbook value is kept.  The transition shown is the composite entry (H I/He II)
   or the strongest component (largest g_k*A_ki); dl/class then refer to the adopted pool value or, when
   the Handbook value is kept, to that shown row.
7. H I / He II: when the adopted ASD row is undesignated (ASD's unclassified blend rows), the transition
   fields are taken from the nearest composite "n - n'" row.

Also collects ASD Ritz/wavenumber pairs for the dispersion-formula identification (notes).

Output: data/handbook_asd.json (one record per Handbook line: cls, adopt 'asd'|'handbook', the chosen ASD
value asd_wl_str / asd_unc / asd_kind, dl_mA = closest ASD value - Handbook, transition fields of the shown
row, asd_file 'main'|'vac1850'; blend, blend_kind, blend_spread_mA, blend_centroid for multiplets; mult_n =
number of multiplet components near a line that was judged to single out one component).
"""
import bisect, csv, io, json, math, os, re, sys
from collections import Counter, defaultdict

ROOT = os.path.join(os.path.dirname(__file__), '..')
H = json.load(open(os.path.join(ROOT, 'data', 'handbook.json')))
HYDROGENIC = {('H', 'I'), ('He', 'II')}


def clean(v):
    v = (v or '').strip()
    m = re.match(r'^="(.*)"$', v)
    if m:
        v = m.group(1)
    return v.strip()


def num(v):
    v = clean(v)
    m = re.match(r'^[\[(]?\s*([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)', v)
    return float(m.group(1)) if m else None


def dec_of(s):
    """decimals of a wavelength string, ignoring non-digit suffixes such as ASD's '+' or '?'"""
    s = re.sub(r'[^\d.]', '', s or '')
    return len(s.split('.')[1]) if '.' in s else 0


def load_asd_file(fn):
    if not os.path.exists(fn):
        return None
    rows = []
    txt = open(fn, encoding='utf-8', errors='replace').read()
    for r in csv.DictReader(io.StringIO(txt)):
        r = {k: clean(v) for k, v in r.items() if k}
        obs_key = next((k for k in r if k.startswith('obs_wl')), None)
        ritz_key = next((k for k in r if k.startswith('ritz_wl')), None)
        obs, ritz = num(r.get(obs_key)), num(r.get(ritz_key))
        if obs is None and ritz is None:
            continue
        rows.append({
            'obs': obs, 'obs_str': r.get(obs_key, '') if obs is not None else '', 'unc_obs': num(r.get('unc_obs_wl')),
            'ritz': ritz, 'ritz_str': r.get(ritz_key, '') if ritz is not None else '', 'unc_ritz': num(r.get('unc_ritz_wl')),
            'wn': num(r.get('wn(cm-1)')), 'wn_str': r.get('wn(cm-1)', ''),
            'intens': r.get('intens', ''), 'inum': num(r.get('intens')),
            'Aki': num(r.get('Aki(s^-1)')), 'acc': r.get('Acc', ''),
            'Ei': r.get('Ei(cm-1)', ''), 'Ek': r.get('Ek(cm-1)', ''),
            'conf_i': r.get('conf_i', ''), 'term_i': r.get('term_i', ''), 'J_i': r.get('J_i', ''),
            'conf_k': r.get('conf_k', ''), 'term_k': r.get('term_k', ''), 'J_k': r.get('J_k', ''),
            'type': r.get('Type', ''), 'line_ref': r.get('line_ref', ''), 'tp_ref': r.get('tp_ref', ''),
        })
    return rows


def load_asd(el, ion):
    return load_asd_file(os.path.join(ROOT, 'data', 'raw', 'nist_asd', f'{el}_{ion}.csv'))


def load_asd_vac1850(el, ion):
    return load_asd_file(os.path.join(ROOT, 'data', 'raw', 'nist_asd', f'{el}_{ion}_vac1850.csv'))


def _n_pr72(wv):
    s2 = (1e4 / wv) ** 2
    return 1 + 1e-8 * (8060.51 + 2480990 / (132.274 - s2) + 17455.7 / (39.32957 - s2))


def air_to_vac_pr72(wa):
    v = wa
    for _ in range(30):
        nv = wa * _n_pr72(v)
        if abs(nv - v) < 1e-12:
            break
        v = nv
    return v


def to_vac_if_needed(w, ws):
    """Fallback only (no <El>_<ion>_vac1850.csv): ASD show_av=4 gives AIR from 1850 A upward;
    return the vacuum value (same number of decimals) for 1850 <= w < 2000, else unchanged."""
    if w is None or not (1850 <= w < 2000):
        return w, ws, False
    d = dec_of(ws)
    v = air_to_vac_pr72(w)
    return v, f'{v:.{d}f}', True


def is_e1(a):
    return a['type'] in ('', None)


def designated(a):
    return bool(a['term_i']) and bool(a['term_k'])


def n_of(conf):
    """principal quantum number of a hydrogenic configuration ('3d' -> 3, '3' -> 3)"""
    m = re.match(r'^(\d+)', conf or '')
    return int(m.group(1)) if m else None


def is_composite(a, hyd):
    """ASD composite entry: hydrogenic 'n - n'' rows (bare n, no term) or rows listing several J values"""
    if hyd:
        return bool(re.fullmatch(r'\d+', a['conf_i'] or '')) and bool(re.fullmatch(r'\d+', a['conf_k'] or ''))
    return designated(a) and (',' in a['J_i'] or ',' in a['J_k'])


def gk(a):
    j = a['J_k']
    m = re.fullmatch(r'(\d+)(?:/(\d+))?', j or '')
    if not m:
        return 1.0
    J = int(m.group(1)) / (int(m.group(2)) if m.group(2) else 1)
    return 2 * J + 1


def pos(a, conv):
    """component position: Ritz value if given, else observed"""
    w = a['ritz'] if a['ritz'] is not None else a['obs']
    if conv:
        w = to_vac_if_needed(w, '')[0]
    return w


def implied_unc(v):
    """stated ASD uncertainty, else 2.5 units of the last digit"""
    return v[3] if v[3] is not None else 2.5 * 10.0 ** (-v[4])


def main():
    lines = H['lines']
    by_spec = defaultdict(list)
    for i, l in enumerate(lines):
        by_spec[(l['el'], l['ion'])].append(i)
    cls, kinds = Counter(), Counter()
    disp = []          # (air, vac_from_wn) pairs for formula identification
    out = []
    missing_spec = []
    for (el, ion), idx in sorted(by_spec.items()):
        asd_main = load_asd(el, ion)
        asd_vac = load_asd_vac1850(el, ion)
        hyd = (el, ion) in HYDROGENIC
        hb_med = defaultdict(list)            # Handbook lines of this spectrum per medium, for "nearest line"
        for i in idx:
            hb_med[lines[i]['medium']].append((lines[i]['wl'], i))
        for v_ in hb_med.values():
            v_.sort()
        if asd_main is None:
            missing_spec.append(f'{el} {ion}')
        for i in idx:
            l = lines[i]
            rec = {'i': i, 'cls': 'none'}
            vac_window = l['medium'] == 'vacuum' and 1849 <= l['wl'] < 2001
            if vac_window and asd_vac is not None:
                asd, conv, rec_src = asd_vac, False, 'vac1850'
            else:
                asd, conv, rec_src = asd_main, (l['medium'] == 'vacuum'), 'main'
            if not asd:
                cls['none'] += 1; out.append(rec); continue
            d = l['dec']
            qh = 10.0 ** (-d)
            lim = max(5 * qh, 0.02)
            win_ok = max(3 * qh, 0.003)
            hbw = [w_ for w_, _ in hb_med[l['medium']]]
            hbi = [i_ for _, i_ in hb_med[l['medium']]]

            def nearest_hb(p):
                k = bisect.bisect_left(hbw, p)
                best = min((j for j in (k - 1, k) if 0 <= j < len(hbw)), key=lambda j: abs(hbw[j] - p))
                return hbi[best]

            cands = []
            used_conv = False
            for a in asd:
                vals = []
                for kind in ('obs', 'ritz'):
                    w = a[kind]
                    if w is None:
                        continue
                    ws = a[kind + '_str']
                    if conv:
                        w, ws, c_ = to_vac_if_needed(w, ws)
                        used_conv |= c_
                    vals.append((kind, w, ws, a['unc_' + kind], dec_of(ws), w - l['wl']))
                if not vals:
                    continue
                dmin = min(abs(v[5]) for v in vals)
                if dmin <= lim:
                    cands.append((dmin, a, vals))
            if not cands:
                cls['none'] += 1; out.append(rec); continue
            # ---- best candidate: nearest, near-ties resolved by E1, A_ki, intensity
            d0 = min(c_[0] for c_ in cands)
            tie = [c_ for c_ in cands if c_[0] <= d0 + min(0.5 * qh, 5e-4)]
            tie.sort(key=lambda c_: (not is_e1(c_[1]), -(c_[1]['Aki'] or 0), -(c_[1]['inum'] or 0), c_[0]))
            _, a, vals = tie[0]

            def classify(vals_):
                v0 = min(vals_, key=lambda v: abs(v[5]))
                q = max(qh, 10.0 ** (-v0[4]))
                if abs(v0[5]) <= 0.5 * q + 1e-4:
                    return 'exact', v0
                if abs(v0[5]) <= max(3 * q, 0.003):
                    return 'close', v0
                return 'loose', v0

            def eligible(v):
                return v[3] is not None and (v[4] >= d or v[3] <= 3 * qh)

            def adopt_from(pool, c_):
                """pool: candidate ASD values consistent with the Handbook line -> (value, 'asd'|'handbook')"""
                ok = [v for v in pool if abs(v[5]) <= win_ok]
                stated = [v for v in ok if eligible(v)]
                if c_ != 'loose' and stated:
                    return min(stated, key=lambda v: (v[3], abs(v[5]))), 'asd'
                finer = [v for v in ok if v[3] is None and v[4] > d]
                if c_ != 'loose' and finer:
                    return max(finer, key=lambda v: (v[4], -abs(v[5]))), 'asd'
                ref = min(ok or pool, key=lambda v: abs(v[5]))
                return ref, 'handbook'

            c, v0 = classify(vals)
            blend = None          # (n_components, kind, centroid, spread)
            n_mult = 0            # multiplet components nearest to this line (diagnostic)
            row = a               # ASD row whose transition is reported
            # ---- unresolved multiplets
            if hyd or (is_e1(a) and designated(a)):
                if hyd:
                    ref = a
                    if not designated(a):          # composite or undesignated best row: use the composite's n - n'
                        comps_ = [c_ for c_ in cands if is_composite(c_[1], True)]
                        ref = min(comps_, key=lambda c_: c_[0])[1] if comps_ else a
                    key = lambda b: (n_of(b['conf_i']), n_of(b['conf_k']))
                    S = [c_ for c_ in cands if is_e1(c_[1]) and designated(c_[1]) and key(c_[1]) == key(ref)
                         and key(ref)[0] is not None]
                    comp = [c_ for c_ in cands if is_composite(c_[1], True) and key(c_[1]) == key(ref)] if d <= 2 else []
                else:
                    key = lambda b: (b['conf_i'], b['term_i'], b['conf_k'], b['term_k'])
                    S = [c_ for c_ in cands if is_e1(c_[1]) and designated(c_[1]) and key(c_[1]) == key(a)
                         and not is_composite(c_[1], False)]
                    comp = [c_ for c_ in cands if is_composite(c_[1], False) and key(c_[1]) == key(a)]
                # components the Handbook lists as separate lines belong to those lines
                S = [c_ for c_ in S if nearest_hb(pos(c_[1], conv)) == i]
                n_mult = len(S)
                if len(S) >= 2:
                    P = [pos(b, conv) for _, b, _ in S]
                    spread = max(P) - min(P)
                    useA = all(b['Aki'] for _, b, _ in S)
                    wts = [gk(b) * (b['Aki'] if useA else (b['inum'] or 1)) for _, b, _ in S]
                    cen = sum(w_ * p for w_, p in zip(wts, P)) / sum(wts)
                    p_near = min(abs(p - l['wl']) for p in P)
                    k_s = max(range(len(S)), key=lambda k: wts[k])
                    b_s = S[k_s][1]
                    v_s = next(v for v in S[k_s][2] if v[0] == ('ritz' if b_s['ritz'] is not None else 'obs'))
                    cnt = Counter(c_[1]['obs_str'] for c_ in S if c_[1]['obs_str'])
                    kind_bl = None
                    if hyd and d <= 2:
                        if spread > 0:
                            kind_bl = 'n-n\''             # coarse H I / He II value = whole n - n' group
                    elif spread > 0:
                        # (a) one ASD observed value shared by >= 2 components (ASD's own blend measurement)
                        for ostr, n_ in cnt.items():
                            if n_ < 2:
                                continue
                            ov = next(v for c_ in S if c_[1]['obs_str'] == ostr for v in c_[2] if v[0] == 'obs')
                            if abs(ov[5]) <= max(3 * max(qh, 10.0 ** (-ov[4])), 0.003) and abs(ov[5]) <= p_near + 1e-9:
                                kind_bl = 'shared-obs'
                        # (b) the blend matters (centroid off the strongest component by more than that component's
                        #     uncertainty) and the Handbook value is the blend: either it is the rounding of both a
                        #     component and the centroid (ambiguous -> blend), or it is closer to the centroid than
                        #     to any component and is the rounding of the centroid or of an ASD observed value that
                        #     was measured on the blend (closer to the centroid than to its own Ritz value)
                        r_c = abs(cen - l['wl'])
                        if not kind_bl and abs(cen - P[k_s]) > implied_unc(v_s):
                            hb_is_blend_obs = any(
                                v[0] == 'obs' and b['ritz'] is not None
                                and abs(v[5]) <= 0.5 * max(qh, 10.0 ** (-v[4])) + 1e-4
                                and abs(v[1] - cen) < abs(v[1] - pos(b, conv))
                                for _, b, vv in S for v in vv)
                            obs_strs = [b['obs_str'] for _, b, _ in S]
                            resolved_in_asd = all(obs_strs) and len(set(obs_strs)) == len(S)
                            if p_near <= qh + 1e-9 and r_c <= qh + 1e-9 and not resolved_in_asd:
                                kind_bl = 'centroid'
                            elif r_c < p_near and (r_c <= 0.5 * qh + 1e-4 or hb_is_blend_obs):
                                kind_bl = 'centroid'
                    if kind_bl:
                        # blend values: shared observed values, observed values measured on the blend, composites
                        pool, seen = [], set()
                        for _, b, vv in S:
                            for v in vv:
                                if v[0] != 'obs':
                                    continue
                                on_blend = b['ritz'] is not None and abs(v[1] - cen) < abs(v[1] - pos(b, conv))
                                if (cnt[b['obs_str']] >= 2 or on_blend) and v[2] not in seen:
                                    pool.append(v); seen.add(v[2])
                        pool += [v for _, b, vv in comp for v in vv]
                        # reported transition: the composite entry (H I/He II <= 2 decimals), else the strongest component
                        row = min(comp, key=lambda c_: c_[0])[1] if (hyd and comp) else b_s
                        chosen, adopt = None, 'handbook'
                        if pool:
                            c, v0 = classify(pool)
                            chosen, adopt = adopt_from(pool, c)
                        if adopt == 'handbook':
                            # keep the Handbook value; class and offset refer to the reported row
                            c, v0 = classify([v for _, b, vv in cands if b is row for v in vv])
                            chosen = v0
                        blend = (len(S), kind_bl, cen, spread)
            if blend is None:
                if c == 'loose' and len(cands) == 1:
                    prec = [v for v in vals if v[3] is not None and v[3] <= qh]
                    if prec:
                        c = 'loose-unique'
                        chosen, adopt = min(prec, key=lambda v: v[3]), 'asd'
                if c != 'loose-unique':
                    chosen, adopt = adopt_from(vals, c)
                if hyd and not designated(row):
                    # ASD's unclassified H I blend rows: show the transition of the nearest composite n - n' row
                    comps_ = [c_ for c_ in cands if is_composite(c_[1], True)]
                    if comps_:
                        row = min(comps_, key=lambda c_: c_[0])[1]
            kind, w, ws, u, da, dl = chosen
            rec.update(cls=c, dl_mA=round(v0[5] * 1000, 4), adopt=adopt, asd_kind=kind, asd_wl=w,
                       asd_wl_str=ws, asd_unc=u, asd_dec=da,
                       asd_obs=row['obs_str'], asd_unc_obs=row['unc_obs'],
                       asd_ritz=row['ritz_str'], asd_unc_ritz=row['unc_ritz'], asd_wn=row['wn_str'],
                       asd_int=row['intens'], Aki=row['Aki'], acc=row['acc'], Ei=row['Ei'], Ek=row['Ek'],
                       lower=' '.join(x for x in (row['conf_i'], row['term_i'], row['J_i']) if x),
                       upper=' '.join(x for x in (row['conf_k'], row['term_k'], row['J_k']) if x),
                       type=row['type'], asd_line_ref=row['line_ref'], asd_tp_ref=row['tp_ref'],
                       n_cand=len(cands), asd_file=rec_src)
            if n_mult >= 2 and not blend:
                rec['mult_n'] = n_mult      # resolved: Handbook value singles out one component
            if blend:
                rec.update(blend=blend[0], blend_kind=blend[1], blend_spread_mA=round(blend[3] * 1000, 3),
                           blend_centroid=round(blend[2], 7))
                kinds[blend[1]] += 1
            if used_conv:
                rec['asd_conv'] = 1      # ASD air value converted to vacuum (P&R 1972; fallback)
            if a['ritz'] and a['wn'] and a['ritz'] > 2000 and a['unc_ritz'] and a['unc_ritz'] < 2e-4:
                disp.append((a['ritz'], 1e8 / a['wn'], a['ritz_str'], a['wn_str']))
            cls[rec['cls']] += 1
            out.append(rec)
    print('match classes', dict(cls))
    print('blends', sum(kinds.values()), dict(kinds))
    print('adoption', dict(Counter((m.get('adopt'), m.get('asd_kind')) for m in out if m['cls'] != 'none')))
    print('spectra missing in ASD:', missing_spec)
    json.dump({'matches': out, 'disp_sample': disp[:20000], 'missing_spectra': missing_spec},
              open(os.path.join(ROOT, 'data', 'handbook_asd.json'), 'w'))


if __name__ == '__main__':
    main()
