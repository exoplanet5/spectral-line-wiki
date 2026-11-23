/* app_core.js -- page bootstrap, global state, datasets, element-line section, explorer, finding list */
(function () {
  'use strict';
  const { $, $$, esc, h, store, fmtWl, fmtSrc, fmtSig, grade, gradeBadge, DataTable, SpectrumView } = UI;
  const S = SpecCore;
  const APP = window.APP = { datasets: {}, state: { med: store.get('med', 'air') } };
  const J = id => { const e = document.getElementById(id); return e ? JSON.parse(e.textContent) : null; };
  APP.J = J;

  // ---------------------------------------------------------------- theme / nav
  function initTheme() {
    const qt = location.search.match(/theme=(light|dark)/);          // test hook, e.g. ?selftest&theme=light
    const t = qt ? qt[1] : store.get('theme', null);
    if (t) document.documentElement.dataset.theme = t;
    $('#themeBtn').addEventListener('click', () => {
      const cur = document.documentElement.dataset.theme ||
        (matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light');
      const nx = cur === 'dark' ? 'light' : 'dark';
      document.documentElement.dataset.theme = nx; store.set('theme', nx);
      window.dispatchEvent(new Event('themechange'));
    });
    $('#navBtn').addEventListener('click', () => document.body.classList.toggle('navopen'));
    $$('nav.toc a').forEach(a => a.addEventListener('click', () => document.body.classList.remove('navopen')));
    // scrollspy
    const links = $$('nav.toc a[href^="#"]');
    const map = new Map(links.map(a => [a.getAttribute('href').slice(1), a]));
    const io = new IntersectionObserver(es => {
      for (const e of es) if (e.isIntersecting) {
        links.forEach(a => a.classList.remove('active'));
        const a = map.get(e.target.id); if (a) a.classList.add('active');
      }
    }, { rootMargin: '-60px 0px -75% 0px' });
    map.forEach((a, id) => { const s = document.getElementById(id); if (s) io.observe(s); });
  }

  // ---------------------------------------------------------------- medium handling
  /** convert a wavelength between display media (only >= 2000 A differs) */
  APP.convMed = function (w, from, to) {
    if (from === to || !isFinite(w)) return w;
    if (to === 'vac') return w >= 1999.35 ? S.airToVac(w, 'peckreeder1972') : w;
    return w >= 2000 ? S.vacToAir(w, 'peckreeder1972') : w;
  };
  APP.setMed = function (m) {
    const old = APP.state.med;
    if (old && old !== m && APP.booted) {
      for (const v of (window.UI_VIEWS || [])) {
        v.range = v.range.map(w => APP.convMed(w, old, m));
        v.markers = (v.markers || []).map(k => Object.assign({}, k, { w: APP.convMed(k.w, old, m) }));
      }
      const q = document.getElementById('qWl');
      if (q && isFinite(parseFloat(q.value))) q.value = APP.convMed(parseFloat(q.value), old, m).toFixed(3);
    }
    APP.state.med = m; store.set('med', m);
    $$('[data-med]').forEach(b => b.classList.toggle('on', b.dataset.med === m));
    window.dispatchEvent(new Event('medchange'));
  };
  /** wavelength of a normalised row in the current display medium (number) */
  APP.wl = r => (APP.state.med === 'air' && r.wAir !== null && r.wAir !== undefined) ? r.wAir : r.wVac;
  APP.medLabel = () => APP.state.med === 'air' ? 'air (vacuum < 2000 Å)' : 'vacuum';

  /** format the wavelength of row r in medium m ('a' or 'v'), fading digits beyond the source precision */
  APP.fmtRowWl = function (r, m, showDec = 3) {
    const v = m === 'a' ? r.wAir : r.wVac;
    if (v === null || v === undefined) return '<span class="muted">—</span>';
    const native = r.nat === m;
    const valid = r.dec;                       // decimals supported by the native source value
    if (native && r.natStr) return fmtSrc(r.natStr, Math.max(showDec, 4));
    // the source itself gives this medium too (e.g. NIST ASD air AND vacuum): show it verbatim
    if (m === 'a' && r.sa) return fmtSrc(r.sa, Math.max(showDec, 4));
    if (m === 'v' && r.sv) return fmtSrc(r.sv, Math.max(showDec, 4));
    const dec = Math.max(showDec, Math.min(valid, 5));
    return fmtWl(v, dec, valid, Math.max(showDec, 4));
  };

  // ---------------------------------------------------------------- element data
  function loadElements() {
    const E = J('data-elements');
    const cols = E.cols, dicts = E.dicts || {};
    const idx = Object.fromEntries(cols.map((c, i) => [c, i]));
    const rows = E.rows.map(a => {
      const r = {};
      for (const c of cols) { const v = a[idx[c]]; if (v !== null && v !== undefined) r[c] = dicts[c] ? dicts[c][v] : v; }
      r.sp = r.el + ' ' + r.ion;
      r.nat = r.med;
      r.natStr = r.w;
      r.dec = r.w.includes('.') ? r.w.split('.')[1].length : 0;
      const w = parseFloat(r.w);
      r.wAir = r.med === 'a' ? w : null;
      r.wVac = r.med === 'v' ? w : r.vac;
      r.air = r.wAir; r.vac = r.wVac;
      r.wn = 1e8 / r.wVac;
      if (r.aw === undefined && r.src && r.src.startsWith('ASD')) r.aw = r.w;
      r.g = grade(r.u);
      return r;
    });
    APP.datasets.elements = { rows, elements: E.elements, stats: E.stats };
    APP.elBySym = Object.fromEntries(E.elements.map(e => [e.symbol, e]));
    return rows;
  }

  // periodic-table layout: [Z, row, col]
  function ptPos(Z) {
    const main = { 1: [1, 1], 2: [1, 18], 3: [2, 1], 4: [2, 2], 5: [2, 13], 6: [2, 14], 7: [2, 15], 8: [2, 16], 9: [2, 17], 10: [2, 18],
      11: [3, 1], 12: [3, 2], 13: [3, 13], 14: [3, 14], 15: [3, 15], 16: [3, 16], 17: [3, 17], 18: [3, 18] };
    if (main[Z]) return main[Z];
    if (Z >= 19 && Z <= 36) return [4, Z - 18];
    if (Z >= 37 && Z <= 54) return [5, Z - 36];
    if (Z === 55 || Z === 56) return [6, Z - 54];
    if (Z >= 57 && Z <= 71) return [9, Z - 57 + 3];
    if (Z >= 72 && Z <= 86) return [6, Z - 72 + 4];
    if (Z === 87 || Z === 88) return [7, Z - 86];
    if (Z >= 89 && Z <= 103) return [10, Z - 89 + 3];
    if (Z >= 104 && Z <= 118) return [7, Z - 104 + 4];
    return null;
  }

  // ---- spectroscopic notation helpers -------------------------------------------------
  const CORES = [['[Rn]', '1s22s22p63s23p63d104s24p64d105s25p64f145d106s26p6'], ['[Xe]', '1s22s22p63s23p63d104s24p64d105s25p6'],
    ['[Kr]', '1s22s22p63s23p63d104s24p6'], ['[Ar]', '1s22s22p63s23p6'], ['[Ne]', '1s22s22p6'], ['[He]', '1s2']];
  /** ground state from the structured Handbook fields {ground, conf, term, core, val, note} -> HTML */
  function fmtGround(sp) {
    if (typeof sp === 'string' || !sp.val) {
      const g = typeof sp === 'string' ? sp : sp.ground;
      return esc(g);
    }
    const sub = sp.val.split('.').map(t => {
      const m = t.match(/^(\d+)([a-z])(\d*)(?:\((.*)\))?$/);
      if (!m) return esc(t);
      return `${m[1]}${m[2]}${m[3] ? `<sup>${m[3]}</sup>` : ''}${m[4] ? `(${fmtTerm(m[4])})` : ''}`;
    }).join('');
    return `${esc(sp.core || '')}${sub} ${fmtTerm(sp.term)}${sp.note ? ` <span class="flag" title="${esc(sp.note)}">note</span>` : ''}`;
  }
  /** term "5D4", "2P*3/2", "6D9/2" -> HTML with superscript multiplicity, ° for odd parity, subscript J */
  function fmtTerm(t) {
    const m = String(t ?? '').trim().match(/^(\d+)([A-Z])(\*|°)?\s*([\d/]+)?$/);
    if (!m) return esc(String(t ?? '').replace(/\*/g, '°'));
    return `<sup>${m[1]}</sup>${m[2]}${m[3] ? '°' : ''}${m[4] ? `<sub>${m[4]}</sub>` : ''}`;
  }
  /** ASD level "3d7.(4F).4s a 5F 5" / "3s.3p 3P* 2" -> pretty HTML */
  function fmtLevel(L) {
    if (!L) return '';
    const parts = String(L).split(' ');
    const out = parts.map((p, i) => {
      if (i === 0) return esc(p).replace(/([spdfghi])(\d+)/g, '$1<sup>$2</sup>').replace(/\((\d)([A-Z])(\*?)(<sup>)?/g, (x, a, b, c) => `(<sup>${a}</sup>${b}${c ? '°' : ''}`).replace(/\./g, '<span class="muted">.</span>');
      const t = p.match(/^(\d+)([A-Z])(\*?)$/);
      if (t) return `<sup>${t[1]}</sup>${t[2]}${t[3] ? '°' : ''}`;
      if (/^[\d/]+$/.test(p) && i === parts.length - 1) return `<sub>J=${p}</sub>`;
      return esc(p).replace(/\*/g, '°');
    });
    return out.join(' ');
  }
  APP.fmtLevel = fmtLevel; APP.fmtGround = fmtGround;

  function lineInfoHtml(r) {
    const lw = APP.state.med === 'air' && r.wAir ? r.wAir : r.wVac;
    return `<b>${esc(r.sp)}</b> &nbsp;<b>${lw.toFixed(Math.max(3, Math.min(r.dec, 4)))}</b> Å ${APP.state.med === 'air' && r.wAir ? 'air' : 'vac'}` +
      `<br>±${UI.fmtU(r.u)} mÅ (${r.us ? 'stated' : 'implied: 2.5 units of the last digit'}) · ${esc(r.src)}${r.bl ? ` · unresolved multiplet (${r.bl} components)` : ''}` +
      (r.int ? `<br>Handbook intensity ${r.int}${r.fl ? ' ' + esc(r.fl) : ''}` : '') +
      (r.lo ? `<br>${fmtLevel(r.lo)} → ${fmtLevel(r.up)}` : '') + (r.A ? `<br>A<sub>ki</sub> = ${fmtSig(r.A, 3)} s⁻¹ ${esc(r.acc || '')}` : '');
  }
  APP.lineInfoHtml = lineInfoHtml;

  function elementColumns() {
    return [
      { key: 'air', label: 'λ<sub>air</sub> (Å)', cls: 'num', title: 'standard-air wavelength (λ ≥ 2000 Å). Native values verbatim; converted values faded beyond the source precision', html: r => APP.fmtRowWl(r, 'a'), sort: r => r.wAir ?? r.wVac },
      { key: 'vac', label: 'λ<sub>vac</sub> (Å)', cls: 'num', title: 'vacuum wavelength (converted with Peck & Reeder 1972 where native is air)', html: r => APP.fmtRowWl(r, 'v'), sort: r => r.wVac },
      { key: 'u', label: 'σ<sub>λ</sub> (mÅ)', cls: 'num', title: 'uncertainty: stated by ASD, or implied = 2.5 units of the last published digit (NIST ASD: implied uncertainties are generally 2.5–25 units)', html: r => `${UI.fmtU(r.u)} ${gradeBadge(r.g, r.us ? 'stated uncertainty' : 'implied: 2.5 units of the last digit')}${r.dx !== undefined ? ` <span class="flag" style="color:var(--warn)" title="the other source (Handbook vs ASD) differs by ${r.dx} mÅ — more than twice the combined uncertainty">⚠</span>` : ''}`, sort: r => r.u },
      { key: 'wn', label: 'σ̃ (cm⁻¹)', cls: 'num', title: 'vacuum wavenumber 10⁸/λ_vac', html: r => r.wn.toFixed(r.dec >= 4 ? 4 : 3), sort: r => r.wn },
      { key: 'sp', label: 'Spectrum', html: r => `<b>${esc(r.sp)}</b>`, sort: r => r.Z * 10 + (r.ion === 'II' ? 1 : 0) },
      { key: 'int', label: 'Int', cls: 'num', title: 'Handbook relative intensity (strongest line of each spectrum = 1000). P = persistent line', html: r => `${r.int ?? ''}${r.fl ? ` <span class="flag">${esc(r.fl)}</span>` : ''}`, sort: r => r.int ?? -1 },
      { key: 'src', label: 'Adopted from', html: r => `<span class="small">${esc(r.src)}</span>`, sort: r => r.src },
      { key: 'hb', label: 'Handbook λ', cls: 'num', title: 'Handbook value verbatim (native medium)', html: r => fmtSrc(r.hb, 4), sort: r => parseFloat(r.hb) },
      { key: 'dl', label: 'Δ<sub>ASD−HB</sub> mÅ', cls: 'num', title: 'ASD minus Handbook (closest ASD value); match class e = exact (rounding), c = close (re-measured), u = loose but unique (ASD adopted), l = loose (tentative, Handbook kept); mult = unresolved multiplet (for mult rows kept at the Handbook value, Δ is the offset of the strongest component)', html: r => (r.m ? `${r.dl > 0 ? '+' : ''}${(+r.dl).toFixed(Math.abs(r.dl) < 1 ? 2 : 1)} <span class="flag">${r.m}</span>` : '<span class="flag">no ASD</span>') + (r.bl ? ` <span class="flag" title="unresolved multiplet of ${r.bl} fine-structure components; the Handbook/observed blend value is kept">mult</span>` : ''), sort: r => r.m ? Math.abs(r.dl) : 1e9 },
      { key: 'lo', label: 'Transition (lower → upper)', html: r => r.lo ? `<span class="small">${fmtLevel(r.lo)} → ${fmtLevel(r.up)}</span>` : (r.lv ? `<span class="small muted">${esc(r.lv.join(' → '))}</span>` : ''), sort: r => r.lo || '' },
      { key: 'Ek', label: 'E<sub>i</sub> – E<sub>k</sub> (cm⁻¹)', cls: 'num', html: r => r.Ei !== undefined && r.Ei !== '' ? `<span class="small">${esc(r.Ei)} – ${esc(r.Ek)}</span>` : '', sort: r => parseFloat(r.Ei) },
      { key: 'A', label: 'A<sub>ki</sub> (s⁻¹)', cls: 'num', html: r => r.A ? `${fmtSig(r.A, 3)} <span class="flag">${esc(r.acc || '')}</span>` : (r.hA ? `${esc(r.hA)}e8 <span class="flag">HB</span>` : ''), sort: r => r.A ?? (r.hA ? +r.hA * 1e8 : -1) },
      { key: 'ref', label: 'Ref', html: r => `<span class="small">${esc(r.ref || '')}${r.lr ? ' · ' + esc(r.lr) : ''}</span>`, sort: r => r.ref },
    ];
  }
  const csvCols = [
    ['spectrum', r => r.sp], ['lambda_air_A', r => r.wAir ? r.wAir.toFixed(5) : ''], ['lambda_vac_A', r => r.wVac.toFixed(5)],
    ['native_medium', r => r.nat === 'a' ? 'air' : 'vacuum'], ['native_value_verbatim', r => r.w], ['uncert_A', r => r.u], ['uncert_stated', r => r.us],
    ['wavenumber_cm-1', r => r.wn.toFixed(4)], ['multiplet_components', r => r.bl || ''], ['handbook_intensity', r => r.int], ['flags', r => r.fl], ['adopted_from', r => r.src], ['handbook_lambda', r => r.hb],
    ['asd_minus_handbook_mA', r => r.dl ?? ''], ['match', r => r.m], ['other_source_minus_adopted_mA_if_gt_2sigma', r => r.dx ?? ''], ['lower', r => r.lo || ''], ['upper', r => r.up || ''], ['Ei_cm-1', r => r.Ei || ''],
    ['Ek_cm-1', r => r.Ek || ''], ['Aki_s-1', r => r.A || ''], ['Aki_acc', r => r.acc || ''], ['handbook_ref', r => r.ref], ['asd_line_ref', r => r.lr || ''],
  ];

  function initElements() {
    const rows = loadElements();
    const els = APP.datasets.elements.elements;
    const counts = {};
    for (const r of rows) counts[r.el] = (counts[r.el] || 0) + 1;
    // periodic table
    const pt = $('#ptable');
    const NAMES = Object.fromEntries(els.map(e => [e.Z, e]));
    for (let Z = 1; Z <= 103; Z++) {
      const p = ptPos(Z); if (!p) continue;
      const e = NAMES[Z];
      const b = h('button', { class: 'el' + (e ? '' : ' none'), style: `grid-row:${p[0]};grid-column:${p[1]}`, title: e ? `${e.name} — ${counts[e.symbol] || 0} lines` : 'not in Handbook' },
        h('span', { class: 'z' }, String(Z)), h('span', { class: 's' }, e ? e.symbol : '·'), h('span', { class: 'n' }, e ? String(counts[e.symbol] || 0) : ''));
      if (e) b.addEventListener('click', () => selectElement(e.symbol, true));
      b.dataset.sym = e ? e.symbol : '';
      pt.append(b);
    }
    pt.append(h('div', { class: 'lbl', style: 'grid-row:9;grid-column:1/3' }, 'Lanthanides'));
    pt.append(h('div', { class: 'lbl', style: 'grid-row:10;grid-column:1/3' }, 'Actinides'));
    pt.append(h('div', { style: 'grid-row:8;grid-column:1;height:6px' }));

    // controls
    const st = { sym: store.get('el', 'Fe'), ion: 'all', minInt: 0, onlyP: false, grade: 'all', all: false };
    APP.elState = st;
    const ionSel = h('select', { 'aria-label': 'ionisation stage' }, h('option', { value: 'all' }, 'I + II'), h('option', { value: 'I' }, 'I (neutral)'), h('option', { value: 'II' }, 'II (singly ionised)'));
    const minInt = h('input', { type: 'number', min: 0, max: 1000, step: 10, value: 0, title: 'minimum Handbook intensity' });
    const onlyP = h('input', { type: 'checkbox' });
    const gSel = h('select', { 'aria-label': 'precision grade' }, h('option', { value: 'all' }, 'any precision'), h('option', { value: 'A' }, 'A: ≤ 1 mÅ'), h('option', { value: 'AB' }, 'A–B: ≤ 5 mÅ'), h('option', { value: 'CD' }, 'C–D: > 5 mÅ'));
    const allBox = h('input', { type: 'checkbox' });
    const wlo = h('input', { type: 'number', step: 'any', placeholder: 'λ min', title: 'minimum wavelength (Å, display medium)' });
    const whi = h('input', { type: 'number', step: 'any', placeholder: 'λ max', title: 'maximum wavelength (Å, display medium)' });
    const filt = r => (st.all || r.el === st.sym) && (st.ion === 'all' || r.ion === st.ion) && ((r.int || 0) >= st.minInt) &&
      (!st.onlyP || r.P || (r.fl || '').includes('P')) &&
      (st.grade === 'all' || (st.grade === 'A' ? r.g === 'A' : st.grade === 'AB' ? (r.g === 'A' || r.g === 'B') : (r.g === 'C' || r.g === 'D'))) &&
      (!wlo.value || APP.wl(r) >= +wlo.value) && (!whi.value || APP.wl(r) <= +whi.value);
    const dt = new DataTable('#elTable', {
      columns: elementColumns(), rows, pageSize: 150, sort: { key: 'vac', desc: false },
      searchKeys: ['sp', 'lo', 'up', 'ref', 'src', 'w'], searchPlaceholder: 'filter: term, config, ref…',
      extraFilters: filt,
      rowClass: r => r.P ? 'pl' : '',
      toolbar: bar => {
        const lab = (t, el) => h('label', { class: 'small' }, t, ' ', el);
        bar.append(lab('Stage', ionSel), lab('Int ≥', minInt), lab('', gSel), wlo, whi, h('label', { class: 'small' }, onlyP, ' persistent only'), h('label', { class: 'small', title: 'search every element at once' }, allBox, ' all elements'));
      },
      csv: { name: 'element_lines.csv', cols: csvCols },
      onRow: r => { APP.elSpec.setRange(APP.wl(r) - 3, APP.wl(r) + 3); APP.elSpec.markers = [{ w: APP.wl(r) }]; APP.elSpec.draw(); },
    });
    const upd = () => { dt.page = 0; dt.refresh(); drawElSpec(); };
    ionSel.addEventListener('change', () => { st.ion = ionSel.value; upd(); });
    minInt.addEventListener('input', () => { st.minInt = +minInt.value || 0; upd(); });
    onlyP.addEventListener('change', () => { st.onlyP = onlyP.checked; upd(); });
    gSel.addEventListener('change', () => { st.grade = gSel.value; upd(); });
    allBox.addEventListener('change', () => { st.all = allBox.checked; upd(); });
    wlo.addEventListener('change', upd); whi.addEventListener('change', upd);
    APP.elTable = dt;

    // element spectrum strip
    APP.elSpec = new SpectrumView('#elSpec', { height: 190, range: [2000, 10000], min: 150, max: 130000, medium: () => APP.state.med });
    function drawElSpec() {
      const sel = rows.filter(r => r.el === st.sym && ((r.int || 0) >= st.minInt));
      const mk = ion => {
        const L = sel.filter(r => r.ion === ion).sort((a, b) => APP.wl(a) - APP.wl(b));
        return {
          name: `${st.sym} ${ion}`, sub: `${L.length} lines`, color: 'rainbow',
          lines: L.map(r => APP.wl(r)), str: L.map(r => r.int ? Math.max(0.08, Math.log10(r.int + 1) / 3) : 0.1),
          info: k => lineInfoHtml(L[k]), label: k => APP.wl(L[k]).toFixed(3),
        };
      };
      const tr = [mk('I'), mk('II')].filter(t => t.lines.length);
      APP.elSpec.o.height = 52 + tr.length * 66;
      APP.elSpec.tracks = tr; APP.elSpec.resize();
    }
    APP.drawElSpec = drawElSpec;

    function selectElement(sym, scroll) {
      st.sym = sym; store.set('el', sym);
      $$('#ptable .el').forEach(b => b.classList.toggle('on', b.dataset.sym === sym));
      const e = APP.elBySym[sym];
      const L = rows.filter(r => r.el === sym);
      const gA = L.filter(r => r.g === 'A').length;
      const sp = Object.entries(e.spectra || {}).map(([ion, s]) =>
        `<span class="chip">${esc(sym)} ${ion}: ground ${fmtGround(s)} · IE ${esc(s.ie_eV)} eV (${esc(s.ie_cm)} cm⁻¹)</span>`).join('');
      $('#elHead').innerHTML = `<span class="big">${esc(e.name)} (${esc(sym)})</span><span>Z = ${e.Z}</span><span>A<sub>r</sub> = ${esc(e.weight || '')}</span>` +
        `<span>${L.length} lines · ${gA} at ≤ 1 mÅ (${(100 * gA / Math.max(1, L.length)).toFixed(0)} %)</span><div class="chips">${sp}</div>`;
      const ext = L.map(r => APP.wl(r));
      if (ext.length) {
        const lo = Math.min(...ext), hi = Math.max(...ext);
        APP.elSpec.range = [Math.max(150, lo - (hi - lo) * 0.02), hi + (hi - lo) * 0.02];
      }
      APP.elSpec.markers = [];
      dt.page = 0; dt.refresh(); drawElSpec();
      if (scroll) $('#elHead').scrollIntoView({ block: 'nearest' });
    }
    APP.selectElement = selectElement;
    selectElement(st.sym, false);
    window.addEventListener('medchange', () => { dt.refresh(); drawElSpec(); });
    window.addEventListener('themechange', () => APP.elSpec.draw());

    // global precision stats
    const s = APP.datasets.elements.stats;
    const tot = rows.length;
    $('#elStats') && ($('#elStats').innerHTML = [
      ['Handbook lines', tot.toLocaleString()],
      ['≤ 1 mÅ, Handbook as published', `${(s['before:A'] || 0).toLocaleString()} (${(100 * (s['before:A'] || 0) / tot).toFixed(0)} %)`],
      ['≤ 1 mÅ after ASD upgrade', `${(s['after:A'] || 0).toLocaleString()} (${(100 * (s['after:A'] || 0) / tot).toFixed(0)} %)`],
      ['still only 0.01 Å (no finer value anywhere)', `${rows.filter(r => !r.us && r.dec <= 2).length.toLocaleString()}`],
      ['adopted from ASD', `${((s['src:ASD Ritz'] || 0) + (s['src:ASD obs'] || 0)).toLocaleString()}`],
    ].map(([l, v]) => `<div class="kpi"><div class="v">${v}</div><div class="l">${l}</div></div>`).join(''));
  }

  // ---------------------------------------------------------------- registry of all line datasets (for explorer / finder)
  APP.lineSets = [];   // {key, name, color, rows (normalised: wAir,wVac,nat,natStr,dec,label,str,info), mode}
  APP.registerLineSet = function (ls) { APP.lineSets.push(ls); };

  function initExplorer() {
    const view = new SpectrumView('#explorer', { height: 60, range: [5880, 5900], min: 150, max: 130000, medium: () => APP.state.med, onRange: r => syncInputs(r) });
    APP.explorer = view;
    const lo = $('#exLo'), hi = $('#exHi');
    function syncInputs(r) { const d = r[1] - r[0] < 0.1 ? 4 : 3; lo.value = r[0].toFixed(d); hi.value = r[1].toFixed(d); }
    const checks = $('#exSets');
    const on = store.get('exsets', null);
    function build() {
      const tracks = [];
      for (const ls of APP.lineSets) {
        if (!ls.cb) {
          ls.cb = h('input', { type: 'checkbox' }); ls.cb.checked = on ? on.includes(ls.key) : ls.defaultOn !== false;
          ls.cb.addEventListener('change', () => { store.set('exsets', APP.lineSets.filter(x => x.cb.checked).map(x => x.key)); build(); });
          checks.append(h('label', { class: 'small', style: 'white-space:nowrap' }, ls.cb, ' ', h('i', { style: `display:inline-block;width:10px;height:10px;border-radius:2px;background:${ls.color === 'rainbow' ? 'linear-gradient(90deg,#7a3cff,#00a0ff,#2ec27e,#ffd400,#ff3b30)' : ls.color}` }), ' ', ls.name));
        }
        if (!ls.cb.checked) continue;
        const R = ls.rows.slice().sort((a, b) => APP.wl(a) - APP.wl(b));
        const f = ls.lab ? 1 + (APP.state.zShift || 0) : 1;
        tracks.push({
          name: (ls.short || ls.name) + (f !== 1 ? ' (shifted)' : ''), sub: `${R.length.toLocaleString()} lines`, color: ls.color, mode: ls.mode,
          lines: R.map(r => APP.wl(r) * f), str: R.map(r => r.str ?? 0.6),
          info: k => ls.info(R[k]) + (f !== 1 ? `<br><b>shifted: ${(APP.wl(R[k]) * f).toFixed(4)} Å</b> (z = ${(f - 1).toExponential(4)})` : ''),
          label: k => ls.label ? ls.label(R[k]) : APP.wl(R[k]).toFixed(3),
        });
      }
      view.o.height = 52 + tracks.length * 56;
      view.tracks = tracks; view.resize();
    }
    APP.buildExplorer = build;
    $('#exGo').addEventListener('click', () => { if (lo.value.trim() && hi.value.trim()) view.setRange(parseFloat(lo.value), parseFloat(hi.value)); });
    $$('[data-exr]').forEach(b => b.addEventListener('click', () => { const [a, c] = b.dataset.exr.split(',').map(Number); view.setRange(a, c); }));
    const sv = $('#exShiftV'), sz = $('#exShiftZ');
    const applyShift = src => {
      if (src === 'v') sz.value = ((parseFloat(sv.value) || 0) / S.C).toPrecision(8);
      else sv.value = ((parseFloat(sz.value) || 0) * S.C).toFixed(3);
      APP.state.zShift = parseFloat(sz.value) || 0;
      $('#exShiftMsg').textContent = APP.state.zShift ? `atomic lab lines × (1 + z) = × ${(1 + APP.state.zShift).toFixed(9)} (optical convention z = v/c)` : '';
      build();
    };
    sv.addEventListener('change', () => applyShift('v')); sz.addEventListener('change', () => applyShift('z'));
    $('#exZin').addEventListener('click', () => { const [a, b] = view.range, m = (a + b) / 2, s = (b - a) / 4; view.setRange(m - s, m + s); });
    $('#exZout').addEventListener('click', () => { const [a, b] = view.range, m = (a + b) / 2, s = (b - a); view.setRange(m - s, m + s); });
    window.addEventListener('medchange', build);
    window.addEventListener('themechange', () => view.draw());
    syncInputs(view.range);
    return build;
  }

  // ---------------------------------------------------------------- finding list (wavelength search across all sets)
  function initFinder() {
    const inp = $('#qWl'), tolIn = $('#qTol'), unit = $('#qUnit');
    const res = $('#finderOut');
    const dt = new DataTable(res, {
      columns: [
        { key: 'set', label: 'Catalogue', html: r => `<span class="small">${esc(r.set)}</span>` },
        { key: 'id', label: 'Species / identification', html: r => `<b>${esc(r.id)}</b>` },
        { key: 'wAir', label: 'λ<sub>air</sub> (Å)', cls: 'num', html: r => APP.fmtRowWl(r.r, 'a'), sort: r => r.r.wAir ?? r.r.wVac },
        { key: 'wVac', label: 'λ<sub>vac</sub> (Å)', cls: 'num', html: r => APP.fmtRowWl(r.r, 'v'), sort: r => r.r.wVac },
        { key: 'd', label: 'Δλ (mÅ)', cls: 'num', html: r => (r.d * 1000).toFixed(1), sort: r => Math.abs(r.d) },
        { key: 'dv', label: 'Δv (km/s)', cls: 'num', html: r => (r.d / r.w0 * S.C).toFixed(3), sort: r => Math.abs(r.d) },
        { key: 'u', label: 'σ<sub>λ</sub> (mÅ)', cls: 'num', html: r => r.r.u !== undefined && r.r.u !== null ? (r.r.u * 1000).toFixed(r.r.u < 0.001 ? 2 : 1) : '', sort: r => r.r.u },
        { key: 'x', label: 'Notes', html: r => `<span class="small">${r.note || ''}</span>` },
      ], rows: [], pageSize: 60, sort: { key: 'd', desc: false }, tall: false,
      onRow: r => { APP.explorer.setRange(r.w0 - 2, r.w0 + 2); APP.explorer.markers = [{ w: r.w0 }]; APP.explorer.draw(); $('#explorer').scrollIntoView({ block: 'center' }); },
    });
    function run(wTxt, keepView) {
      const w0 = parseFloat(wTxt ?? inp.value);
      if (!isFinite(w0)) { $('#finderMsg').textContent = 'enter a wavelength in Å'; return; }
      const tol = Math.abs(parseFloat(tolIn.value));
      if (!(tol > 0)) { $('#finderMsg').textContent = 'enter a tolerance > 0'; return; }
      const tolA = unit.value === 'kms' ? w0 * tol / S.C : unit.value === 'mA' ? tol / 1000 : tol;
      const out = [];
      for (const ls of APP.lineSets) for (const r of ls.rows) {
        const w = APP.wl(r); const d = w - w0;
        if (Math.abs(d) <= tolA) out.push({ set: ls.name, id: ls.ident(r), r, d, w0, note: ls.note ? ls.note(r) : '' });
      }
      dt.setRows(out);
      $('#finderMsg').textContent = `${out.length} line(s) within ±${tolA >= 0.1 ? tolA.toFixed(3) : (tolA * 1000).toFixed(1) + ' m'}Å of ${w0} Å (${APP.medLabel()})`;
      APP.explorer.markers = [{ w: w0 }];
      if (keepView) APP.explorer.draw(); else APP.explorer.setRange(w0 - Math.max(tolA * 1.5, 0.5), w0 + Math.max(tolA * 1.5, 0.5));
    }
    $('#qGo').addEventListener('click', () => run());
    inp.addEventListener('keydown', e => { if (e.key === 'Enter') run(); });
    const top = $('#topQ');
    top.addEventListener('keydown', e => { if (e.key === 'Enter') { inp.value = top.value; run(top.value); $('#finder').scrollIntoView(); } });
    window.addEventListener('medchange', () => { if (inp.value) run(undefined, true); });
    APP.find = run;
  }

  APP.boot = function () {
    initTheme();
    $$('[data-med]').forEach(b => b.addEventListener('click', () => APP.setMed(b.dataset.med)));
    APP.setMed(APP.state.med);
    initElements();
    APP.registerLineSet({
      key: 'el', lab: true, name: 'NIST atomic lines (Handbook + ASD)', short: 'Atomic I/II', color: 'rainbow', rows: APP.datasets.elements.rows.map(r => Object.assign(r, { str: r.int ? Math.max(0.08, Math.log10(r.int + 1) / 3) : 0.1 })),
      info: lineInfoHtml, ident: r => r.sp, label: r => r.sp, note: r => `Int ${r.int ?? ''} ${esc(r.fl || '')} · ${esc(r.src)}`, defaultOn: true,
    });
    for (const f of (APP.initSections || [])) { try { f(); } catch (e) { console.error(e); } }
    initExplorer();
    APP.buildExplorer();
    initFinder();
    APP.find();                       // demo query (Na D2) so the cross-catalogue view is populated on load
    APP.booted = true;
  };
  // self-test hook for headless checks:  index.html?selftest
  const ERR = [];
  window.addEventListener('error', e => ERR.push(String(e.message) + ' @' + (e.lineno || '')));
  window.addEventListener('unhandledrejection', e => ERR.push('rejection: ' + e.reason));
  const _ce = console.error; console.error = (...a) => { ERR.push(a.map(x => x && x.stack ? x.stack : String(x)).join(' ')); _ce.apply(console, a); };
  document.addEventListener('DOMContentLoaded', () => {
    try { APP.boot(); } catch (e) { ERR.push('boot: ' + (e.stack || e)); }
    if (location.search.includes('selftest')) setTimeout(() => {
      const rep = { errors: ERR, sets: APP.lineSets.map(l => [l.key, l.rows.length]),
        tables: $$('.dt').map(d => [d.id, d.querySelectorAll('tbody tr').length, (d.querySelector('.count') || {}).textContent]),
        bary: ($('#bcOut') || {}).textContent?.slice(0, 300), frames: ($('#frOut') || {}).textContent?.slice(0, 300),
        explorerTracks: APP.explorer ? APP.explorer.tracks.length : -1 };
      const pre = document.createElement('pre'); pre.id = 'selftest'; pre.textContent = JSON.stringify(rep); document.body.append(pre);
      const sc = location.search.match(/scroll=([\w-]+)/);
      if (sc && document.getElementById(sc[1])) {      // headless screenshots: hide everything before the target
        const tgt = document.getElementById(sc[1]); const sec = tgt.closest('section') || tgt;
        for (const s of $$('main section')) { if (s === sec) break; s.style.display = 'none'; }
        if (tgt !== sec) for (const ch of Array.from(sec.children)) { if (ch === tgt || ch.contains(tgt)) break; ch.style.display = 'none'; }
        window.dispatchEvent(new Event('resize'));
      }
    }, 1500);
  });
})();
