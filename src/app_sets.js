/* app_sets.js -- generic catalogue sets (airglow, Fraunhofer, solar, telluric) + static reference tables.
 * Each data block  data-<key>  = {meta:{name, short, color, mode, host, extra:[{k,label,num,title,fmt}] , sort}, cols, rows}
 * normalised columns: a (air A|null), v (vac A), n ('a'|'v' native medium), s (native string), d (valid decimals),
 *                     u (uncertainty A|null), us (1 if stated), id (identification), st (0..1 plot strength) + extras
 */
(function () {
  'use strict';
  const { $, esc, h, fmtSig, grade, gradeBadge, DataTable, SpectrumView } = UI;
  const APP = window.APP;
  (APP.initSections = APP.initSections || []).push(initSets);

  const SET_KEYS = ['ag_atomic', 'ag_palace', 'ag_optical', 'ag_nir', 'fh_classic', 'solar', 'telluric'];

  function loadSet(key) {
    const D = APP.J('data-' + key);
    if (!D) return null;
    const m = D.meta, dicts = D.dicts || {};
    const idx = Object.fromEntries(D.cols.map((c, i) => [c, i]));
    const rows = D.rows.map(a => {
      const r = {};
      for (const c of D.cols) { const v = a[idx[c]]; if (v !== null && v !== undefined) r[c] = dicts[c] ? dicts[c][v] : v; }
      r.nat = r.n;
      r.natStr = r.s ?? (r.n === 'a' ? r.sa : r.sv);
      r.dec = r.natStr && r.natStr.includes('.') ? r.natStr.split('.')[1].replace(/[^\d]/g, '').length : 0;
      r.wAir = r.sa ? parseFloat(r.sa) : (r.a ?? null);
      r.wVac = r.sv ? parseFloat(r.sv) : r.v;
      if (r.wAir === null && r.wVac >= 2000 && r.n === 'v') r.wAir = SpecCore.vacToAir(r.wVac, 'peckreeder1972');
      if (r.wVac === undefined && r.wAir !== null) r.wVac = SpecCore.airToVac(r.wAir, 'peckreeder1972');
      if (!r.id && m.idCols) r.id = m.idCols.map(k => r[k]).filter(Boolean).join(' ');
      if (!r.lab && m.labKey) r.lab = r[m.labKey];
      r.str = r.st ?? (m.stKey ? Math.max(0.04, Math.min(1, r[m.stKey] || 0)) : 0.6);
      if (m.calc === 'telluric' && r.wAir && r.dps !== undefined) r.aobs = r.wAir + r.dps / 1000 * r.wAir / r.wVac;
      r.g = r.u !== null && r.u !== undefined ? grade(r.u) : '';
      return r;
    });
    return { meta: m, rows };
  }

  function fmtExtra(x, r) {
    const v = r[x.k];
    if (v === null || v === undefined || v === '') return '';
    if (x.fmt === 'sig') return fmtSig(v, x.n || 3);
    if (x.fmt === 'fix') return (+v).toFixed(x.n ?? 2);
    if (x.fmt === 'html') return String(v);
    const t = String(v);
    if (t.length > 40) {               // long free text: give it width instead of height; truncate very long notes
      const shown = t.length > 260 ? t.slice(0, 250) + '…' : t;
      return `<span class="lt"${t.length > 260 ? ` title="${esc(t)}"` : ''}>${esc(shown)}</span>`;
    }
    return esc(t);
  }

  const SUP = { 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
  const SUB = { 0: '₀', 1: '₁', 2: '₂', 3: '₃', 4: '₄', 5: '₅', 6: '₆', 7: '₇', 8: '₈', 9: '₉' };
  /** chemical formatting of leading molecule tokens: 16O18O -> ¹⁶O¹⁸O, H2O -> H₂O, CO2 -> CO₂ */
  function chem(id) {
    const toks = String(id).split(' ');
    for (let i = 0; i < Math.min(2, toks.length); i++) {
      const t = toks[i];
      if (!/^(\d{1,2}[A-Z][a-z]?\d*|[A-Z][a-z]?\d*)+(\(\d+\))?$/.test(t) || /^[A-Z][a-z]?$/.test(t) || !/[A-Z]/.test(t)) continue;
      const m = t.match(/^(.*?)(\(\d+\))?$/);
      // two-digit numbers are mass numbers (superscript), single digits are atom counts (subscript)
      toks[i] = m[1].replace(/\d+/g, n => [...n].map(d => (n.length >= 2 ? SUP : SUB)[d]).join('')) + (m[2] || '');
    }
    return toks.join(' ');
  }
  APP.chem = chem;

  function setColumns(meta) {
    const cols = [
      { key: 'a', label: 'λ<sub>air</sub> (Å)', cls: 'num', html: r => APP.fmtRowWl(r, 'a'), sort: r => r.wAir ?? r.wVac },
      { key: 'v', label: 'λ<sub>vac</sub> (Å)', cls: 'num', html: r => APP.fmtRowWl(r, 'v'), sort: r => r.wVac },
      { key: 'u', label: 'σ<sub>λ</sub> (mÅ)', cls: 'num', title: 'stated uncertainty (or bound) of the position, in mÅ', html: r => r.u !== null && r.u !== undefined ? `${UI.fmtU(r.u)} ${gradeBadge(r.g, r.us ? 'stated by the source (or a source-code bound)' : 'implied: 2.5 units of the last published digit')}` : '<span class="muted" title="no position uncertainty given by the source">—</span>', sort: r => r.u ?? 9 },
      { key: 'id', label: meta.idLabel || 'Identification', html: r => `<b>${esc(meta.chem ? chem(r.id) : r.id)}</b>`, sort: r => r.id },
    ];
    for (const x of (meta.extra || [])) cols.push({ key: x.k, label: x.label, cls: x.num ? 'num' : '', title: x.title || '', html: r => fmtExtra(x, r),
      sort: r => { const v = r[x.k]; return (x.num && typeof v === 'string') ? (isFinite(parseFloat(v)) ? parseFloat(v) : null) : v; } });
    return cols;
  }

  function initSets() {
    APP.sets = {};
    for (const key of SET_KEYS) {
      const S = loadSet(key);
      if (!S) continue;
      APP.sets[key] = S;
      const m = S.meta;
      if (m.host && $(m.host)) {
        const cols = setColumns(m);
        const csv = [['lambda_air_A', r => r.wAir ?? ''], ['lambda_vac_A', r => r.wVac], ['native_medium', r => r.nat === 'a' ? 'air' : 'vacuum'],
          ['native_value_verbatim', r => r.natStr], ['uncert_A', r => r.u ?? ''], ['identification', r => r.id]]
          .concat((m.extra || []).map(x => [`${x.k} [${x.label.replace(/<[^>]+>/g, '')}]`, r => (x.k === 'aobs' && typeof r.aobs === 'number') ? r.aobs.toFixed(5) : (r[x.k] ?? '')]));
        const tb = new DataTable(m.host, {
          columns: cols, rows: S.rows, pageSize: m.pageSize || 100, sort: { key: m.sortKey || 'v', desc: !!m.sortDesc },
          searchKeys: ['id'].concat((m.extra || []).filter(x => !x.num).map(x => x.k)), csv: { name: key + '.csv', cols: csv },
          toolbar: m.filters ? (bar, dt) => addFilters(bar, dt, m) : null,
          extraFilters: r => (!S.filt || S.filt(r)),
          onRow: r => { const w = APP.wl(r); APP.explorer.setRange(w - 1.5, w + 1.5); APP.explorer.markers = [{ w }]; APP.explorer.draw(); },
          tall: true,
        });
        S.table = tb;
        S.setFilt = f => { S.filt = f; tb.page = 0; tb.refresh(); };
        window.addEventListener('medchange', () => tb.refresh());
      }
      if (m.register !== false) {
        APP.registerLineSet({
          key, name: m.name, short: m.short, color: m.color, mode: m.mode, rows: S.rows, defaultOn: m.defaultOn !== false,
          ident: r => r.id, label: r => r.lab || r.id,
          note: r => (m.note ? m.note(r) : (m.extra || []).slice(0, 2).map(x => r[x.k] !== undefined && r[x.k] !== null && r[x.k] !== '' ? `${x.label.replace(/<[^>]+>/g, '')}: ${fmtExtra(x, r)}` : '').filter(Boolean).join(' · ')),
          info: r => `<b>${esc(r.id)}</b><br><b>${APP.wl(r).toFixed(Math.min(4, Math.max(3, r.dec)))}</b> Å ${APP.state.med === 'air' && r.wAir ? 'air' : 'vac'}` +
            (r.u ? ` ± ${(r.u * 1000).toFixed(1)} mÅ` : '') + `<br><span class="muted">${esc(m.name)}</span>` +
            (m.extra || []).slice(0, 4).map(x => r[x.k] !== undefined && r[x.k] !== null && r[x.k] !== '' ? `<br>${x.label}: ${fmtExtra(x, r)}` : '').join(''),
        });
      }
    }
    // static reference tables
    for (const el of document.querySelectorAll('script[type="application/json"][id^="data-st_"]')) {
      const key = el.id.slice(5);
      const T = APP.J(el.id);
      if (!T || !T.host || !$(T.host)) continue;
      const host = $(T.host);
      const idx = Object.fromEntries(T.cols.map((c, i) => [c.k, i]));
      new DataTable(host, {
        columns: T.cols.map(c => ({ key: c.k, label: c.label, cls: c.num ? 'num' : '', title: c.title || '', html: r => c.html ? String(r[c.k] ?? '') : esc(r[c.k] ?? ''), sort: r => r[c.k] })),
        rows: T.rows.map(a => Object.fromEntries(T.cols.map((c, i) => [c.k, a[i]]))), pageSize: 200, tall: T.tall ?? false,
        searchKeys: T.search ? T.cols.filter(c => !c.num).map(c => c.k) : [], csv: { name: key + '.csv', cols: T.cols.map(c => [`${c.k} [${String(c.label).replace(/<[^>]+>/g, '')}]`, r => (r[c.k] ?? '').toString().replace(/<[^>]+>/g, '')]) },
      });
      if (T.caption) host.append(h('p', { class: 'small muted', html: T.caption }));
    }
    // airglow and Fraunhofer overview strips
    APP.extraViews = [];
    const mkView = (sel, keys, range) => {
      if (!$(sel)) return;
      const v = new SpectrumView(sel, { height: 60, range, min: 900, max: 26000, medium: () => APP.state.med });
      const draw = () => {
        const tr = [];
        for (const k of keys) {
          const S = APP.sets[k]; if (!S) continue;
          const R = S.rows.slice().sort((a, b) => APP.wl(a) - APP.wl(b));
          tr.push({ name: S.meta.short || S.meta.name, sub: `${R.length.toLocaleString()} lines`, color: S.meta.color, mode: S.meta.mode,
            lines: R.map(APP.wl), str: R.map(r => r.str ?? 0.6), info: i => APP.lineSets.find(l => l.key === k).info(R[i]), label: i => R[i].lab || R[i].id });
        }
        v.o.height = 52 + tr.length * 60; v.tracks = tr; v.resize();
      };
      draw(); window.addEventListener('medchange', draw); window.addEventListener('themechange', () => v.draw());
    };
    mkView('#agSpec', ['ag_atomic', 'ag_palace', 'ag_nir', 'ag_optical'], [5500, 6400]);
    mkView('#fhSpec', ['fh_classic', 'solar', 'telluric'], [6270, 6320]);

    // overview KPIs and contents list, from what is actually loaded
    const all = APP.lineSets.reduce((s, l) => s + l.rows.length, 0);   // includes the element set
    const kp = [[all.toLocaleString(), 'catalogue lines on this page']];
    const el = APP.datasets.elements.rows;
    kp.push([el.length.toLocaleString(), 'atomic I/II lines, Z = 1–99']);
    for (const k of SET_KEYS) if (APP.sets[k]) kp.push([APP.sets[k].rows.length.toLocaleString(), APP.sets[k].meta.name]);
    $('#kpis').innerHTML = kp.map(([v, l]) => `<div class="kpi"><div class="v">${v}</div><div class="l">${esc(l)}</div></div>`).join('');
    const ANCH = { ag_atomic: '#ag-atomic', ag_palace: '#ag-optical', ag_optical: '#ag-optical', ag_nir: '#ag-nir', ag_oh1: '#ag-nir', fh_classic: '#fh-classic', solar: '#fh-solar', telluric: '#fh-telluric' };
    const items = [['#why', 'Why 10⁻³ Å: velocity equivalents, spectrograph centroid precision, 45 physical effects, verdict'],
      ['#elements', `${el.length.toLocaleString()} element lines (NIST Handbook, upgraded with NIST ASD), with transitions and A-values`]];
    for (const k of SET_KEYS) if (APP.sets[k]) items.push([ANCH[k] || '#', `${APP.sets[k].rows.length.toLocaleString()} — ${APP.sets[k].meta.name}`]);
    items.push(['#ag-bands', 'Molecular airglow/aurora band systems (heads and origins)']);
    items.push(['#cal-airvac', 'Air ↔ vacuum converter comparing eight standard formulas plus the full Ciddor model (T, P, humidity, CO₂)']);
    items.push(['#cal-bary', 'Barycentric correction calculator (JPL DE440s, Wright & Eastman 2014), validated against astropy to < 1 mm/s; forward model of an observed wavelength']);
    items.push(['#cal-lg', 'LSR, Galactocentric, Local Group and CMB frame corrections']);
    items.push(['#cal-cosmo', 'Redshift composition and redshift-drift (Sandage–Loeb) calculators']);
    items.push(['#cal-sys', 'Line-shift checklist']);
    $('#contentsList').innerHTML = items.map(([a, t]) => `<li><a href="${a}">${esc(t)}</a></li>`).join('');
  }

  function addFilters(bar, dt, m) {
    const S = () => APP.sets[Object.keys(APP.sets).find(k => APP.sets[k].meta === m)];
    for (const f of m.filters) {
      const sel = h('select', { 'aria-label': f.label }, h('option', { value: '' }, f.type === 'min' ? f.label + ' any' : f.label + ': all'), ...f.values.map(v => h('option', { value: v }, f.type === 'min' ? `${f.label} ${v}` : v)));
      sel.addEventListener('change', () => {
        const s = S(); s.fv = s.fv || {}; s.fv[f.k] = { v: sel.value, min: f.type === 'min' };
        s.setFilt(r => Object.entries(s.fv).every(([k, c]) => !c.v || (c.min ? (r[k] ?? -Infinity) >= +c.v : String(r[k]) === c.v)));
      });
      bar.append(sel);
    }
  }
})();
