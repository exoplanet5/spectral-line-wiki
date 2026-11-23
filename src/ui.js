/* ui.js -- generic UI components: helpers, wavelength formatting, DataTable, SpectrumView, LineChart */
(function (root) {
  'use strict';
  const $ = (s, el = document) => el.querySelector(s);
  const $$ = (s, el = document) => Array.from(el.querySelectorAll(s));
  const esc = s => String(s ?? '').replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  const h = (tag, attrs = {}, ...kids) => {
    const e = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs || {})) {
      if (k === 'class') e.className = v; else if (k === 'html') e.innerHTML = v;
      else if (k.startsWith('on')) e.addEventListener(k.slice(2), v); else if (v !== null && v !== undefined && v !== false) e.setAttribute(k, v);
    }
    for (const k of kids.flat()) if (k !== null && k !== undefined) e.append(k instanceof Node ? k : document.createTextNode(k));
    return e;
  };
  const store = {
    get(k, d) { try { const v = localStorage.getItem('slw:' + k); return v === null ? d : JSON.parse(v); } catch (e) { return d; } },
    set(k, v) { try { localStorage.setItem('slw:' + k, JSON.stringify(v)); } catch (e) { /* ignore */ } },
  };

  /** Format a wavelength. `dec` = number of decimals to show, `valid` = decimals supported by the source
   * (digits beyond `valid` are faded).  `pad` = total decimals to align to (invisible padding). */
  function fmtWl(v, dec = 3, valid = dec, pad = dec) {
    if (v === null || v === undefined || !isFinite(v)) return '';
    let s = v.toFixed(dec);
    let out;
    if (valid < dec) {
      const cut = s.length - (dec - valid);
      out = esc(s.slice(0, cut)) + '<span class="xd">' + s.slice(cut) + '</span>';
    } else out = esc(s);
    if (pad > dec) out += '<span class="pad">' + '0'.repeat(pad - dec) + '</span>';
    return out;
  }
  /** Show a source wavelength string verbatim, padded to `pad` decimals */
  function fmtSrc(str, pad = 4) {
    if (!str) return '';
    const d = str.includes('.') ? str.split('.')[1].length : 0;
    return esc(str) + (pad > d ? '<span class="pad">' + (d === 0 ? '.' : '') + '0'.repeat(pad - d) + '</span>' : '');
  }
  function fmtNum(v, d = 3) { return v === null || v === undefined || !isFinite(v) ? '' : v.toFixed(d); }
  function fmtSig(v, n = 3) {
    if (v === null || v === undefined || !isFinite(v)) return '';
    if (v === 0) return '0';
    const a = Math.abs(v);
    if (a >= 1e-3 && a < 1e5) return v.toPrecision(n).replace(/\.?0+$/, m => m.includes('.') ? '' : m);
    return v.toExponential(n - 1);
  }
  /** uncertainty in A -> compact mA string */
  function fmtU(uA) {
    if (uA === null || uA === undefined || !isFinite(uA)) return '';
    const m = uA * 1000;
    if (m === 0) return '0';
    if (m < 0.01) return m.toExponential(0).replace('e-', '×10⁻').replace(/⁻(\d)/g, (x, d) => '⁻' + '⁰¹²³⁴⁵⁶⁷⁸⁹'[d]).replace(/(\d)$/, d => '⁰¹²³⁴⁵⁶⁷⁸⁹'[d]);
    if (m < 0.995) return m.toFixed(2);
    if (m < 9.95) return m.toFixed(1);
    return m.toFixed(0);
  }
  function grade(uncA) {        // precision grade of a wavelength from its uncertainty / implied precision (A)
    if (uncA === null || uncA === undefined) return '';
    if (uncA <= 0.001 + 1e-9) return 'A';
    if (uncA <= 0.005 + 1e-9) return 'B';
    if (uncA <= 0.010 + 1e-9) return 'C';
    return 'D';
  }
  const gradeBadge = (g, title) => g ? `<span class="badge g-${g}" title="${esc(title || '')}">${g}</span>` : '';

  /** approximate visible colour of a (vacuum/air) wavelength in A -> css rgb */
  function wlColor(wA, alpha = 1) {
    const w = wA / 10;   // nm
    let r = 0, g = 0, b = 0;
    if (w >= 380 && w < 440) { r = -(w - 440) / 60; b = 1; }
    else if (w < 490 && w >= 440) { g = (w - 440) / 50; b = 1; }
    else if (w < 510 && w >= 490) { g = 1; b = -(w - 510) / 20; }
    else if (w < 580 && w >= 510) { r = (w - 510) / 70; g = 1; }
    else if (w < 645 && w >= 580) { r = 1; g = -(w - 645) / 65; }
    else if (w <= 780 && w >= 645) { r = 1; }
    let f = 0;
    if (w >= 380 && w < 420) f = 0.3 + 0.7 * (w - 380) / 40;
    else if (w >= 420 && w < 701) f = 1;
    else if (w >= 701 && w <= 780) f = 0.3 + 0.7 * (780 - w) / 80;
    const gm = x => Math.round(255 * Math.pow(x * f, 0.8));
    if (w < 380 || w > 780) return null;
    return `rgba(${gm(r)},${gm(g)},${gm(b)},${alpha})`;
  }
  const cssVar = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
  function isLightTheme() {
    const m = cssVar('--plot-bg').match(/^#?([0-9a-f]{2})([0-9a-f]{2})([0-9a-f]{2})$/i);
    return m ? (parseInt(m[1], 16) + parseInt(m[2], 16) + parseInt(m[3], 16)) / 3 > 128 : false;
  }
  /** stick colour for a wavelength, darkened on light backgrounds for contrast */
  function stickColor(wA, alpha = 0.95) {
    const c = wlColor(wA, alpha); if (!c) return null;
    if (!isLightTheme()) return c;
    return c.replace(/rgba\((\d+),(\d+),(\d+),/, (x, r, g, b) => `rgba(${Math.round(r * 0.55)},${Math.round(g * 0.55)},${Math.round(b * 0.55)},`);
  }

  function download(name, text, type = 'text/csv') {
    const a = document.createElement('a');
    a.href = URL.createObjectURL(new Blob([text], { type }));
    a.download = name; document.body.append(a); a.click();
    setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 500);
  }

  // ======================================================================= DataTable
  /**
   * new DataTable(host, {columns, rows, pageSize, sort:{key,desc}, searchKeys, extraFilters:(row)=>bool,
   *                      toolbar:(barEl)=>void, csv:{name, cols:[[label, row=>value]]}, rowClass, onRow})
   * column: {key, label, title, cls, html:(row)=>string, sort:(row)=>number|string}
   */
  class DataTable {
    constructor(host, opt) {
      this.host = typeof host === 'string' ? $(host) : host;
      this.o = Object.assign({ pageSize: 100, searchKeys: [] }, opt);
      this.rows = opt.rows || [];
      this.view = this.rows;
      this.page = 0;
      this.sort = opt.sort || null;
      this.q = '';
      this.build();
      this.refresh();
    }
    build() {
      const o = this.o;
      this.host.classList.add('dt');
      this.bar = h('div', { class: 'bar' });
      if (o.searchKeys.length) {
        this.search = h('input', { type: 'search', placeholder: o.searchPlaceholder || 'filter text…', 'aria-label': 'filter' });
        this.search.addEventListener('input', () => { this.q = this.search.value.trim().toLowerCase(); this.page = 0; this.refresh(); });
        this.bar.append(this.search);
      }
      if (o.toolbar) o.toolbar(this.bar, this);
      this.countEl = h('span', { class: 'count' });
      this.bar.append(this.countEl);
      if (o.csv) this.bar.append(h('button', { class: 'small', onclick: () => this.exportCsv(), title: 'download the filtered rows as CSV' }, 'CSV'));
      this.scroll = h('div', { class: 'tscroll' + (o.tall === false ? '' : ' tall') });
      this.table = h('table', { class: 't' });
      this.thead = h('thead'); this.tbody = h('tbody');
      this.table.append(this.thead, this.tbody); this.scroll.append(this.table);
      const tr = h('tr');
      for (const c of o.columns) {
        const th = h('th', { class: c.cls || '', title: c.title || '', html: c.label, tabindex: 0, 'aria-sort': 'none' });
        th.addEventListener('keydown', e => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); th.click(); } });
        th.addEventListener('click', () => {
          if (this.sort && this.sort.key === c.key) this.sort.desc = !this.sort.desc; else this.sort = { key: c.key, desc: false };
          this.page = 0; this.refresh();
        });
        c.th = th; tr.append(th);
      }
      this.thead.append(tr);
      this.pager = h('div', { class: 'pager' });
      this.host.append(this.bar, this.scroll, this.pager);
      if (o.onRow) this.tbody.addEventListener('click', e => {
        const r = e.target.closest('tr'); if (!r || e.target.closest('a')) return;
        const row = this.pageRows[+r.dataset.i]; if (row) o.onRow(row, r, this);
      });
    }
    setRows(rows) { this.rows = rows; this.page = 0; this.refresh(); }
    refresh() {
      const o = this.o, q = this.q;
      let v = this.rows;
      if (o.extraFilters) v = v.filter(o.extraFilters);
      if (q) {
        const terms = q.split(/\s+/);
        v = v.filter(r => {
          const s = (r._s || (r._s = o.searchKeys.map(k => typeof k === 'function' ? k(r) : r[k]).join(' ').toLowerCase()));
          return terms.every(t => s.includes(t));
        });
      }
      if (this.sort) {
        const c = o.columns.find(c => c.key === this.sort.key);
        const f = c && (c.sort || (r => r[c.key]));
        if (f) {
          const d = this.sort.desc ? -1 : 1;
          v = v.slice().sort((a, b) => {
            const x = f(a), y = f(b);
            if (x === y) return 0;
            if (x === null || x === undefined || x === '') return 1;
            if (y === null || y === undefined || y === '') return -1;
            return (x < y ? -1 : 1) * d;
          });
        }
        for (const cc of o.columns) cc.th.classList.toggle('sorted', cc.key === this.sort.key), cc.th.classList.toggle('desc', cc.key === this.sort.key && this.sort.desc);
      }
      this.view = v;
      this.render();
    }
    render() {
      const o = this.o, ps = o.pageSize, n = this.view.length;
      const pages = Math.max(1, Math.ceil(n / ps));
      this.page = Math.min(this.page, pages - 1);
      const a = this.page * ps;
      this.pageRows = this.view.slice(a, a + ps);
      const html = [];
      this.pageRows.forEach((r, i) => {
        const cls = o.rowClass ? o.rowClass(r) : '';
        html.push(`<tr data-i="${i}"${cls ? ` class="${cls}"` : ''}>`);
        for (const c of o.columns) html.push(`<td class="${c.cls || ''}">${c.html ? c.html(r) : esc(r[c.key])}</td>`);
        html.push('</tr>');
      });
      this.tbody.innerHTML = html.join('');
      this.countEl.textContent = `${n.toLocaleString()} rows` + (n !== this.rows.length ? ` of ${this.rows.length.toLocaleString()}` : '');
      this.pager.innerHTML = '';
      if (pages > 1) {
        const mk = (t, p, dis) => h('button', { class: 'small', disabled: dis || null, onclick: () => { this.page = p; this.render(); this.scroll.scrollTop = 0; } }, t);
        this.pager.append(mk('«', 0, this.page === 0), mk('‹', this.page - 1, this.page === 0),
          h('span', {}, ` page ${this.page + 1} / ${pages} `), mk('›', this.page + 1, this.page >= pages - 1), mk('»', pages - 1, this.page >= pages - 1));
      }
    }
    exportCsv() {
      const cols = this.o.csv.cols;
      const q = v => { v = v === null || v === undefined ? '' : String(v); return /[",\n]/.test(v) ? '"' + v.replace(/"/g, '""') + '"' : v; };
      const lines = [cols.map(c => q(c[0])).join(',')];
      for (const r of this.view) lines.push(cols.map(c => q(c[1](r))).join(','));
      download(this.o.csv.name, lines.join('\n') + '\n');
    }
  }

  // ======================================================================= SpectrumView
  /**
   * Multi-track stick-spectrum viewer with zoom (wheel / pinch / buttons) and pan (drag).
   * tracks: [{name, color, lines:Float64Array|[] of wl, str:[] (0..1), info:(i)=>html, mode:'emission'|'absorption'}]
   */
  class SpectrumView {
    constructor(host, opt) {
      this.host = typeof host === 'string' ? $(host) : host;
      this.o = Object.assign({ height: 260, min: 900, max: 25000, rainbow: true }, opt);
      this.tracks = [];
      this.range = opt.range || [3800, 7000];
      this.host.classList.add('plot');
      this.cv = h('canvas', { role: 'img', 'aria-label': 'stick spectrum; click the plot, then scroll to zoom and drag to pan', tabindex: 0 }); this.tip = h('div', { class: 'tip' });
      (root.UI_VIEWS = root.UI_VIEWS || []).push(this);
      this.host.append(this.cv, this.tip);
      this.ctx = this.cv.getContext('2d');
      this.markers = [];
      const ro = new ResizeObserver(() => this.resize()); ro.observe(this.host);
      this.bind();
    }
    setTracks(t) { this.tracks = t; this.draw(); }
    setRange(a, b) {
      if (!isFinite(a) || !isFinite(b)) return;
      if (a > b) [a, b] = [b, a];
      const span = Math.max(b - a, 0.004);
      const mid = (a + b) / 2;
      a = mid - span / 2; b = mid + span / 2;
      if (a < this.o.min) { b += this.o.min - a; a = this.o.min; }
      if (b > this.o.max) { a -= b - this.o.max; b = this.o.max; }
      this.range = [Math.max(a, this.o.min), Math.min(b, this.o.max)];
      this.draw();
      if (this.o.onRange) this.o.onRange(this.range);
    }
    resize() {
      const dpr = window.devicePixelRatio || 1;
      const w = this.host.clientWidth;
      const nT = Math.max(1, this.tracks.length);
      const H = this.o.height || (40 + nT * 52);
      this.W = w; this.H = H;
      this.cv.width = Math.round(w * dpr); this.cv.height = Math.round(H * dpr);
      this.cv.style.height = H + 'px';
      this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      this.draw();
    }
    x(w) { const [a, b] = this.range; return this.L + (w - a) / (b - a) * (this.W - this.L - this.R); }
    wAt(px) { const [a, b] = this.range; return a + (px - this.L) / (this.W - this.L - this.R) * (b - a); }
    bind() {
      const cv = this.cv;
      cv.addEventListener('wheel', e => {
        // let the page scroll unless the plot was clicked (active) or a modifier / pinch (ctrlKey) is used
        if (!this.active && !e.ctrlKey && !e.metaKey && !e.altKey && !e.shiftKey) { this.hint = true; this.draw(); return; }
        e.preventDefault();
        const r = cv.getBoundingClientRect(), px = e.clientX - r.left;
        const w0 = this.wAt(px), f = Math.exp((e.deltaY || 0) * 0.0018);
        const [a, b] = this.range;
        if (e.shiftKey || Math.abs(e.deltaX) > Math.abs(e.deltaY)) {
          const d = (e.deltaX || e.deltaY) * (b - a) / this.W; this.setRange(a + d, b + d); return;
        }
        this.setRange(w0 - (w0 - a) * f, w0 + (b - w0) * f);
      }, { passive: false });
      let drag = null;
      cv.addEventListener('pointerdown', e => { this.active = true; this.hint = false; this.host.classList.add('active'); drag = { x: e.clientX, r: this.range.slice() }; cv.setPointerCapture(e.pointerId); });
      cv.addEventListener('pointermove', e => {
        const r = cv.getBoundingClientRect();
        if (drag) {
          const d = (e.clientX - drag.x) / (this.W - this.L - this.R) * (drag.r[1] - drag.r[0]);
          this.setRange(drag.r[0] - d, drag.r[1] - d); this.tip.style.display = 'none'; return;
        }
        this.hover(e.clientX - r.left, e.clientY - r.top);
      });
      cv.addEventListener('pointerup', () => { drag = null; });
      cv.addEventListener('pointerleave', () => { this.tip.style.display = 'none'; this.hot = null; this.active = false; this.hint = false; this.host.classList.remove('active'); this.draw(); });
      cv.addEventListener('dblclick', e => {
        const r = cv.getBoundingClientRect(), w0 = this.wAt(e.clientX - r.left), [a, b] = this.range;
        this.setRange(w0 - (b - a) / 6, w0 + (b - a) / 6);
      });
    }
    trackGeom() {
      const top = 26, bottom = this.H - 34, n = Math.max(1, this.tracks.length);
      const th = (bottom - top) / n;
      return this.tracks.map((t, i) => ({ y0: top + i * th, y1: top + (i + 1) * th - 4 }));
    }
    visible(t) {
      const [a, b] = this.range, L = t.lines;
      // binary search bounds (lines sorted)
      let lo = 0, hi = L.length;
      while (lo < hi) { const m = (lo + hi) >> 1; if (L[m] < a) lo = m + 1; else hi = m; }
      let lo2 = lo, hi2 = L.length;
      while (lo2 < hi2) { const m = (lo2 + hi2) >> 1; if (L[m] <= b) lo2 = m + 1; else hi2 = m; }
      return [lo, lo2];
    }
    draw() {
      if (!this.W) return;
      const c = this.ctx, W = this.W, H = this.H;
      this.L = W < 560 ? 74 : 112; this.R = W < 560 ? 6 : 12;
      const bg = cssVar('--plot-bg'), grid = cssVar('--plot-grid'), ink = cssVar('--plot-ink');
      c.clearRect(0, 0, W, H); c.fillStyle = bg; c.fillRect(0, 0, W, H);
      const [a, b] = this.range;
      // rainbow band at top
      if (this.o.rainbow) {
        for (let px = this.L; px < W - this.R; px += 2) {
          const col = wlColor(this.wAt(px), 0.9);
          c.fillStyle = col || (this.wAt(px) < 3800 ? 'rgba(120,80,170,.25)' : 'rgba(150,40,40,.22)');
          c.fillRect(px, 6, 2, 10);
        }
      }
      // ticks
      const span = b - a, target = (W - this.L) / 110;
      const raw = span / target, p10 = Math.pow(10, Math.floor(Math.log10(raw)));
      const step = [1, 2, 5, 10].map(m => m * p10).find(s => s >= raw) || 10 * p10;
      const dec = Math.max(0, -Math.floor(Math.log10(step) + 1e-9));
      c.font = '11px ' + cssVar('--mono'); c.textAlign = 'center'; c.fillStyle = ink;
      const geo = this.trackGeom();
      for (let t = Math.ceil(a / step) * step; t <= b; t += step) {
        const x = this.x(t);
        c.strokeStyle = grid; c.lineWidth = 1; c.beginPath(); c.moveTo(x + .5, 20); c.lineTo(x + .5, H - 30); c.stroke();
        c.fillText(t.toFixed(dec), x, H - 16);
      }
      const med = this.o.medium ? this.o.medium() : 'air';
      c.textAlign = 'right'; c.fillText(`λ (Å, ${med === 'air' ? (a < 2000 ? (b <= 2000 ? 'vacuum' : 'vacuum < 2000 Å, air above') : 'air') : 'vacuum'})`, W - this.R, H - 3);
      if (this.hint) { c.textAlign = 'left'; c.fillStyle = cssVar('--accent'); c.fillText('click the plot, then scroll to zoom (or hold ⌘/Ctrl)', this.L + 4, H - 3); c.fillStyle = ink; }
      // tracks
      this.tracks.forEach((t, i) => {
        const g = geo[i];
        c.fillStyle = ink; c.textAlign = 'left'; c.font = '600 11.5px ' + cssVar('--sans');
        c.fillText(t.name, 8, (g.y0 + g.y1) / 2 + 4);
        if (t.sub) { c.font = '10px ' + cssVar('--sans'); c.fillStyle = cssVar('--muted'); c.fillText(t.sub, 8, (g.y0 + g.y1) / 2 + 17); }
        c.strokeStyle = grid; c.beginPath(); c.moveTo(this.L, g.y1 + 2.5); c.lineTo(W - this.R, g.y1 + 2.5); c.stroke();
        const [i0, i1] = this.visible(t);
        const n = i1 - i0;
        const absorb = t.mode === 'absorption';
        if (absorb) { c.fillStyle = cssVar('--panel'); c.fillRect(this.L, g.y0, W - this.L - this.R, g.y1 - g.y0); }
        const lw = n > 4000 ? 1 : 1.4;
        for (let k = i0; k < i1; k++) {
          const x = this.x(t.lines[k]);
          const s = t.str ? t.str[k] : 1;
          const hgt = Math.max(2, s * (g.y1 - g.y0));
          c.strokeStyle = t.color === 'rainbow' ? (stickColor(t.lines[k], .95) || ink) : t.color;
          c.globalAlpha = absorb ? Math.min(1, 0.25 + s) : (n > 3000 ? 0.7 : 1);
          c.lineWidth = lw;
          c.beginPath();
          if (absorb) { c.moveTo(x, g.y0); c.lineTo(x, g.y0 + hgt); }
          else { c.moveTo(x, g.y1); c.lineTo(x, g.y1 - hgt); }
          c.stroke();
        }
        c.globalAlpha = 1;
        // labels when zoomed in
        if (t.label && n > 0 && n < 60) {
          c.font = '10px ' + cssVar('--mono'); c.fillStyle = ink; c.textAlign = 'left';
          let lastEnd = -1e9;
          for (let k = i0; k < i1; k++) {
            const x = this.x(t.lines[k]);
            if (x < lastEnd) continue;
            let s = t.label(k); if (!s) continue;
            s = String(s); if (s.length > 18) s = s.slice(0, 17) + '…';
            c.fillText(s, x + 2, g.y0 + 9);
            lastEnd = x + c.measureText(s).width + 8;
          }
        }
      });
      // markers (e.g. search target)
      for (const m of this.markers) {
        const x = this.x(m.w); if (x < this.L || x > W - this.R) continue;
        c.strokeStyle = m.color || cssVar('--bad'); c.setLineDash([4, 3]); c.beginPath(); c.moveTo(x, 18); c.lineTo(x, H - 30); c.stroke(); c.setLineDash([]);
      }
      if (this.hot) {
        const x = this.x(this.hot.w); c.strokeStyle = cssVar('--accent'); c.lineWidth = 1; c.beginPath(); c.moveTo(x, 18); c.lineTo(x, H - 30); c.stroke();
      }
    }
    hover(px, py) {
      const geo = this.trackGeom();
      const ti = geo.findIndex(g => py >= g.y0 - 2 && py <= g.y1 + 3);
      if (ti < 0 || px < this.L) { this.tip.style.display = 'none'; return; }
      const t = this.tracks[ti]; const [i0, i1] = this.visible(t);
      let best = -1, bd = 7;
      for (let k = i0; k < i1; k++) { const d = Math.abs(this.x(t.lines[k]) - px); if (d < bd) { bd = d; best = k; } }
      if (best < 0) { this.tip.style.display = 'none'; this.hot = null; this.draw(); return; }
      this.hot = { w: t.lines[best] }; this.draw();
      this.tip.innerHTML = t.info ? t.info(best) : t.lines[best].toFixed(3);
      this.tip.style.display = 'block';
      const tw = this.tip.offsetWidth;
      this.tip.style.left = Math.min(this.W - tw - 6, Math.max(4, px + 12)) + 'px';
      this.tip.style.top = Math.max(4, py - 10) + 'px';
    }
  }

  // ======================================================================= LineChart (SVG)
  /** lineChart(host, {series:[{name,color,pts:[[x,y]],dash}], xlab, ylab, xlog, ylog, xr, yr, hlines, height}) */
  function lineChart(host, o) {
    host = typeof host === 'string' ? $(host) : host;
    const W = host.clientWidth || 700, H = o.height || 300, m = { l: 62, r: 16, t: 14, b: 42 };
    const allx = o.series.flatMap(s => s.pts.map(p => p[0])), ally = o.series.flatMap(s => s.pts.map(p => p[1]));
    const tx = v => o.xlog ? Math.log10(v) : v, ty = v => o.ylog ? Math.log10(v) : v;
    let [x0, x1] = o.xr || [Math.min(...allx), Math.max(...allx)];
    let [y0, y1] = o.yr || [Math.min(...ally), Math.max(...ally)];
    if (!o.yr) { const pad = (ty(y1) - ty(y0)) * 0.06 || 1; if (!o.ylog) { y0 -= pad; y1 += pad; } }
    const X = v => m.l + (tx(v) - tx(x0)) / (tx(x1) - tx(x0)) * (W - m.l - m.r);
    const Y = v => H - m.b - (ty(v) - ty(y0)) / (ty(y1) - ty(y0)) * (H - m.t - m.b);
    const ticks = (a, b, log) => {
      if (log) {
        const out = [], dec = Math.log10(b / a);
        const mults = dec <= 1.5 ? [1, 2, 3, 5] : dec <= 3.5 ? [1, 2, 5] : [1];
        for (let e = Math.floor(Math.log10(a)); e <= Math.ceil(Math.log10(b)); e++) for (const m of mults) { const v = m * Math.pow(10, e); if (v >= a * 0.999 && v <= b * 1.001) out.push(v); }
        return out;
      }
      const raw = (b - a) / 6, p = Math.pow(10, Math.floor(Math.log10(raw))), st = [1, 2, 2.5, 5, 10].map(k => k * p).find(s => s >= raw);
      const out = []; for (let v = Math.ceil(a / st) * st; v <= b + st * 1e-9; v += st) out.push(+v.toPrecision(12)); return out;
    };
    const fmt = v => (Math.abs(v) < 1e-2 && v !== 0) || Math.abs(v) >= 1e7 ? v.toExponential(0).replace('e+', 'e')
      : (Math.abs(v) >= 1e4 ? Math.round(v).toLocaleString('en-US').replace(/,/g, '\u2009') : String(+v.toPrecision(6)));
    let s = `<svg class="chart" width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" role="img" aria-label="${esc(o.aria || o.ylab || 'chart')}">`;
    s += '<g class="grid">';
    for (const v of ticks(x0, x1, o.xlog)) s += `<line x1="${X(v)}" x2="${X(v)}" y1="${m.t}" y2="${H - m.b}"/>`;
    for (const v of ticks(y0, y1, o.ylog)) s += `<line x1="${m.l}" x2="${W - m.r}" y1="${Y(v)}" y2="${Y(v)}"/>`;
    s += '</g>';
    for (const v of ticks(x0, x1, o.xlog)) s += `<text x="${X(v)}" y="${H - m.b + 16}" text-anchor="middle">${fmt(v)}</text>`;
    for (const v of ticks(y0, y1, o.ylog)) s += `<text x="${m.l - 6}" y="${Y(v) + 4}" text-anchor="end">${fmt(v)}</text>`;
    s += `<text x="${(m.l + W - m.r) / 2}" y="${H - 6}" text-anchor="middle">${esc(o.xlab || '')}</text>`;
    s += `<text transform="translate(14 ${(m.t + H - m.b) / 2}) rotate(-90)" text-anchor="middle">${esc(o.ylab || '')}</text>`;
    for (const hl of (o.hlines || [])) {
      const left = hl.pos === 'left', below = hl.below;
      s += `<line x1="${m.l}" x2="${W - m.r}" y1="${Y(hl.y)}" y2="${Y(hl.y)}" stroke="${hl.color}" stroke-dasharray="5 4"/><text x="${left ? m.l + 6 : W - m.r - 4}" y="${Y(hl.y) + (below ? 13 : -4)}" text-anchor="${left ? 'start' : 'end'}" style="fill:${hl.color}">${esc(hl.label || '')}</text>`;
    }
    for (const se of o.series) {
      const pts = se.pts.filter(p => isFinite(p[0]) && isFinite(p[1]) && (!o.ylog || p[1] > 0));
      if (se.points) { for (const p of pts) s += `<circle cx="${X(p[0])}" cy="${Y(p[1])}" r="${se.r || 3.5}" fill="${se.color}"><title>${esc(se.name)}: ${esc(p[2] || '')}</title></circle>`; continue; }
      s += `<path d="${pts.map((p, i) => (i ? 'L' : 'M') + X(p[0]).toFixed(1) + ' ' + Y(p[1]).toFixed(1)).join('')}" fill="none" stroke="${se.color}" stroke-width="${se.width || 2}"${se.dash ? ` stroke-dasharray="${se.dash}"` : ''}/>`;
    }
    s += '</svg>';
    host.innerHTML = s;
    if (o.legend !== false) {
      const lg = h('div', { class: 'legend' });
      for (const se of o.series) lg.append(h('span', { html: `<i style="background:${se.color}"></i>${esc(se.name)}` }));
      host.append(lg);
    }
  }

  root.UI = { $, $$, esc, h, store, fmtWl, fmtSrc, fmtNum, fmtSig, fmtU, grade, gradeBadge, wlColor, cssVar, download, DataTable, SpectrumView, lineChart };
})(typeof window !== 'undefined' ? window : globalThis);
