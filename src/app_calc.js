/* app_calc.js -- calculators: precision, air/vacuum, barycentric, frames, cosmology */
(function () {
  'use strict';
  const { $, $$, esc, h, store, fmtSig, lineChart, cssVar } = UI;
  const S = SpecCore;
  const APP = window.APP;
  (APP.initSections = APP.initSections || []).push(initCalcs);

  const SERIES_COLORS = ['#1f77b4', '#d62728', '#2ca02c', '#9467bd', '#ff7f0e', '#17becf', '#8c564b', '#e377c2'];
  const FORMULAS = [
    ['ciddor1996', 'Ciddor 1996 (IAG 1999)'], ['peckreeder1972', 'Peck & Reeder 1972 eq. 3 (NIST ASD)'], ['peckreeder1972eq2', 'Peck & Reeder 1972 eq. 2 (specutils, PyAstronomy)'], ['edlen1966', 'Edlén 1966'],
    ['birchdowns1994', 'Birch & Downs 1994'], ['morton2000', 'Morton 2000 (= Birch & Downs 1994; VALD)'], ['edlen1953', 'Edlén 1953 (IAU 1957)'], ['greisen2006', 'Greisen et al. 2006 eq. 65 (FITS WCS; 0 °C convention!)'],
  ];
  APP.FORMULAS = FORMULAS;

  function parseAngle(str, hours) {     // "10 20 30.5", "10:20:30", "155.2"
    str = String(str).trim();
    if (/^[-+]?\d+(\.\d*)?$/.test(str)) return parseFloat(str);
    const m = str.replace(/[hdms°'"]/g, ' ').trim().split(/[\s:]+/).map(Number);
    if (m.some(x => !isFinite(x))) return NaN;
    const sgn = /^\s*-/.test(str) ? -1 : 1;
    const v = Math.abs(m[0]) + (m[1] || 0) / 60 + (m[2] || 0) / 3600;
    return sgn * v * (hours ? 15 : 1);
  }
  APP.parseAngle = parseAngle;

  // ------------------------------------------------------------------ 1. precision calculator
  function initPrecision() {
    const wl = $('#pcWl'), dl = $('#pcDl'), dv = $('#pcDv'), out = $('#pcOut');
    function upd(src) {
      const w = +wl.value;
      if (src === 'v') dl.value = (+dv.value / S.C * w * 1000).toPrecision(6);
      else dv.value = (+dl.value / 1000 / w * S.C).toPrecision(6);
      const dA = +dl.value / 1000;
      const nu = S.C_M / (w * 1e-10);                       // Hz
      const dnu = nu * dA / w;
      const sig = 1e8 / w, dsig = sig * dA / w;
      out.innerHTML = `Δλ = ${(+dl.value).toPrecision(4)} mÅ at ${w} Å ⇔ Δv = ${(+dv.value * 1000).toPrecision(4)} m/s ⇔ Δν = ${(dnu / 1e6).toPrecision(4)} MHz ⇔ Δσ̃ = ${dsig.toExponential(3)} cm⁻¹ ⇔ fractional ${(dA / w).toExponential(3)}`;
    }
    wl.addEventListener('input', () => upd('l')); dl.addEventListener('input', () => upd('l')); dv.addEventListener('input', () => upd('v'));
    upd('l');
    // chart: velocity equivalent of 1 mA, and of 10 mA
    const pts1 = [], pts10 = [], pts01 = [];
    for (let w = 900; w <= 25000; w *= 1.03) { pts1.push([w, 0.001 / w * S.C_M]); pts10.push([w, 0.01 / w * S.C_M]); pts01.push([w, 0.0001 / w * S.C_M]); }
    const draw = () => lineChart('#pcChart', {
      series: [{ name: '10 mÅ (0.01 Å: 43 % of Handbook lines)', color: '#d62728', pts: pts10 }, { name: '1 mÅ (this wiki’s target)', color: '#1f77b4', pts: pts1 }, { name: '0.1 mÅ', color: '#2ca02c', pts: pts01 }],
      xlab: 'wavelength (Å)', ylab: 'velocity equivalent (m/s)', xlog: true, ylog: true, height: 320, yr: [0.5, 5000],
      hlines: [{ y: 636, color: '#d9730d', label: 'solar gravitational redshift 636 m/s' }, { y: 465, color: '#8e5cc4', label: 'Earth rotation 465 m/s', pos: 'left', below: true }, { y: 1, color: '#8c564b', label: '1 m/s RV (EPRV goal: 0.1 m/s)', pos: 'left' }],
    });
    draw(); window.addEventListener('resize', draw); window.addEventListener('themechange', draw);
  }

  // ------------------------------------------------------------------ 2. air / vacuum converter
  function initAirVac() {
    const wl = $('#avWl'), dir = $('#avDir'), T = $('#avT'), P = $('#avP'), RH = $('#avRH'), CO2 = $('#avCO2'), out = $('#avOut');
    function upd() {
      const w = +wl.value; if (!(w > 0)) return;
      const toAir = dir.value === 'v2a';
      const rows = FORMULAS.map(([k, name]) => {
        const r = toAir ? S.vacToAir(w, k) : S.airToVac(w, k);
        return [name, r, k];
      });
      const ref = rows[0][1];
      let html = `<table class="t"><thead><tr><th>Formula (standard air)</th><th class="num">${toAir ? 'λ<sub>air</sub>' : 'λ<sub>vac</sub>'} (Å)</th><th class="num">Δ vs Ciddor (mÅ)</th><th class="num">Δv (m/s)</th><th class="num">n − 1</th></tr></thead><tbody>`;
      for (const [name, r, k] of rows) {
        const vac = toAir ? w : r;
        html += `<tr><td>${esc(name)}</td><td class="num">${r.toFixed(5)}</td><td class="num">${((r - ref) * 1000).toFixed(4)}</td><td class="num">${((r - ref) / w * S.C_M).toFixed(2)}</td><td class="num">${S.AIR[k](1e4 / vac).toExponential(6)}</td></tr>`;
      }
      if (w < 2000) html += `<tr><td colspan="5" class="small muted">Below 2000 Å wavelengths are quoted in vacuum (IAU/NIST convention); no air formula is validated below 1850 Å, and the formulas diverge by several mÅ there.</td></tr>`;
      // Ciddor full model at ambient conditions (air->vac: iterate with the ambient index)
      let vacW = w;
      if (!toAir) {
        for (let i = 0; i < 30; i++) { const nv = w * (1 + S.ciddorFull(vacW, +T.value, +P.value, +RH.value / 100, +CO2.value)); if (Math.abs(nv - vacW) < 1e-12) { vacW = nv; break; } vacW = nv; }
      }
      const nAmb = S.ciddorFull(vacW, +T.value, +P.value, +RH.value / 100, +CO2.value);
      const nStd = S.ciddorFull(vacW, 15, 101325, 0, 450);
      const ambAir = vacW / (1 + nAmb);
      html += `<tr><td><b>Ciddor 1996 full model</b> at T = ${T.value} °C, P = ${P.value} Pa, RH = ${RH.value} %, CO₂ = ${CO2.value} ppm</td><td class="num">${(toAir ? ambAir : vacW).toFixed(5)}</td><td class="num">${(((toAir ? ambAir : vacW) - ref) * 1000).toFixed(3)}</td><td class="num">${(((toAir ? ambAir : vacW) - ref) / w * S.C_M).toFixed(1)}</td><td class="num">${nAmb.toExponential(6)}</td></tr>`;
      html += `<tr><td class="small muted" colspan="5">Full-model check at standard conditions: n−1 = ${nStd.toExponential(7)} vs standard-air formula ${S.AIR.ciddor1996(1e4 / vacW).toExponential(7)}. ` +
        `Vacuum−air difference: ${((vacW - vacW / (1 + nStd))).toFixed(4)} Å (${((vacW - vacW / (1 + nStd)) / vacW * S.C).toFixed(2)} km/s).</td></tr>`;
      html += '</tbody></table>';
      out.innerHTML = html;
    }
    [wl, dir, T, P, RH, CO2].forEach(e => e.addEventListener('input', upd));
    $$('[data-avp]').forEach(b => b.addEventListener('click', () => { const [t, p, rh] = b.dataset.avp.split(',').map(Number); T.value = t; P.value = p; RH.value = rh; upd(); }));
    upd();
    // comparison chart (differences vs Ciddor in mA, vac->air)
    const draw = () => {
      const series = FORMULAS.slice(1).filter(([k]) => k !== 'greisen2006').map(([k, name], i) => {
        const pts = [];
        for (let w = 2000; w <= 25000; w *= 1.02) pts.push([w, (S.vacToAir(w, k) - S.vacToAir(w, 'ciddor1996')) * 1000]);
        return { name: name + ' − Ciddor', color: SERIES_COLORS[i + 1], pts };
      });
      lineChart('#avChart', { series, xlab: 'vacuum wavelength (Å)', ylab: 'λ_air difference (mÅ)', xlog: true, height: 320, yr: [-1.2, 1.2], hlines: [{ y: 1, color: '#999', label: '+1 mÅ' }, { y: -1, color: '#999', label: '−1 mÅ', below: true }] });
      const pts = [];
      for (let w = 2000; w <= 25000; w *= 1.02) pts.push([w, w - S.vacToAir(w, 'ciddor1996')]);
      lineChart('#avChart2', { series: [{ name: 'λ_vac − λ_air (Ciddor 1996)', color: '#1f77b4', pts }], xlab: 'vacuum wavelength (Å)', ylab: 'λ_vac − λ_air (Å)', xlog: true, height: 240 });
    };
    draw(); window.addEventListener('resize', draw); window.addEventListener('themechange', draw);
  }

  // ------------------------------------------------------------------ 3. barycentric correction
  const SITES = [
    ['Mauna Kea — Keck', -155.47833, 19.82833, 4160], ['Mauna Kea — Subaru', -155.47610, 19.82552, 4139], ['Cerro Paranal — VLT (ESPRESSO)', -70.40417, -24.62722, 2635],
    ['La Silla — ESO 3.6 m (HARPS)', -70.73169, -29.26097, 2400], ['Roque de los Muchachos — TNG (HARPS-N)', -17.88917, 28.75400, 2387],
    ['Kitt Peak — WIYN (NEID)', -111.60028, 31.95833, 2091], ['Kitt Peak — Mayall (DESI)', -111.59972, 31.96400, 2120], ['Apache Point — SDSS 2.5 m', -105.82028, 32.78028, 2788],
    ['Las Campanas — Magellan', -70.69264, -29.01458, 2380], ['Cerro Pachón — Gemini S / Rubin', -70.73667, -30.24074, 2722], ['McDonald — HET', -104.01472, 30.68139, 2026],
    ['Calar Alto (CARMENES)', -2.54625, 37.22361, 2168], ['Haute-Provence (SOPHIE)', 5.71333, 43.93083, 650], ['Xinglong — LAMOST', 117.57500, 40.39583, 900],
    ['Lijiang — Gaomeigu 2.4 m', 100.03017, 26.69500, 3200], ['Lenghu — Saishiteng (MUST site)', 93.89580, 38.60680, 4200], ['Siding Spring — AAT', 149.06611, -31.27333, 1164],
    ['Sutherland — SALT', 20.81069, -32.37592, 1798], ['Geocentre (no rotation; geoid Earth potential)', 0, 90, -6356752.3],
  ];
  APP.SITES = SITES;
  function initBary() {
    const ra = $('#bcRa'), dec = $('#bcDec'), site = $('#bcSite'), lon = $('#bcLon'), lat = $('#bcLat'), hh = $('#bcH'), t = $('#bcT'), dut = $('#bcDut1'), out = $('#bcOut');
    SITES.forEach(([n], i) => site.append(h('option', { value: i }, n)));
    site.value = store.get('site', 0);
    const fill = () => { const s = SITES[+site.value]; lon.value = s[1]; lat.value = s[2]; hh.value = s[3]; store.set('site', +site.value); };
    site.addEventListener('change', () => { fill(); upd(); });
    fill();
    if (!t.value) t.value = new Date().toISOString().slice(0, 19);
    $('#bcNow').addEventListener('click', () => { t.value = new Date().toISOString().slice(0, 19); upd(); });
    function upd() {
      try {
        const a = parseAngle(ra.value, /[:\sh]/.test(ra.value.trim())), d = parseAngle(dec.value, false);
        const r = S.baryCorrection({ raDeg: a, decDeg: d, lonDeg: +lon.value, latDeg: +lat.value, hM: +hh.value, isoUtc: t.value, dut1: +dut.value || 0 });
        const g = S.icrsToGal(a, d);
        const ms = v => (v * 1000).toFixed(3);
        out.innerHTML = `<table class="t"><tbody>
          <tr><td><b>Barycentric correction v<sub>B</sub> = c z<sub>B</sub></b> (Wright &amp; Eastman 2014 eq. 28)</td><td class="num"><b>${ms(r.vcorr_kms)}</b> m/s</td></tr>
          <tr><td>· Earth orbital velocity (geocentre wrt SSB, DE440s) projected on target</td><td class="num">${ms(r.v_orbit_kms)} m/s</td></tr>
          <tr><td>· Earth rotation (observer wrt geocentre) projected on target</td><td class="num">${ms(r.v_rot_kms)} m/s</td></tr>
          <tr><td>· relativistic terms: γ − 1 (time dilation) and Σ GM/(rc) (gravitational redshift at observer)</td><td class="num">${ms(r.gamma_term_kms)} / ${ms(r.grav_kms)} m/s</td></tr>
          <tr><td>Heliocentric correction (legacy; Sun instead of SSB, non-relativistic)</td><td class="num">${ms(r.v_helio_kms)} m/s</td></tr>
          <tr><td>Heliocentric − barycentric</td><td class="num">${((r.v_helio_kms - r.vcorr_kms) * 1000).toFixed(3)} m/s</td></tr>
          <tr><td>Rømer delay (observer → SSB, along target) / BJD<sub>TDB</sub></td><td class="num">${r.roemer_s.toFixed(3)} s / ${r.bjd_tdb.toFixed(7)}</td></tr>
          <tr><td>Time scales: TAI−UTC / JD<sub>UTC</sub> / JD<sub>TDB</sub></td><td class="num">${r.ts.taiMinusUtc} s / ${r.ts.utc.toFixed(6)} / ${r.ts.tdb.toFixed(6)}</td></tr>
          <tr><td>Target Galactic coordinates</td><td class="num">l = ${g.l.toFixed(4)}°, b = ${g.b.toFixed(4)}°</td></tr>
          <tr><td>Shift of a 5000 Å line due to v<sub>B</sub></td><td class="num">${(5000 * r.vcorr_kms / S.C).toFixed(4)} Å</td></tr>
          </tbody></table>
          <p class="small muted">Apply as 1 + z<sub>true</sub> = (1 + z<sub>meas</sub>)(1 + z<sub>B</sub>), i.e. v<sub>true</sub> ≈ v<sub>meas</sub> + v<sub>B</sub> + v<sub>meas</sub>v<sub>B</sub>/c. Validated against astropy 7.2 <code>radial_velocity_correction</code> with DE440s: rms 0.26 mm/s, max 0.68 mm/s over 300 random cases (1990–2030).</p>`;
        APP.lastTarget = { ra: a, dec: d, l: g.l, b: g.b, vb: r.vcorr_kms, zb: r.zb, iso: t.value };
        if (APP.updForward) APP.updForward();
        drawYear(a, d);
        if (APP.updFrames) APP.updFrames();
      } catch (e) { out.innerHTML = `<p class="note warn">${esc(e.message)}</p>`; }
    }
    function drawYear(a, d) {
      const y0 = +t.value.slice(0, 4);
      const pts = [], ptsR = [];
      for (let k = 0; k <= 365; k += 1) {
        const dt = new Date(Date.UTC(y0, 0, 1) + k * 86400e3).toISOString().slice(0, 19);
        try {
          const r = S.baryCorrection({ raDeg: a, decDeg: d, lonDeg: +lon.value, latDeg: +lat.value, hM: +hh.value, isoUtc: dt, dut1: 0 });
          pts.push([k, r.vcorr_kms]); ptsR.push([k, r.v_rot_kms * 50]);
        } catch (e) { /* outside span */ }
      }
      lineChart('#bcChart', { series: [{ name: `v_B at 00:00 UTC each day of ${y0} (km/s)`, color: '#1f77b4', pts }, { name: 'rotation term × 50 (km/s)', color: '#ff7f0e', pts: ptsR, width: 1.2 }], xlab: `day of ${y0}`, ylab: 'km/s', height: 260 });
    }
    [ra, dec, lon, lat, hh, t, dut].forEach(e => e.addEventListener('change', upd));
    upd();
    window.addEventListener('themechange', upd);
  }

  // ------------------------------------------------------------------ 4. kinematic frames (LSR, GSR, LG, CMB)
  function initFrames() {
    const F = APP.J('data-frames');
    APP.frames = F;
    const vin = $('#frV'), out = $('#frOut');
    function upd() {
      const T = APP.lastTarget; if (!T) return;
      const vB = vin.value === '' ? 0 : +vin.value;
      let html = `<table class="t"><thead><tr><th>Frame</th><th>Definition (solar motion)</th><th class="num">Correction Δv (km/s)</th><th class="num">v in frame (km/s)</th><th class="num">z in frame</th></tr></thead><tbody>`;
      const zB = vB / S.C;
      for (const f of F.frames) {
        let dv;
        if (f.uvw) dv = S.uvwProjection(f.uvw, T.l, T.b);
        else if (f.apex) dv = S.apexProjection(f.v, f.apex[0], f.apex[1], T.l, T.b);
        else if (f.formula === 'yts77') { const L = T.l * Math.PI / 180, B = T.b * Math.PI / 180; dv = -79 * Math.cos(L) * Math.cos(B) + 296 * Math.sin(L) * Math.cos(B) - 36 * Math.sin(B); }
        else if (f.formula === 'iau76') { dv = 300 * Math.sin(T.l * Math.PI / 180) * Math.cos(T.b * Math.PI / 180); }
        f._dv = dv;
        // exact special-relativistic boost for a direction given in the barycentric (catalogue) frame:
        // 1+z_frame = (1+z_bary) / [gamma (1 - beta cos(theta))]   (to first order (1+z_bary)(1 + v cos(theta)/c))
        const vmag = f.uvw ? Math.hypot(...f.uvw) : (f.v || Math.abs(dv));
        const gam = 1 / Math.sqrt(1 - (vmag / S.C) ** 2);
        const z = (1 + zB) / (gam * (1 - dv / S.C)) - 1;
        html += `<tr><td><b>${esc(f.name)}</b></td><td class="small">${f.desc}</td><td class="num">${dv >= 0 ? '+' : ''}${dv.toFixed(2)}</td><td class="num">${(vB + dv).toFixed(2)}</td><td class="num">${z.toExponential(6)}</td></tr>`;
      }
      html += '</tbody></table>';
      out.innerHTML = html + `<p class="small muted">Target l = ${T.l.toFixed(3)}°, b = ${T.b.toFixed(3)}° (from the barycentric calculator above). Input v is the barycentric (≈ heliocentric) radial velocity; each frame adds the projection of the Sun's velocity relative to that frame onto the line of sight. Redshifts are composed exactly, 1 + z<sub>frame</sub> = (1 + z<sub>bary</sub>) / [γ(1 − β cos θ)] with θ the catalogue (barycentric) angle. Relative to the first-order (1 + z<sub>bary</sub>)(1 + β cos θ) the second-order terms β²(cos²θ − ½) reach ±7.6×10⁻⁷ (±0.23 km/s) for the CMB dipole; using γ(1 + β cos θ) with the catalogue angle would err by β² sin²θ ≤ 1.5×10⁻⁶ (0.46 km/s).</p>`;
    }
    vin.addEventListener('input', upd);
    APP.updFrames = upd;
    upd();
  }

  // ------------------------------------------------------------------ 5. cosmology: redshift drift & composition
  function initCosmo() {
    const H0 = $('#cdH0'), Om = $('#cdOm'), OL = $('#cdOL'), z = $('#cdZ'), yr = $('#cdYr'), out = $('#cdOut');
    const E = (zz, om, ol) => Math.sqrt(om * (1 + zz) ** 3 + (1 - om - ol) * (1 + zz) ** 2 + ol);
    function upd() {
      const h0 = +H0.value, om = +Om.value, ol = +OL.value, zz = +z.value, t = +yr.value;
      const H0s = h0 / 3.0856775814913673e19;             // s^-1
      const sec = t * 365.25 * 86400;
      const zdot = H0s * ((1 + zz) - E(zz, om, ol));        // per second
      const dz = zdot * sec, dv = S.C_M * dz / (1 + zz);
      const lam = 1215.67 * (1 + zz);
      out.innerHTML = `ż = ${(zdot * 365.25 * 86400).toExponential(3)} yr⁻¹ · over ${t} yr: Δz = ${dz.toExponential(3)}, Δv = ${(dv * 100).toFixed(3)} cm/s, Δλ(Lyα at ${lam.toFixed(1)} Å) = ${(lam * dv / S.C_M).toExponential(3)} Å`;
      const pts = [];
      for (let q = 0.01; q <= 5; q += 0.02) pts.push([q, S.C_M * H0s * ((1 + q) - E(q, om, ol)) / (1 + q) * t * 365.25 * 86400 * 100]);
      lineChart('#cdChart', { series: [{ name: `Δv over ${t} yr (cm/s), H0=${h0}, Ωm=${om}, ΩΛ=${ol}`, color: '#1f77b4', pts }], xlab: 'redshift z', ylab: 'Δv (cm/s)', height: 250, hlines: [{ y: 0, color: '#999', label: '' }] });
    }
    [H0, Om, OL, z, yr].forEach(e => e.addEventListener('input', upd));
    upd();
    window.addEventListener('themechange', upd);
    // composition calculator
    const zo = $('#czObs'), vb = $('#czVb'), vp = $('#czCmb'), oc = $('#czOut');
    function upd2() {
      const zObs = +zo.value, zB = +vb.value / S.C, dCMB = +vp.value / S.C;
      const zBary = (1 + zObs) * (1 + zB) - 1;
      const zC = (1 + zBary) * (1 + dCMB) - 1;
      const zAdd = zObs + zB + dCMB;
      oc.innerHTML = `z<sub>bary</sub> = ${zBary.toFixed(8)} · z<sub>CMB</sub> = ${zC.toFixed(8)} · naive additive z = ${zAdd.toFixed(8)} → error ${((zAdd - zC)).toExponential(2)} (${((zAdd - zC) * S.C).toFixed(3)} km/s, ${(((zAdd - zC) / (1 + zC)) * S.C * 1000).toFixed(1)} m/s in the rest frame)`;
    }
    [zo, vb, vp].forEach(e => e.addEventListener('input', upd2));
    upd2();
  }

  // ------------------------------------------------------------------ 6. forward model (eq. 7.1 end to end)
  function frameDv(f, T) {
    if (!f) return { dv: 0, vmag: 0 };
    let dv;
    if (f.uvw) dv = S.uvwProjection(f.uvw, T.l, T.b);
    else if (f.apex) dv = S.apexProjection(f.v, f.apex[0], f.apex[1], T.l, T.b);
    else if (f.formula === 'yts77') { const L = T.l * Math.PI / 180, B = T.b * Math.PI / 180; dv = -79 * Math.cos(L) * Math.cos(B) + 296 * Math.sin(L) * Math.cos(B) - 36 * Math.sin(B); }
    else if (f.formula === 'iau76') dv = 300 * Math.sin(T.l * Math.PI / 180) * Math.cos(T.b * Math.PI / 180);
    return { dv, vmag: f.uvw ? Math.hypot(...f.uvw) : (f.v || Math.abs(dv)) };
  }
  function initForward() {
    const lam = $('#fmLam'), v = $('#fmV'), z = $('#fmZ'), fr = $('#fmFrame'), gr = $('#fmGrav'), med = $('#fmMed'), out = $('#fmOut');
    if (!lam) return;
    fr.append(h('option', { value: '' }, 'Solar-System barycentre (no frame boost)'));
    for (const f of APP.frames.frames) fr.append(h('option', { value: f.key }, f.name));
    function upd() {
      const T = APP.lastTarget; if (!T) return;
      const L0 = parseFloat(lam.value); if (!(L0 > 0)) { out.innerHTML = '<p class="note warn">enter a laboratory vacuum wavelength</p>'; return; }
      const zIn = parseFloat(z.value) || 0, vIn = parseFloat(v.value) || 0;
      const zF = zIn !== 0 ? zIn : vIn / S.C;
      const f = APP.frames.frames.find(x => x.key === fr.value);
      const { dv, vmag } = frameDv(f, T);
      const gam = 1 / Math.sqrt(1 - (vmag / S.C) ** 2);
      const zSrc = f ? (1 + zF) * gam * (1 - dv / S.C) - 1 : zF;          // to barycentric (inverse of the §7.7 boost)
      const zG = (parseFloat(gr.value) || 0) / S.C_M;
      const L1 = L0 * (1 + zSrc), L2 = L1 * (1 + zG), L3 = L2 / (1 + T.zb);
      let n = 0, mlab = 'vacuum';
      if (med.value === 'std') { n = S.AIR.peckreeder1972(1e4 / L3); mlab = 'standard air (P&R 1972)'; }
      if (med.value === 'amb') { n = S.ciddorFull(L3, +$('#avT').value, +$('#avP').value, +$('#avRH').value / 100, +$('#avCO2').value); mlab = `air at ${$('#avT').value} °C, ${$('#avP').value} Pa, ${$('#avRH').value} % RH`; }
      const L4 = L3 / (1 + n);
      const row = (name, fac, Lp, Ln) => `<tr><td>${name}</td><td class="num">${fac}</td><td class="num">${Ln.toFixed(5)}</td><td class="num">${((Ln - Lp) * 1000).toFixed(2)}</td><td class="num">${((Ln - Lp) / Lp * S.C).toFixed(4)}</td></tr>`;
      out.innerHTML = `<table class="t"><thead><tr><th>Step</th><th class="num">factor</th><th class="num">λ (Å)</th><th class="num">Δλ (mÅ)</th><th class="num">≈ km/s</th></tr></thead><tbody>
        <tr><td>laboratory, vacuum</td><td class="num"></td><td class="num">${L0.toFixed(5)}</td><td></td><td></td></tr>
        ${row(`source motion${f ? ` (z = ${zF.toExponential(5)} in ${esc(f.name)} → barycentric z = ${zSrc.toExponential(6)})` : ` (z = ${zSrc.toExponential(6)}, barycentric)`}`, `×${(1 + zSrc).toFixed(9)}`, L0, L1)}
        ${row('source gravitational redshift', `×${(1 + zG).toFixed(12)}`, L1, L2)}
        ${row(`observer motion: ÷(1 + z<sub>B</sub>), v<sub>B</sub> = ${(T.zb * S.C).toFixed(4)} km/s at ${esc(T.iso)} UTC`, `÷${(1 + T.zb).toFixed(9)}`, L2, L3)}
        ${row(`medium: ${esc(mlab)}`, n ? `÷${(1 + n).toFixed(9)}` : '×1', L3, L4)}
        <tr><td><b>expected observed wavelength</b></td><td></td><td class="num"><b>${L4.toFixed(5)}</b></td><td class="num">${((L4 - L0) * 1000).toFixed(1)}</td><td></td></tr>
        </tbody></table><p class="small muted">Frame boost: 1 + z<sub>bary</sub> = (1 + z<sub>frame</sub>) γ (1 − β cos θ) with θ the catalogue angle to the frame apex (inverse of §7.7). Not included: the source's convective blueshift and pressure shifts, instrumental zero points.</p>`;
    }
    [lam, v, z, fr, gr, med].forEach(e => e.addEventListener('input', upd));
    ['#avT', '#avP', '#avRH', '#avCO2'].forEach(id => $(id).addEventListener('input', upd));
    APP.updForward = upd; upd();
  }

  function initCalcs() {
    APP.ephem = APP.J('data-ephem'); S.setEphemeris(APP.ephem);
    initPrecision(); initAirVac(); initFrames(); initBary(); initForward(); initCosmo();
  }
})();
