/* spectro_core.js -- computational core of the Spectral Line Wiki.
 * Runs in the browser (window.SpecCore) and in node (module.exports) so that it can be
 * regression-tested against astropy / ERFA.
 *
 * Contents
 *   Time scales      UTC -> TAI -> TT -> TDB, Julian dates, ERA
 *   Ephemeris        Chebyshev evaluation of the DE440s fit (Earth, Sun; ICRS, SSB)
 *   Earth rotation   WGS84 geodetic -> ITRS -> GCRS via CIO/ERA with IAU 2006/2000A leading terms
 *   RV correction    Wright & Eastman (2014) eq. 28 as implemented in astropy
 *   Frames           ICRS <-> Galactic; LSRK, LSRD, GSR, LG, CMB projections
 *   Air/vacuum       standard-air dispersion formulas (Edlen 1953/1966, Peck & Reeder 1972,
 *                    Birch & Downs 1994, Ciddor 1996 incl. T/P/RH/CO2, Morton 2000, Greisen 2006,
 *                    VALD/Piskunov inverse) with exact iterative inverses
 *   Doppler          optical / radio / relativistic conventions, redshift composition
 */
(function (root) {
  'use strict';
  const C = 299792.458;               // km/s (exact)
  const C_M = 299792458.0;            // m/s
  const DAY = 86400.0;
  const AS2R = Math.PI / 180 / 3600;
  const D2R = Math.PI / 180;
  const GM_SUN = 1.32712440041e20;    // m^3 s^-2 (IAU 2015 nominal, TDB-compatible; astropy GM_sun)
  const GM_EARTH = 3.986004e14;       // astropy consts.GM_earth (IAU 2015 nominal)
  const GM_JUP = 1.2668653e17;
  const GM_MOON = 6.6743e-11 * 7.34767309e22;
  const AU_KM = 149597870.7;

  // ------------------------------------------------------------------ time
  // TAI-UTC (s) from IERS Bulletin C; [MJD of start, value]
  const LEAP = [[41317, 10], [41499, 11], [41683, 12], [42048, 13], [42413, 14], [42778, 15],
    [43144, 16], [43509, 17], [43874, 18], [44239, 19], [44786, 20], [45151, 21], [45516, 22],
    [46247, 23], [47161, 24], [47892, 25], [48257, 26], [48804, 27], [49169, 28], [49534, 29],
    [50083, 30], [50630, 31], [51179, 32], [53736, 33], [54832, 34], [56109, 35], [57204, 36],
    [57754, 37]];
  function taiMinusUtc(mjdUtc) {
    let v = 10;
    for (const [m, s] of LEAP) if (mjdUtc >= m) v = s;
    return v;
  }
  // Calendar date (UTC, proleptic Gregorian) -> JD
  function jdFromDate(Y, M, D, h = 0, mi = 0, s = 0) {
    if (M <= 2) { Y -= 1; M += 12; }
    const A = Math.floor(Y / 100), B = 2 - A + Math.floor(A / 4);
    const jd0 = Math.floor(365.25 * (Y + 4716)) + Math.floor(30.6001 * (M + 1)) + D + B - 1524.5;
    return [jd0, (h + mi / 60 + s / 3600) / 24];     // two-part JD for precision
  }
  function parseISO(str) {
    const m = String(str).trim().match(/^(\d{4})-(\d{2})-(\d{2})(?:[T ](\d{2}):(\d{2})(?::(\d{2}(?:\.\d*)?))?)?Z?$/);
    if (!m) throw new Error('bad ISO time: ' + str);
    return jdFromDate(+m[1], +m[2], +m[3], +(m[4] || 0), +(m[5] || 0), +(m[6] || 0));
  }
  // TDB - TT (s), Fairhead & Bretagnon (1990) leading terms (USNO Circ. 179 eq. 2.6), ~10 us
  function tdbMinusTt(jdTT) {
    const T = (jdTT - 2451545.0) / 36525;
    return 0.001657 * Math.sin(628.3076 * T + 6.2401) + 0.000022 * Math.sin(575.3385 * T + 4.2970)
      + 0.000014 * Math.sin(1256.6152 * T + 6.1969) + 0.000005 * Math.sin(606.9777 * T + 4.0212)
      + 0.000005 * Math.sin(52.9691 * T + 0.4444) + 0.000002 * Math.sin(21.3299 * T + 5.5431)
      + 0.000010 * T * Math.sin(628.3076 * T + 4.2490);
  }
  /** UTC two-part JD -> {utc, tt, tdb, ut1} (single-number JDs plus two-part tt) */
  function timeScales(jdUtc2, dut1 = 0) {
    const [a, b] = jdUtc2;
    const mjd = a + b - 2400000.5;
    const dat = taiMinusUtc(mjd);
    const ttB = b + (dat + 32.184) / DAY;
    const tt = a + ttB;
    const tdb = tt + tdbMinusTt(tt) / DAY;
    return { utc: a + b, tt, tdb, tdb2: [a, ttB + tdbMinusTt(tt) / DAY], ut1: a + b + dut1 / DAY, taiMinusUtc: dat };
  }
  function era(jdUt1) {            // Earth rotation angle (rad), IERS 2010 eq. 5.15
    const Du = jdUt1 - 2451545.0;
    const f = jdUt1 % 1.0;          // (JD - 2451545.0) mod 1 == JD mod 1, as in ERFA eraEra00
    let t = 2 * Math.PI * ((f + 0.7790572732640 + 0.00273781191135448 * Du) % 1.0);
    if (t < 0) t += 2 * Math.PI;
    return t;
  }

  // ------------------------------------------------------------------ ephemeris
  let EPH = null;
  function setEphemeris(e) { EPH = e; }
  function chebEval(block, jdTdb) {
    const { seg_days: seg, deg, scale, coef } = block;
    const n = deg + 1;
    let i = Math.floor((jdTdb - EPH.jd0_tdb) / seg);
    if (i < 0 || i >= coef.length) throw new Error('date outside ephemeris span (1990-2066)');
    const x = 2 * (jdTdb - (EPH.jd0_tdb + i * seg)) / seg - 1;
    const c = coef[i];
    const out = [0, 0, 0];
    for (let k = 0; k < 3; k++) {     // Clenshaw
      let b1 = 0, b2 = 0;
      for (let j = n - 1; j >= 1; j--) { const t = 2 * x * b1 - b2 + c[k * n + j]; b2 = b1; b1 = t; }
      out[k] = (x * b1 - b2 + c[k * n]) * scale;
    }
    return out;
  }
  const earthPos = jd => chebEval(EPH.earth_pos, jd);   // km
  const earthVel = jd => chebEval(EPH.earth_vel, jd);   // km/s
  const sunPos = jd => chebEval(EPH.sun_pos, jd);
  const sunVel = jd => chebEval(EPH.sun_vel, jd);

  // ------------------------------------------------------------------ vectors
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const scl = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
  const norm = a => Math.hypot(a[0], a[1], a[2]);
  const mxv = (M, v) => [dot(M[0], v), dot(M[1], v), dot(M[2], v)];
  const mxm = (A, B) => A.map(r => [0, 1, 2].map(j => r[0] * B[0][j] + r[1] * B[1][j] + r[2] * B[2][j]));
  const tr = A => [0, 1, 2].map(i => [A[0][i], A[1][i], A[2][i]]);
  function R1(t) { const c = Math.cos(t), s = Math.sin(t); return [[1, 0, 0], [0, c, s], [0, -s, c]]; }
  function R2(t) { const c = Math.cos(t), s = Math.sin(t); return [[c, 0, -s], [0, 1, 0], [s, 0, c]]; }
  function R3(t) { const c = Math.cos(t), s = Math.sin(t); return [[c, s, 0], [-s, c, 0], [0, 0, 1]]; }
  function unitRaDec(raDeg, decDeg) {
    const a = raDeg * D2R, d = decDeg * D2R;
    return [Math.cos(d) * Math.cos(a), Math.cos(d) * Math.sin(a), Math.sin(d)];
  }

  // ------------------------------------------------------------------ Earth orientation
  // CIP X, Y (rad): IAU 2006 precession polynomial + largest IAU 2000A nutation terms.
  // Truncation error ~0.1 arcsec => < 0.3 mm/s on the rotational velocity.
  function cipXY(jdTT) {
    const T = (jdTT - 2451545.0) / 36525;
    const Om = (450160.398036 - 6962890.5431 * T) * AS2R;                // Moon node
    const F = (335779.526232 + 1739527262.8478 * T) * AS2R;
    const D = (1072260.703692 + 1602961601.2090 * T) * AS2R;
    const lp = (1287104.79305 + 129596581.0481 * T) * AS2R;
    const l = (485868.249036 + 1717915923.2178 * T) * AS2R;
    // nutation in longitude/obliquity (arcsec), leading IAU 2000A terms
    const dpsi = -17.2064161 * Math.sin(Om) - 1.3170907 * Math.sin(2 * F - 2 * D + 2 * Om)
      - 0.2276413 * Math.sin(2 * F + 2 * Om) + 0.2074554 * Math.sin(2 * Om) + 0.1475877 * Math.sin(lp)
      + 0.0711159 * Math.sin(l) - 0.0516821 * Math.sin(lp + 2 * F - 2 * D + 2 * Om)
      - 0.0387298 * Math.sin(2 * F + Om) - 0.0301461 * Math.sin(l + 2 * F + 2 * Om);
    const deps = 9.2052331 * Math.cos(Om) + 0.5730336 * Math.cos(2 * F - 2 * D + 2 * Om)
      + 0.0978459 * Math.cos(2 * F + 2 * Om) - 0.0897492 * Math.cos(2 * Om)
      + 0.0224386 * Math.cos(lp + 2 * F - 2 * D + 2 * Om) + 0.0200728 * Math.cos(2 * F + Om)
      + 0.0129025 * Math.cos(l + 2 * F + 2 * Om);
    const eps0 = 84381.406 * AS2R;
    const X = (-0.016617 + 2004.191898 * T - 0.4297829 * T * T - 0.19861834 * T ** 3) * AS2R
      + dpsi * Math.sin(eps0) * AS2R;
    const Y = (-0.006951 - 0.025896 * T - 22.4072747 * T * T + 0.00190059 * T ** 3) * AS2R
      + deps * AS2R;
    return [X, Y];
  }
  // GCRS <- CIRS matrix Q (IERS Conventions 2010 eq. 5.10), s ~ -XY/2
  function cioQ(jdTT) {
    const [X, Y] = cipXY(jdTT);
    const s = -X * Y / 2 + (-0.0000094) * AS2R;   // leading s+XY/2 constant term is 94 uas
    const E = Math.atan2(Y, X), d = Math.atan(Math.sqrt((X * X + Y * Y) / (1 - X * X - Y * Y)));
    // Q = R3(-E) R2(-d) R3(E) R3(s)
    return mxm(mxm(mxm(R3(-E), R2(-d)), R3(E)), R3(s));
  }
  const OMEGA_E = 2 * Math.PI * 1.00273781191135448 / DAY;   // rad/s (dERA/dt)
  function geodeticToItrs(lonDeg, latDeg, hM) {           // WGS84 -> km
    const a = 6378.137, f = 1 / 298.257223563, e2 = f * (2 - f);
    const la = latDeg * D2R, lo = lonDeg * D2R, h = hM / 1000;
    const N = a / Math.sqrt(1 - e2 * Math.sin(la) ** 2);
    return [(N + h) * Math.cos(la) * Math.cos(lo), (N + h) * Math.cos(la) * Math.sin(lo), (N * (1 - e2) + h) * Math.sin(la)];
  }
  /** observer geocentric position (km) and velocity (km/s) in GCRS */
  function observerGcrs(lonDeg, latDeg, hM, ts) {
    const r = geodeticToItrs(lonDeg, latDeg, hM);         // polar motion neglected (<1 mm/s)
    const th = era(ts.ut1);
    const M = mxm(cioQ(ts.tt), R3(-th));                  // GCRS <- TIRS
    const v_tirs = [-OMEGA_E * r[1], OMEGA_E * r[0], 0];
    return { pos: mxv(M, r), vel: mxv(M, v_tirs) };
  }

  // ------------------------------------------------------------------ barycentric correction
  /**
   * Barycentric RV correction, Wright & Eastman (2014) eq. 28 as in astropy
   * SkyCoord.radial_velocity_correction (no Shapiro delay, no stellar-motion term).
   * Returns velocities in km/s (optical convention vcorr = c*zb) plus diagnostics.
   */
  function baryCorrection({ raDeg, decDeg, lonDeg, latDeg, hM, isoUtc, dut1 = 0 }) {
    const ts = timeScales(parseISO(isoUtc), dut1);
    const n = unitRaDec(raDeg, decDeg);
    const vE = earthVel(ts.tdb), pE = earthPos(ts.tdb);
    const obs = observerGcrs(lonDeg, latDeg, hM, ts);
    const vObs = add(vE, obs.vel);
    const beta = scl(vObs, 1 / C);
    const gamma = 1 / Math.sqrt(1 - dot(beta, beta));
    const pS = sunPos(ts.tdb);
    const rSun = norm(sub(pS, pE)) * 1e3;                 // m (from geocentre, as astropy)
    // geocentre option: use the geoid-level Earth potential instead of diverging (L_G c = 0.209 m/s)
    const rEarth = Math.max(norm(obs.pos) * 1e3, 6356752.3);
    const rJup = 5.2 * AU_KM * 1e3;                        // mean; term is <1e-3 m/s
    const rMoon = 384400e3;
    const grOverC = -(GM_JUP / rJup + GM_MOON / rMoon + GM_EARTH / rEarth + GM_SUN / rSun) / (C_M * C_M);
    const zb = gamma * (1 + dot(beta, n)) / (1 + grOverC) - 1;
    const vHelio = add(sub(vE, sunVel(ts.tdb)), obs.vel);
    const pObs = add(pE, obs.pos);
    const roemer = dot(pObs, n) / C;                       // s
    return {
      vcorr_kms: zb * C, zb,
      v_orbit_kms: dot(vE, n), v_rot_kms: dot(obs.vel, n),
      v_helio_kms: dot(vHelio, n),
      grav_kms: -grOverC * C, gamma_term_kms: (gamma - 1) * C,
      bjd_tdb: ts.tdb + roemer / DAY, roemer_s: roemer,
      ts,
    };
  }

  // ------------------------------------------------------------------ Galactic & cosmic frames
  // ICRS -> Galactic rotation matrix (Hipparcos Vol.1, eq. 1.5.11)
  const A_G = [[-0.0548755604162154, -0.8734370902348850, -0.4838350155487132],
               [0.4941094278755837, -0.4448296299600112, 0.7469822444972189],
               [-0.8676661490190047, -0.1980763734312015, 0.4559837761750669]];
  function icrsToGal(raDeg, decDeg) {
    const g = mxv(A_G, unitRaDec(raDeg, decDeg));
    let l = Math.atan2(g[1], g[0]) / D2R; if (l < 0) l += 360;
    return { l, b: Math.asin(Math.max(-1, Math.min(1, g[2]))) / D2R, xyz: g };
  }
  function galToIcrs(lDeg, bDeg) {
    const g = [Math.cos(bDeg * D2R) * Math.cos(lDeg * D2R), Math.cos(bDeg * D2R) * Math.sin(lDeg * D2R), Math.sin(bDeg * D2R)];
    const e = mxv(tr(A_G), g);
    let ra = Math.atan2(e[1], e[0]) / D2R; if (ra < 0) ra += 360;
    return { ra, dec: Math.asin(e[2]) / D2R };
  }
  /** projection of an apex velocity (|v| toward (l,b)) on direction (l0,b0) */
  function apexProjection(v, lA, bA, l, b) {
    return v * (Math.cos(b * D2R) * Math.cos(bA * D2R) * Math.cos((l - lA) * D2R) + Math.sin(b * D2R) * Math.sin(bA * D2R));
  }
  function uvwProjection([U, V, W], l, b) {   // U toward GC, V toward rotation, W toward NGP
    return U * Math.cos(b * D2R) * Math.cos(l * D2R) + V * Math.cos(b * D2R) * Math.sin(l * D2R) + W * Math.sin(b * D2R);
  }

  // ------------------------------------------------------------------ air / vacuum
  // sigma = 1e4 / lambda_vac[A]  (vacuum wavenumber in um^-1). Each returns n-1 of standard air.
  const AIR = {
    edlen1953: s => 6.4328e-5 + 2.94981e-2 / (146 - s * s) + 2.5540e-4 / (41 - s * s),
    edlen1966: s => 1e-8 * (8342.13 + 2406030 / (130 - s * s) + 15997 / (38.9 - s * s)),
    peckreeder1972: s => 1e-8 * (8060.51 + 2480990 / (132.274 - s * s) + 17455.7 / (39.32957 - s * s)),   // eq. 3 (5 parameters; NIST ASD)
    peckreeder1972eq2: s => 1e-8 * (5791817 / (238.0185 - s * s) + 167909 / (57.362 - s * s)),            // eq. 2 (4 parameters; specutils/PyAstronomy)
    birchdowns1994: s => 1e-8 * (8342.54 + 2406147 / (130 - s * s) + 15998 / (38.9 - s * s)),
    morton2000: s => 8.34254e-5 + 2.406147e-2 / (130 - s * s) + 1.5998e-4 / (38.9 - s * s),
    ciddor1996: s => 0.05792105 / (238.0185 - s * s) + 0.00167917 / (57.362 - s * s),
  };
  // Greisen et al. 2006 (FITS WCS III) eq. 65 (= IUGG/IAG 1999 phase-refractivity formula): a 0 C, 101325 Pa, dry, 375 ppm CO2 convention, NOT
  // 15 C standard air (-77 mA at 5000 A).  Evaluated at the vacuum wavelength as Greisen et al. recommend
  // (lambda_a = lambda/n(lambda)); the exact inverse of eq. 65 (n at the air wavelength) differs by < 0.08 mA.
  AIR.greisen2006 = s => 1e-6 * (287.6155 + 1.62887 * s * s + 0.01360 * s ** 4);

  function vacToAir(wlVacA, formula = 'ciddor1996') {
    const f = AIR[formula]; if (!f) throw new Error('unknown formula ' + formula);
    const s = 1e4 / wlVacA;
    return wlVacA / (1 + f(s));
  }
  function airToVac(wlAirA, formula = 'ciddor1996') {   // exact inverse by fixed-point iteration
    if (formula === 'greisen2006') return wlAirA * (1 + AIR.greisen2006(1e4 / wlAirA));   // eq. 64: n at the AIR wavelength
    let v = wlAirA;
    for (let i = 0; i < 20; i++) {
      const nv = wlAirA * (1 + AIR[formula](1e4 / v));
      if (Math.abs(nv - v) < 1e-12) { v = nv; break; }
      v = nv;
    }
    return v;
  }
  // VALD3 (N. Piskunov) closed-form air->vac and IAU/Morton vac->air
  function airToVacVALD(wlAirA) {
    const s = 1e4 / wlAirA;
    const n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s * s) + 0.0001599740894897 / (38.92568793293 - s * s);
    return wlAirA * n;
  }

  // Ciddor (1996) full model: T [C], P [Pa], RH [0-1], xCO2 [ppm]; returns n-1 at vacuum lambda
  function ciddorFull(wlVacA, T = 15, P = 101325, RH = 0, xc = 450) {
    const s2 = (1e4 / wlVacA) ** 2;
    const TK = T + 273.15;
    // reference refractivities
    const ras = 1e-8 * (5792105 / (238.0185 - s2) + 167917 / (57.362 - s2));
    const rws = 1.022e-8 * (295.235 + 2.6422 * s2 - 0.032380 * s2 * s2 + 0.004028 * s2 ** 3);
    const raxs = ras * (1 + 0.534e-6 * (xc - 450));
    const Ma = 1e-3 * (28.9635 + 12.011e-6 * (xc - 400));
    const Mw = 0.018015, R = 8.314510;
    // saturation vapour pressure (IAPWS / Ciddor eq.)
    const A = 1.2378847e-5, B = -1.9121316e-2, Cc = 33.93711047, D = -6.3431645e3;
    // Ciddor (1996) App. A (Davis 1992) over water; App. C (Marti & Mauersberger 1993) over ice below 0 C
    const svp = T >= 0 ? Math.exp(A * TK * TK + B * TK + Cc + D / TK) : Math.pow(10, -2663.5 / TK + 12.537);
    const fenh = 1.00062 + 3.14e-8 * P + 5.6e-7 * T * T;
    const xw = RH * fenh * svp / P;
    function Z(p, t, x) {
      const tc = t - 273.15;
      const a0 = 1.58123e-6, a1 = -2.9331e-8, a2 = 1.1043e-10, b0 = 5.707e-6, b1 = -2.051e-8, c0 = 1.9898e-4, c1 = -2.376e-6, d = 1.83e-11, e = -0.765e-8;
      return 1 - (p / t) * (a0 + a1 * tc + a2 * tc * tc + (b0 + b1 * tc) * x + (c0 + c1 * tc) * x * x) + (p / t) ** 2 * (d + e * x * x);
    }
    const Za = Z(101325, 288.15, 0);
    const Zw = Z(1333, 293.15, 1);
    const rhoaxs = 101325 * Ma / (Za * R * 288.15);
    const rhows = 1333 * Mw / (Zw * R * 293.15);
    const Zm = Z(P, TK, xw);
    const rhoa = P * Ma * (1 - xw) / (Zm * R * TK);
    const rhow = P * Mw * xw / (Zm * R * TK);
    return (rhoa / rhoaxs) * raxs + (rhow / rhows) * rws;
  }

  // ------------------------------------------------------------------ Doppler
  const doppler = {
    zFromOptical: v => v / C,
    zFromRelativistic: v => Math.sqrt((1 + v / C) / (1 - v / C)) - 1,
    radioToZ: v => 1 / (1 - v / C) - 1,
    relativisticFromZ: z => C * ((1 + z) ** 2 - 1) / ((1 + z) ** 2 + 1),
    compose: (...zs) => zs.reduce((p, z) => p * (1 + z), 1) - 1,
  };

  const SpecCore = {
    C, C_M, timeScales, parseISO, jdFromDate, taiMinusUtc, tdbMinusTt, era, cipXY, cioQ,
    setEphemeris, earthPos, earthVel, sunPos, sunVel, observerGcrs, geodeticToItrs,
    baryCorrection, icrsToGal, galToIcrs, apexProjection, uvwProjection, unitRaDec,
    AIR, vacToAir, airToVac, airToVacVALD, ciddorFull, doppler,
  };
  if (typeof module !== 'undefined' && module.exports) module.exports = SpecCore;
  else root.SpecCore = SpecCore;
})(typeof window !== 'undefined' ? window : globalThis);
