const S = require('../src/spectro_core.js');
S.setEphemeris(require('../data/ephem_cheb.json'));
const E = require('../data/ephem_cheb.json');
// 1) ephemeris test vectors
let e1 = 0, e2 = 0;
for (const t of E.test_vectors) {
  const v = S.earthVel(t.jd), p = S.earthPos(t.jd);
  for (let k = 0; k < 3; k++) { e1 = Math.max(e1, Math.abs(v[k] - t.earth_vel[k])); e2 = Math.max(e2, Math.abs(p[k] - t.earth_pos[k])); }
}
console.log(`ephemeris: max |dv| = ${(e1 * 1e6).toFixed(3)} mm/s, max |dr| = ${e2.toFixed(3)} km`);
// 2) astropy barycentric cases
const cases = require('./bary_cases.json');
let worst = 0, worstH = 0, worstB = 0, worstT = 0, sq = 0;
for (const c of cases) {
  const r = S.baryCorrection({ raDeg: c.ra, decDeg: c.dec, lonDeg: c.lon, latDeg: c.lat, hM: c.h, isoUtc: c.iso.replace(/\.\d+$/, m => m), dut1: c.dut1 });
  const d = (r.vcorr_kms - c.vcorr) * 1e6;      // mm/s
  sq += d * d; worst = Math.max(worst, Math.abs(d));
  worstH = Math.max(worstH, Math.abs(r.v_helio_kms - c.vhelio) * 1e6);
  worstT = Math.max(worstT, Math.abs(r.ts.tdb - c.jd_tdb) * 86400e6);
  worstB = Math.max(worstB, Math.abs(r.bjd_tdb - c.bjd) * 86400);
}
console.log(`barycentric vs astropy (${cases.length} cases): rms ${Math.sqrt(sq / cases.length).toFixed(3)} mm/s, max ${worst.toFixed(3)} mm/s`);
console.log(`heliocentric max diff ${worstH.toFixed(2)} mm/s;  TDB max diff ${worstT.toFixed(2)} us;  BJD max diff ${(worstB * 1e3).toFixed(4)} ms`);
