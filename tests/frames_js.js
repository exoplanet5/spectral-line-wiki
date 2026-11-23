const S = require('../src/spectro_core.js');
const [F, ra, dec] = JSON.parse(require('fs').readFileSync(0, 'utf8'));
const out = [];
for (let i = 0; i < ra.length; i++) {
  const g = S.icrsToGal(ra[i], dec[i]); const o = {};
  for (const k of ['lsrk', 'lsrd', 'lsr', 'gsr']) { const f = F[k]; o[k] = f.uvw ? S.uvwProjection(f.uvw, g.l, g.b) : S.apexProjection(f.v, f.apex[0], f.apex[1], g.l, g.b); }
  const a = S.galToIcrs(g.l, g.b); o.rt = Math.hypot(a.ra - ra[i], a.dec - dec[i]);
  out.push(o);
}
console.log(JSON.stringify(out));
