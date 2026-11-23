const S = require('../src/spectro_core.js');
const wl = [2000, 2500, 3000, 4000, 5000, 6563, 8000, 10000, 15000, 22000];
const out = {};
for (const k of Object.keys(S.AIR)) out[k] = wl.map(w => [S.vacToAir(w, k), S.airToVac(S.vacToAir(w, k), k)]);
out.ciddorFull = [[5000, 15, 101325, 0, 450], [5000, 0, 61500, 0.2, 420], [6563, 20, 101325, 0.5, 450], [8000, 12, 74300, 0.15, 420], [5000, -5, 60000, 0.3, 420]].map(a => [a, S.ciddorFull(...a)]);
console.log(JSON.stringify({ wl, out }));
