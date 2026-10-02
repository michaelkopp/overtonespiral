// Single-frame frequency accuracy: reassignment vs quadratic interpolation, against the Cramer-Rao bound.
// usage: node accuracy.js > accuracy.json
const C = require('../core.js');
const sr = 48000, fa = new C.FrameAnalyzer(sr, {}), N = fa.N, binHz = fa.binHz;
let seed = 12345;
const uni = () => { seed = (seed * 16807) % 2147483647; return seed / 2147483647; };
const gauss = () => { const u = Math.max(1e-12, uni()), v = uni(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v); };
const cents = (a, b) => 1200 * Math.log2(a / b);
const x = new Float64Array(N);
const A = 0.1, F0 = 1000;
const out = { N, sr, binHz, snr: [], reass: [], qifft: [], crb: [], bias: { off: [], reass: [], qifft: [] } };

function estimate(f, ph, sigma, tau) {
  for (let n = 0; n < N; n++) x[n] = A * Math.exp(-n / sr / tau) * Math.cos(2 * Math.PI * f * n / sr + ph) + sigma * gauss();
  fa.frame(x, 0);
  // largest bin near f
  let k = Math.round(f / binHz);
  for (const j of [k - 1, k + 1]) if (fa.db[j] > fa.db[k]) k = j;
  const a = fa.db[k - 1], b = fa.db[k], c = fa.db[k + 1], den = a - 2 * b + c;
  const q = (k + (den < 0 ? 0.5 * (a - c) / den : 0)) * binHz;
  return { r: fa.ifr[k], q };
}

// RMS error vs SNR (per-sample SNR = A^2 / (2 sigma^2)), stationary sinusoid, random offset and phase
for (let s = -20; s <= 60; s += 5) {
  const snr = Math.pow(10, s / 10), sigma = A / Math.sqrt(2 * snr);
  let er = 0, eq = 0; const M = 300;
  for (let m = 0; m < M; m++) {
    const f = F0 + uni() * binHz, ph = 2 * Math.PI * uni();
    const e = estimate(f, ph, sigma, Infinity);
    er += cents(e.r, f) ** 2; eq += cents(e.q, f) ** 2;
  }
  // Rife & Boorstyn (1974): var(omega) >= 12 / (SNR N (N^2 - 1)), omega in rad/sample
  const sdHz = sr / (2 * Math.PI) * Math.sqrt(12 / (snr * N * (N * N - 1)));
  out.snr.push(s); out.reass.push(Math.sqrt(er / M)); out.qifft.push(Math.sqrt(eq / M)); out.crb.push(1200 / Math.LN2 * sdHz / F0);
}
// noise-free error vs fractional bin offset, stationary and decaying (tau = 2 s)
for (let o = 0; o <= 0.5001; o += 0.025) {
  const f = (Math.round(F0 / binHz) + o) * binHz;
  let wr = 0, wq = 0;
  for (const ph of [0, 0.7, 1.9, 3.1, 4.4, 5.6]) { const e = estimate(f, ph, 0, Infinity); wr = Math.max(wr, Math.abs(cents(e.r, f))); wq = Math.max(wq, Math.abs(cents(e.q, f))); }
  out.bias.off.push(o); out.bias.reass.push(wr); out.bias.qifft.push(wq);
}
let wd = 0;
for (let o = 0; o <= 0.5001; o += 0.05) { const f = (Math.round(F0 / binHz) + o) * binHz; for (const ph of [0, 1.3, 2.6, 4.0]) { const e = estimate(f, ph, 0, 2); wd = Math.max(wd, Math.abs(cents(e.r, f))); } }
out.decayWorst = wd;
process.stdout.write(JSON.stringify(out));
