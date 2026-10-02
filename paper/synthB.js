// Case B of the paper: a synthetic ring-out with the parameters of the recorded bowl (close splits, comparable SNR).
// usage: node synthB.js [noise standard deviation, default 0.003 as in Table 2; the paper also quotes 0.0015 and 0]
const C = require('../core.js');
const sr = 48000, dur = 30, n = sr * dur, x = new Float32Array(n);
let seed = 7; const uni = () => { seed = (seed * 16807) % 2147483647; return seed / 2147483647; };
const gauss = () => Math.sqrt(-2 * Math.log(Math.max(1e-12, uni()))) * Math.cos(2 * Math.PI * uni());
const tS = 3.0, sigma = +(process.argv[2] || 0.003);
const tau = dbs => 20 / Math.LN10 / dbs;                 // decay in dB/s -> amplitude time constant
const parts = [[238.87, 0.05, tau(1.1)], [237.76, 0.05 * 0.127, tau(1.1)], [618.46, 0.07, tau(0.8)], [619.42, 0.07 * 0.188, tau(0.8)],
               [1143.27, 0.05, tau(4.8)], [1142.42, 0.05 * 0.55, tau(4.8)], [1799.24, 0.013, tau(27)], [1802.22, 0.009, tau(27)]];
for (let i = 0; i < n; i++) {
  const t = i / sr; let v = sigma * gauss();
  if (t >= tS) { const u = t - tS; for (const [f, a, tc] of parts) v += a * Math.exp(-u / tc) * Math.sin(2 * Math.PI * f * u + f); if (u < 0.003) v += 0.3 * (uni() - 0.5); }
  x[i] = v;
}
(async () => {
  const opts = { threshold: 0.7, dynRange: 60, maxWindow: 5, minDur: 1, frame: { fmin: 25, fmax: 12500 } };
  const ctx = await C.analyzeSignal(x, sr, opts);
  const c = ctx.candidates.find(c => c.score > 0);
  const res = C.analyzeWindow(ctx.frames, c.a, c.bTrim, C.windowOpts(ctx, opts));
  C.finishWindow(ctx, res, c.onset != null ? Math.round(c.onset * sr) : null);
  console.log('sigma', sigma, 'window', (c.a * ctx.hopSec).toFixed(2), (c.bTrim * ctx.hopSec + ctx.winSec).toFixed(2));
  const cents = (a, b) => 1200 * Math.log2(a / b);
  for (const t of res.tones.filter(t => t.conf >= 0.7 && t.rel >= -60 && t.dur >= 1).sort((p, q) => p.f - q.f)) {
    const truth = parts.reduce((b, p) => Math.abs(p[0] - t.f) < Math.abs(b[0] - t.f) ? p : b);
    const pt = t.partner ? parts.reduce((b, p) => Math.abs(p[0] - t.partner.f) < Math.abs(b[0] - t.partner.f) ? p : b) : null;
    console.log([t.f.toFixed(4), 'true', truth[0], 'err ct', cents(t.f, truth[0]).toFixed(4), 'track err ct', cents(t.fTrack || t.f, truth[0]).toFixed(3),
      'u', t.unc.toFixed(3), 'snr', t.snrMed.toFixed(1), 'conf', t.conf.toFixed(2), 'T', (t.T || 0).toFixed(2),
      t.partner ? `partner ${t.partner.f.toFixed(4)} (true ${pt[0]}, err ${(t.partner.f - pt[0]).toFixed(4)} Hz) rel ${t.partner.rel.toFixed(1)} dB beat ${t.beat.toFixed(4)}` : (t.beat ? `beat ${t.beat.toFixed(3)} ${t.beatSrc}` : '')].join(' '));
  }
})();
