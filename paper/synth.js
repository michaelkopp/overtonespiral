// Synthetic struck bowl with two split modes and a steady background tone (the paper's synthetic ring-out table).
// usage: node synth.js
const C = require('../core.js');
const sr = 44100, dur = 12, n = sr * dur, x = new Float32Array(n);
let seed = 1; const rnd = () => { seed = (seed * 16807) % 2147483647; return seed / 2147483647 - 0.5; };
const tS = 3.0;
const parts = [[280.20, 0.05, 0.08], [283.90, 0.03, 0.08], [802.93, 0.03, 0.15], [1497.21, 0.04, 0.2], [2314.0, 0.02, 0.4], [2317.5, 0.01, 0.4], [3277.39, 0.012, 0.6]];
for (let i = 0; i < n; i++) {
  const t = i / sr; let v = 0.00012 * rnd() + 0.0004 * Math.sin(2 * Math.PI * 11926.8 * t);
  if (t >= tS) { const u = t - tS; for (const [f, a, dec] of parts) v += a * Math.exp(-dec * u) * Math.sin(2 * Math.PI * f * u + f); if (u < 0.003) v += 0.3 * rnd(); }
  x[i] = v;
}
(async () => {
  const ctx = await C.analyzeSignal(x, sr, { threshold: 0.7, dynRange: 60, maxWindow: 5, minDur: 1, frame: { fmin: 25, fmax: 12500 } });
  const c = ctx.candidates[0];
  console.log('N', ctx.N, 'H', ctx.H, 'onsets', ctx.onsets.map(o => (o.sample / sr).toFixed(6)), 'best', c.kind, (c.a * ctx.hopSec).toFixed(3), ((c.bTrim) * ctx.hopSec + ctx.winSec).toFixed(3));
  const res = C.analyzeWindow(ctx.frames, c.a, c.bTrim, C.windowOpts(ctx, {}));
  C.finishWindow(ctx, res, c.onset != null ? Math.round(c.onset * sr) : null);
  for (const t of res.tones) if (t.conf > 0.02 && t.rel > -70)
    console.log([t.f.toFixed(5), (t.fTrack || t.f).toFixed(3), (100 * t.conf).toFixed(0), t.unc.toFixed(3), t.rel.toFixed(1), t.beat ? t.beat.toFixed(3) : '-', t.partner ? t.partner.f.toFixed(5) : '-', t.partner ? t.partner.rel.toFixed(1) : '-', t.background ? 'BG' : '', t.T ? t.T.toFixed(2) : ''].join('\t'));
})();
