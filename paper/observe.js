// Real-time display (Observe mode): the measurements of Section 9.5 (Table 3), from synthetic signals.
// usage: node observe.js > observe.json      (needs nothing but Node.js)
// Compares the unweighted reassigned spectrum shown frame by frame (FrameAnalyzer, 170 ms frames, cell power averaged)
// with the real-time spectrum of Observe mode (LiveSpectrum), with and without its two weights.
'use strict';
const C = require('../core.js');
const SR = 48000, CPO = 600, F0 = 440 * Math.pow(2, (24 - 69) / 12), NC = 9 * CPO, RANGE = 60;
const cellOf = f => Math.log2(f / F0) * CPO;

let seed = 7; const rnd = () => { seed = (seed * 1664525 + 1013904223) >>> 0; return seed / 4294967296; };
const gauss = () => Math.sqrt(-2 * Math.log(rnd() + 1e-12)) * Math.cos(2 * Math.PI * rnd());
// partials {f, a, on, off}; white noise of standard deviation `noise`; every `click` s a tap, a 12 ms burst of
// noise decaying with a 2 ms time constant, of peak standard deviation 0.25
function synth(dur, parts, noise, click) {
  const n = Math.round(dur * SR), x = new Float32Array(n);
  for (const p of parts) {
    let ph = rnd() * 2 * Math.PI; const on = Math.round((p.on || 0) * SR), off = Math.round((p.off != null ? p.off : dur) * SR);
    for (let i = on; i < Math.min(n, off); i++) { ph += 2 * Math.PI * p.f / SR; x[i] += p.a * Math.min(1, (i - on) / 96, (off - i) / 96) * Math.cos(ph); }
  }
  if (noise) for (let i = 0; i < n; i++) x[i] += noise * gauss();
  if (click) for (let t = click; t < dur; t += click) { const i0 = Math.round(t * SR); for (let i = 0; i < 0.012 * SR; i++) x[i0 + i] += 0.25 * gauss() * Math.exp(-i / (0.002 * SR)); }
  return x;
}

// A synthetic scene with several sources (no recording needed): a struck bowl (four partials, three of them
// split pairs, decaying), a sung melody with vibrato and 60 ms glides between its notes, taps, and pink noise.
const TAPS = [1.2, 2.3, 3.9, 4.8], NOTES = [[293.66, 0.5], [329.63, 1.6], [392.0, 2.6], [349.23, 3.5], [329.63, 4.4]], MEL_END = 5.4;
function scene(dur = 6, notes = NOTES) {
  const n = Math.round(dur * SR), x = new Float32Array(n);
  seed = 11;
  // bowl, struck at 0: [f (Hz), amplitude, decay (dB/s), partner offset (Hz), partner level (dB)]
  for (const [f, a, dec, df, pl] of [[238.87, 0.12, 2, -1.11, -18], [618.46, 0.2, 1.5, 0.95, -14], [1143.27, 0.13, 3, -0.84, -5], [1799.24, 0.05, 20, 2.95, -1]]) {
    for (const [ff, aa] of [[f, a], [f + df, a * Math.pow(10, pl / 20)]]) {
      const ph0 = rnd() * 2 * Math.PI, k = Math.log(10) / 20 * dec;
      for (let i = 0; i < n; i++) { const t = i / SR; x[i] += aa * Math.exp(-k * t) * Math.min(1, i / 48) * Math.cos(2 * Math.PI * ff * t + ph0); }
    }
  }
  // voice: 20 harmonics, two formants, vibrato 5.5 Hz +-0.8 %
  const ph = new Float64Array(21).map(() => rnd() * 2 * Math.PI);
  for (let i = Math.round(0.5 * SR); i < Math.round(MEL_END * SR); i++) {
    const t = i / SR; let j = 0; while (j < notes.length - 1 && t >= notes[j + 1][1]) j++;
    let f = notes[j][0];
    if (j > 0 && t - notes[j][1] < 0.06) f = notes[j - 1][0] * Math.pow(notes[j][0] / notes[j - 1][0], (t - notes[j][1]) / 0.06);
    const vib = 1 + 0.008 * Math.sin(2 * Math.PI * 5.5 * (t - 0.5)) * Math.min(1, (t - 0.5) / 0.4);
    const env = Math.min(1, (t - 0.5) / 0.05, (MEL_END - t) / 0.08) * 0.05;
    let v = 0;
    for (let k = 1; k <= 20; k++) {
      const fk = k * f * vib; if (fk > 9000) break; ph[k] += 2 * Math.PI * fk / SR;
      v += (Math.exp(-Math.pow(Math.log(fk / 650), 2) / 0.5) + 0.5 * Math.exp(-Math.pow(Math.log(fk / 2600), 2) / 0.15) + 0.05) / Math.sqrt(k) * Math.sin(ph[k]);
    }
    x[i] += env * v;
  }
  for (const tt of TAPS) { const i0 = Math.round(tt * SR); for (let i = 0; i < 0.012 * SR; i++) x[i0 + i] += 0.25 * gauss() * Math.exp(-i / (0.002 * SR)); }
  let y = 0; for (let i = 0; i < n; i++) { y = 0.98 * y + 0.02 * gauss(); x[i] += 0.01 * y + 0.0008 * gauss(); }
  return x;
}

// display levels v in [0, 1] per cell, as the app maps them, sampled every 20 ms; t = time of the newest sample
function run(m, x) {
  const out = []; let nextT = 0;
  const take = (t, pow, top) => {
    if (t + 1e-9 < nextT) return; nextT += 0.02;
    const v = new Float32Array(NC), lo = top - RANGE;
    for (let c = 0; c < NC; c++) { const u = (10 * Math.log10(pow[c] + 1e-24) - lo) / RANGE; v[c] = u < 0 ? 0 : u > 1 ? 1 : u; }
    out.push({ t, v });
  };
  const visTop = (pow) => { let mx = -300; for (let c = Math.ceil(cellOf(65.41)); c < Math.floor(cellOf(65.41 * 128)); c++) { const d = 10 * Math.log10(pow[c] + 1e-24); if (d > mx) mx = d; } return mx; };
  if (m.kind === 'measure') {
    const fa = new C.FrameAnalyzer(SR, { N: 8192, hop: 512, fmin: 25, fmax: 12500 });
    const live = new Float64Array(NC), one = new Float64Array(NC); let top = -200;
    for (let s = -8192; s + 8192 <= x.length; s += 512) {
      fa.frame(x, s); one.fill(0); fa.accumulateCells(one, { f0: F0, cellCents: 2, range: 100 });
      for (let i = 0; i < NC; i++) live[i] = 0.6 * live[i] + 0.4 * one[i];
      top = Math.max(visTop(live), top - 3 * 512 / SR, -50);        // the app: maximum on view, falling 3 dB/s
      take((s + 8192) / SR, live, top);
    }
  } else {
    const ls = new C.LiveSpectrum(SR, Object.assign({ fmin: 25, fmax: 12500, f0: F0, cellCents: 2, nCells: NC, floorRange: RANGE + 25 }, m.opts));
    ls.setView(Math.ceil(cellOf(65.41)), Math.floor(cellOf(65.41 * 128)));
    const pow = new Float64Array(NC);
    for (let s = -ls.N; s + ls.N <= x.length; s += ls.H) { ls.frame(x, s); take((s + ls.N) / SR, ls.power(pow), ls.top); m.last = ls; }
  }
  return out;
}
const mean = a => a.reduce((u, v) => u + v, 0) / a.length;
const sd = a => { const m = mean(a); return Math.sqrt(mean(a.map(v => (v - m) * (v - m)))); };
function background(S, partsF, t0, t1, thr = 0.03) {
  const cs = []; for (let c = Math.round(cellOf(80)); c <= Math.round(cellOf(4500)); c++) if (partsF.every(f => Math.abs(c - cellOf(f)) > 30)) cs.push(c);
  let vis = 0, n = 0; for (const s of S) if (s.t >= t0 && s.t <= t1) for (const c of cs) { n++; if (s.v[c] > thr) vis++; }
  return vis / n;
}
const peakV = (s, f) => { const c0 = Math.round(cellOf(f)); let m = 0; for (let c = c0 - 15; c <= c0 + 15; c++) m = Math.max(m, s.v[c]); return m; };

const methods = [
  { name: 'unweighted power average', kind: 'measure' },
  { name: 'Observe, no weights', opts: { coherence: false, snrGate: false } },
  { name: 'Observe, level weight only', opts: { coherence: false } },
  { name: 'Observe', opts: {} },
  { name: 'Observe, fast smoothing', opts: { attack: 0.008, release: 0.03 } },
  { name: 'Observe, steady smoothing', opts: { attack: 0.04, release: 0.17 } },
  { name: 'Observe, short frames', opts: { frame: 0.04 } },
  { name: 'Observe, long frames', opts: { frame: 0.16 } },
  { name: 'Observe, long frames and steady smoothing (the default of the implementation)', opts: { frame: 0.16, attack: 0.04, release: 0.17 } },
];
const PARTS = [196.4, 541.7, 1043.2, 1694.6, 2486.3, 3412.8].map((f, i) => ({ f, a: [1, 0.5, 0.25, 0.1, 0.04, 0.02][i] }));
const sigs = {
  quiet: synth(3, PARTS, 0.05), noisy: synth(3, PARTS, 0.2),
  step: synth(3, [{ f: 440, a: 0.5 }, { f: 660, a: 0.5, on: 1, off: 2 }], 0.003),
  clicks: synth(3, [{ f: 440, a: 0.1 }], 0.003, 0.5),
  noisyLong: synth(8, PARTS, 0.2),
};
if (require.main !== module) { module.exports = { synth, scene, run, SR, NC, F0, cellOf, TAPS, NOTES }; return; }
const res = [];
const SCENE = scene(), VIEW = [Math.ceil(cellOf(65.41)), Math.floor(cellOf(65.41 * 128))];
const visible = (S, a, b) => { const v = S.filter(s => s.t >= a && s.t < b); let n = 0; for (const s of v) for (let c = VIEW[0]; c <= VIEW[1]; c++) if (s.v[c] > 0.03) n++; return n / v.length * 2; };   // cents lit
for (const m of methods) {
  const r = { method: m.name };
  for (const k of ['quiet', 'noisy']) {
    const S = run(m, sigs[k]);
    r['bg_' + k] = background(S, PARTS.map(p => p.f), 1, 2.9);
    r['flicker_' + k] = sd(S.filter(s => s.t >= 1 && s.t <= 2.9).map(s => peakV(s, 3412.8))) * RANGE;   // weakest partial, -34 dB
  }
  const St = run(m, sigs.step), vss = mean(St.filter(s => s.t > 1.5 && s.t < 1.95).map(s => peakV(s, 660)));
  const when = (cond, t0) => { const p = St.find(s => s.t > t0 && cond(peakV(s, 660))); return p ? +(1000 * (p.t - t0)).toFixed(0) : null; };
  r.full_ms = when(v => v >= 0.9 * vss, 1); r.gone_ms = when(v => v <= 0.1 * vss, 2);
  r.bg_clicks = background(run(m, sigs.clicks), [440], 0.5, 2.9);
  const Ss = run(m, SCENE), win = 0.15, avg = a => a.reduce((u, v) => u + v, 0) / a.length;
  r.scene_calm = avg([[0.9, 1.15], [2.0, 2.25], [3.0, 3.4], [5.6, 5.95]].map(([a, b]) => visible(Ss, a, b)));
  r.scene_taps = avg(TAPS.map(t => visible(Ss, t, t + win)));
  r.scene_glides = avg(NOTES.slice(1).map(([f, t]) => visible(Ss, t, t + win)));
  if (m.kind !== 'measure') {
    // read-out: the steady list of the stronger noise, renewed every second from 2 to 8 s
    const ls = new C.LiveSpectrum(SR, Object.assign({ fmin: 25, fmax: 12500, f0: F0, cellCents: 2, nCells: NC, floorRange: RANGE + 25 }, m.opts));
    ls.setView(Math.ceil(cellOf(65.41)), Math.floor(cellOf(65.41 * 128)));
    const L = new C.SteadyList(), x = sigs.noisyLong;
    const err = PARTS.map(() => []); let next = 2;
    for (let s = -ls.N; s + ls.N <= x.length; s += ls.H) {
      ls.frame(x, s); const at = ls.frames * ls.H / SR; L.step(ls, at);
      if (at < next - 1e-9) continue; next += 1;
      const P = L.rows(at);
      PARTS.forEach((p, i) => {
        let best = null; for (const q of P) { const e = 1200 * Math.log2(q.f / p.f); if (Math.abs(e) < 20 && (best == null || Math.abs(e) < Math.abs(best))) best = e; }
        if (best != null) err[i].push(best);
      });
    }
    r.readout_rms_cents = err.map(e => (e.length ? Math.sqrt(mean(e.map(v => v * v))) : null));
    r.readout_found = err.map(e => e.length);
  }
  res.push(r);
  console.error(m.name.padEnd(28), `bg ${(100 * r.bg_quiet).toFixed(1)}% / ${(100 * r.bg_noisy).toFixed(1)}%`, `flicker ${r.flicker_quiet.toFixed(1)} / ${r.flicker_noisy.toFixed(1)} dB`,
    `full ${r.full_ms} gone ${r.gone_ms} ms`, `clicks ${(100 * r.bg_clicks).toFixed(1)}%`, `scene lit ${r.scene_calm.toFixed(0)}/${r.scene_taps.toFixed(0)}/${r.scene_glides.toFixed(0)} cents`, r.readout_rms_cents ? 'readout ' + r.readout_rms_cents.map(v => v == null ? '-' : v.toFixed(2)).join('/') + ' c, found ' + r.readout_found.join('/') : '');
}
// cost of one Observe frame (4 FFTs of 4096 points and the cell update)
{
  const ls = new C.LiveSpectrum(SR, { f0: F0, cellCents: 2, nCells: NC }), x = sigs.noisy;
  for (let k = 0; k < 200; k++) ls.frame(x, k * 256);
  const t = process.hrtime.bigint(); for (let k = 0; k < 2000; k++) ls.frame(x, (k * 256) % 100000);
  const ms = Number(process.hrtime.bigint() - t) / 1e6 / 2000;
  console.error(`Observe frame: ${ms.toFixed(3)} ms (hop ${(1000 * ls.H / SR).toFixed(2)} ms)`);
}
// Stability of the list (default display): entries that join or leave between successive lists, per second,
// for the peaks of the display read four times a second and for the steady list renewed every second.
const lists = {};
for (const [name, x] of [['melody', SCENE], ['held note', scene(6, [[329.63, 0.5]])], ['noise only', synth(6, [], 0.001)]]) {
  const ls = new C.LiveSpectrum(SR, { fmin: 25, fmax: 12500, f0: F0, cellCents: 2, nCells: NC, floorRange: RANGE + 25 });
  ls.setView(Math.ceil(cellOf(65.41)), Math.floor(cellOf(65.41 * 128)));
  const L = new C.SteadyList(), raw = [], steady = []; let nextR = 0;
  for (let s = -ls.N; s + ls.N <= x.length; s += ls.H) {
    ls.frame(x, s); const at = ls.frames * ls.H / SR;
    if (L.step(ls, at) && at >= 2) raw.push(ls.peaks(60, 40).map(p => p.f));
    if (at >= nextR) { nextR = at + 1; if (at >= 2) steady.push(L.rows(at).map(r => r.f)); }
  }
  const near = (u, V) => V.some(v => Math.abs(1200 * Math.log2(u / v)) < 25);
  const changes = (A, dt) => { let c = 0; for (let i = 1; i < A.length; i++) c += A[i - 1].filter(u => !near(u, A[i])).length + A[i].filter(u => !near(u, A[i - 1])).length; return A.length > 1 ? c / ((A.length - 1) * dt) : 0; };
  const len = A => { const n = A.map(a => a.length); return n.length ? [Math.min(...n), Math.max(...n)] : [0, 0]; };
  lists[name] = { raw_changes_per_s: +changes(raw, 0.25).toFixed(1), raw_rows: len(raw), steady_changes_per_s: +changes(steady, 1).toFixed(1), steady_rows: len(steady) };
  console.error(`list, ${name}: peaks read 4/s ${lists[name].raw_changes_per_s} changes/s, ${lists[name].raw_rows.join('-')} rows; steady list ${lists[name].steady_changes_per_s} changes/s, ${lists[name].steady_rows.join('-')} rows`);
}
console.log(JSON.stringify({ methods: res, lists }, null, 1));
