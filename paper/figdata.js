// Data for the paper's figures and tables, computed with the published analysis core.
// usage: node figdata.js [recording.f32] > figdata.json
// The recording is raw 32-bit float, mono, 48 kHz (not included in the repository).
const C = require('../core.js');
const fs = require('fs');

const out = {};
const buf = fs.readFileSync(process.argv[2] || __dirname + '/../bowl_mono.f32');
const x = new Float32Array(buf.buffer, buf.byteOffset, buf.length / 4);
const sr = 48000;

// same biquad as the core's onset detector, to plot both band envelopes
function hp(fc) {
  const w0 = 2 * Math.PI * fc / sr, cw = Math.cos(w0), al = Math.sin(w0) / Math.SQRT2, a0 = 1 + al;
  return { b0: (1 + cw) / 2 / a0, b1: -(1 + cw) / a0, b2: (1 + cw) / 2 / a0, a1: -2 * cw / a0, a2: (1 - al) / a0, x1: 0, x2: 0, y1: 0, y2: 0 };
}
function bandEnv(fc, B) {
  const g = hp(fc), e = []; let acc = 0, c = 0;
  for (let i = 0; i < x.length; i++) {
    const v = x[i], y = g.b0 * v + g.b1 * g.x1 + g.b2 * g.x2 - g.a1 * g.y1 - g.a2 * g.y2;
    g.x2 = g.x1; g.x1 = v; g.y2 = g.y1; g.y1 = y; acc += y * y;
    if (++c === B) { e.push(10 * Math.log10(acc / B + 1e-20)); acc = 0; c = 0; }
  }
  return e;
}

(async () => {
  const opts = { threshold: 0.7, dynRange: 60, maxWindow: 5, minDur: 1, frame: { fmin: 25, fmax: 12500 } };   // as the app (anaOpts)
  const ctx = await C.analyzeSignal(x, sr, opts);
  const H = ctx.H, N = ctx.N, hop = ctx.hopSec;
  out.rec = { sr, N, H, hopSec: hop, winSec: ctx.winSec, binHz: ctx.binHz, dur: x.length / sr, nFrames: ctx.nFrames };
  const B = ctx.envBlock;
  const e0 = bandEnv(150, B), e1 = bandEnv(1000, B);
  // decimate the 5 ms envelopes by 4 for plotting
  out.env = { dt: 4 * B / sr, low: [], high: [] };
  for (let i = 0; i + 4 <= e0.length; i += 4) {
    out.env.low.push(10 * Math.log10((Math.pow(10, e0[i] / 10) + Math.pow(10, e0[i + 1] / 10) + Math.pow(10, e0[i + 2] / 10) + Math.pow(10, e0[i + 3] / 10)) / 4));
    out.env.high.push(10 * Math.log10((Math.pow(10, e1[i] / 10) + Math.pow(10, e1[i + 1] / 10) + Math.pow(10, e1[i + 2] / 10) + Math.pow(10, e1[i + 3] / 10)) / 4));
  }
  out.onsets = ctx.onsets.map(o => ({ t: o.sample / sr, jump: o.jump, peak: o.peak }));
  out.candidates = ctx.candidates.map(c => ({ kind: c.kind, t0: c.a * hop, t1: (c.bTrim || c.b) * hop + ctx.winSec, score: c.score, count: c.count }));
  const best = ctx.candidates.find(c => c.score > 0);
  // the app re-analyses the chosen window [a, bTrim)
  const res = C.analyzeWindow(ctx.frames, best.a, best.bTrim, C.windowOpts(ctx, opts));
  C.finishWindow(ctx, res, Math.round(best.onset * sr));
  out.window = { t0: best.a * hop, t1: best.bTrim * hop + ctx.winSec, onset: best.onset, a: best.a, b: best.bTrim };
  const listed = res.tones.filter(t => t.conf >= 0.7 && t.rel >= -60 && t.dur >= 1).sort((p, q) => p.f - q.f);
  const f1 = listed[0].f;
  const row = t => ({ f: t.f, fTrack: t.fTrack || t.f, note: C.noteInfo(t.f).str, ratio: t.f / f1, rel: t.rel, conf: t.conf, unc: t.unc,
    sigmaR: t.sigmaR, snrMed: t.snrMed, dur: t.dur, beat: t.beat || null, beatSrc: t.beatSrc || null,
    partner: t.partner ? t.partner.f : null, partnerRel: t.partner ? t.partner.rel : null, T: t.T || null,
    parts: t.parts, coreFirst: t.coreFirst * hop, coreLast: t.coreLast * hop });
  out.listed = listed.map(row);
  out.rejected = res.tones.filter(t => !(t.conf >= 0.7 && t.rel >= -60 && t.dur >= 1) && t.rel >= -60 && t.conf >= 0.05).map(row);
  out.allCount = res.tones.length;
  // level and frequency series of the listed partials
  out.series = listed.map(t => ({ f: t.f, t: t.series.fi.map(i => i * hop + ctx.winSec / 2), d: t.series.d.slice(), fr: t.series.f.slice() }));
  // zoom spectra as used by the refinement step
  out.zoom = [];
  const winS0 = res.a * H, winS1 = res.b * H + N;
  for (const t of listed) {
    if (!t.refined) continue;
    const s0 = Math.max(winS0, t.coreFirst * H);
    const s1 = Math.min(winS1, t.coreLast * H + N, s0 + 5 * sr);
    const T = (s1 - s0) / sr;
    let near = Infinity; for (const u of res.tones.filter(u => u.conf >= 0.2 && !u.background)) if (u !== t) near = Math.min(near, Math.abs(u.f - t.f));
    const Bz = Math.max(4 / T, Math.min(12, 0.45 * near));
    const z = C.zoomSpectrum(x, s0, s1, t.fTrack, sr, Bz, Math.min(0.05, 0.2 / T));
    out.zoom.push({ f: t.f, fTrack: t.fTrack, partner: t.partner ? t.partner.f : null, beat: t.beat, T, B: Bz, df: z.df, db: Array.from(z.db) });
  }
  // one frame around the 618 Hz partial: magnitude and reassigned frequency of every bin
  const fa = new C.FrameAnalyzer(sr, opts.frame);
  const fi = best.a + Math.round(1.5 / hop);
  const fr = fa.frame(x, fi * H);
  const kc = Math.round(listed[1].f / fa.binHz);
  out.frame = { t: fr.t, binHz: fa.binHz, k: [], db: [], ifr: [], floor: [] };
  for (let k = kc - 14; k <= kc + 14; k++) { out.frame.k.push(k); out.frame.db.push(fa.db[k]); out.frame.ifr.push(fa.ifr[k]); out.frame.floor.push(fa.floor[k]); }
  out.frame.peaks = []; for (let j = 0; j < fr.peaks.length; j += 3) if (Math.abs(fr.peaks[j] - kc * fa.binHz) < 15 * fa.binHz) out.frame.peaks.push([fr.peaks[j], fr.peaks[j + 1], fr.peaks[j + 2]]);
  process.stdout.write(JSON.stringify(out));
})();
