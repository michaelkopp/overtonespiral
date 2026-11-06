/* Overtone Spiral — analysis core (no DOM). Works in browser and Node.
   Authors: Claude (Anthropic) and Michael Kopp.
   Copyright (C) 2026 Michael Kopp. SPDX-License-Identifier: GPL-3.0-or-later */
(function (root) {
  'use strict';

  const SHARP = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B'];
  const FLAT = ['C', 'Db', 'D', 'Eb', 'E', 'F', 'Gb', 'G', 'Ab', 'A', 'Bb', 'B'];
  const TWO_PI = 2 * Math.PI;

  /* ---------- note naming ---------- */
  function noteInfo(f, a4 = 440, flats = false) {
    const midi = 69 + 12 * Math.log2(f / a4);
    let n = Math.round(midi);
    let c = Math.round((midi - n) * 100);
    if (c === 50) { c = -50; n += 1; }
    const pc = ((n % 12) + 12) % 12;
    const name = (flats ? FLAT : SHARP)[pc];
    const octave = Math.floor(n / 12) - 1;
    const sign = c < 0 ? '-' : '+';
    return { midi, n, pc, name, octave, cents: (midi - Math.round(midi)) * 100, centsInt: c,
             str: name + octave + sign + Math.abs(c) };
  }

  /* ---------- real FFT (N real -> N/2+1 complex) via N/2 complex radix-2 ---------- */
  class RealFFT {
    constructor(N) {
      this.N = N; const M = N >> 1; this.M = M;
      const bits = Math.round(Math.log2(M));
      this.rev = new Uint32Array(M);
      for (let i = 0; i < M; i++) { let r = 0, x = i; for (let b = 0; b < bits; b++) { r = (r << 1) | (x & 1); x >>= 1; } this.rev[i] = r; }
      this.cosM = new Float64Array(M >> 1); this.sinM = new Float64Array(M >> 1);
      for (let k = 0; k < (M >> 1); k++) { this.cosM[k] = Math.cos(TWO_PI * k / M); this.sinM[k] = Math.sin(TWO_PI * k / M); }
      this.cosN = new Float64Array(M + 1); this.sinN = new Float64Array(M + 1);
      for (let k = 0; k <= M; k++) { this.cosN[k] = Math.cos(TWO_PI * k / N); this.sinN[k] = Math.sin(TWO_PI * k / N); }
      this.zr = new Float64Array(M); this.zi = new Float64Array(M);
    }
    forward(x, re, im) {
      const M = this.M, zr = this.zr, zi = this.zi, rev = this.rev, cM = this.cosM, sM = this.sinM;
      for (let n = 0; n < M; n++) { const r = rev[n]; zr[r] = x[2 * n]; zi[r] = x[2 * n + 1]; }
      for (let size = 2; size <= M; size <<= 1) {
        const half = size >> 1, step = M / size;
        for (let i = 0; i < M; i += size) {
          for (let j = 0, k = 0; j < half; j++, k += step) {
            const wr = cM[k], wi = -sM[k];
            const a = i + j, b = a + half;
            const tr = wr * zr[b] - wi * zi[b], ti = wr * zi[b] + wi * zr[b];
            zr[b] = zr[a] - tr; zi[b] = zi[a] - ti; zr[a] += tr; zi[a] += ti;
          }
        }
      }
      const cN = this.cosN, sN = this.sinN;
      for (let k = 0; k <= M; k++) {
        const ka = k === M ? 0 : k, kb = k === 0 ? 0 : M - k;
        const a = zr[ka], b = zi[ka], c = zr[kb], d = zi[kb];
        const fer = (a + c) * 0.5, fei = (b - d) * 0.5, fOr = (b + d) * 0.5, fOi = -(a - c) * 0.5;
        const cs = cN[k], sn = sN[k];
        re[k] = fer + cs * fOr + sn * fOi;
        im[k] = fei + cs * fOi - sn * fOr;
      }
    }
  }

  function blackmanHarris(N) {
    const w = new Float64Array(N);
    for (let n = 0; n < N; n++) {
      const x = TWO_PI * n / N;
      w[n] = 0.35875 - 0.48829 * Math.cos(x) + 0.14128 * Math.cos(2 * x) - 0.01168 * Math.cos(3 * x);
    }
    return w;
  }

  function nextPow2(v) { let p = 1; while (p < v) p <<= 1; return p; }

  // Time derivative of the Blackman-Harris window, per sample.
  function blackmanHarrisDeriv(N) {
    const d = new Float64Array(N), c = TWO_PI / N;
    for (let n = 0; n < N; n++) {
      const x = c * n;
      d[n] = c * (0.48829 * Math.sin(x) - 2 * 0.14128 * Math.sin(2 * x) + 3 * 0.01168 * Math.sin(3 * x));
    }
    return d;
  }

  /* ---------- frame analyser ---------- */
  // For one frame: magnitude (dBFS-calibrated for sinusoids), the reassigned frequency of every bin,
  // a local noise floor, and the spectral peaks.
  // Reassigned frequency (Kodera, Gendrin & de Villedary 1978; Auger & Flandrin 1995): with X_h the
  // spectrum through window h and X_dh the spectrum through its time derivative, a component inside
  // the main lobe of bin k sits at  k - N/(2 pi) * Im(X_dh conj(X_h)) / |X_h|^2  bins.
  // Everything is computed within the single frame.
  class FrameAnalyzer {
    constructor(sr, opts = {}) {
      this.sr = sr;
      this.N = opts.N || nextPow2(Math.round(sr * 0.17));
      this.H = opts.hop || (this.N >> 4);
      this.fmin = opts.fmin || 30;
      this.fmax = Math.min(opts.fmax || 12500, sr * 0.45);
      this.snrMin = opts.snrMin != null ? opts.snrMin : 8;
      this.dynMax = opts.dynMax || 90;
      this.absMin = opts.absMin != null ? opts.absMin : -125;
      this.maxPeaks = opts.maxPeaks || 80;
      this.cohMax = opts.cohMax || 0.45;          // max lobe incoherence (bins) for weak peaks
      this.cohSnr = opts.cohSnr || 20;            // peaks above this SNR (dB) skip the coherence test
      const N = this.N, K = N / 2 + 1;
      this.K = K;
      this.binHz = sr / N;
      this.win = blackmanHarris(N);
      this.dwin = blackmanHarrisDeriv(N);
      let s = 0, s2 = 0; for (let i = 0; i < N; i++) { s += this.win[i]; s2 += this.win[i] * this.win[i]; }
      this.sumW = s; this.sumW2 = s2;
      this.ampScale = 2 / s;                           // |X| * ampScale = sinusoid amplitude
      this.fft = new RealFFT(N);
      this.buf = new Float64Array(N);
      this.re = new Float64Array(K); this.im = new Float64Array(K);
      this.dre = new Float64Array(K); this.dim = new Float64Array(K);
      this.db = new Float32Array(K); this.mag2 = new Float32Array(K);
      this.ifr = new Float32Array(K);
      this.floor = new Float32Array(K);
      this.kmin = Math.max(3, Math.ceil(this.fmin / this.binHz));
      this.kmax = Math.min(K - 3, Math.floor(this.fmax / this.binHz));
      this.reScale = N / TWO_PI;                       // radians per sample -> bins
      // noise-floor bands (1/3 octave, at least 32 bins wide)
      this.bands = [];
      let k0 = this.kmin;
      while (k0 < this.kmax) {
        const w = Math.max(32, Math.round(k0 * (Math.pow(2, 1 / 3) - 1)));
        const k1 = Math.min(this.kmax + 1, k0 + w);
        this.bands.push([k0, k1, (k0 + k1 - 1) / 2]);
        k0 = k1;
      }
      let maxW = 0; for (const b of this.bands) maxW = Math.max(maxW, b[1] - b[0]);
      this.tmp = new Float32Array(maxW);
      this.bandMed = new Float32Array(this.bands.length);
    }

    // FFT of x[start .. start+N) (zero outside) through window w into (re, im)
    spectrum(x, start, w, re, im) {
      const N = this.N, b = this.buf, L = x.length;
      for (let n = 0; n < N; n++) { const j = start + n; b[n] = (j >= 0 && j < L) ? x[j] * w[n] : 0; }
      this.fft.forward(b, re, im);
    }

    // Analyse the frame whose window starts at sample `start`.
    frame(x, start) {
      this.lastX = x; this.lastStart = start;
      this.spectrum(x, start, this.win, this.re, this.im);
      this.spectrum(x, start, this.dwin, this.dre, this.dim);
      return this.analyseCurrent(start);
    }

    analyseCurrent(start) {
      const { re, im, dre, dim, db, mag2, ifr, kmin, kmax, K } = this;
      const sc = this.ampScale, binHz = this.binHz, rs = this.reScale;
      let fmax = -300, esum = 0;
      for (let k = 0; k < K; k++) {
        const m2 = re[k] * re[k] + im[k] * im[k];
        mag2[k] = m2;
        const d = 10 * Math.log10(m2 * sc * sc + 1e-30);
        db[k] = d;
        if (k >= kmin && k <= kmax) {
          esum += m2;
          if (d > fmax) fmax = d;
          ifr[k] = m2 > 0 ? (k - rs * (dim[k] * re[k] - dre[k] * im[k]) / m2) * binHz : k * binHz;
        }
      }
      // noise floor: band medians, linearly interpolated
      const bands = this.bands, tmp = this.tmp, med = this.bandMed;
      for (let i = 0; i < bands.length; i++) {
        const [a, b] = bands[i]; const n = b - a;
        const t = tmp.subarray(0, n);
        for (let k = a; k < b; k++) t[k - a] = db[k];
        t.sort();
        med[i] = t[n >> 1];
      }
      // lossy codecs empty whole bands; keep the floor from collapsing there
      let minLow = Infinity;
      for (let i = 0; i < bands.length; i++) if (bands[i][2] * binHz < 5000 && med[i] < minLow) minLow = med[i];
      if (minLow < Infinity) for (let i = 0; i < bands.length; i++) if (med[i] < minLow - 12) med[i] = minLow - 12;
      const fl = this.floor;
      let bi = 0;
      for (let k = kmin; k <= kmax; k++) {
        while (bi < bands.length - 1 && bands[bi + 1][2] <= k) bi++;
        if (k <= bands[0][2]) fl[k] = med[0];
        else if (bi >= bands.length - 1) fl[k] = med[bands.length - 1];
        else { const c0 = bands[bi][2], c1 = bands[bi + 1][2]; const u = (k - c0) / (c1 - c0); fl[k] = med[bi] * (1 - u) + med[bi + 1] * u; }
      }
      // peaks
      const peaks = [];
      const lo = Math.max(fmax - this.dynMax, this.absMin);
      for (let k = kmin + 2; k <= kmax - 2; k++) {
        const v = db[k];
        if (v < lo || v <= db[k - 1] || v < db[k + 1] || v < db[k - 2] || v < db[k + 2]) continue;
        if (v - fl[k] < this.snrMin) continue;
        // Bins inside a sinusoid's main lobe all reassign to the same frequency, noise peaks
        // scatter (Auger & Flandrin 1995).  Applied to weak peaks only: split modes (beating
        // doublets) of strong partials legitimately reassign to different frequencies.
        if (v - fl[k] < this.cohSnr) {
          const coh = Math.max(Math.abs(ifr[k - 1] - ifr[k]), Math.abs(ifr[k + 1] - ifr[k])) / binHz;
          if (coh > this.cohMax) continue;
        }
        const a = db[k - 1], c = db[k + 1], den = a - 2 * v + c;
        const d = den < 0 ? Math.max(-0.5, Math.min(0.5, 0.5 * (a - c) / den)) : 0;
        const fq = (k + d) * binHz;                    // quadratic interpolation of the log magnitude
        const amp = v - 0.25 * (a - c) * d;
        const fr = ifr[k];
        const f = Math.abs(fr - fq) < 0.6 * binHz ? fr : fq;
        peaks.push(f, amp, amp - fl[k]);
      }
      // keep the strongest maxPeaks
      let P = peaks;
      if (peaks.length / 3 > this.maxPeaks) {
        const idx = []; for (let i = 0; i < peaks.length; i += 3) idx.push(i);
        idx.sort((u, v) => peaks[v + 1] - peaks[u + 1]);
        P = []; for (let j = 0; j < this.maxPeaks; j++) { const i = idx[j]; P.push(peaks[i], peaks[i + 1], peaks[i + 2]); }
      }
      const energyDb = 10 * Math.log10(esum * sc * sc / 2 + 1e-30);
      return { start, t: (start + this.N / 2) / this.sr, peaks: Float32Array.from(P), energy: energyDb, max: fmax };
    }

    // Reassigned power spectrum of the current frame, collected into log-frequency cells: each bin's
    // power is moved to its reassigned frequency; bins whose estimate lies outside the window's main
    // lobe are discarded.
    accumulateCells(cells, cfg, weight = 1) {
      const { mag2, ifr, db, kmin, kmax, binHz } = this;
      const lim = 3.5 * binHz;
      const f0 = cfg.f0, invCell = 1200 / cfg.cellCents, nC = cells.length;
      let fm = -300; for (let k = kmin; k <= kmax; k++) if (db[k] > fm) fm = db[k];
      const lo = fm - (cfg.range || 90) - 12;
      const cal = weight * 4 / (this.N * this.sumW2);   // sum of lobe power -> amplitude^2
      for (let k = kmin; k <= kmax; k++) {
        if (db[k] < lo) continue;
        const f = ifr[k];
        if (!(f > 0) || Math.abs(f - k * binHz) > lim) continue;
        const c = Math.round(Math.log2(f / f0) * invCell);
        if (c < 0 || c >= nC) continue;
        cells[c] += mag2[k] * cal;
      }
    }

    // The same for the spiral of a measured window, with the two weights of the real-time display (LiveSpectrum):
    // the lobe coherence exp(-(d f_hat / d f / 0.3)^2), from two more spectra of the current frame, and the
    // level above the local noise floor (smoothstep from 3 to 15 dB). The measurement does not use it.
    accumulateCellsWeighted(cells, cfg, weight = 1) {
      const N = this.N, K = this.K;
      if (!this.twin) {
        this.twin = new Float64Array(N); this.tdwin = new Float64Array(N);
        for (let n = 0; n < N; n++) { const t = n - N / 2; this.twin[n] = t * this.win[n]; this.tdwin[n] = t * this.dwin[n]; }
        this.tre = new Float64Array(K); this.tim = new Float64Array(K); this.tdre = new Float64Array(K); this.tdim = new Float64Array(K);
      }
      this.spectrum(this.lastX, this.lastStart, this.twin, this.tre, this.tim);
      this.spectrum(this.lastX, this.lastStart, this.tdwin, this.tdre, this.tdim);
      const { re, im, dre, dim, tre, tim, tdre, tdim, mag2, ifr, db, floor, kmin, kmax, binHz } = this;
      const lim = 3.5 * binHz, f0 = cfg.f0, inv = 1200 / cfg.cellCents, nC = cells.length;
      const cal = weight * 4 / (N * this.sumW2);
      for (let k = kmin; k <= kmax; k++) {
        const p = mag2[k]; if (!(p > 0)) continue;
        const snr = db[k] - floor[k] - 1.6;        // re the band's mean noise power (the floor is the median of dB)
        if (snr <= 3) continue;
        const f = ifr[k]; if (!(f > 0) || Math.abs(f - k * binHz) > lim) continue;
        const xr = re[k], xi = im[k];
        const dr = (dre[k] * xr + dim[k] * xi) / p, di = (dim[k] * xr - dre[k] * xi) / p;
        const tr = (tre[k] * xr + tim[k] * xi) / p, ti = (tim[k] * xr - tre[k] * xi) / p;
        const m = 1 + (tdre[k] * xr + tdim[k] * xi) / p - (dr * tr - di * ti);
        const w = Math.exp(-m * m / 0.09) * smoothstep(snr, 3, 15);
        if (w < 1e-3) continue;
        const pos = Math.log2(f / f0) * inv, e = p * w * cal, j0 = Math.floor(pos);
        for (let j = j0 - 1; j <= j0 + 2; j++) { const t = 1 - Math.abs(j - pos) / 2; if (t > 0 && j >= 0 && j < nC) cells[j] += 0.5 * t * e; }
      }
    }
  }

  /* ---------- real-time spectrum (Observe mode) ---------- */
  // A display for watching sound as it happens, built from single-frame quantities and power only.
  // Every frame (hop about 5 ms):
  //  1. Reassigned frequency of each bin, as in FrameAnalyzer.
  //  2. Lobe coherence, the derivative of the reassigned frequency with respect to the bin frequency,
  //     d f_hat / d f = 1 + Re( X_tdh / X_h - X_dh X_th / X_h^2 ),
  //     from two more spectra through the time-weighted windows t h(t) and t h'(t). It is 0 for a
  //     stationary sinusoid (all bins of its main lobe reassign to one frequency), 1 for an impulse,
  //     and scattered for noise. It is the mixed partial derivative of the STFT phase (plus one),
  //     with which Nelson (2001, 2002) told line components from impulses; the closed form is that
  //     of Auger, Chassande-Mottin & Flandrin (2012).
  //  3. Each bin's power, weighted by exp(-(d f_hat / d f / 0.3)^2) and by its level above the local
  //     noise floor (smoothstep from 3 to 15 dB), is moved to its reassigned frequency on the
  //     display's log-frequency cells, spread linearly over the neighbouring cells.
  //  4. The cell levels, in dB, follow these deposits with a short attack and a slightly longer
  //     release (first-order, as in a level meter). Averaging in dB favours what is present frame
  //     after frame; sporadic noise deposits barely rise above the floor.
  // No quantity is carried from one frame to the next except the smoothed cell levels.
  class LiveSpectrum {
    constructor(sr, opts = {}) {
      this.sr = sr;
      const N = this.N = opts.N || nextPow2(Math.round(sr * (opts.frame || 0.085)));
      this.H = opts.hop || Math.min(N >> 2, nextPow2(Math.round(sr * 0.004)));   // 256 samples at 44.1/48 kHz
      const K = this.K = N / 2 + 1;
      this.binHz = sr / N;
      this.kmin = Math.max(4, Math.ceil((opts.fmin || 25) / this.binHz));
      this.kmax = Math.min(K - 5, Math.floor(Math.min(opts.fmax || 12500, 0.45 * sr) / this.binHz));
      this.win = blackmanHarris(N); this.dwin = blackmanHarrisDeriv(N);
      this.twin = new Float64Array(N); this.tdwin = new Float64Array(N);
      let s2 = 0;
      for (let n = 0; n < N; n++) { const t = n - N / 2; this.twin[n] = t * this.win[n]; this.tdwin[n] = t * this.dwin[n]; s2 += this.win[n] * this.win[n]; }
      this.cal = 4 / (N * s2);                       // summed lobe power -> sinusoid amplitude^2
      this.fft = new RealFFT(N);
      this.buf = new Float64Array(N);
      const A = () => new Float64Array(K);
      this.Xr = A(); this.Xi = A(); this.Dr = A(); this.Di = A(); this.Tr = A(); this.Ti = A(); this.TDr = A(); this.TDi = A();
      this.P = A(); this.fl = A();
      this.sigma = opts.sigma || 0.3;
      this.useCoherence = opts.coherence !== false; this.useSnr = opts.snrGate !== false;   // both on; off only for evaluation
      // noise-floor bands (1/3 octave, at least 32 bins)
      this.bands = [];
      for (let k0 = this.kmin; k0 <= this.kmax;) {
        const w = Math.max(32, Math.round(k0 * (Math.pow(2, 1 / 3) - 1))), k1 = Math.min(this.kmax + 1, k0 + w);
        this.bands.push([k0, k1, (k0 + k1 - 1) / 2]); k0 = k1;
      }
      this.tmp = new Float64Array(N); this.bandFl = new Float64Array(this.bands.length);
      // display cells
      this.f0 = opts.f0; this.cellCents = opts.cellCents || 2; this.nC = opts.nCells;
      this.inv = 1200 / this.cellCents;
      this.c0 = Math.max(0, Math.floor(Math.log2(this.kmin * this.binHz / this.f0) * this.inv) - 2);
      this.c1 = Math.min(this.nC - 1, Math.ceil(Math.log2(this.kmax * this.binHz / this.f0) * this.inv) + 2);
      this.v0 = this.c0; this.v1 = this.c1;          // cells on view: the top level and the peaks are taken there
      this.dep = new Float64Array(this.nC); this.depN = new Float64Array(this.nC); this.depS = new Float64Array(this.nC);
      this.L = new Float32Array(this.nC);            // smoothed level, dB (amplitude^2 re full scale)
      this.num = new Float64Array(this.nC); this.den = new Float64Array(this.nC);   // for frequency read-out
      this.snrN = new Float64Array(this.nC);         // deposit-weighted level above the noise floor (dB), for the read-out
      this.floorRange = opts.floorRange || 85;       // levels are held no lower than top - floorRange
      this.topMin = opts.topMin != null ? opts.topMin : -50;
      this.setTimes(opts.attack || 0.02, opts.release || 0.07);
      this.reset();
    }
    setTimes(attack, release) {
      const dt = this.H / this.sr;
      this.attack = attack; this.release = release;
      this.aA = 1 - Math.exp(-dt / attack); this.aR = 1 - Math.exp(-dt / release);
    }
    setView(c0, c1) { this.v0 = Math.max(this.c0, c0); this.v1 = Math.min(this.c1, c1); }
    reset() { this.top = -300; this.L.fill(-300); this.num.fill(0); this.den.fill(0); this.snrN.fill(0); this.frames = 0; }
    spectrum(x, start, w, re, im) {
      const N = this.N, b = this.buf, L = x.length;
      for (let n = 0; n < N; n++) { const j = start + n; b[n] = (j >= 0 && j < L) ? x[j] * w[n] : 0; }
      this.fft.forward(b, re, im);
    }
    // Analyse the frame x[start .. start+N) and update the display levels.
    frame(x, start) {
      this.spectrum(x, start, this.win, this.Xr, this.Xi);
      this.spectrum(x, start, this.dwin, this.Dr, this.Di);
      this.spectrum(x, start, this.twin, this.Tr, this.Ti);
      this.spectrum(x, start, this.tdwin, this.TDr, this.TDi);
      const { Xr, Xi, Dr, Di, Tr, Ti, TDr, TDi, P, fl, kmin, kmax, bands } = this;
      for (let k = kmin; k <= kmax; k++) P[k] = Xr[k] * Xr[k] + Xi[k] * Xi[k];
      // local noise floor: median power in each band (median of an exponential variate = ln 2 x mean)
      const tmp = this.tmp, bf = this.bandFl;
      for (let i = 0; i < bands.length; i++) {
        const [a, b] = bands[i], t = tmp.subarray(0, b - a);
        for (let k = a; k < b; k++) t[k - a] = P[k];
        t.sort(); bf[i] = Math.log(t[(b - a) >> 1] / Math.LN2 + 1e-40);
      }
      for (let k = kmin, i = 0; k <= kmax; k++) {
        while (i < bands.length - 1 && bands[i + 1][2] <= k) i++;
        let v;
        if (k <= bands[0][2] || i >= bands.length - 1) v = bf[Math.min(i, bands.length - 1)];
        else { const u = (k - bands[i][2]) / (bands[i + 1][2] - bands[i][2]); v = bf[i] * (1 - u) + bf[i + 1] * u; }
        fl[k] = Math.exp(v);
      }
      // deposit weighted power at the reassigned frequencies
      const dep = this.dep, depN = this.depN, depS = this.depS, c0 = this.c0, c1 = this.c1;
      for (let c = c0; c <= c1; c++) { dep[c] = 0; depN[c] = 0; depS[c] = 0; }
      const rs = this.N / TWO_PI, bin = this.binHz, lim = 3.5 * bin, is2 = 1 / (this.sigma * this.sigma), cal = this.cal;
      for (let k = kmin; k <= kmax; k++) {
        const p = P[k]; if (!(p > 0)) continue;
        const snr = this.useSnr ? 10 * Math.log10(p / (fl[k] + 1e-40)) : 99;
        if (snr <= 3) continue;
        const xr = Xr[k], xi = Xi[k];
        const dr = (Dr[k] * xr + Di[k] * xi) / p, di = (Di[k] * xr - Dr[k] * xi) / p;       // X_dh / X_h
        const f = (k - rs * di) * bin;
        if (!(f > 0) || Math.abs(f - k * bin) > lim) continue;
        const tr = (Tr[k] * xr + Ti[k] * xi) / p, ti = (Ti[k] * xr - Tr[k] * xi) / p;       // X_th / X_h
        const m = 1 + (TDr[k] * xr + TDi[k] * xi) / p - (dr * tr - di * ti);              // d f_hat / d f
        const w = (this.useCoherence ? Math.exp(-m * m * is2) : 1) * smoothstep(snr, 3, 15);
        if (w < 1e-3) continue;
        const pos = Math.log2(f / this.f0) * this.inv, e = p * w * cal;
        // triangular spread over the neighbouring cells, centred on the continuous position
        const j0 = Math.floor(pos);
        for (let j = j0 - 1; j <= j0 + 2; j++) {
          const t = 1 - Math.abs(j - pos) / 2; if (t <= 0 || j < c0 || j > c1) continue;
          dep[j] += 0.5 * t * e; depN[j] += 0.5 * t * e * pos; depS[j] += 0.5 * t * e * Math.min(snr, 60);
        }
      }
      // smoothed levels in dB, held no lower than a floor below the tracked top
      let mx = -300;
      for (let c = this.v0; c <= this.v1; c++) if (dep[c] > 0) { const d = 10 * Math.log10(dep[c]); if (d > mx) mx = d; }
      this.top = Math.max(mx, this.top - 3 * this.H / this.sr, this.topMin);
      const floor = this.top - this.floorRange, L = this.L, aA = this.aA, aR = this.aR, num = this.num, den = this.den, snrN = this.snrN;
      for (let c = c0; c <= c1; c++) {
        const x = dep[c] > 0 ? Math.max(floor, 10 * Math.log10(dep[c])) : floor;
        const y = L[c] < floor - 20 ? floor - 20 : L[c];
        L[c] = y + (x > y ? aA : aR) * (x - y);
        num[c] += aR * (depN[c] - num[c]); den[c] += aR * (dep[c] - den[c]); snrN[c] += aR * (depS[c] - snrN[c]);
      }
      this.frames++;
    }
    // display power (amplitude^2) per cell
    power(out) {
      out.fill(0);
      for (let c = this.c0; c <= this.c1; c++) out[c] = Math.pow(10, this.L[c] / 10);
      return out;
    }
    // the strongest peaks of the display: [{f, level (dB), rel (dB re strongest), snr (dB above the noise floor)}],
    // sorted by frequency
    peaks(dynDb = 60, maxN = 40, minSep = 6) {
      const L = this.L, c0 = Math.max(this.c0 + minSep, this.v0), c1 = Math.min(this.c1 - minSep, this.v1), floor = this.top - this.floorRange;
      let mx = -300; for (let c = c0; c <= c1; c++) if (L[c] > mx) mx = L[c];
      const lo = Math.max(mx - dynDb, floor + 12), out = [];
      for (let c = c0; c <= c1; c++) {
        const v = L[c]; if (v < lo) continue;
        let isMax = true;
        for (let j = 1; j <= minSep && isMax; j++) if (L[c - j] >= v || L[c + j] > v) isMax = false;
        if (!isMax) continue;
        let n = 0, d = 0, q = 0; for (let j = c - 4; j <= c + 4; j++) { n += this.num[j]; d += this.den[j]; q += this.snrN[j]; }
        const pos = d > 0 ? n / d : c;
        out.push({ f: this.f0 * Math.pow(2, (Math.abs(pos - c) < 4 ? pos : c) / this.inv), level: v, rel: v - mx, snr: d > 0 ? q / d : 0 });
      }
      out.sort((a, b) => b.level - a.level);
      return out.slice(0, maxN).sort((a, b) => a.f - b.f);
    }
  }

  /* ---------- steady list of the partials sounding now (Observe mode) ---------- */
  // The display follows the sound within a tenth of a second; a list of numbers that changed as fast could
  // not be read. The peaks of the display are read four times a second and followed from one read-out to
  // the next by frequency (within 25 cents). A peak joins the list once it has been found in 3 of the last 4
  // read-outs (about 0.75 s) with a median level within `range` dB of the strongest, or of `refMin` dB (re full
  // scale) when nothing stronger sounds, and at least `minSnr` dB above the local noise floor; the last two keep
  // background noise out of the list when nothing else sounds. A listed one is marked as fading after 3
  // read-outs without it and leaves after 6 (1.5 s), or when it falls `margin` dB below either limit.
  // Frequency and level shown are medians over the last 1.5 s. At most `max` rows.
  class SteadyList {
    constructor(opts = {}) {
      this.every = opts.every || 0.25; this.join = opts.join || 3; this.of = opts.of || 4;
      this.fadeAfter = opts.fadeAfter || 3; this.leaveAfter = opts.leaveAfter || 6;
      this.range = opts.range || 40; this.margin = opts.margin != null ? opts.margin : 6;
      this.tol = opts.tol || 25; this.max = opts.max || 16; this.avg = opts.avg || 1.5;
      this.minSnr = opts.minSnr != null ? opts.minSnr : 12;
      this.refMin = opts.refMin != null ? opts.refMin : -50;
      this.reset();
    }
    reset() { this.tracks = []; this.next = 0; }
    // call after each analysed frame with the time of its newest sample (s); reads the display when due
    step(ls, at) {
      if (at + 1e-9 < this.next) return false;
      this.next = (this.next > 0 && at - this.next < this.every ? this.next : at) + this.every;
      this.update(ls.peaks(this.range + this.margin, 48), at);
      return true;
    }
    update(P, at) {
      const T = this.tracks;
      for (const t of T) t.hit = null;
      let ref = this.refMin; for (const p of P) if (p.level > ref) ref = p.level;   // the level the range is counted from
      // strongest peaks first, each to the nearest free track within the tolerance
      for (const p of P.slice().sort((a, b) => b.level - a.level)) {
        let best = null, bd = this.tol;
        for (const t of T) { if (t.hit) continue; const d = Math.abs(1200 * Math.log2(p.f / t.fRef)); if (d < bd) { bd = d; best = t; } }
        if (best) best.hit = p;
        else T.push({ fRef: p.f, hist: [], seen: [], miss: 0, listed: false, hit: p });
      }
      for (const t of T) {
        t.seen.push(!!t.hit); if (t.seen.length > this.of) t.seen.shift();
        if (t.hit) { t.fRef = t.hit.f; t.hist.push([at, t.hit.f, t.hit.rel, t.hit.snr, t.hit.level - ref]); t.miss = 0; } else t.miss++;
        while (t.hist.length && t.hist[0][0] < at - 2 * this.avg) t.hist.shift();
        const lev = this.level(t, at, 4), snr = this.level(t, at, 3);
        if (!t.listed) { if (t.seen.filter(Boolean).length >= this.join && lev >= -this.range && snr >= this.minSnr) t.listed = true; }
        else if (t.miss >= this.leaveAfter || lev < -this.range - this.margin || snr < this.minSnr - this.margin) t.gone = true;
      }
      // tracks that have run into each other (a gliding or wide-vibrato partial) are merged
      T.sort((a, b) => a.fRef - b.fRef);
      for (let i = 1; i < T.length; i++) {
        const a = T[i - 1], b = T[i];
        if (a.gone || b.gone || Math.abs(1200 * Math.log2(b.fRef / a.fRef)) >= 0.6 * this.tol) continue;
        const keep = (a.listed && !b.listed) || (a.listed === b.listed && a.hist.length >= b.hist.length) ? a : b;
        (keep === a ? b : a).gone = true;
      }
      this.tracks = T.filter(t => !t.gone && (t.listed || t.miss < this.of));
    }
    level(t, at, k = 2) {                     // median over the last `avg` seconds of the relative level (k = 2), the
                                              // level above the noise floor (k = 3) or the level re the reference (k = 4)
      const h = t.hist.filter(q => q[0] >= at - this.avg);
      return h.length ? median(h.map(q => q[k])) : t.hist.length ? t.hist[t.hist.length - 1][k] : -Infinity;
    }
    // the rows to show: {f, rel, fading}, sorted by frequency. When more than `max` are listed, the strongest
    // are shown, rows already on show counting 3 dB stronger, so that near-equal partials do not take turns,
    // and fading rows 6 dB weaker, so that they give way to partials that sound now.
    rows(at) {
      const R = [];
      for (const t of this.tracks) {
        if (!t.listed) continue;
        let h = t.hist.filter(q => q[0] >= at - this.avg); if (!h.length) h = t.hist.slice(-1);
        if (!h.length) continue;
        R.push({ t, f: median(h.map(q => q[1])), rel: median(h.map(q => q[2])), fading: t.miss >= this.fadeAfter });
      }
      const prio = r => r.rel + (r.fading ? -6 : r.t.shown ? 3 : 0);
      R.sort((a, b) => prio(b) - prio(a));
      for (const t of this.tracks) t.shown = false;
      const out = R.slice(0, this.max);
      for (const r of out) r.t.shown = true;
      return out.sort((a, b) => a.f - b.f).map(r => ({ f: r.f, rel: r.rel, fading: r.fading }));
    }
  }

  /* ---------- onset detector (time domain, block energy in two high-passed bands) ---------- */
  // Band 0 (>150 Hz) follows the overall level and catches new tonal sound; band 1 (>1 kHz)
  // isolates the broadband click of a strike, which is masked in band 0 when the instrument
  // is already sounding.
  function hpBiquad(fc, sr) {
    const w0 = TWO_PI * fc / sr, cw = Math.cos(w0), al = Math.sin(w0) / Math.SQRT2, a0 = 1 + al;
    return { b0: (1 + cw) / 2 / a0, b1: -(1 + cw) / a0, b2: (1 + cw) / 2 / a0, a1: -2 * cw / a0, a2: (1 - al) / a0,
             x1: 0, x2: 0, y1: 0, y2: 0, acc: 0, hist: new Float64Array(64) };
  }
  class OnsetDetector {
    constructor(sr, opts = {}) {
      this.sr = sr;
      this.block = Math.max(64, Math.round(sr * 0.005));
      this.jump = [opts.jumpLow || 10, opts.jumpHigh || 12];
      this.gateDb = opts.gateDb != null ? opts.gateDb : -75;
      this.refractory = Math.round((opts.refractory || 0.3) * sr / this.block);
      this.f = [hpBiquad(opts.hpLow || 150, sr), hpBiquad(opts.hpHigh || 1000, sr)];
      this.cnt = 0; this.blockIdx = 0;
      this.last = -1e9;
      this.envelope = [];         // band-0 block energies (dB)
      this.pending = null;
    }
    // returns completed onsets {sample, jump, level, peak}; `peak` = max band-0 level in the 100 ms after
    push(x, offset = 0, n = x.length) {
      const out = [], F = this.f, B = this.block;
      for (let i = offset; i < offset + n; i++) {
        const v = x[i];
        for (let q = 0; q < 2; q++) {
          const g = F[q];
          const y = g.b0 * v + g.b1 * g.x1 + g.b2 * g.x2 - g.a1 * g.y1 - g.a2 * g.y2;
          g.x2 = g.x1; g.x1 = v; g.y2 = g.y1; g.y1 = y; g.acc += y * y;
        }
        if (++this.cnt === B) {
          this.cnt = 0;
          const bi = this.blockIdx++;
          let fire = false, jumpMax = -99;
          let E0 = 0;
          for (let q = 0; q < 2; q++) {
            const g = F[q], e = g.acc / B; g.acc = 0;
            const E = 10 * Math.log10(e + 1e-20);
            if (q === 0) { E0 = E; this.envelope.push(E); }
            const Hh = g.hist, L = Hh.length;
            let s = 0, c = 0;
            for (let j = 4; j <= 24; j++) { const idx = bi - j; if (idx < 0) break; s += Hh[idx % L]; c++; }
            Hh[bi % L] = e;
            if (c >= 8) {
              const jump = E - 10 * Math.log10(s / c + 1e-20);
              if (jump >= this.jump[q]) fire = true;
              jumpMax = Math.max(jumpMax, jump);
            }
          }
          if (this.pending) {
            this.pending.peak = Math.max(this.pending.peak, E0);
            if (bi - this.pending.block >= 20) { out.push(this.pending); this.pending = null; }
          }
          if (fire && E0 >= this.gateDb && bi - this.last >= this.refractory) {
            this.last = bi;
            if (this.pending) out.push(this.pending);
            this.pending = { sample: Math.max(0, bi - 1) * B, block: bi, jump: jumpMax, level: E0, peak: E0 };
          }
        }
      }
      return out;
    }
    flush() { const p = this.pending; this.pending = null; return p ? [p] : []; }
  }

  /* ---------- partial tracking ---------- */
  function trackFrames(frames, a, b, opt) {
    const hopSec = opt.hopSec, binHz = opt.binHz;
    const maxGap = Math.max(2, Math.round((opt.maxGapSec || 0.12) / hopSec));
    const tolC = opt.tolCents || 25;
    const active = [], done = [];
    for (let i = a; i < b; i++) {
      const fr = frames[i]; if (!fr) continue;
      const P = fr.peaks, n = P.length / 3;
      const order = []; for (let j = 0; j < n; j++) order.push(j);
      order.sort((u, v) => P[3 * v + 1] - P[3 * u + 1]);
      for (const tr of active) tr.used = false;
      for (const j of order) {
        const f = P[3 * j], d = P[3 * j + 1], s = P[3 * j + 2];
        const tolHz = Math.max(0.5 * binHz, f * (Math.pow(2, tolC / 1200) - 1));
        let best = null, bd = Infinity;
        for (const tr of active) {
          if (tr.used) continue;
          const dd = Math.abs(f - tr.fRef);
          if (dd < tolHz && dd < bd) { bd = dd; best = tr; }
        }
        if (best) {
          best.fi.push(i); best.f.push(f); best.d.push(d); best.s.push(s);
          best.fRef += 0.3 * (f - best.fRef); best.last = i; best.used = true;
        } else {
          active.push({ fi: [i], f: [f], d: [d], s: [s], fRef: f, last: i, used: true });
        }
      }
      for (let q = active.length - 1; q >= 0; q--) {
        if (i - active[q].last > maxGap) { done.push(active[q]); active.splice(q, 1); }
      }
    }
    for (const tr of active) done.push(tr);
    return done;
  }

  function wmedian(vals, w) {
    const idx = vals.map((_, i) => i).sort((p, q) => vals[p] - vals[q]);
    let tot = 0; for (const x of w) tot += x;
    let acc = 0;
    for (const i of idx) { acc += w[i]; if (acc >= tot / 2) return vals[i]; }
    return vals[idx[idx.length - 1]];
  }
  function median(a) { const s = Float64Array.from(a).sort(); const n = s.length; return n ? (n % 2 ? s[(n - 1) / 2] : 0.5 * (s[n / 2 - 1] + s[n / 2])) : NaN; }
  function smoothstep(x, a, b) { const t = Math.max(0, Math.min(1, (x - a) / (b - a))); return t * t * (3 - 2 * t); }

  // Statistics + confidence for a set of points belonging to one tone.
  function toneStats(pts, opt) {
    const hopSec = opt.hopSec, winSec = opt.winSec;
    const n = pts.f.length;
    const w = pts.s.map(s => Math.pow(10, Math.min(s, 45) / 10));
    const fm = wmedian(pts.f, w);
    // cents deviations, robust outlier rejection
    const dc = pts.f.map(f => 1200 * Math.log2(f / fm));
    const absd = dc.map(Math.abs);
    const mad = Math.max(0.3, median(absd));
    const lim = Math.max(6 * 1.4826 * mad, 3);
    let sw = 0, sf = 0, nIn = 0, wAll = 0;
    for (let i = 0; i < n; i++) { wAll += w[i]; if (absd[i] <= lim) { sw += w[i]; sf += w[i] * dc[i]; nIn++; } }
    const mc = sf / sw;
    let sv = 0;
    for (let i = 0; i < n; i++) if (absd[i] <= lim) sv += w[i] * (dc[i] - mc) * (dc[i] - mc);
    const sigma = Math.sqrt(sv / sw);
    // robust spread: weighted median absolute deviation.  Beating doublets throw the frame
    // frequency off at the beat minima, where the level (and hence the weight) is low.
    const sigmaR = 1.4826 * Math.max(0.03, wmedian(absd, w));
    const limF = Math.max(3 * sigmaR, 1.5);
    let swF = 0, sfF = 0;
    for (let i = 0; i < n; i++) if (absd[i] <= limF) { swF += w[i]; sfF += w[i] * dc[i]; }
    const f = fm * Math.pow(2, (swF > 0 ? sfF / swF : 0) / 1200);
    const first = pts.fi[0], last = pts.fi[n - 1];
    // core interval: shortest run of frames holding 90 % of the (SNR) weight
    let wt = 0; for (const x of w) wt += x;
    let i0 = 0, acc = 0, bestLen = Infinity, c0 = 0, c1 = n - 1;
    for (let i1 = 0; i1 < n; i1++) {
      acc += w[i1];
      while (acc - w[i0] >= 0.9 * wt) { acc -= w[i0]; i0++; }
      if (acc >= 0.9 * wt) { const len = pts.fi[i1] - pts.fi[i0] + 1; if (len < bestLen) { bestLen = len; c0 = i0; c1 = i1; } }
    }
    const coreN = c1 - c0 + 1, coreSpan = pts.fi[c1] - pts.fi[c0] + 1;
    const dur = n * hopSec;
    const presence = coreN / coreSpan;
    const dmax = Math.max(...pts.d);
    const snrMed = median(pts.s);
    const snrMax = Math.max(...pts.s);
    const nEff = Math.max(1, dur / winSec);
    const unc = Math.sqrt(sigmaR * sigmaR + 0.05) / Math.sqrt(nEff);
    // confidence components (0..1)
    const cSnr = smoothstep(snrMed, 6, 24);
    const cDur = 1 - Math.exp(-dur / (opt.tauDur || 0.2));
    const cPres = Math.sqrt(Math.min(1, presence));
    const cStab = 1 / (1 + Math.pow(sigmaR / (opt.stabCents || 6), 2));
    const cInl = sw / wAll;                       // weight fraction of consistent frames
    const parts = { snr: cSnr, dur: cDur, pres: cPres, stab: cStab, inl: cInl, bg: 1 };
    return { f, sigma, sigmaR, unc, dur, presence, dmax, snrMed, snrMax, conf: confOf(parts), first, last, n, parts,
             coreFirst: pts.fi[c0], coreLast: pts.fi[c1] };
  }
  function confOf(p) { return p.snr * p.dur * p.pres * p.stab * Math.sqrt(p.inl) * p.bg; }

  // Merge track fragments with the same frequency, compute stats, filter artefacts.
  function analyzeWindow(frames, a, b, opt) {
    const tracks = trackFrames(frames, a, b, opt).filter(t => t.f.length >= 3);
    // per-track weighted median frequency
    for (const t of tracks) { t.fm = median(t.f); t.dm = Math.max(...t.d); }
    tracks.sort((p, q) => p.fm - q.fm);
    const groups = [];
    for (const t of tracks) {
      // two peaks less than 2 bins apart cannot coexist in one frame: such tracks are fragments
      const g = groups.find(g => Math.abs(t.fm - g.fm) < 2 * opt.binHz ||
        (Math.abs(1200 * Math.log2(t.fm / g.fm)) < (opt.mergeCents || 12) && !overlaps(g, t)));
      if (g) { g.items.push(t); g.fm = (g.fm * g.n + t.fm * t.f.length) / (g.n + t.f.length); g.n += t.f.length; g.fi0 = Math.min(g.fi0, t.fi[0]); g.fi1 = Math.max(g.fi1, t.fi[t.fi.length - 1]); }
      else groups.push({ items: [t], fm: t.fm, n: t.f.length, fi0: t.fi[0], fi1: t.fi[t.fi.length - 1] });
    }
    function overlaps(g, t) {
      const t0 = t.fi[0], t1 = t.fi[t.fi.length - 1];
      let ov = 0;
      for (const u of g.items) { const u0 = u.fi[0], u1 = u.fi[u.fi.length - 1]; ov += Math.max(0, Math.min(t1, u1) - Math.max(t0, u0)); }
      return ov > 0.25 * (t1 - t0 + 1);
    }
    let tones = [];
    for (const g of groups) {
      const pts = { fi: [], f: [], d: [], s: [] };
      const items = g.items.slice().sort((p, q) => p.fi[0] - q.fi[0]);
      if (items.length === 1) { const t = items[0]; pts.fi.push(...t.fi); pts.f.push(...t.f); pts.d.push(...t.d); pts.s.push(...t.s); }
      else {
        const byFrame = new Map();
        for (const t of items) for (let k = 0; k < t.fi.length; k++) {
          const q = byFrame.get(t.fi[k]);
          if (!q || t.d[k] > q[1]) byFrame.set(t.fi[k], [t.f[k], t.d[k], t.s[k]]);
        }
        for (const fi of [...byFrame.keys()].sort((u, v) => u - v)) { const q = byFrame.get(fi); pts.fi.push(fi); pts.f.push(q[0]); pts.d.push(q[1]); pts.s.push(q[2]); }
      }
      if (pts.f.length * opt.hopSec < (opt.minDur || 0.08)) continue;
      const st = toneStats(pts, opt);
      st.series = pts;
      tones.push(st);
    }
    if (!tones.length) return { a, b, tones: [], score: 0, ref: -Infinity };
    // reference level = strongest tone peak
    // suppress weak neighbours within the main lobe of a much stronger tone (window side lobes, codec artefacts)
    tones.sort((p, q) => q.dmax - p.dmax);
    const kept = [];
    const lobeHz = 4.5 * opt.binHz;
    for (const t of tones) {
      const strong = kept.find(k => Math.abs(k.f - t.f) < lobeHz && k.dmax - t.dmax > 25);
      if (!strong) kept.push(t);
    }
    tones = kept;
    // the strongest *confident* tone sets the reference level
    let ref = -Infinity;
    for (const t of tones) if (t.conf >= 0.3 && t.dmax > ref) ref = t.dmax;
    if (ref === -Infinity) ref = tones[0].dmax;
    for (const t of tones) t.rel = t.dmax - ref;
    tones.sort((p, q) => p.f - q.f);
    return { a, b, tones, ref };
  }

  /* ---------- high-resolution zoom spectrum (doublets, beating) ---------- */
  // Blackman-Harris weighting over [s0, s1), heterodyne of f0 to 0 Hz, decimation to ~400 Hz with a
  // triangular (boxcar squared) kernel, then a direct DFT on a fine grid over +-B Hz.
  function zoomSpectrum(x, s0, s1, f0, sr, B, df) {
    s0 = Math.max(0, Math.floor(s0)); s1 = Math.min(x.length, Math.floor(s1));
    const L = s1 - s0; if (L < sr * 0.3) return null;
    const D = Math.max(1, Math.round(sr / 400));
    const M = Math.floor((L - 2 * D) / D) + 1; if (M < 16) return null;
    const zr = new Float64Array(M), zi = new Float64Array(M);
    const w0 = TWO_PI * f0 / sr, er = Math.cos(w0), ei = -Math.sin(w0);
    const u0 = TWO_PI / L, ur = Math.cos(u0), ui = Math.sin(u0);
    let cr = 1, ci = 0, wr = 1, wi = 0;
    for (let n = 0; n < L; n++) {
      const c1 = wr, c2 = 2 * c1 * c1 - 1, c3 = 4 * c1 * c1 * c1 - 3 * c1;
      const v = x[s0 + n] * (0.35875 - 0.48829 * c1 + 0.14128 * c2 - 0.01168 * c3);
      const yr = v * cr, yi = v * ci;
      const m1 = (n / D) | 0, k1 = n - m1 * D;
      if (m1 < M) { const t = k1 + 1; zr[m1] += t * yr; zi[m1] += t * yi; }
      if (m1 >= 1 && m1 - 1 < M) { const t = D - k1; zr[m1 - 1] += t * yr; zi[m1 - 1] += t * yi; }
      let q = cr * er - ci * ei; ci = cr * ei + ci * er; cr = q;
      q = wr * ur - wi * ui; wi = wr * ui + wi * ur; wr = q;
      if ((n & 1023) === 1023) { let g = 1 / Math.hypot(cr, ci); cr *= g; ci *= g; g = 1 / Math.hypot(wr, wi); wr *= g; wi *= g; }
    }
    const P = Math.floor(2 * B / df) + 1, db = new Float32Array(P);
    const dt = D / sr;
    for (let p = 0; p < P; p++) {
      const off = -B + p * df, a = -TWO_PI * off * dt, ar = Math.cos(a), ai = Math.sin(a);
      let rr = 1, ri = 0, sr_ = 0, si = 0;
      for (let m = 0; m < M; m++) {
        sr_ += zr[m] * rr - zi[m] * ri; si += zr[m] * ri + zi[m] * rr;
        const q = rr * ar - ri * ai; ri = rr * ai + ri * ar; rr = q;
      }
      db[p] = 10 * Math.log10(sr_ * sr_ + si * si + 1e-30);
    }
    return { f0, B, df, db, T: L / sr };
  }
  function zoomPeaks(z) {
    const { db, df, B } = z, out = [];
    for (let p = 1; p < db.length - 1; p++) {
      if (db[p] > db[p - 1] && db[p] >= db[p + 1]) {
        const a = db[p - 1], b = db[p], c = db[p + 1], den = a - 2 * b + c;
        const d = den < 0 ? 0.5 * (a - c) / den : 0;
        out.push({ off: -B + (p + d) * df, db: b - 0.25 * (a - c) * d });
      }
    }
    return out.sort((u, v) => v.db - u.db);
  }

  // Refine the frequency of each long-lived tone on a zoomed long-window spectrum and detect a
  // beating partner (split mode).  x: samples, `offset`: sample index in x of frame 0's start.
  function refineTones(x, sr, res, opt) {
    const H = opt.H, N = opt.N, offset = opt.offset || 0;
    const winS0 = res.a * H - offset, winS1 = res.b * H + N - offset;
    const cands = res.tones.filter(t => t.conf >= 0.2 && !t.background);
    for (const t of cands) {
      const s0 = Math.max(winS0, t.coreFirst * H - offset);
      const s1 = Math.min(winS1, t.coreLast * H + N - offset, s0 + Math.round((opt.maxSec || 5) * sr));
      const T = (s1 - s0) / sr;
      if (T < 0.8) continue;
      const f0 = t.fTrack != null ? t.fTrack : t.f;          // the frequency tracked from the frames (kept if called again)
      let near = Infinity;
      for (const u of cands) if (u !== t) near = Math.min(near, Math.abs(u.f - f0));
      const B = Math.max(4 / T, Math.min(12, 0.45 * near));
      const z = zoomSpectrum(x, s0, s1, f0, sr, B, Math.min(0.05, 0.2 / T));
      if (!z) continue;
      const pk = zoomPeaks(z); if (!pk.length) continue;
      const main = pk[0];
      // an unresolved pair is tracked at its weighted mean, up to about one frame bin from either mode
      if (Math.abs(main.off) > Math.max(1.0, f0 * (Math.pow(2, 8 / 1200) - 1), 0.8 * (opt.binHz || 0))) continue;
      t.fTrack = f0; t.f = f0 + main.off; t.refined = true; t.T = T;
      const sepMin = 3 / T, depth = Math.min(25, Math.max(10, t.snrMed - 10));
      const partner = pk.find(q => q !== main && Math.abs(q.off - main.off) >= sepMin && q.db >= main.db - depth &&
                                     Math.abs(q.off) < B - 2 * z.df);
      if (partner) {
        t.partner = { f: t.fTrack + partner.off, rel: partner.db - main.db };
        t.beat = Math.abs(partner.off - main.off);
        // a resolved doublet explains frame-to-frame frequency wobble
        t.parts.stab = Math.max(t.parts.stab, 0.95);
        t.conf = confOf(t.parts);
        t.beatSrc = 'spectrum';
      }
    }
    // Beats too fast-decaying to separate spectrally still show as a periodic amplitude envelope.
    for (const t of res.tones) if (!t.beat && t.conf >= 0.2 && !t.background) envelopeBeat(t, opt.hopSec || H / sr);
    // A partial without a resolved partner may be an unresolved pair, measured as a blend of its two modes: the zoom
    // spectrum's shift of the frequency, or half the beat rate seen in the envelope, whichever is larger, is added to
    // the uncertainty.
    for (const t of res.tones) {
      if (t.partner) continue;
      const shift = t.refined ? Math.abs(1200 * Math.log2(t.f / t.fTrack)) : 0;
      const half = t.beat ? 1200 * Math.log2((t.f + t.beat / 2) / t.f) : 0;
      const e = Math.max(shift, half);
      if (e > 0) { const u0 = t.uncFrames != null ? t.uncFrames : t.unc; t.uncFrames = u0; t.unc = Math.sqrt(u0 * u0 + e * e); }
    }
  }

  function envelopeBeat(t, hopSec) {
    const s = t.series, i0 = t.first, i1 = t.last, n = i1 - i0 + 1;    // the whole track: a fast decay leaves a short core
    if (n * hopSec < 0.6) return;
    const env = new Float64Array(n).fill(NaN);
    for (let k = 0; k < s.fi.length; k++) { const j = s.fi[k] - i0; if (j >= 0 && j < n) env[j] = s.d[k]; }
    let last = -1;                                  // fill gaps linearly
    for (let j = 0; j < n; j++) if (!isNaN(env[j])) {
      if (last >= 0 && j - last > 1) for (let q = last + 1; q < j; q++) env[q] = env[last] + (env[j] - env[last]) * (q - last) / (j - last);
      last = j;
    }
    let first = 0; while (first < n && isNaN(env[first])) first++;
    for (let j = 0; j < first; j++) env[j] = env[first];
    for (let j = last + 1; j < n; j++) env[j] = env[last];
    // detrend (linear in dB = exponential decay)
    let mx = 0, my = 0; for (let j = 0; j < n; j++) { mx += j; my += env[j]; } mx /= n; my /= n;
    let sxy = 0, sxx = 0; for (let j = 0; j < n; j++) { sxy += (j - mx) * (env[j] - my); sxx += (j - mx) * (j - mx); }
    const sl = sxy / sxx, r = new Float64Array(n);
    let e0 = 0; for (let j = 0; j < n; j++) { r[j] = env[j] - my - sl * (j - mx); e0 += r[j] * r[j]; }
    const sorted = Float64Array.from(r).sort(), swing = sorted[Math.floor(0.95 * (n - 1))] - sorted[Math.floor(0.05 * (n - 1))];
    if (swing < 2 || e0 <= 0) return;
    const L0 = Math.max(2, Math.round(0.15 / hopSec)), L1 = Math.min(Math.round(2.5 / hopSec), Math.floor(n / 2));
    const ac = new Float64Array(L1 + 2);
    for (let L = L0 - 1; L <= L1 + 1; L++) { let a = 0; for (let j = 0; j + L < n; j++) a += r[j] * r[j + L]; ac[L] = a / e0 * n / (n - L); }
    let bestL = -1, bestV = 0;
    for (let L = L0; L <= L1; L++) if (ac[L] > ac[L - 1] && ac[L] >= ac[L + 1] && ac[L] > bestV) { bestV = ac[L]; bestL = L; }
    if (bestL < 0 || bestV < 0.5 || n < 2 * bestL) return;
    const a = ac[bestL - 1], b = ac[bestL], c = ac[bestL + 1], den = a - 2 * b + c;
    const Lr = bestL + (den < 0 ? 0.5 * (a - c) / den : 0);
    t.beat = 1 / (Lr * hopSec); t.beatSrc = 'envelope'; t.beatSwing = swing;
  }

  // Background (room hum, electronics) is not part of the instrument.  A partial is background if it
  // was already sounding in the reference frames at about the same level and does not decay:
  //  mode 'attack': reference = the 0.5 s before the attack, compared with the first 0.4 s after it;
  //  mode 'quiet' : reference = the quietest frames of the signal, compared with the tone's median level.
  function markBackground(frames, base, ref, res, opt, mode) {
    if (ref.length < 5) return;
    const hopSec = opt.hopSec;
    for (const t of res.tones) {
      if (t.background) continue;
      const tolHz = Math.max(1.5 * opt.binHz, t.f * (Math.pow(2, 20 / 1200) - 1));
      const lv = []; let tot = 0;
      for (const i of ref) {
        const fr = frames[i - base]; if (!fr) continue; tot++;
        const P = fr.peaks; let m = -Infinity;
        for (let j = 0; j < P.length; j += 3) if (Math.abs(P[j] - t.f) < tolHz && P[j + 1] > m) m = P[j + 1];
        if (m > -Infinity) lv.push(m);
      }
      if (!tot || lv.length < 0.5 * tot) continue;
      const ser = t.series;
      let n = 0, sx = 0, sy = 0, sxx = 0, sxy = 0;
      for (let k = 0; k < ser.fi.length; k++) { const xx = ser.fi[k] * hopSec, yy = ser.d[k]; n++; sx += xx; sy += yy; sxx += xx * xx; sxy += xx * yy; }
      const slope = n > 2 ? (n * sxy - sx * sy) / (n * sxx - sx * sx) : 0;
      let rise;
      if (mode === 'attack') {
        const lim = res.a + Math.round(0.4 / hopSec);
        let post = -Infinity; for (let k = 0; k < ser.fi.length && ser.fi[k] <= lim; k++) post = Math.max(post, ser.d[k]);
        rise = post - Math.max(...lv);
        if (rise < 4 && Math.abs(slope) < 0.3) t.background = true;
      } else {
        rise = median(ser.d) - median(lv);
        if (rise < 6 && Math.abs(slope) < 0.3) t.background = true;
      }
      t.slope = slope; t.rise = rise;
      if (t.background) { t.parts.bg = 0.05; t.conf = confOf(t.parts); }
    }
  }

  // Frames in the quietest 10 % of a frame list, if they are clearly quieter than `level`.
  function quietFrames(frames, base, a, b, level) {
    const idx = [];
    for (let i = a; i < b; i++) if (frames[i - base]) idx.push(i);
    if (idx.length < 50) return [];
    const en = idx.map(i => frames[i - base].energy).sort((u, v) => u - v);
    const p10 = en[Math.floor(0.1 * en.length)];
    if (level - p10 < 10) return [];
    const out = idx.filter(i => frames[i - base].energy <= p10);
    const step = Math.max(1, Math.floor(out.length / 200));
    return out.filter((_, k) => k % step === 0);
  }
  function windowLevel(frames, base, a, b) {
    const e = []; for (let i = a; i < b; i++) { const fr = frames[i - base]; if (fr) e.push(fr.energy); }
    return e.length ? median(e) : -Infinity;
  }

  function scoreWindow(res, thr, dyn, minDur = 0) {
    let s = 0, n = 0;
    for (const t of res.tones) {
      // measurement quality first: a few partials at high SNR beat many weak ones
      if (t.conf >= thr && t.rel >= -dyn && t.dur >= minDur) { const q = smoothstep(t.snrMed, 8, 45); s += t.conf * q * q; n++; }
    }
    return { score: s, count: n };
  }

  /* ---------- offline analysis of a whole signal ---------- */
  // Calls progress(frac) periodically; returns a promise.
  async function analyzeSignal(x, sr, opts = {}, progress) {
    const fa = new FrameAnalyzer(sr, opts.frame || {});
    const H = fa.H, N = fa.N;
    const nFrames = Math.max(0, Math.floor((x.length - N) / H) + 1);
    const frames = new Array(nFrames);
    const energy = new Float32Array(nFrames);
    let tLast = Date.now();
    for (let i = 0; i < nFrames; i++) {
      frames[i] = fa.frame(x, i * H);
      energy[i] = frames[i].energy;
      if (progress && (i & 63) === 0 && Date.now() - tLast > 40) {
        progress(i / nFrames);
        await new Promise(r => setTimeout(r, 0));
        tLast = Date.now();
      }
    }
    const od = new OnsetDetector(sr, opts.onset || {});
    let onsets = od.push(x, 0, x.length).concat(od.flush());
    const env = od.envelope, envBlock = od.block;
    const ctx = { sr, H, N, hopSec: H / sr, winSec: N / sr, binHz: fa.binHz, frames, nFrames, energy, onsets, env, envBlock, fa, x };
    ctx.candidates = findCandidates(ctx, opts);
    return ctx;
  }

  function windowOpts(ctx, opts) {
    return Object.assign({ hopSec: ctx.hopSec, winSec: ctx.winSec, binHz: ctx.binHz }, opts.track || {});
  }

  // Background check (after an attack) and high-resolution refinement for a window of a whole signal.
  function finishWindow(ctx, res, onsetSample) {
    const wopt = { hopSec: ctx.hopSec, binHz: ctx.binHz };
    if (onsetSample != null) {
      const preB = Math.floor((onsetSample - ctx.N) / ctx.H), preA = Math.max(0, preB - Math.round(0.5 / ctx.hopSec));
      const ref = []; for (let i = preA; i <= preB; i++) ref.push(i);
      markBackground(ctx.frames, 0, ref, res, wopt, 'attack');
    }
    const q = quietFrames(ctx.frames, 0, 0, ctx.nFrames, windowLevel(ctx.frames, 0, res.a, res.b));
    markBackground(ctx.frames, 0, q, res, wopt, 'quiet');
    if (ctx.x) refineTones(ctx.x, ctx.sr, res, { H: ctx.H, N: ctx.N, offset: 0, hopSec: ctx.hopSec, binHz: ctx.binHz });
    return res;
  }

  function findCandidates(ctx, opts = {}) {
    const thr = opts.threshold != null ? opts.threshold : 0.7;
    const dyn = opts.dynRange || 60, minDur = opts.minDur || 0;
    const D = Math.round((opts.maxWindow || 5) / ctx.hopSec);
    const wopt = windowOpts(ctx, opts);
    const settle = Math.round((opts.settle || 0.015) * ctx.sr);
    const cands = [];
    const strongOn = ctx.onsets;
    for (let oi = 0; oi < strongOn.length; oi++) {
      const o = strongOn[oi];
      // first frame whose window starts after the impact
      const a = Math.ceil((o.sample + settle) / ctx.H);
      if (a >= ctx.nFrames - 5) continue;
      let b = Math.min(ctx.nFrames, a + D);
      // a later, comparably strong excitation ends the window
      for (let oj = oi + 1; oj < strongOn.length; oj++) {
        const p = strongOn[oj];
        if (p.sample <= o.sample) continue;
        const fb = Math.floor(p.sample / ctx.H) - Math.ceil(ctx.N / ctx.H);
        if (fb >= b) break;
        if (p.peak >= o.peak - 3) { b = Math.max(a + 5, fb); break; }
      }
      const res = analyzeWindow(ctx.frames, a, b, wopt);
      finishWindow(ctx, res, o.sample);
      const sc = scoreWindow(res, thr, dyn, minDur);
      cands.push(Object.assign(res, { kind: 'onset', onset: o.sample / ctx.sr, score: sc.score * 1.15, count: sc.count }));
    }
    // sliding windows for sustained sounds
    const step = Math.max(1, Math.round(D / 2));
    const onF = strongOn.map(o => o.sample / ctx.H);
    for (let a = 0; a + Math.min(D, ctx.nFrames) <= ctx.nFrames; a += step) {
      const b = Math.min(ctx.nFrames, a + D);
      // a window that contains an attack mixes two sounds; the attack's own window covers it
      if (onF.some(f => f > a - ctx.N / ctx.H && f < b)) { if (b === ctx.nFrames) break; continue; }
      const res = analyzeWindow(ctx.frames, a, b, wopt);
      finishWindow(ctx, res, null);
      const sc = scoreWindow(res, thr, dyn, minDur);
      cands.push(Object.assign(res, { kind: 'sustain', score: sc.score, count: sc.count }));
      if (b === ctx.nFrames) break;
    }
    // trim each window's end to the last frame where a confident tone is still present
    for (const c of cands) {
      let lastF = c.a;
      for (const t of c.tones) if (t.conf >= thr && t.rel >= -dyn && t.dur >= minDur) lastF = Math.max(lastF, t.last + 1);
      c.bTrim = Math.max(Math.min(c.b, lastF), Math.min(c.b, c.a + Math.round(0.5 / ctx.hopSec)));
    }
    cands.sort((p, q) => q.score - p.score);
    return cands;
  }

  const api = { noteInfo, RealFFT, blackmanHarris, FrameAnalyzer, LiveSpectrum, SteadyList, OnsetDetector, trackFrames, analyzeWindow, finishWindow,
                zoomSpectrum, zoomPeaks, refineTones, markBackground, quietFrames, windowLevel, confOf,
                scoreWindow, analyzeSignal, findCandidates, windowOpts, toneStats, nextPow2, median, smoothstep };
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.OSCore = api;
})(typeof self !== 'undefined' ? self : this);
