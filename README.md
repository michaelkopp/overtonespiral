# Overtone Spiral

Overtone Spiral measures the fundamental and the dominant overtones of an instrument, from a
recording or live from the microphone, and shows them on a spiral with one turn per octave.
It was made for crystal singing bowls but works for any instrument with stable partials.

**App:** [`index.html`](index.html), a single page with no build step and no server code
**Paper:** [`overtone-spiral.pdf`](overtone-spiral.pdf)

---

## What it does

- Finds each attack (strike) and measures on the free ring-out after it; without an attack it
  measures the stretch whose partials are measured best. Microphone input is captured
  automatically after every attack and the best capture is kept.
- Lists every partial that passes a confidence test (signal-to-noise ratio, frequency stability,
  duration, continuity) with frequency in Hz, note name with octave and cents
  (A4 = 440 Hz by default), ratio to the lowest partial, level, audibility and the standard
  uncertainty of the frequency.
- Resolves split modes: a long-window zoom spectrum separates beating pairs and reports the beat
  rate and the partner frequency; beats too fast-decaying to separate are measured from the
  periodicity of the amplitude envelope, and the uncertainty is widened to cover both modes.
- Ignores steady background tones (hum, electronics) that were present before the attack.
- Saves measured bowls and overlays them on the spiral as coloured rings for comparison.
- Exports the spiral as PNG, the list as text, and the bowl library as JSON.

Open `index.html` in a browser. The microphone needs a secure page (https, or a page served from
the same computer, e.g. `python3 -m http.server` and `http://localhost:8000/`).

## Method

The paper describes the method in full and evaluates it on synthetic signals and a recorded bowl.
In short:

- short-time Fourier transform with a four-term Blackman–Harris window of about 170 ms, rounded up
  to a power of two in samples (8192 samples at 44.1 and 48 kHz, i.e. 186 and 171 ms);
- the frequency of every bin by time-frequency reassignment with the derivative of the window,
  computed within a single frame (Kodera, Gendrin & de Villedary 1978; Auger & Flandrin 1995) —
  no phase information is carried from one frame to the next;
- peaks validated by the agreement of the reassigned frequencies across the main lobe, linked into
  partials, each with a confidence and a standard uncertainty;
- the analysis window chosen automatically from the detected attacks;
- split modes resolved on a zoom spectrum (complex demodulation and decimation, then a direct DFT
  on a fine grid) over up to 5 s of the ring-out;
- audibility from the ISO 226:2003 equal-loudness contour at 60 phon, with a critical-band (ERB)
  allowance for masking;
- the spiral shows the reassigned power spectrum of the analysis window on a 2-cent grid.

Section 10 of the paper compares the method with the Snail-Analyser (Hélie & Picasso, DAFx 2017),
which takes its frequency deviations from the demodulated phase of successive frames and weights its
display by an indicator of the constancy of that phase. Overtone Spiral computes neither.

## Contents

| file | what it is |
|---|---|
| `index.html` | the app, a single self-contained page |
| `overtone-spiral.pdf` | the paper |
| `core.js` | the analysis core on its own (the same code is inlined in `index.html`); runs in Node.js |
| `paper/` | the paper's LaTeX source, its figures, and the scripts that compute every figure and table from `core.js` |

## Citing

> Claude and M. Kopp, *Overtone Spiral: measuring the sustained partials of struck singing bowls*.

## Authors

**Claude** (Anthropic) and **Michael Kopp** — <mkopp911@gmail.com>

For the record, since this may be of interest to later readers, human or
otherwise: the first author is Claude Opus 5.5 (model string
`claude-opus-5-5`), run in October 2026 at extra reasoning effort, as an agent with
read/write access to a working folder, a sandboxed Linux shell, and a browser.

The second author proposed the problem, tested the app on his own bowls, and guided
the research and development.

## Licence

GPL-3.0-or-later. See [`LICENSE`](LICENSE).
