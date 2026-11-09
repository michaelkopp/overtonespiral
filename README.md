# Overtone Spiral

Overtone Spiral measures the fundamental and the dominant overtones of an instrument, from a
recording or live from the microphone. It shows the partials on a spiral with one turn per octave
and as a list that can be exported. It was made for crystal singing bowls but works for any
instrument with stable partials.

It has two modes, chosen at the top of the page: **Measure** measures the partials of a bowl after
each strike, with a confidence and an uncertainty for each; **Observe** shows the spectrum in real
time, for watching one or several instruments as they sound.

**App:** <https://overtonespiral.tystnadsklangen.se/>
**Paper:** [`overtone-spiral.pdf`](overtone-spiral.pdf)

---

## How to use

1. Click **Measure** at the top of the page. The app opens in the mode it was last used in; in
   *Observe* it shows the spectrum in real time but measures nothing (see below).
2. Click **Microphone** and allow the browser to use your microphone.
3. Strike your singing bowl.
4. Wait a few seconds while the app records and analyses the ring-out. A line under the spiral
   shows the progress and says when the result is ready; the list of partials (the fundamental
   and the overtones) then appears. Every further strike is analysed too, and the best
   measurement is kept.
5. Click **Microphone** again to switch it off.
6. Inspect the overtones on the spiral or in the list. **Export list** saves the list as a text
   file, **Export image** saves the spiral as a picture.
7. Click **Save current** under *Saved bowls* to keep the measurement for later use, for example
   to compare it with another bowl: saved bowls are drawn on the spiral as coloured rings, so
   matching overtones show up in the same place. Saved bowls stay in this browser on this
   device; *Export* there makes a file copy that *Import* reads back.

Instead of using the microphone, you can analyse a recording with **Open audio…** or by dropping
the file onto the page.

### On a phone, also offline

On an iPhone, open the app in Safari, tap **Share → Add to Home Screen**, and open it once from the
home screen while online. From then on it starts without a network connection: a service worker keeps
the page (its fonts are part of it) and its icons in the phone's cache and serves them from there,
without using the network. A new version is picked up when the app is opened with a connection; if it
is still running in the background, close it in the app switcher first. (Android: Chrome's menu →
*Add to Home screen* or *Install app*.) The home-screen app keeps its own saved bowls, separate from
Safari's; use *Export* and *Import* under *Saved bowls* to move them across.

If the phone interrupts the microphone (the screen goes dark, a call comes in), the app reconnects it
when you return. Should the phone not allow that without a touch, the status line asks you to tap
anywhere.

The version is shown at the foot of the page, with whether this copy is available offline.
*Settings → About this copy* says more and keeps a short event log (page loads, interruptions of the
microphone, and every request the app made to the network), which *Copy log* puts on the clipboard.

### Observing in real time

1. Click **Observe** at the top of the page.
2. Click **Microphone** and play: the spiral shows what sounds now. The list *Sounding now* names
   the partials that have sounded for about ¾ s, with their frequencies and levels; it is renewed
   once a second, and a partial that stops is dimmed before it leaves. Nothing is captured or
   measured in this mode, so attacks do not matter and several instruments can play at once.
3. **Hold** (or the space bar) freezes the picture; press it again to continue.
4. With a recording, Observe shows the spectrum at the cursor on the timeline: click or drag there
   to look at any moment, press the space bar to play from it, or **▶ All** to play from the start.
5. **Save current** keeps the partials sounding at that moment as a set of rings, like a measured
   bowl. *Settings → Observe* sets the frame length (long by default; shorter follows faster
   changes, longer separates close partials) and how quickly the display follows the sound
   (steady by default).

## What it does

- Finds each attack (strike) and measures on the free ring-out after it; without an attack it
  measures the stretch whose partials are measured best. Microphone input is captured
  automatically after every attack and the best capture is kept.
- Measures with fixed parameters, so that a measurement is always made the same way; the settings
  for the list (confidence threshold, weakest and shortest partial listed) only choose which of the
  measured partials are shown.
- Lists every partial that passes a confidence test (signal-to-noise ratio, frequency stability,
  duration, continuity) with frequency in Hz, note name with octave and cents
  (A4 = 440 Hz by default), ratio to the lowest partial, level, audibility and the standard
  uncertainty of the frequency.
- Resolves split modes: a long-window zoom spectrum separates beating pairs and reports the beat
  rate and the partner frequency; beats too fast-decaying to separate are measured from the
  periodicity of the amplitude envelope, and the uncertainty is widened to cover both modes.
- Ignores steady background tones (hum, electronics) that were present before the attack.
- Saves measured bowls and overlays them on the spiral as coloured rings for comparison.
- Shows, in Observe mode, a real-time spiral that keeps noise and the clicks of strikes and taps
  off the display, with a list of the partials sounding at the moment.
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
- the spiral shows the reassigned power spectrum of the analysis window on a 2-cent grid, weighted
  as in the real-time display below so that noise and transients stay dark (the measurement itself
  does not use these weights); while the microphone is on, it shows the real-time display;
- the real-time display (Observe) moves the power of each bin of 85 ms frames, every 5.3 ms, to its
  reassigned frequency, weighted by its level above the local noise floor and by the agreement of
  the bins of its lobe (the derivative of the reassigned frequency with respect to frequency, 0 for
  a sinusoid and 1 for an impulse; Nelson 2001, 2002), and smooths the levels in decibels with a
  short attack and release, as a level meter does.

Section 11 of the paper compares the method with the Snail-Analyser (Hélie & Picasso, DAFx 2017),
which takes its frequency deviations from the demodulated phase of successive frames and weights its
display by an indicator of the constancy of that phase. Overtone Spiral computes neither, in either
mode: its weights are computed within single frames, and only levels are carried from one frame to
the next.

## Contents

| file | what it is |
|---|---|
| `index.html` | the app, a single self-contained page |
| `sw.js`, `manifest.webmanifest`, `icons/` | what makes the page installable and usable offline (service worker, web app manifest, icons) |
| `overtone-spiral.pdf` | the paper |
| `core.js` | the analysis core on its own (the same code is inlined in `index.html`); runs in Node.js |
| `paper/` | the paper's LaTeX source, its figures, and the scripts that compute every figure and table from `core.js` |

## Citing

> Claude and M. Kopp, *Overtone Spiral: measuring the sustained partials of struck singing bowls*.

## Authors

**Claude** (Anthropic) and **Michael Kopp** — <mkopp911@gmail.com>

For the record, since this may be of interest to later readers, human or
otherwise: the first author is Claude Opus 5.5 (model string
`claude-opus-5-5`), run in October and November 2026 at extra reasoning effort, as an agent with
read/write access to a working folder, a sandboxed Linux shell, and a browser.

The second author proposed the problem, tested the app on his own bowls, and guided
the research and development.

## Licence

GPL-3.0-or-later. See [`LICENSE`](LICENSE).

The page embeds subsets of three fonts under the SIL Open Font License 1.1: Spectral (Copyright 2017 The
Spectral Project Authors), Hanken Grotesk (Copyright 2021 The Hanken Grotesk Project Authors) and Azeret Mono
(Copyright 2021 The Azeret Project Authors).
