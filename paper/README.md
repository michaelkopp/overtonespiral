# Paper source

`overtone-spiral.tex` is the source of [`../overtone-spiral.pdf`](../overtone-spiral.pdf). It builds with
pdfLaTeX (`latexmk -pdf overtone-spiral.tex`) from the figure files in this folder.

The measurements in Section 10 come from the analysis core, [`../core.js`](../core.js), through these scripts
(the timings of Section 10.5 apart, which depend on the machine):

| script | produces |
|---|---|
| `accuracy.js` | `accuracy.json`: single-frame accuracy against the Cramér–Rao bound (Section 10.1, Fig. 3) |
| `synth.js` | Table 2, case A: the widely split synthetic ring-out (Section 10.2), printed to the console |
| `synthB.js` | Table 2, case B: a synthetic ring-out with the recorded bowl's parameters; the argument is the noise level |
| `figdata.js` | `figdata.json`: the recorded bowl (Section 10.3, Figs. 2, 4, 5 and Table 1) |
| `plots.py` | `fig_frame.pdf`, `fig_accuracy.pdf`, `fig_timeline.pdf`, `fig_zoom.pdf` from the two JSON files, and the decay rates and beat swings quoted in Section 10.3 |
| `spiralfig.js` | `fig_spiral.png` and the text export of Table 1, drawn by the app itself |
| `observe.js` | `observe.json`: the real-time display on synthetic signals and a synthetic scene (Section 10.4, Table 3) |
| `observefig.js` | `fig_observe.png`: that scene on the spiral, drawn by the app itself (Fig. 6) |

```
node accuracy.js > accuracy.json
node synth.js
node synthB.js 0.003                           # also 0.0015 and 0, quoted in Section 10.2
node figdata.js bowl.f32 > figdata.json        # raw 32-bit float, mono, 48 kHz
python3 plots.py                               # needs numpy and matplotlib
python3 -m http.server 8765 --directory ..     # in another terminal
node spiralfig.js fig_spiral.png bowl.wav      # needs Playwright
node observe.js > observe.json
node observefig.js fig_observe.png             # needs Playwright and the server above
```

The bowl recording is not part of the repository, but `figdata.json` is, so `plots.py` runs without it;
`accuracy.js`, `synth.js`, `synthB.js` and `observe.js` need nothing but Node.js.
