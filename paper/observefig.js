// Fig. 7: the synthetic scene of observe.js on the spiral, as the app draws it (light theme), at two moments,
// for the unweighted reassigned power averaged frame by frame, Observe without the lobe-coherence weight, and Observe.
// usage: node observefig.js fig_observe.png [page-url]
// Serve the app first, e.g. python3 -m http.server 8765 in the repository root; needs Playwright.
'use strict';
const { chromium } = require('playwright');
const fs = require('fs');
const O = require('./observe.js');
const FONTS = {
  'termes-r.otf': '/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyretermes-regular.otf',
  'termes-b.otf': '/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyretermes-bold.otf',
  'mono.ttf': '/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf',
};
const TIMES = [2.34, 2.64];            // 40 ms after the tap at 2.3 s; 40 ms into the glide from E4 to G4 at 2.6 s
const METHODS = [
  { name: 'unweighted power average', kind: 'measure' },
  { name: 'Observe, level weight only', opts: { coherence: false } },
  { name: 'Observe', opts: {} },
];
(async () => {
  const out = process.argv[2] || 'fig_observe.png', url = process.argv[3] || 'http://localhost:8765/index.html';
  const x = O.scene();
  // display power and top at each moment, from the same runs as Table 3
  const panels = [];
  const runs = METHODS.map(m => O.run(m, x.subarray(0, Math.round(3 * O.SR))));
  for (const t of TIMES) for (const [i, m] of METHODS.entries()) {
    const s = runs[i].find(q => q.t >= t - 1e-9);
    // run() returns display levels v; hand them to the app as powers over a 60 dB range below 0 dB
    panels.push({ label: m.name, t: +s.t.toFixed(3), pow: Array.from(s.v, v => (v > 0 ? Math.pow(10, (60 * v - 60) / 10) : 0)) });
  }
  const b = await chromium.launch();
  const p = await (await b.newContext({ viewport: { width: 1300, height: 900 } })).newPage();
  const local = Object.values(FONTS).every(f => fs.existsSync(f));
  if (local) await p.route('**/__fonts/*', r => { const n = r.request().url().split('/').pop(); r.fulfill({ body: fs.readFileSync(FONTS[n]), contentType: 'font/otf' }); });
  await p.goto(url);
  if (local) await p.addStyleTag({ content: `
    @font-face { font-family: Spectral; src: url(/__fonts/termes-r.otf); font-weight: 400 500; }
    @font-face { font-family: Spectral; src: url(/__fonts/termes-b.otf); font-weight: 600 700; }
    @font-face { font-family: "Azeret Mono"; src: url(/__fonts/mono.ttf); font-weight: 400 600; }` });
  await p.evaluate(() => Promise.all(['500 20px Spectral', '600 20px Spectral', '400 20px "Azeret Mono"'].map(f => document.fonts.load(f))));
  await p.waitForFunction(() => window.__os);
  const data = await p.evaluate((panels) => {
    const { S, toHeights, drawSpiral, NCELLS } = window.__os;
    Object.assign(S, { low: 2, oct: 7, range: 60, lobe: 22, a4: 440 });
    const W = 900, scale = W / 700, cols = 3, rows = panels.length / cols;
    const cv = document.createElement('canvas'); cv.width = cols * W; cv.height = rows * W;
    const g = cv.getContext('2d');
    const tmp = document.createElement('canvas'); tmp.width = W; tmp.height = W; const tg = tmp.getContext('2d');
    const h = new Float32Array(NCELLS);
    panels.forEach((q, i) => {
      toHeights(Float64Array.from(q.pow), h, 60, 0);
      drawSpiral(tg, W, { heights: h, tones: [], theme: 'light', scale, overlays: [] });
      g.drawImage(tmp, (i % cols) * W, Math.floor(i / cols) * W);
    });
    return cv.toDataURL('image/png');
  }, panels);
  fs.writeFileSync(out, Buffer.from(data.split(',')[1], 'base64'));
  await b.close();
  console.log(JSON.stringify(panels.map(q => [q.label, q.t])));
})();
