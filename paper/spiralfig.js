// Render the spiral of the recorded bowl, as the app draws it, for the paper (light theme).
// usage: node spiralfig.js out.png recording.wav [page-url]
// Serve the app first, e.g. python3 -m http.server 8765 in the repository root; needs Playwright.
const { chromium } = require('playwright');
const fs = require('fs');
const FONTS = {
  'termes-r.otf': '/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyretermes-regular.otf',
  'termes-b.otf': '/usr/share/texmf/fonts/opentype/public/tex-gyre/texgyretermes-bold.otf',
  'mono.ttf': '/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf',
};
(async () => {
  const out = process.argv[2] || 'fig_spiral.png', wav = process.argv[3], url = process.argv[4] || 'http://localhost:8765/index.html';
  const b = await chromium.launch();
  const p = await (await b.newContext({ viewport: { width: 1300, height: 900 } })).newPage();
  // Where the page's web fonts cannot be loaded, substitute local ones (the paper used TeX Gyre Termes and DejaVu Sans Mono).
  const local = Object.values(FONTS).every(f => fs.existsSync(f));
  if (local) await p.route('**/__fonts/*', r => { const n = r.request().url().split('/').pop(); r.fulfill({ body: fs.readFileSync(FONTS[n]), contentType: 'font/otf' }); });
  await p.goto(url);
  if (local) await p.addStyleTag({ content: `
    @font-face { font-family: Spectral; src: url(/__fonts/termes-r.otf); font-weight: 400 500; }
    @font-face { font-family: Spectral; src: url(/__fonts/termes-b.otf); font-weight: 600 700; }
    @font-face { font-family: "Azeret Mono"; src: url(/__fonts/mono.ttf); font-weight: 400 600; }` });
  await p.evaluate(() => Promise.all(['500 20px Spectral', '600 20px Spectral', '400 20px "Azeret Mono"', '600 20px "Azeret Mono"'].map(f => document.fonts.load(f))));
  await p.setInputFiles('#fileIn', wav);
  await p.waitForFunction(() => window.__os && window.__os.st.cells && window.__os.st.sel && !window.__os.st.busy, null, { timeout: 90000 });
  await p.waitForTimeout(500);
  const data = await p.evaluate(() => {
    const W = 1800, scale = W / 700;
    const cv = document.createElement('canvas'); cv.width = W; cv.height = W;
    const g = cv.getContext('2d');
    const { st, S, listTones, noteOf, toHeights, displayPower, drawSpiral, buildText, NCELLS } = window.__os;
    const tones = listTones(st.result);
    const h = new Float32Array(NCELLS);
    toHeights(displayPower(st.cells, st.result), h, S.range, null);
    drawSpiral(g, W, { heights: h, tones: tones.map(t => ({ f: t.f, n: t.n, label: noteOf(t.f).str })), theme: 'light', scale, overlays: [] });
    return { url: cv.toDataURL('image/png'), tones: tones.map(t => [t.n, noteOf(t.f).str, t.f.toFixed(2), t.aud && t.aud.cat.key, (100 * t.conf).toFixed(0)]),
             sel: [st.sel.a * st.ctx.hopSec, st.sel.b * st.ctx.hopSec + st.ctx.winSec], text: buildText() };
  });
  fs.writeFileSync(out, Buffer.from(data.url.split(',')[1], 'base64'));
  fs.writeFileSync(out.replace(/\.png$/, '.txt'), data.text);
  console.log(JSON.stringify({ tones: data.tones, sel: data.sel }));
  await b.close();
})();
