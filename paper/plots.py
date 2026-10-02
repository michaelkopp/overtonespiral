"""Figures for the paper, from figdata.json and accuracy.json (see figdata.js, accuracy.js)."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['Liberation Serif'], 'mathtext.fontset': 'stix',
    'font.size': 8, 'axes.labelsize': 8, 'axes.titlesize': 8, 'legend.fontsize': 7, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.minor.width': 0.4, 'ytick.minor.width': 0.4,
    'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
    'lines.linewidth': 0.9, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02, 'pdf.fonttype': 42,
})
W = 3.3
D = json.load(open('figdata.json'))
A = json.load(open('accuracy.json'))
COL = ['#1f4e79', '#b03a2e', '#6c3483', '#b7950b']

def tag(ax, s):
    ax.set_title(s, loc='left', pad=2, fontsize=8)

# ---------- Fig. 4: timeline ----------
fig, (a1, a2) = plt.subplots(2, 1, figsize=(W, 3.1), gridspec_kw={'height_ratios': [1, 1.15], 'hspace': 0.42})
env = D['env']; t = np.arange(len(env['low'])) * env['dt']
a1.plot(t, env['low'], color='0.15', lw=0.6, label=r'$>150$ Hz')
a1.plot(t, env['high'], color='#b03a2e', lw=0.6, label=r'$>1$ kHz')
w = D['window']
a1.axvspan(w['t0'], w['t1'], color='#1f4e79', alpha=0.13, lw=0)
for o in D['onsets']:
    a1.axvline(o['t'], color='#1f4e79', lw=0.7, ls=(0, (3, 2)))
a1.set_xlim(0, D['rec']['dur']); a1.set_ylim(-122, -18)
a1.set_xlabel('time (s)'); a1.set_ylabel('block energy (dB)')
a1.legend(loc='lower right', frameon=False, handlelength=1.6, borderaxespad=0.3, ncol=2, columnspacing=1.2)
tag(a1, '(a)')
for k, s in enumerate(D['series']):
    tt = np.array(s['t']); dd = np.array(s['d'])
    ref = max(max(q['d']) for q in D['series'])
    a2.plot(tt, dd - ref, color=COL[k], lw=0.7, label=D['listed'][k]['note'])
a2.set_xlim(w['t0'] - 0.05, w['t1'] + 0.05)
a2.set_ylim(-64, 3)
a2.legend(loc='lower right', frameon=False, handlelength=1.6, ncol=2, columnspacing=1.2, borderaxespad=0.3)
a2.set_xlabel('time (s)'); a2.set_ylabel('partial level (dB)')
tag(a2, '(b)')
fig.savefig('fig_timeline.pdf')
plt.close(fig)

# ---------- Fig. 2: one frame, reassigned frequencies ----------
F = D['frame']; bh = F['binHz']
sl = slice(2, len(F['k']) - 2)
F = {k: (np.array(v)[sl] if isinstance(v, list) and len(v) == len(D['frame']['k']) else v) for k, v in F.items()}
fk = F['k'] * bh
fig, (b1, b2) = plt.subplots(2, 1, figsize=(W, 2.9), sharex=True, gridspec_kw={'height_ratios': [1, 1.1], 'hspace': 0.22})
b1.plot(fk, F['db'], 'o', ms=2.2, color='0.1', zorder=3)
b1.plot(fk, F['db'], color='0.6', lw=0.5)
b1.plot(fk, F['floor'], color='#b03a2e', lw=0.8, ls=(0, (4, 2)), label='noise floor')
b1.set_ylabel('level (dBFS)')
b1.legend(loc='upper right', frameon=False, handlelength=2)
tag(b1, '(a)')
ifr = np.array(F['ifr'])
pk = D['listed'][1]['f']
b2.plot(fk, fk, color='0.6', lw=0.6, ls=(0, (2, 2)), label='bin frequency')
fp = F['peaks'][0][0]
coh = np.abs(ifr - fp) < 0.45 * bh
b2.axhline(fp, color='#1f4e79', lw=0.6)
b2.plot(fk[coh], ifr[coh], 'o', ms=2.6, color='0.1', label='reassigned, coherent', zorder=3)
b2.plot(fk[~coh], ifr[~coh], 'o', ms=2.6, mfc='white', mec='0.35', mew=0.6, label='reassigned, other bins', zorder=3)
b2.set_ylim(fk[0] - 3, fk[-1] + 3)
b2.set_xlim(fk[0] - 3, fk[-1] + 3)
b2.set_xlabel('bin frequency $k f_s/N$ (Hz)'); b2.set_ylabel('frequency (Hz)')
b2.legend(loc='lower right', frameon=False, handlelength=2)
b2.text(fk[0] + 1, fp + 3, f'{fp:.2f} Hz', color='#1f4e79', fontsize=6.8, va='bottom')
tag(b2, '(b)')
fig.savefig('fig_frame.pdf')
plt.close(fig)

# ---------- Fig. 3: single-frame accuracy ----------
fig, ax = plt.subplots(figsize=(W, 2.1))
snr = np.array(A['snr'])
ax.semilogy(snr, A['reass'], color='0.1', marker='o', ms=2.4, label='reassigned frequency')
ax.semilogy(snr, A['qifft'], color='#b03a2e', marker='s', ms=2.2, ls=(0, (4, 2)), label='quadratic interpolation')
ax.semilogy(snr, A['crb'], color='#1f4e79', ls=(0, (1, 1.5)), label='Cramér–Rao bound')
ax.set_xlabel('SNR per sample (dB)'); ax.set_ylabel('RMS error (cents)')
ax.set_xlim(-21, 61); ax.set_ylim(3e-5, 4)
ax.legend(loc='upper right', frameon=False, handlelength=2.4)
fig.savefig('fig_accuracy.pdf')
plt.close(fig)

# ---------- Fig. 5: zoom spectra ----------
Z = [z for z in D['zoom'] if z['partner']]
fig, axs = plt.subplots(len(Z), 1, figsize=(W, 0.95 * len(Z) + 0.45), sharex=True, gridspec_kw={'hspace': 0.32})
for ax, z, lab in zip(axs, Z, ['(a)', '(b)', '(c)']):
    off = -z['B'] + np.arange(len(z['db'])) * z['df']
    db = np.array(z['db']); db -= db.max()
    ax.plot(off, db, color='0.1', lw=0.7)
    ax.axvline(0, color='0.55', lw=0.6, ls=(0, (2, 2)))
    m = z['f'] - z['fTrack']; p = z['partner'] - z['fTrack']
    ax.plot([m], [0], 'v', color='#1f4e79', ms=3.5)
    pi = np.argmin(np.abs(off - p)); ax.plot([p], [db[pi]], 'v', color='#b03a2e', ms=3.5)
    ax.set_xlim(-3, 3); ax.set_ylim(-62, 6)
    ax.set_yticks([-60, -40, -20, 0])
    note = next(q['note'] for q in D['listed'] if abs(q['f'] - z['f']) < 0.01)
    ax.set_title(f"{note}: {z['f']:.2f} Hz, partner {z['partner']:.2f} Hz, beat {z['beat']:.2f} Hz", loc='right', pad=2, fontsize=6.8)
    tag(ax, lab)
axs[len(Z) // 2].set_ylabel('level (dB)')
axs[-1].set_xlabel('frequency offset from the tracked partial (Hz)')
fig.savefig('fig_zoom.pdf')
plt.close(fig)
# numbers quoted in Section 9.3: decay rate and beat swing of each listed partial over its core interval
# (frames whose start lies in [coreFirst, coreLast]; series times are frame centres),
# and the swing implied by the partner level, 20 log10((1 + r) / (1 - r))
for q, l in zip(D['series'], D['listed']):
    t = np.array(q['t']); y = np.array(q['d'])
    ts = t - D['rec']['winSec'] / 2
    m = (ts >= l['coreFirst'] - 1e-6) & (ts <= l['coreLast'] + 1e-6)
    c = np.polyfit(t[m], y[m], 1); r = y[m] - np.polyval(c, t[m])
    pr = l['partnerRel']
    implied = 20 * np.log10((1 + 10 ** (pr / 20)) / (1 - 10 ** (pr / 20))) if pr is not None else float('nan')
    print(f"{l['note']:8s} decay {-c[0]:5.1f} dB/s  swing (5-95 %) {np.percentile(r, 95) - np.percentile(r, 5):4.1f} dB  implied by partner {implied:4.1f} dB")
