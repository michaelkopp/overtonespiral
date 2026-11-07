/* Overtone Spiral — service worker: keeps the app usable without a network connection, for example when it
   has been added to the home screen of a phone. The app is one HTML page; it is kept in the cache together with
   its manifest, icons and web fonts. Pages are fetched from the network when it answers within a few seconds
   (so that updates arrive), otherwise taken from the cache.
   Authors: Claude (Anthropic) and Michael Kopp. Copyright (C) 2026 Michael Kopp. SPDX-License-Identifier: GPL-3.0-or-later */
'use strict';
const CACHE = 'overtone-spiral-v1';
const CORE = ['./', 'index.html', 'manifest.webmanifest', 'icons/apple-touch-icon.png', 'icons/icon-192.png', 'icons/icon-512.png',
  'icons/icon-maskable-192.png', 'icons/icon-maskable-512.png', 'icons/favicon-64.png'];
const FONT_CSS = 'https://fonts.googleapis.com/css2?family=Azeret+Mono:wght@400;500&family=Hanken+Grotesk:wght@400;500;600&family=Spectral:ital,wght@0,500;0,600;1,400;1,500&display=swap';
const FONT_HOSTS = ['fonts.googleapis.com', 'fonts.gstatic.com'];

// the web fonts: the stylesheet, then every font file it names (optional: the app falls back to system fonts)
async function cacheFonts(cache) {
  const res = await fetch(FONT_CSS, { mode: 'cors' });
  if (!res.ok) return;
  await cache.put(FONT_CSS, res.clone());
  const urls = [...new Set((await res.text()).match(/https:\/\/fonts\.gstatic\.com\/[^)'"\s]+/g) || [])];
  await Promise.all(urls.map(u => fetch(u, { mode: 'cors' }).then(r => (r.ok ? cache.put(u, r) : null)).catch(() => null)));
}

self.addEventListener('install', event => {
  event.waitUntil((async () => {
    const cache = await caches.open(CACHE);
    await cache.addAll(CORE.map(u => new Request(u, { cache: 'reload' })));
    try { await cacheFonts(cache); } catch (e) { /* offline or blocked: system fonts are used */ }
    await self.skipWaiting();
  })());
});

self.addEventListener('activate', event => {
  event.waitUntil((async () => {
    for (const key of await caches.keys()) if (key !== CACHE) await caches.delete(key);
    await self.clients.claim();
  })());
});

// pages: network first (at most 4 s), then the cached copy of the app
async function page(request) {
  const cache = await caches.open(CACHE);
  try {
    const res = await Promise.race([fetch(request), new Promise((_, no) => setTimeout(() => no(new Error('timeout')), 4000))]);
    if (res && res.ok) { cache.put('index.html', res.clone()); return res; }
    throw new Error('status ' + (res && res.status));
  } catch (e) {
    return (await cache.match(request, { ignoreSearch: true })) || (await cache.match('index.html')) || (await cache.match('./')) || Response.error();
  }
}
// everything else the app uses (its icons, manifest, fonts): the cached copy at once, refreshed in the background
async function asset(request, event) {
  const cache = await caches.open(CACHE);
  const hit = await cache.match(request, { ignoreSearch: new URL(request.url).origin === location.origin });
  const update = fetch(request).then(res => { if (res && (res.ok || res.type === 'opaque')) cache.put(request, res.clone()); return res; }).catch(() => null);
  if (hit) { event.waitUntil(update); return hit; }
  return (await update) || Response.error();
}

self.addEventListener('fetch', event => {
  const req = event.request;
  if (req.method !== 'GET') return;
  const url = new URL(req.url);
  if (req.mode === 'navigate' && url.origin === location.origin) { event.respondWith(page(req)); return; }
  const mine = url.origin === location.origin && CORE.some(c => url.pathname.endsWith(c.replace('./', '/')) || (c === './' && url.pathname.endsWith('/')));
  if (mine || FONT_HOSTS.includes(url.hostname)) event.respondWith(asset(req, event));
});
