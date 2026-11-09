/* Overtone Spiral — service worker: keeps the app usable without a network connection, for example when it
   has been added to the home screen of a phone. The app is one HTML page (its fonts included); it is kept in
   the cache together with its manifest and icons, and everything is answered from the cache without touching
   the network (a phone in flight mode then has no reason to warn about the missing connection). A new version of the app
   comes with a new version of this file, which the browser looks for by itself when the app is opened online.
   Authors: Claude (Anthropic) and Michael Kopp. Copyright (C) 2026 Michael Kopp. SPDX-License-Identifier: GPL-3.0-or-later */
'use strict';
const VERSION = 'd15f8b6807c2';          // the build id, set by the build from the sources of everything kept here
const CACHE = 'overtone-spiral-' + VERSION;
const CORE = ['./', 'index.html', 'manifest.webmanifest', 'icons/apple-touch-icon.png', 'icons/icon-192.png', 'icons/icon-512.png',
  'icons/icon-maskable-192.png', 'icons/icon-maskable-512.png', 'icons/favicon-64.png'];
const SCOPE = new URL(self.registration.scope).pathname;
// notes for the app's event log (Settings → About this copy)
const say = msg => self.clients.matchAll({ includeUncontrolled: true, type: 'window' }).then(cs => cs.forEach(c => c.postMessage({ osLog: msg }))).catch(() => {});
const short = u => { const x = new URL(u); return (x.origin === location.origin ? '' : x.host) + x.pathname.slice(0, 60); };

self.addEventListener('install', event => {
  event.waitUntil((async () => {
    say(`installing build ${VERSION.slice(0, 6)}`);
    const cache = await caches.open(CACHE);
    await cache.addAll(CORE.map(u => new Request(u, { cache: 'reload' })));
    say(`build ${VERSION.slice(0, 6)} installed, ${(await cache.keys()).length} files kept`);
    await self.skipWaiting();
  })());
});

self.addEventListener('activate', event => {
  event.waitUntil((async () => {
    for (const key of await caches.keys()) if (key.startsWith('overtone-spiral-') && key !== CACHE) await caches.delete(key);
    await self.clients.claim();
  })());
});

// From the cache. The network only for a file of the app that is not there (the first visit, before the cache
// is filled), and not while the device reports that it is offline.
async function cached(request, key) {
  const cache = await caches.open(CACHE);
  const url = new URL(request.url), same = url.origin === location.origin;
  const hit = await cache.match(key || request, { ignoreSearch: same, ignoreVary: true });
  if (hit) return hit;
  if (!same || self.navigator.onLine === false) {
    say(`not in the offline copy, not fetched: ${short(request.url)}`);
    return request.destination === 'style' ? new Response('', { headers: { 'Content-Type': 'text/css' } }) : Response.error();
  }
  say(`network request: ${short(request.url)}`);
  try {
    const res = await fetch(request);
    if (res && res.ok && !key) cache.put(request, res.clone());
    return res;
  } catch (e) { say(`network request failed: ${short(request.url)}`); throw e; }
}

self.addEventListener('fetch', event => {
  const req = event.request;
  if (req.method !== 'GET') return;
  const url = new URL(req.url), same = url.origin === location.origin;
  // the app itself (other pages on the site, such as the paper, are left to the network)
  if (req.mode === 'navigate') {
    if (same && (url.pathname === SCOPE || url.pathname === SCOPE + 'index.html')) event.respondWith(cached(req, 'index.html'));
    return;
  }
  const mine = same && CORE.some(c => c !== './' && url.pathname === SCOPE + c);
  if (mine) event.respondWith(cached(req));
});
