/*
 * Similarity collapse for the heat equation (posts/pde-symmetries.md).
 *
 *   <div class="collapse-fig"></div>
 *
 * Left: temperature profiles theta = erfc(x / 2 sqrt(alpha t)) in a slab of concrete at
 * several times. Right: the same data plotted against x / 2 sqrt(alpha t), where it all
 * lands on one curve.
 */
(function () {
  'use strict';

  const ALPHA = 1e-6;                 // m^2/s, roughly concrete
  const HOURS = [0.5, 2, 8, 24];      // the faint reference profiles
  const T_MIN = 0.25;
  const T_MAX = 24;
  const X_MAX = 30;                   // cm
  const ETA_MAX = 2.5;
  const COLORS = ['#b3b2d6', '#8f8cc0', '#6d6a9c', '#4b4878'];
  const NOW = '#2f78c4';
  const MONO = '"Courier Prime", "Courier New", monospace';
  const SERIF = '"Source Serif 4", Georgia, serif';
  const REDUCED_MOTION = window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches;

  const STYLE = `
    .collapse-fig canvas { display: block; width: 100%; border: 1px solid var(--line, #e4e3ec);
      border-radius: 4px; background: #fff; }
    .collapse-controls { display: flex; flex-wrap: wrap; align-items: center; gap: 0.5em 1.1em; margin-top: 0.6em;
      color: var(--muted, #6c6b80); font-family: var(--mono, monospace); font-size: 0.85rem; }
    .collapse-controls label { display: inline-flex; align-items: center; gap: 0.5em; }
    .collapse-controls input[type="range"] { width: 10em; accent-color: var(--theorem, #6d6a9c); }
    .collapse-controls button { padding: 0.3em 0.75em; border: 1px solid var(--line, #e4e3ec); border-radius: 4px;
      background: #fff; color: var(--ink, #333); font: inherit; cursor: pointer; }
    .collapse-legend { display: inline-flex; flex-wrap: wrap; gap: 0.2em 0.9em; }
    .collapse-legend span::before { content: ""; display: inline-block; width: 0.8em; height: 0.8em;
      margin-right: 0.35em; vertical-align: -0.05em; background: var(--c); }
    .collapse-readout { min-height: 1.5em; margin-top: 0.45em; color: var(--ink, #333);
      font-family: var(--mono, monospace); font-size: 0.85rem; }
    figure:has(> .collapse-fig) figcaption { text-align: left; }
  `;

  // erf, Abramowitz & Stegun 7.1.26 (error below 1.5e-7)
  function erf(z) {
    const s = Math.sign(z);
    const x = Math.abs(z);
    const t = 1 / (1 + 0.3275911 * x);
    const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x);
    return s * y;
  }
  const theta = (xCm, hours) => 1 - erf(eta(xCm, hours));
  const eta = (xCm, hours) => (xCm / 100) / (2 * Math.sqrt(ALPHA * hours * 3600));

  function el(tag, className, text) {
    const e = document.createElement(tag);
    if (className) e.className = className;
    if (text) e.textContent = text;
    return e;
  }

  function figure(root) {
    const canvas = el('canvas');
    canvas.setAttribute('role', 'img');
    canvas.setAttribute('aria-label', 'Temperature profiles at different times spread out with depth, but plotted against x over 2 root alpha t they all fall on one curve.');
    const controls = el('div', 'collapse-controls');
    const readout = el('div', 'collapse-readout');
    readout.setAttribute('aria-live', 'polite');

    let s = REDUCED_MOTION ? 1 : 0;   // position of the time slider, 0..1 on a log scale
    let playing = !REDUCED_MOTION;
    let pause = 0;
    const hoursAt = (p) => T_MIN * Math.pow(T_MAX / T_MIN, p);

    const lab = el('label', '', 'time ');
    const slider = el('input');
    slider.type = 'range';
    slider.min = '0';
    slider.max = '1';
    slider.step = '0.001';
    slider.value = String(s);
    slider.addEventListener('input', () => {
      s = +slider.value;
      playing = false;
      play.textContent = 'play';
      draw();
    });
    lab.appendChild(slider);
    const play = el('button', '', playing ? 'pause' : 'play');
    play.type = 'button';
    play.addEventListener('click', () => {
      playing = !playing;
      if (playing && s >= 1) s = 0;
      play.textContent = playing ? 'pause' : 'play';
    });
    const legend = el('span', 'collapse-legend');
    HOURS.forEach((h, i) => {
      const item = el('span', '', h + ' h');
      item.style.setProperty('--c', COLORS[i]);
      legend.appendChild(item);
    });
    const nowItem = el('span', '', 'now');
    nowItem.style.setProperty('--c', NOW);
    legend.appendChild(nowItem);
    controls.append(lab, play, legend);
    root.append(canvas, controls, readout);

    let ctx = null;
    let W = 0;
    let H = 0;
    let panels = [];

    function layout() {
      W = root.clientWidth;
      const wide = W >= 560;
      const ph = Math.round(Math.max(210, Math.min(300, (wide ? W / 2 : W) * 0.62)));
      H = wide ? ph : 2 * ph;
      const pw = wide ? W / 2 : W;
      panels = [0, 1].map((i) => {
        const ox = wide ? i * pw : 0;
        const oy = wide ? 0 : i * ph;
        return { x: ox + 52, y: oy + 30, w: pw - 52 - 18, h: ph - 30 - 46 };
      });
      const dpr = Math.min(window.devicePixelRatio || 1, 2);
      canvas.width = Math.round(W * dpr);
      canvas.height = Math.round(H * dpr);
      canvas.style.height = H + 'px';
      ctx = canvas.getContext('2d');
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      draw();
    }

    function axes(p, xMax, xTicks, xLabel, title) {
      ctx.strokeStyle = '#d6d5e2';
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(p.x, p.y);
      ctx.lineTo(p.x, p.y + p.h);
      ctx.lineTo(p.x + p.w, p.y + p.h);
      ctx.stroke();
      ctx.font = '11px ' + MONO;
      ctx.fillStyle = '#6c6b80';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'top';
      xTicks.forEach((v) => {
        const X = p.x + (v / xMax) * p.w;
        ctx.fillText(String(v), X, p.y + p.h + 5);
        ctx.beginPath();
        ctx.moveTo(X, p.y + p.h);
        ctx.lineTo(X, p.y + p.h + 3);
        ctx.stroke();
      });
      ctx.fillText(xLabel, p.x + p.w / 2, p.y + p.h + 22);
      ctx.textAlign = 'right';
      ctx.textBaseline = 'middle';
      [0, 0.5, 1].forEach((v) => ctx.fillText(String(v), p.x - 6, p.y + p.h - v * p.h));
      ctx.save();
      ctx.translate(p.x - 36, p.y + p.h / 2);
      ctx.rotate(-Math.PI / 2);
      ctx.textAlign = 'center';
      ctx.font = 'italic 14px ' + SERIF;
      ctx.fillText('θ', 0, 0);
      ctx.restore();
      ctx.textAlign = 'left';
      ctx.textBaseline = 'alphabetic';
      ctx.font = '12px ' + MONO;
      ctx.fillStyle = '#2a2a31';
      ctx.fillText(title, p.x, p.y - 12);
    }

    function curve(p, xMax, f, color, width, from, to) {
      ctx.beginPath();
      for (let i = 0; i <= 160; i++) {
        const v = from + (to - from) * (i / 160);
        const X = p.x + (v / xMax) * p.w;
        const Y = p.y + p.h - f(v) * p.h;
        if (i === 0) ctx.moveTo(X, Y);
        else ctx.lineTo(X, Y);
      }
      ctx.strokeStyle = color;
      ctx.lineWidth = width;
      ctx.stroke();
    }

    function marker(p, xMax, v, th, color) {
      ctx.beginPath();
      ctx.arc(p.x + (v / xMax) * p.w, p.y + p.h - th * p.h, 3.2, 0, 2 * Math.PI);
      ctx.fillStyle = color;
      ctx.fill();
    }

    function draw() {
      if (!ctx) return;
      ctx.clearRect(0, 0, W, H);
      const [L, R] = panels;
      const now = hoursAt(s);
      axes(L, X_MAX, [0, 10, 20, 30], 'depth x (cm)', 'profiles at different times');
      axes(R, ETA_MAX, [0, 0.5, 1, 1.5, 2, 2.5], 'x / 2√(αt)', 'the same data, rescaled');

      // The single similarity profile, under everything on the right.
      curve(R, ETA_MAX, (e) => 1 - erf(e), '#e4e3ec', 6, 0, ETA_MAX);

      const series = HOURS.map((h, i) => ({ h, color: COLORS[i], width: 1.5 })).concat([{ h: now, color: NOW, width: 2.5 }]);
      series.forEach(({ h, color, width }) => {
        curve(L, X_MAX, (xc) => theta(xc, h), color, width, 0, X_MAX);
        for (let xc = 0; xc <= X_MAX; xc += 2.5) {
          const th = theta(xc, h);
          marker(L, X_MAX, xc, th, color);
          const e = eta(xc, h);
          if (e <= ETA_MAX) marker(R, ETA_MAX, e, th, color);
        }
      });

      const half = 0.4769 * 2 * Math.sqrt(ALPHA * now * 3600) * 100; // depth where theta = 1/2, in cm
      readout.textContent = 't = ' + now.toFixed(1) + ' h: the temperature rise is half the surface’s at a depth of ' + half.toFixed(1) + ' cm';
    }

    let width = 0;
    const check = () => {
      if (root.clientWidth && root.clientWidth !== width) {
        width = root.clientWidth;
        layout();
      }
    };
    if ('ResizeObserver' in window) new ResizeObserver(check).observe(root);
    else window.addEventListener('resize', check);
    check();

    let prev = null;
    function frame(time) {
      if (!canvas.isConnected) return; // the post was re-rendered
      const dt = prev === null ? 0 : Math.min(0.05, (time - prev) / 1000);
      prev = time;
      if (playing) {
        if (s < 1) {
          s = Math.min(1, s + dt / 7);
        } else if ((pause += dt) > 1.5) {
          s = 0;
          pause = 0;
        }
        slider.value = String(s);
        draw();
      }
      requestAnimationFrame(frame);
    }
    requestAnimationFrame(frame);
  }

  if (!document.getElementById('collapse-fig-style')) {
    const style = el('style');
    style.id = 'collapse-fig-style';
    style.textContent = STYLE;
    document.head.appendChild(style);
  }
  document.querySelectorAll('.collapse-fig:not([data-ready])').forEach((root) => {
    root.setAttribute('data-ready', '');
    figure(root);
  });
})();
