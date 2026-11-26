/*
 * Interactive figures for "Chow's Theorem, Physics, and Control".
 *
 *   <div class="chow-fig" data-fig="loop"></div>  the loop X, Y, -X, -Y lifted to a distribution in R^3
 *   <div class="chow-fig" data-fig="gas"></div>   two gases in thermal contact run around a piston cycle
 *
 * Plain canvas, no libraries. Drag a 3D view to rotate it.
 */
(function () {
  'use strict';

  const COLOR = {
    ink: '#2a2a31',
    muted: '#6c6b80',
    faint: '#c9c8d6',
    grid: '#eeedf3',
    plane: 'rgba(109, 106, 156, 0.13)',
    planeEdge: 'rgba(109, 106, 156, 0.42)',
    planeNow: 'rgba(47, 120, 196, 0.20)',
    planeNowEdge: 'rgba(47, 120, 196, 0.9)',
    path: '#2f78c4',
    shadow: 'rgba(108, 107, 128, 0.6)',
    area: 'rgba(47, 120, 196, 0.10)',
    gap: '#c2553f',
    leaf: 'rgba(127, 148, 85, 0.55)',
  };
  const MONO = '"Courier Prime", "Courier New", monospace';
  const SERIF = '"Source Serif 4", Georgia, serif';
  const REDUCED_MOTION = window.matchMedia && window.matchMedia('(prefers-reduced-motion: reduce)').matches;

  const STYLE = `
    .chow-fig canvas { display: block; width: 100%; border: 1px solid var(--line, #e4e3ec); border-radius: 4px;
      background: #fff; cursor: grab; touch-action: pan-y; }
    .chow-fig canvas:active { cursor: grabbing; }
    .chow-controls { display: flex; flex-wrap: wrap; align-items: center; gap: 0.5em 1.1em; margin-top: 0.6em;
      color: var(--muted, #6c6b80); font-family: var(--mono, monospace); font-size: 0.85rem; }
    .chow-seg { display: inline-flex; overflow: hidden; border: 1px solid var(--line, #e4e3ec); border-radius: 4px; }
    .chow-seg button, .chow-btn { padding: 0.3em 0.75em; border: 0; background: #fff; color: var(--ink, #333);
      font: inherit; cursor: pointer; }
    .chow-seg button + button { border-left: 1px solid var(--line, #e4e3ec); }
    .chow-seg button[aria-pressed="true"] { background: var(--theorem, #6d6a9c); color: #fff; }
    .chow-btn { border: 1px solid var(--line, #e4e3ec); border-radius: 4px; }
    .chow-btn:hover, .chow-seg button:hover { color: var(--link, #2f78c4); }
    .chow-seg button[aria-pressed="true"]:hover { color: #fff; }
    .chow-controls label { display: inline-flex; align-items: center; gap: 0.5em; }
    .chow-controls input[type="range"] { width: 8.5em; accent-color: var(--theorem, #6d6a9c); }
    figure:has(> .chow-fig) figcaption { text-align: left; }
    .chow-readout { min-height: 1.5em; margin-top: 0.45em; color: var(--ink, #333);
      font-family: var(--mono, monospace); font-size: 0.85rem; }
  `;

  // ---------------------------------------------------------------------------
  // Helpers
  // ---------------------------------------------------------------------------

  const clamp = (x, a, b) => Math.min(b, Math.max(a, x));
  const lerp = (a, b, f) => a + (b - a) * f;

  function el(tag, className, text) {
    const e = document.createElement(tag);
    if (className) e.className = className;
    if (text) e.textContent = text;
    return e;
  }

  function segmented(labels, onPick) {
    const wrap = el('span', 'chow-seg');
    const buttons = labels.map((label, i) => {
      const b = el('button', '', label);
      b.type = 'button';
      b.setAttribute('aria-pressed', String(i === 0));
      b.addEventListener('click', () => {
        buttons.forEach((x, j) => x.setAttribute('aria-pressed', String(j === i)));
        onPick(i);
      });
      wrap.appendChild(b);
      return b;
    });
    return wrap;
  }

  function button(label, onClick) {
    const b = el('button', 'chow-btn', label);
    b.type = 'button';
    b.addEventListener('click', onClick);
    return b;
  }

  function setupCanvas(canvas, w, h) {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    canvas.width = Math.round(w * dpr);
    canvas.height = Math.round(h * dpr);
    canvas.style.height = h + 'px';
    const ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return ctx;
  }

  /** Orthographic camera: turn about the vertical axis (yaw), then look down (pitch). */
  function project(view, p) {
    const cy = Math.cos(view.yaw), sy = Math.sin(view.yaw);
    const x1 = p[0] * cy - p[1] * sy;
    const y1 = p[0] * sy + p[1] * cy;
    const cp = Math.cos(view.pitch), sp = Math.sin(view.pitch);
    return [view.cx + view.s * x1, view.cy - view.s * (p[2] * cp + y1 * sp), y1 * cp - p[2] * sp];
  }

  function enableDrag(canvas, view, redraw) {
    let last = null;
    canvas.addEventListener('pointerdown', (e) => {
      last = [e.clientX, e.clientY];
      canvas.setPointerCapture(e.pointerId);
    });
    canvas.addEventListener('pointermove', (e) => {
      if (!last) return;
      view.yaw += (e.clientX - last[0]) * 0.01;
      view.pitch = clamp(view.pitch + (e.clientY - last[1]) * 0.006, 0.08, 1.3);
      last = [e.clientX, e.clientY];
      redraw();
    });
    const stop = () => { last = null; };
    canvas.addEventListener('pointerup', stop);
    canvas.addEventListener('pointercancel', stop);
  }

  /** Calls step(dt) every frame while the canvas is on the page and on screen. */
  function runLoop(canvas, step) {
    let visible = true;
    if ('IntersectionObserver' in window) {
      new IntersectionObserver((entries) => { visible = entries[0].isIntersecting; }).observe(canvas);
    }
    let prev = null;
    function frame(now) {
      if (!canvas.isConnected) return; // the post was re-rendered
      const dt = prev === null ? 0 : Math.min(0.05, (now - prev) / 1000);
      prev = now;
      if (visible) step(dt);
      requestAnimationFrame(frame);
    }
    requestAnimationFrame(frame);
  }

  function onResize(root, fn) {
    let width = 0;
    const check = () => {
      if (root.clientWidth && root.clientWidth !== width) {
        width = root.clientWidth;
        fn();
      }
    };
    if ('ResizeObserver' in window) new ResizeObserver(check).observe(root);
    else window.addEventListener('resize', check);
    check();
  }

  function stroke(ctx, pts, color, width, dash) {
    if (pts.length < 2) return;
    ctx.beginPath();
    ctx.moveTo(pts[0][0], pts[0][1]);
    for (let i = 1; i < pts.length; i++) ctx.lineTo(pts[i][0], pts[i][1]);
    ctx.strokeStyle = color;
    ctx.lineWidth = width || 1;
    ctx.setLineDash(dash || []);
    ctx.lineJoin = 'round';
    ctx.lineCap = 'round';
    ctx.stroke();
    ctx.setLineDash([]);
  }

  function fill(ctx, pts, color, edge) {
    ctx.beginPath();
    ctx.moveTo(pts[0][0], pts[0][1]);
    for (let i = 1; i < pts.length; i++) ctx.lineTo(pts[i][0], pts[i][1]);
    ctx.closePath();
    ctx.fillStyle = color;
    ctx.fill();
    if (edge) {
      ctx.strokeStyle = edge;
      ctx.lineWidth = 1;
      ctx.stroke();
    }
  }

  function dot(ctx, p, r, color, hollow) {
    ctx.beginPath();
    ctx.arc(p[0], p[1], r, 0, 2 * Math.PI);
    if (hollow) {
      ctx.fillStyle = '#fff';
      ctx.fill();
      ctx.strokeStyle = color;
      ctx.lineWidth = 2;
      ctx.stroke();
    } else {
      ctx.fillStyle = color;
      ctx.fill();
    }
  }

  function label(ctx, p, text, opts) {
    const o = opts || {};
    ctx.font = o.font || '13px ' + MONO;
    ctx.fillStyle = o.color || COLOR.muted;
    ctx.textAlign = o.align || 'center';
    ctx.textBaseline = o.baseline || 'middle';
    ctx.fillText(text, p[0] + (o.dx || 0), p[1] + (o.dy || 0));
  }

  /** A small parallelogram: the plane spanned by u and v at p (all in world coordinates). */
  function patch(P, p, u, v, size) {
    const c = (a, b) => P([p[0] + size * (a * u[0] + b * v[0]), p[1] + size * (a * u[1] + b * v[1]), p[2] + size * (a * u[2] + b * v[2])]);
    return [c(-1, -1), c(1, -1), c(1, 1), c(-1, 1)];
  }

  function floorGrid(ctx, P, x0, x1, y0, y1, z, step) {
    for (let x = x0; x <= x1 + 1e-9; x += step) stroke(ctx, [P([x, y0, z]), P([x, y1, z])], COLOR.grid, 1);
    for (let y = y0; y <= y1 + 1e-9; y += step) stroke(ctx, [P([x0, y, z]), P([x1, y, z])], COLOR.grid, 1);
    stroke(ctx, [P([x0, y0, z]), P([x1, y0, z]), P([x1, y1, z]), P([x0, y1, z]), P([x0, y0, z])], COLOR.faint, 1);
  }

  // ---------------------------------------------------------------------------
  // Figure 1: the loop X, Y, -X, -Y lifted to D
  // ---------------------------------------------------------------------------

  const LEG_SECONDS = 1.1;
  const LEG_NAMES = ['flowing along X', 'flowing along Y', 'flowing along −X', 'flowing along −Y'];

  function loopFigure(root) {
    const canvas = el('canvas');
    canvas.setAttribute('role', 'img');
    canvas.setAttribute('aria-label', 'A curve tangent to the distribution goes around a square and ends up above its starting point.');
    const readout = el('div', 'chow-readout');
    readout.setAttribute('aria-live', 'polite');
    const controls = el('div', 'chow-controls');

    let twisted = true;   // [X,Y] = dz, or the integrable comparison
    let t = 1.0;          // side of the square
    let sigma = 0;        // progress: 0..4 legs, then a pause
    let playing = !REDUCED_MOTION;
    if (!playing) sigma = 4;

    controls.appendChild(segmented(['[X,Y] = ∂z', '[X,Y] = 0'], (i) => { twisted = i === 0; restart(); }));
    const sizeLabel = el('label', '', 'loop size t ');
    const slider = el('input');
    slider.type = 'range';
    slider.min = '0.4';
    slider.max = '1.4';
    slider.step = '0.01';
    slider.value = String(t);
    slider.addEventListener('input', () => { t = +slider.value; restart(); });
    sizeLabel.appendChild(slider);
    controls.appendChild(sizeLabel);
    const playButton = button(playing ? 'pause' : 'play', () => {
      playing = !playing;
      if (playing && sigma >= 4) sigma = 0;
      playButton.textContent = playing ? 'pause' : 'play';
      draw();
    });
    controls.appendChild(playButton);
    root.append(canvas, controls, readout);

    function restart() {
      sigma = playing ? 0 : 4;
      draw();
    }

    // Distribution: X = dx -/+ (y/2) dz, Y = dy + (x/2) dz
    const X = (p) => [1, 0, twisted ? -p[1] / 2 : p[1] / 2];
    const Y = (p) => [0, 1, p[0] / 2];

    // The lifted square, starting at (-h, -h, 0).
    function pathAt(s) {
      const h = t / 2;
      const corners = [[-h, -h], [h, -h], [h, h], [-h, h], [-h, -h]];
      const k = Math.min(3, Math.floor(s));
      const f = clamp(s - k, 0, 1);
      const x = lerp(corners[k][0], corners[k + 1][0], f);
      const y = lerp(corners[k][1], corners[k + 1][1], f);
      // twisted: dz = (x dy - y dx)/2 climbs t^2/4 on every leg; flat: z = xy/2 + const
      const z = twisted ? (t * t / 4) * clamp(s, 0, 4) : x * y / 2 - h * h / 2;
      return [x, y, z];
    }

    const view = { yaw: -0.62, pitch: 0.42, s: 1, cx: 0, cy: 0 };
    let ctx = null;
    let W = 0;
    let H = 0;
    const ZC = 0.5; // vertical centre of the scene
    const FLOOR = -1.0;
    const P = (p) => project(view, [p[0], p[1], p[2] - ZC]);

    function layout() {
      W = root.clientWidth;
      H = Math.round(clamp(W * 0.66, 260, 440));
      ctx = setupCanvas(canvas, W, H);
      view.s = Math.min(W * 0.2, H * 0.205);
      view.cx = W / 2;
      view.cy = H / 2 + 6;
      draw();
    }

    function draw() {
      if (!ctx) return;
      ctx.clearRect(0, 0, W, H);
      const s = Math.min(sigma, 4);
      const done = sigma >= 4;
      const h = t / 2;

      floorGrid(ctx, P, -1.5, 1.5, -1.5, 1.5, FLOOR, 0.5);
      label(ctx, P([0, -1.85, FLOOR]), 'x', { font: 'italic 15px ' + SERIF });
      label(ctx, P([1.85, 0, FLOOR]), 'y', { font: 'italic 15px ' + SERIF });
      stroke(ctx, [P([-1.5, -1.5, FLOOR]), P([-1.5, -1.5, 2.1])], COLOR.faint, 1);
      label(ctx, P([-1.5, -1.5, 2.25]), 'z', { font: 'italic 15px ' + SERIF });

      // Trail and its shadow on the floor.
      const trail = [];
      for (let u = 0; u <= s + 1e-9; u += 0.02) trail.push(pathAt(Math.min(u, s)));
      trail.push(pathAt(s));
      const square = [[-h, -h], [h, -h], [h, h], [-h, h]].map((q) => P([q[0], q[1], FLOOR]));
      if (done) fill(ctx, square, COLOR.area);
      stroke(ctx, trail.map((q) => P([q[0], q[1], FLOOR])), COLOR.shadow, 1.5);

      // The surface z = xy/2 + c that the flat distribution's curves stay on.
      if (!twisted) {
        const c = -h * h / 2;
        for (let a = -1.4; a <= 1.4 + 1e-9; a += 0.35) {
          const l1 = [];
          const l2 = [];
          for (let b = -1.4; b <= 1.4 + 1e-9; b += 0.1) {
            l1.push(P([a, b, a * b / 2 + c]));
            l2.push(P([b, a, a * b / 2 + c]));
          }
          stroke(ctx, l1, COLOR.leaf, 1);
          stroke(ctx, l2, COLOR.leaf, 1);
        }
      }

      // The planes D_p on a grid at height 0, back to front.
      const planes = [];
      for (let i = -2; i <= 2; i++) {
        for (let j = -2; j <= 2; j++) {
          const p = [i * 0.6, j * 0.6, 0];
          planes.push({ depth: P(p)[2], pts: patch(P, p, X(p), Y(p), 0.15) });
        }
      }
      planes.sort((a, b) => b.depth - a.depth).forEach((q) => fill(ctx, q.pts, COLOR.plane, COLOR.planeEdge));

      const start = pathAt(0);
      const now = pathAt(s);
      stroke(ctx, [P(start), P([start[0], start[1], FLOOR])], COLOR.faint, 1, [4, 4]);
      stroke(ctx, [P(now), P([now[0], now[1], FLOOR])], COLOR.faint, 1, [4, 4]);
      stroke(ctx, trail.map(P), COLOR.path, 2.5);

      // The plane the curve is allowed to move in right now.
      fill(ctx, patch(P, now, X(now), Y(now), 0.28), COLOR.planeNow, COLOR.planeNowEdge);

      dot(ctx, P(start), 5, COLOR.ink, true);
      if (done && twisted) {
        stroke(ctx, [P(start), P(now)], COLOR.gap, 2.5);
        const mid = P([start[0], start[1], now[2] / 2]);
        label(ctx, mid, 'Δz = t²', { color: COLOR.gap, align: 'left', dx: 10, font: 'bold 13px ' + MONO });
      }
      dot(ctx, P(now), 5.5, COLOR.path);

      label(ctx, [14, 18], done ? (twisted ? 'back above the start' : 'back at the start') : LEG_NAMES[Math.floor(s)], { align: 'left' });

      const dz = now[2] - start[2];
      if (done) {
        readout.textContent = twisted
          ? 'The loop climbed Δz = ' + dz.toFixed(2) + ', the area of the ' + t.toFixed(2) + ' × ' + t.toFixed(2) + ' square.'
          : 'The loop closed: Δz = 0. Curves tangent to this distribution stay on one surface.';
      } else {
        readout.textContent = 'Δz so far: ' + (Math.abs(dz) < 0.005 ? '0.00' : dz.toFixed(2));
      }
    }

    enableDrag(canvas, view, draw);
    onResize(root, layout);
    runLoop(canvas, (dt) => {
      if (!playing) return;
      sigma += dt / LEG_SECONDS;
      if (sigma > 5.6) sigma = 0;
      draw();
    });
  }

  // ---------------------------------------------------------------------------
  // Figure 2: two gases in thermal contact, around a piston cycle
  // ---------------------------------------------------------------------------

  const T0 = 300;           // K
  const A = Math.log(2);    // each piston moves between V and 2V
  const KAPPA = -0.5;       // twist of the hypothetical heat form
  const CYCLES = 3;         // cycles shown for the hypothetical gas before starting over
  const GAS_LEG_SECONDS = 1.0;

  function gasFigure(root) {
    const canvas = el('canvas');
    canvas.setAttribute('role', 'img');
    canvas.setAttribute('aria-label', 'Two cylinders of gas run through a piston cycle, next to the path the state takes in (V1, V2, T) space.');
    const readout = el('div', 'chow-readout');
    readout.setAttribute('aria-live', 'polite');
    const controls = el('div', 'chow-controls');

    let real = true;
    let sigma = 0;        // progress in legs; four legs per cycle
    let playing = !REDUCED_MOTION;
    let trail = [];
    if (!playing) sigma = 4;

    controls.appendChild(segmented(['real gas', 'non-integrable δQ'], (i) => { real = i === 0; restart(); }));
    const playButton = button(playing ? 'pause' : 'play', () => {
      playing = !playing;
      playButton.textContent = playing ? 'pause' : 'play';
      draw();
    });
    controls.appendChild(playButton);
    controls.appendChild(button('restart', () => restart()));
    root.append(canvas, controls, readout);

    function cyclesShown() { return real ? 1 : CYCLES; }

    function restart() {
      sigma = playing ? 0 : 4 * cyclesShown();
      trail = [];
      if (!playing) fullTrail();
      draw();
    }

    // State after s legs, in log coordinates u = ln(V1/V0), v = ln(V2/V0), w = ln(T/T0).
    // Real gas: dw = -(du + dv)/3, so T^3 V1 V2 is constant.
    // Hypothetical gas: dw = -(du + dv)/3 + KAPPA (u dv - v du)/2, which is not integrable.
    function stateAt(s) {
      const n = Math.floor(s / 4);
      const r = s - 4 * n;
      const k = Math.min(3, Math.floor(r));
      const f = clamp(r - k, 0, 1);
      let u = 0, v = 0, twist = 0;
      if (k === 0) { u = f * A; }
      else if (k === 1) { u = A; v = f * A; twist = f * A * A / 2; }
      else if (k === 2) { u = (1 - f) * A; v = A; twist = A * A / 2 + f * A * A / 2; }
      else { v = (1 - f) * A; twist = A * A; }
      let w = -(u + v) / 3;
      if (!real) w += KAPPA * (twist + n * A * A);
      return { V1: Math.exp(u), V2: Math.exp(v), T: T0 * Math.exp(w), u, v, leg: k, cycle: n };
    }

    // Slopes of the plane ker(dQ) at a state, as dT/dV1 and dT/dV2.
    function slopes(st) {
      const a = real ? -1 / 3 : -1 / 3 - KAPPA * st.v / 2;
      const b = real ? -1 / 3 : -1 / 3 + KAPPA * st.u / 2;
      return [st.T * a / st.V1, st.T * b / st.V2];
    }

    const view = { yaw: -0.7, pitch: 0.36, s: 1, cx: 0, cy: 0 };
    let ctx = null;
    let W = 0;
    let H = 0;
    let wide = true;
    let panel = { x: 0, y: 0, w: 0, h: 0 };
    // V in [1, 2] -> [-1, 1]; T in [60, 300] K -> [-1, 1]
    const T_FLOOR = 60;
    const world = (V1, V2, T) => [(V1 - 1.5) * 2, (V2 - 1.5) * 2, (T - 180) / 120];
    const P = (p) => project(view, p);
    const PS = (V1, V2, T) => P(world(V1, V2, T));

    function layout() {
      W = root.clientWidth;
      wide = W >= 560;
      if (wide) {
        H = Math.round(clamp(W * 0.5, 300, 380));
        panel = { x: 0, y: 0, w: Math.round(W * 0.36), h: H };
        view.s = Math.min((W - panel.w) * 0.19, H * 0.27);
        view.cx = panel.w + (W - panel.w) / 2;
        view.cy = H / 2 + 8;
      } else {
        const plotH = Math.round(W * 0.8);
        H = 215 + plotH;
        panel = { x: 0, y: 0, w: W, h: 215 };
        view.s = Math.min(W * 0.19, plotH * 0.27);
        view.cx = W / 2;
        view.cy = 215 + plotH / 2 + 4;
      }
      ctx = setupCanvas(canvas, W, H);
      draw();
    }

    function gasColor(T, alpha) {
      const f = clamp((T - 100) / 200, 0, 1);
      const r = Math.round(lerp(95, 222, f));
      const g = Math.round(lerp(140, 112, f));
      const b = Math.round(lerp(222, 90, f));
      return 'rgba(' + r + ',' + g + ',' + b + ',' + alpha + ')';
    }

    function drawPistons(st) {
      const cw = Math.min(64, panel.w * 0.24);
      const ch = Math.min(150, panel.h * 0.52);
      const cx = panel.x + panel.w / 2;
      const top = panel.y + (wide ? panel.h * 0.22 : 34);
      const bottom = top + ch;
      const left = cx - cw;
      const right = cx + cw;
      const ins = 7; // insulation thickness

      // insulation
      ctx.fillStyle = '#f1f0f6';
      ctx.fillRect(left - ins, top, 2 * cw + 2 * ins, ch + ins);
      ctx.strokeStyle = COLOR.faint;
      ctx.lineWidth = 1;
      ctx.strokeRect(left - ins, top, 2 * cw + 2 * ins, ch + ins);

      [st.V1, st.V2].forEach((V, i) => {
        const x = i === 0 ? left : cx;
        const level = bottom - ch * (V / 2.25);
        ctx.fillStyle = '#fff';
        ctx.fillRect(x, top, cw, ch);
        ctx.fillStyle = gasColor(st.T, 0.6);
        ctx.fillRect(x, level, cw, bottom - level);
        // piston + rod
        ctx.fillStyle = '#8a899c';
        ctx.fillRect(x + 2, level - 7, cw - 4, 7);
        ctx.fillRect(x + cw / 2 - 2, top - 16, 4, level - 7 - (top - 16));
        label(ctx, [x + cw / 2, bottom + ins + 12], 'gas ' + (i + 1), { font: '12px ' + MONO });
      });
      // the wall between the gases conducts heat
      stroke(ctx, [[cx, top], [cx, bottom]], gasColor(st.T, 1), 2, [3, 3]);

      const tText = 'T = ' + st.T.toFixed(1) + ' K';
      label(ctx, [cx, bottom + ins + 34], tText, { font: 'bold 14px ' + MONO, color: COLOR.ink });
      if (!real) label(ctx, [cx, bottom + ins + 52], 'cycle ' + Math.min(st.cycle + 1, CYCLES), { font: '12px ' + MONO });
    }

    function drawPlot(st) {
      const FLOOR = -1;
      floorGrid(ctx, P, -1, 1, -1, 1, FLOOR, 0.5);
      label(ctx, P([0, -1.3, FLOOR]), 'V₁', { font: 'italic 15px ' + SERIF });
      label(ctx, P([1.3, 0, FLOOR]), 'V₂', { font: 'italic 15px ' + SERIF });
      stroke(ctx, [P([-1, -1, FLOOR]), P([-1, -1, 1.15])], COLOR.faint, 1);
      label(ctx, P([-1, -1, 1.3]), 'T', { font: 'italic 15px ' + SERIF });
      [100, 200, 300].forEach((T) => {
        const q = PS(1, 1, T);
        stroke(ctx, [q, [q[0] - 5, q[1]]], COLOR.faint, 1);
        label(ctx, q, T + ' K', { align: 'right', dx: -8, font: '11px ' + MONO });
      });

      if (real) {
        // The adiabat through the starting state: T^3 V1 V2 = T0^3.
        for (let a = 1; a <= 2 + 1e-9; a += 0.25) {
          const l1 = [];
          const l2 = [];
          for (let b = 1; b <= 2 + 1e-9; b += 0.05) {
            l1.push(PS(a, b, T0 * Math.pow(a * b, -1 / 3)));
            l2.push(PS(b, a, T0 * Math.pow(a * b, -1 / 3)));
          }
          stroke(ctx, l1, COLOR.leaf, 1);
          stroke(ctx, l2, COLOR.leaf, 1);
        }
      }

      const start = PS(1, 1, T0);
      stroke(ctx, [start, PS(1, 1, T_FLOOR)], COLOR.faint, 1, [4, 4]);
      const cur = PS(st.V1, st.V2, st.T);
      stroke(ctx, [cur, PS(st.V1, st.V2, T_FLOOR)], COLOR.faint, 1, [4, 4]);
      stroke(ctx, trail.map((q) => PS(q[0], q[1], T_FLOOR)), COLOR.shadow, 1.5);
      stroke(ctx, trail.map((q) => PS(q[0], q[1], q[2])), COLOR.path, 2.5);

      // The plane ker(dQ) at the current state.
      const [a, b] = slopes(st);
      const d = 0.16;
      const corner = (i, j) => PS(st.V1 + i * d, st.V2 + j * d, st.T + a * i * d + b * j * d);
      fill(ctx, [corner(-1, -1), corner(1, -1), corner(1, 1), corner(-1, 1)], COLOR.planeNow, COLOR.planeNowEdge);

      dot(ctx, start, 5, COLOR.ink, true);
      dot(ctx, cur, 5.5, COLOR.path);
    }

    function draw() {
      if (!ctx) return;
      ctx.clearRect(0, 0, W, H);
      const st = stateAt(sigma);
      if (!wide) stroke(ctx, [[16, panel.h], [W - 16, panel.h]], COLOR.grid, 1);
      drawPistons(st);
      drawPlot(st);
      const legNames = ['expanding gas 1', 'expanding gas 2', 'compressing gas 1', 'compressing gas 2'];
      const total = 4 * cyclesShown();
      const label0 = sigma >= total ? 'cycle finished' : legNames[st.leg];
      label(ctx, [wide ? panel.w + 14 : 14, wide ? 18 : panel.h + 18], label0, { align: 'left' });

      const finished = Math.min(Math.floor(sigma / 4 + 1e-9), cyclesShown());
      if (real) {
        readout.textContent = finished >= 1
          ? 'After the cycle: T = ' + stateAt(4).T.toFixed(1) + ' K, back where it started.'
          : 'T = ' + st.T.toFixed(1) + ' K';
      } else {
        readout.textContent = finished >= 1
          ? 'After ' + finished + (finished === 1 ? ' cycle' : ' cycles') + ': T = ' + stateAt(4 * finished).T.toFixed(1) +
            ' K. Same volumes, no heat exchanged, but colder.'
          : 'T = ' + st.T.toFixed(1) + ' K';
      }
    }

    function record() {
      const st = stateAt(sigma);
      trail.push([st.V1, st.V2, st.T]);
    }

    function fullTrail() {
      trail = [];
      for (let s = 0; s <= 4 * cyclesShown() + 1e-9; s += 0.02) {
        const st = stateAt(s);
        trail.push([st.V1, st.V2, st.T]);
      }
    }

    enableDrag(canvas, view, draw);
    onResize(root, () => {
      if (!playing && !trail.length) fullTrail();
      layout();
    });
    runLoop(canvas, (dt) => {
      if (!playing) return;
      const total = 4 * cyclesShown();
      if (sigma < total) {
        sigma = Math.min(total, sigma + dt / GAS_LEG_SECONDS);
        record();
      } else {
        sigma += dt / GAS_LEG_SECONDS;
        if (sigma > total + (real ? 1.2 : 2.4)) {
          sigma = 0;
          trail = [];
        }
      }
      draw();
    });
  }

  // ---------------------------------------------------------------------------

  function init() {
    if (!document.getElementById('chow-fig-style')) {
      const style = el('style');
      style.id = 'chow-fig-style';
      style.textContent = STYLE;
      document.head.appendChild(style);
    }
    document.querySelectorAll('.chow-fig:not([data-ready])').forEach((root) => {
      root.setAttribute('data-ready', '');
      if (root.dataset.fig === 'gas') gasFigure(root);
      else loopFigure(root);
    });
  }

  init();
})();
