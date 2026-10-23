/* Biothing Sound Lab growth: free forms grown by natural rules, tested with a wave simulator.
 *
 * A form is a list of capsules {a, b, r} (a sphere when a == b) in metres, z up.
 * Hanging forms hang from z = 0 downwards; standing forms stand on z = 0.
 *
 * Evaluation: the form is voxelised, projected onto two vertical planes (homogenised along
 * the projection axis), and each projection is placed in a 2D room section simulated with a
 * finite-difference time-domain (FDTD) wave solver. The porous material damps air motion in
 * proportion to its flow resistivity. The decay of the sound energy gives a 2D absorption
 * length, which times the form's extent gives an estimate of its absorption in m2 (Sabine).
 * Good for comparing forms; not a substitute for measurement. */
(function (root) {
  "use strict";

  const C0 = 343.2, RHO0 = 1.204;
  const BANDS = [250, 500, 1000];
  const CELLS_PER_WAVE = 6; // grid cells per wavelength at the highest band
  const lerp = (a, b, u) => a + (b - a) * u;
  const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const mul = (a, k) => [a[0] * k, a[1] * k, a[2] * k];
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const norm = (a) => { const l = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0] / l, a[1] / l, a[2] / l]; };
  const perp = (d) => norm(Math.abs(d[2]) < 0.9 ? cross(d, [0, 0, 1]) : cross(d, [1, 0, 0]));

  function rng(seed) {
    let a = seed >>> 0;
    return () => {
      a = (a + 0x6d2b79f5) >>> 0;
      let t = Math.imul(a ^ (a >>> 15), a | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }
  const seedOf = (g) => Math.floor(g.reduce((acc, x, i) => acc + x * 1e4 * (i + 1), 7)) >>> 0;

  // ------------------------------------------------------------ growth rules
  // Each rule turns genes (numbers in [0, 1]) into capsules at an arbitrary scale.
  // `down` is true for hanging forms (growth points downwards).
  const RULES = {
    tree: {
      name: "Tree and coral", imitates: "trees, branching coral, lungs",
      idea: "Each branch splits into smaller branches: lots of surface and edges for very little material, at many sizes at once.",
      genes: ["branches", "spread", "length", "taper", "levels", "twist", "thickness", "bend"],
      build(g, down) {
        const n = 2 + Math.floor(g[0] * 2.999), spread = lerp(15, 80, g[1]) * Math.PI / 180, lr = lerp(0.5, 0.85, g[2]);
        const rr = lerp(0.45, 0.8, g[3]), twist = lerp(0, Math.PI, g[5]), r0 = lerp(0.02, 0.12, g[6]), bend = lerp(-0.4, 0.4, g[7]);
        let levels = 3 + Math.floor(g[4] * 3.999);
        while (Math.pow(n, levels + 1) > 1400) levels--;
        const axis = down ? [0, 0, -1] : [0, 0, 1], prims = [];
        const grow = (p, d, L, r, lvl, az) => {
          const q = add(p, mul(d, L));
          prims.push({ a: p, b: q, r });
          if (lvl === levels) return;
          const u = perp(d), v = cross(d, u);
          for (let k = 0; k < n; k++) {
            const ang = az + (k * 2 * Math.PI) / n;
            const side = add(mul(u, Math.cos(ang)), mul(v, Math.sin(ang)));
            let nd = norm(add(mul(d, Math.cos(spread)), mul(side, Math.sin(spread))));
            nd = norm(add(nd, mul(axis, bend)));
            grow(q, nd, L * lr, r * rr, lvl + 1, az + twist);
          }
        };
        grow([0, 0, 0], axis, 1, r0, 0, 0);
        return prims;
      },
      describe(g) {
        return `${2 + Math.floor(g[0] * 2.999)} branches at each split, spreading ${Math.round(lerp(15, 80, g[1]))}°, each generation ${Math.round(lerp(0.5, 0.85, g[2]) * 100)}% as long.`;
      },
    },
    snow: {
      name: "Snow crystal", imitates: "snowflakes, frost, fresh snow",
      idea: "Six arms with side branches that grow fastest at their tips. Fresh snow is about 90% air and almost silent.",
      genes: ["plates", "spacing", "side", "subbranches", "thickness", "armTaper"],
      build(g) {
        const plates = 1 + Math.floor(g[0] * 3.999), spacing = lerp(0.08, 0.3, g[1]), side = lerp(0.3, 0.8, g[2]);
        const sub2 = g[3] > 0.5, t = lerp(0.01, 0.05, g[4]), taper = lerp(0.4, 1, g[5]), prims = [];
        for (let k = 0; k < plates; k++) {
          const phi = (k * Math.PI) / plates, h = [Math.cos(phi), Math.sin(phi), 0], z = [0, 0, 1];
          const inPlane = (a) => add(mul(h, Math.cos(a)), mul(z, Math.sin(a)));
          for (let arm = 0; arm < 6; arm++) {
            const a = (arm * Math.PI) / 3 + Math.PI / 6, d = inPlane(a);
            prims.push({ a: [0, 0, 0], b: d, r: t });
            for (let s = spacing; s < 0.95; s += spacing) {
              const base = mul(d, s), len = side * (1 - s), rb = t * lerp(1, taper, s);
              for (const sgn of [-1, 1]) {
                const bd = inPlane(a + sgn * Math.PI / 3), tip = add(base, mul(bd, len));
                prims.push({ a: base, b: tip, r: rb * 0.8 });
                if (sub2) for (let s2 = 0.35; s2 < 0.9; s2 += 0.3) {
                  const b2 = add(base, mul(bd, len * s2)), d2 = inPlane(a + sgn * Math.PI / 3 - sgn * Math.PI / 3);
                  prims.push({ a: b2, b: add(b2, mul(d2, len * (1 - s2) * 0.5)), r: rb * 0.6 });
                }
              }
            }
          }
        }
        return prims;
      },
      describe(g) {
        return `${1 + Math.floor(g[0] * 3.999)} crystal plates crossing each other, side branches every ${Math.round(lerp(0.08, 0.3, g[1]) * 100)}% of an arm${g[3] > 0.5 ? ", with a second level of branches" : ""}.`;
      },
    },
    coral: {
      name: "Folding coral", imitates: "lettuce coral, brain coral, kale leaves",
      idea: "The edge grows faster than the centre, so the sheet ruffles. Ruffles add surface and edges without adding weight.",
      genes: ["discs", "ruffles", "amplitude", "growth", "rings", "thickness", "spacing"],
      build(g, down) {
        const discs = 1 + Math.floor(g[0] * 4.999), m = 3 + Math.floor(g[1] * 9.999), amp = lerp(0, 0.35, g[2]), grow = lerp(1, 3, g[3]);
        const rings = 3 + Math.floor(g[4] * 3.999), t = lerp(0.008, 0.03, g[5]), gap = lerp(0.25, 0.6, g[6]), N = 28, prims = [];
        const sgn = down ? -1 : 1;
        for (let j = 0; j < discs; j++) {
          const zc = sgn * (0.15 + j * gap), pt = (r, th) => [r * Math.cos(th), r * Math.sin(th), zc + amp * Math.pow(r, grow) * Math.sin(m * th + j)];
          for (let k = 1; k <= rings; k++) {
            const r = k / rings;
            for (let i = 0; i < N; i++) prims.push({ a: pt(r, (i * 2 * Math.PI) / N), b: pt(r, ((i + 1) * 2 * Math.PI) / N), r: t });
          }
          for (let s = 0; s < 8; s++) {
            const th = (s * 2 * Math.PI) / 8;
            prims.push({ a: [0, 0, zc], b: pt(1 / rings, th), r: t });
            for (let k = 1; k < rings; k++) prims.push({ a: pt(k / rings, th), b: pt((k + 1) / rings, th), r: t });
          }
        }
        prims.push({ a: [0, 0, 0], b: [0, 0, sgn * (0.15 + (discs - 1) * gap)], r: t * 1.6 });
        return prims;
      },
      describe(g) {
        return `${1 + Math.floor(g[0] * 4.999)} ruffled sheets on a stem, ${3 + Math.floor(g[1] * 9.999)} folds around the edge.`;
      },
    },
    shell: {
      name: "Spiral shell", imitates: "nautilus, snail and turret shells",
      idea: "A tube that grows along a spiral, wider each turn, with chambers inside. Hollow chambers catch the lower pitches.",
      genes: ["turns", "expansion", "pitch", "tube", "rings", "lattice", "chambers", "thickness"],
      build(g, down) {
        const turns = lerp(1.5, 4, g[0]), b = lerp(0.08, 0.25, g[1]), pitch = lerp(0, 0.6, g[2]), k = lerp(0.25, 0.6, g[3]);
        const perTurn = 6 + Math.floor(g[4] * 8.999), nphi = 6 + Math.floor(g[5] * 6.999), septa = g[6], tw = lerp(0.01, 0.04, g[7]);
        const total = Math.max(4, Math.round(turns * perTurn)), end = turns * 2 * Math.PI, rho = (th) => Math.exp(b * th);
        const sgn = down ? -1 : 1, rEnd = rho(end), prims = [];
        const center = (th) => [rho(th) * Math.cos(th), rho(th) * Math.sin(th), sgn * pitch * (rEnd - rho(th)) * 2];
        const rand = rng(seedOf(g)), ringsPts = [];
        for (let i = 0; i <= total; i++) {
          const th = (i / total) * end, c = center(th), T = norm(sub(center(th + 0.01), c));
          let N = norm(sub([Math.cos(th), Math.sin(th), 0], mul(T, dot(T, [Math.cos(th), Math.sin(th), 0]))));
          const B = cross(T, N), rt = k * rho(th);
          const ring = [];
          for (let j = 0; j < nphi; j++) { const f = (j / nphi) * 2 * Math.PI; ring.push(add(c, add(mul(N, rt * Math.cos(f)), mul(B, rt * Math.sin(f))))); }
          ringsPts.push({ c, ring });
        }
        const r = tw * rEnd;
        ringsPts.forEach(({ c, ring }, i) => {
          for (let j = 0; j < nphi; j++) {
            prims.push({ a: ring[j], b: ring[(j + 1) % nphi], r });
            if (i > 0) prims.push({ a: ringsPts[i - 1].ring[j], b: ring[j], r });
          }
          if (i > 0 && i < ringsPts.length - 1 && rand() < septa * 0.5) for (const p of ring) prims.push({ a: c, b: p, r: r * 0.8 });
        });
        return prims;
      },
      describe(g) {
        return `${lerp(1.5, 4, g[0]).toFixed(1)} turns, ${g[2] > 0.3 ? "a tall turret" : g[2] > 0.1 ? "a low cone" : "flat like a nautilus"}, open lattice walls${g[6] > 0.2 ? " with inner chambers" : ""}.`;
      },
    },
    foam: {
      name: "Foam cloud", imitates: "sea foam, frogspawn, clouds of seeds",
      idea: "Many balls of different sizes, packed loosely so sound can reach between them.",
      genes: ["count", "spread", "aspect", "gap", "fill"],
      build(g) {
        const n = 8 + Math.floor(g[0] * 52), spread = g[1], aspect = lerp(0.4, 2.5, g[2]), gap = lerp(0, 0.6, g[3]), fill = lerp(0.15, 0.45, g[4]);
        const rand = rng(seedOf(g)), shape = Array.from({ length: n }, () => 1 - spread * rand() * 0.8);
        const base = Math.cbrt((fill * aspect) / shape.reduce((a, s) => a + s * s * s, 0));
        const balls = shape.map((s) => {
          let p;
          do p = [rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1]; while (dot(p, p) > 1);
          return { p: [p[0], p[1], p[2] * aspect], r: s * base };
        });
        for (let it = 0; it < 60; it++) {
          for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
            const d = sub(balls[j].p, balls[i].p), dist = Math.hypot(...d) || 1e-6, want = (balls[i].r + balls[j].r) * (1 + gap);
            if (dist < want) { const push = mul(d, (want - dist) / dist / 2); balls[i].p = sub(balls[i].p, push); balls[j].p = add(balls[j].p, push); }
          }
          for (const b of balls) { const e = Math.hypot(b.p[0], b.p[1], b.p[2] / aspect); if (e > 1) b.p = mul(b.p, 1 / e); }
        }
        return balls.map((b) => ({ a: b.p, b: b.p, r: b.r }));
      },
      describe(g) {
        return `${8 + Math.floor(g[0] * 52)} balls, ${g[1] > 0.5 ? "very different sizes" : "similar sizes"}, ${lerp(0.4, 2.5, g[2]) > 1.2 ? "tall" : lerp(0.4, 2.5, g[2]) < 0.8 ? "flat" : "round"} cloud.`;
      },
    },
  };

  // Every rule also gets a "material" gene (last) choosing the soft material.
  const geneCount = (rule) => RULES[rule].genes.length + 1;

  function bounds(prims) {
    const min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
    for (const p of prims) for (const q of [p.a, p.b]) for (let i = 0; i < 3; i++) { min[i] = Math.min(min[i], q[i] - p.r); max[i] = Math.max(max[i], q[i] + p.r); }
    return { min, max };
  }

  // Scale and place a form inside its envelope: hanging from z = 0, or standing on z = 0.
  function fit(prims, { width, height, mount }) {
    const { min, max } = bounds(prims);
    const s = Math.min(width / (max[0] - min[0]), width / (max[1] - min[1]), height / (max[2] - min[2]));
    const cx = (min[0] + max[0]) / 2, cy = (min[1] + max[1]) / 2, cz = mount === "hanging" ? max[2] : min[2];
    const tr = (q) => [(q[0] - cx) * s, (q[1] - cy) * s, (q[2] - cz) * s];
    return prims.map((p) => ({ a: tr(p.a), b: tr(p.b), r: p.r * s }));
  }

  function envelope(size, mount) {
    return { width: size, height: mount === "hanging" ? size : Math.min(2.2, size * 1.5), mount };
  }

  function grow(rule, genes, opts) {
    const raw = RULES[rule].build(genes, opts.mount === "hanging");
    return fit(raw, envelope(opts.size, opts.mount));
  }

  // ------------------------------------------------------------- voxels ----
  function segDist(c, a, b) {
    const ab = sub(b, a), t = Math.max(0, Math.min(1, dot(sub(c, a), ab) / (dot(ab, ab) || 1e-12)));
    return Math.hypot(...sub(c, add(a, mul(ab, t))));
  }
  // Solid fraction per voxel (thin parts spread over neighbouring voxels so their volume is kept).
  function voxelize(prims, dx) {
    const { min, max } = bounds(prims);
    const o = min.map((m) => m - dx), n = max.map((m, i) => Math.max(1, Math.ceil((m + dx - o[i]) / dx)));
    const frac = new Float32Array(n[0] * n[1] * n[2]), R = 0.75 * dx, NORM = 0.589;
    for (const p of prims) {
      const lo = [0, 1, 2].map((i) => Math.max(0, Math.floor((Math.min(p.a[i], p.b[i]) - p.r - R - o[i]) / dx)));
      const hi = [0, 1, 2].map((i) => Math.min(n[i] - 1, Math.ceil((Math.max(p.a[i], p.b[i]) + p.r + R - o[i]) / dx)));
      const thin = p.r < dx / 2, sphere = p.a[0] === p.b[0] && p.a[1] === p.b[1] && p.a[2] === p.b[2];
      for (let k = lo[2]; k <= hi[2]; k++) for (let j = lo[1]; j <= hi[1]; j++) for (let i = lo[0]; i <= hi[0]; i++) {
        const c = [o[0] + (i + 0.5) * dx, o[1] + (j + 0.5) * dx, o[2] + (k + 0.5) * dx], d = segDist(c, p.a, p.b);
        let f;
        if (!thin) f = Math.max(0, Math.min(1, (p.r - d) / dx + 0.5));
        else if (sphere) f = d < dx / 2 ? (4 / 3) * Math.PI * p.r ** 3 / dx ** 3 : 0;
        else f = Math.max(0, 1 - d / R) * (Math.PI * p.r * p.r) / (dx * dx) / NORM;
        if (f > 0) frac[i + n[0] * (j + n[1] * k)] += f;
      }
    }
    let vol = 0;
    for (let i = 0; i < frac.length; i++) { if (frac[i] > 1) frac[i] = 1; vol += frac[i]; }
    return { frac, n, o, dx, volume: vol * dx ** 3 };
  }

  // A slice of the voxels perpendicular to x (axis 0) or y (axis 1), at a fraction `at` (0..1)
  // of the form's extent along that axis. `extent` is that extent in metres.
  function slice(vox, axis, at) {
    const { frac, n, dx } = vox, h = axis === 0 ? 1 : 0, nh = n[h], nz = n[2];
    let lo = Infinity, hi = -Infinity;
    for (let k = 0; k < nz; k++) for (let j = 0; j < n[1]; j++) for (let i = 0; i < n[0]; i++) {
      if (!frac[i + n[0] * (j + n[1] * k)]) continue;
      const a = axis === 0 ? i : j;
      lo = Math.min(lo, a); hi = Math.max(hi, a);
    }
    const layer = Math.round(lo + at * (hi - lo)), grid = new Float32Array(nh * nz);
    for (let k = 0; k < nz; k++) for (let b = 0; b < nh; b++) {
      const i = axis === 0 ? layer : b, j = axis === 0 ? b : layer;
      grid[b + nh * k] = frac[i + n[0] * (j + n[1] * k)];
    }
    return { grid, nh, nz, dx, extent: (hi - lo + 1) * dx, x0: vox.o[h], z0: vox.o[2] };
  }
  // Mean solid fraction along an axis (used only for drawing the form in the wave view).
  function project(vox, axis) {
    const { frac, n, dx } = vox, h = axis === 0 ? 1 : 0, nh = n[h], nz = n[2], grid = new Float32Array(nh * nz);
    for (let k = 0; k < nz; k++) for (let j = 0; j < n[1]; j++) for (let i = 0; i < n[0]; i++) {
      const b = axis === 0 ? j : i;
      grid[b + nh * k] = Math.max(grid[b + nh * k], frac[i + n[0] * (j + n[1] * k)]);
    }
    return { grid, nh, nz, dx, extent: n[axis] * dx, x0: vox.o[h], z0: vox.o[2] };
  }
  function downsample(pr, factor) {
    if (factor === 1) return pr;
    const nh = Math.ceil(pr.nh / factor), nz = Math.ceil(pr.nz / factor), grid = new Float32Array(nh * nz);
    for (let k = 0; k < pr.nz; k++) for (let i = 0; i < pr.nh; i++) grid[Math.floor(i / factor) + nh * Math.floor(k / factor)] += pr.grid[i + pr.nh * k] / (factor * factor);
    return { ...pr, grid, nh, nz, dx: pr.dx * factor };
  }

  // ------------------------------------------------------------------ FDTD -
  // 2D room section (width x height, rigid walls), the projected form placed at (cx, z), a
  // tone burst near a floor corner. Returns the energy decay rate and the absorption length.
  function roomFor(mount) { return mount === "hanging" ? { w: 6, h: 4, cx: 3.4 } : { w: 6, h: 3.5, cx: 3.4 }; }

  function makeSim(pr, sigma, freq, room, mount) {
    const dx = pr.dx, nx = Math.round(room.w / dx), nz = Math.round(room.h / dx), dt = (0.68 * dx) / C0;
    const p = new Float32Array(nx * nz), vx = new Float32Array((nx + 1) * nz), vz = new Float32Array(nx * (nz + 1)), s = new Float32Array(nx * nz);
    const ox = Math.round((room.cx + pr.x0) / dx), oz = mount === "hanging" ? Math.round((room.h + pr.z0) / dx) : Math.round(pr.z0 / dx);
    for (let k = 0; k < pr.nz; k++) for (let i = 0; i < pr.nh; i++) {
      const X = ox + i, Z = oz + k;
      if (X >= 0 && X < nx && Z >= 0 && Z < nz) s[X + nx * Z] = Math.max(s[X + nx * Z], pr.grid[i + pr.nh * k]);
    }
    // Two parallel paths in every face (as in the rest of the app): a share s of the air passes
    // through the felt and is slowed by its flow resistance; the rest flows freely around it.
    const base = dt / (RHO0 * dx), pf = [];
    const addFace = (arr, idx, L, R, sf) => { if (sf > 0) pf.push({ arr, idx, L, R, sf }); };
    for (let k = 0; k < nz; k++) for (let i = 1; i < nx; i++) addFace(0, i + (nx + 1) * k, i - 1 + nx * k, i + nx * k, (s[i - 1 + nx * k] + s[i + nx * k]) / 2);
    for (let k = 1; k < nz; k++) for (let i = 0; i < nx; i++) addFace(1, i + nx * k, i + nx * (k - 1), i + nx * k, (s[i + nx * (k - 1)] + s[i + nx * k]) / 2);
    const nf = pf.length, fArr = new Uint8Array(nf), fIdx = new Int32Array(nf), fL = new Int32Array(nf), fR = new Int32Array(nf);
    const fS = new Float32Array(nf), fK1 = new Float32Array(nf), fK2 = new Float32Array(nf), q = new Float32Array(nf);
    pf.forEach((f, j) => {
      const a = (sigma * dt) / (2 * RHO0);
      fArr[j] = f.arr; fIdx[j] = f.idx; fL[j] = f.L; fR[j] = f.R; fS[j] = f.sf; fK1[j] = (1 - a) / (1 + a); fK2[j] = base / (1 + a);
    });
    const src = Math.round(0.4 / dx) + nx * Math.round(0.4 / dx), cycles = 6, tb = cycles / freq, pc = RHO0 * C0 * C0 * dt / dx;
    let step = 0;
    return {
      nx, nz, dx, dt, p, s, get t() { return step * dt; },
      step(nsteps = 1) {
        for (let n = 0; n < nsteps; n++, step++) {
          for (let j = 0; j < nf; j++) q[j] = fK1[j] * q[j] - fK2[j] * (p[fR[j]] - p[fL[j]]);
          for (let k = 0; k < nz; k++) {
            const row = nx * k, rv = (nx + 1) * k;
            for (let i = 1; i < nx; i++) vx[rv + i] -= base * (p[row + i] - p[row + i - 1]);
          }
          for (let k = 1; k < nz; k++) {
            const row = nx * k;
            for (let i = 0; i < nx; i++) vz[row + i] -= base * (p[row + i] - p[row - nx + i]);
          }
          for (let k = 0; k < nz; k++) {
            const row = nx * k, rv = (nx + 1) * k;
            for (let i = 0; i < nx; i++) p[row + i] -= pc * (vx[rv + i + 1] - vx[rv + i] + vz[row + nx + i] - vz[row + i]);
          }
          // The felt share of each porous face carries q instead of the free-air velocity.
          for (let j = 0; j < nf; j++) {
            const d = pc * fS[j] * (q[j] - (fArr[j] ? vz[fIdx[j]] : vx[fIdx[j]]));
            p[fR[j]] += d; p[fL[j]] -= d;
          }
          const t = step * dt;
          if (t < tb) p[src] += Math.sin(2 * Math.PI * freq * t) * Math.sin((Math.PI * t) / tb) ** 2;
        }
      },
      energy() {
        let e = 0, ev = 0;
        for (let i = 0; i < p.length; i++) e += p[i] * p[i];
        for (let i = 0; i < vx.length; i++) ev += vx[i] * vx[i];
        for (let i = 0; i < vz.length; i++) ev += vz[i] * vz[i];
        for (let j = 0; j < nf; j++) { const v = fArr[j] ? vz[fIdx[j]] : vx[fIdx[j]]; ev += fS[j] * (q[j] * q[j] - v * v); }
        return e / (RHO0 * C0 * C0) + RHO0 * ev;
      },
      burstEnd: tb,
    };
  }

  function absorptionLength(pr, sigma, freq, room, mount) {
    const sim = makeSim(pr, sigma, freq, room, mount), t0 = sim.burstEnd + 0.025, t1 = t0 + 0.12;
    const ts = [], es = [], every = 8;
    while (sim.t < t1) {
      sim.step(every);
      if (sim.t >= t0) { ts.push(sim.t); es.push(Math.log(Math.max(sim.energy(), 1e-30))); }
    }
    const n = ts.length, mt = ts.reduce((a, b) => a + b, 0) / n, me = es.reduce((a, b) => a + b, 0) / n;
    let num = 0, den = 0;
    for (let i = 0; i < n; i++) { num += (ts[i] - mt) * (es[i] - me); den += (ts[i] - mt) ** 2; }
    const decay = Math.max(0, -num / den);
    return { decay, length: (Math.PI * room.w * room.h * decay) / C0 };
  }

  // Full evaluation of one genome. material: {sigma, density}.
  function evaluate({ rule, genes, mount, size, material, budget, weights }) {
    const prims = grow(rule, genes, { mount, size });
    const vox = voxelize(prims, C0 / BANDS[BANDS.length - 1] / CELLS_PER_WAVE);
    const mass = vox.volume * material.density;
    // Strip method: absorption per metre in two slices per direction, times the extent.
    const room = roomFor(mount), slices = [];
    for (const axis of [0, 1]) for (const at of [0.3, 0.7]) slices.push({ axis, sl: slice(vox, axis, at) });
    const area = BANDS.map((f, bi) => {
      const factor = 2 ** (BANDS.length - 1 - bi);
      const perAxis = [0, 1].map((axis) => {
        const mine = slices.filter((x) => x.axis === axis);
        const mean = mine.reduce((acc, { sl }) => acc + absorptionLength(downsample(sl, factor), material.sigma, f, room, mount).length, 0) / mine.length;
        return mean * mine[0].sl.extent;
      });
      return (perAxis[0] + perAxis[1]) / 2;
    });
    const w = weights || [1 / 3, 1 / 3, 1 / 3];
    const score = area.reduce((a, x, i) => a + x * w[i], 0);
    const over = mass > budget;
    const fitness = score - (over ? score * (mass / budget - 1) + 0.05 : 0);
    return { rule, genes, mount, size, prims, mass, volume: vox.volume, area, score, fitness, perKg: score / Math.max(mass, 1e-6) };
  }

  // A flat 50 mm panel of the same material and mass on a wall, through the same solver,
  // so forms are always compared like for like.
  function panelBaseline({ material, mass, thickness = 0.05 }) {
    const side = Math.sqrt(mass / (material.density * thickness)), room = roomFor("standing");
    return BANDS.map((f) => {
      const dx = C0 / f / CELLS_PER_WAVE, nh = Math.max(1, Math.round(thickness / dx)), nz = Math.max(1, Math.round(side / dx));
      const fill = Math.min(1, thickness / (nh * dx));
      const pr = { grid: new Float32Array(nh * nz).fill(fill), nh, nz, dx, x0: -room.cx, z0: 1.0, extent: side };
      return absorptionLength(pr, material.sigma, f, room, "standing").length * side;
    });
  }

  // One generation of a genetic algorithm over gene vectors in [0, 1].
  function nextGeneration(pop, fit, rand, { elite = 2, mutation = 0.18 } = {}) {
    const normal = () => Math.sqrt(-2 * Math.log(1 - rand())) * Math.cos(2 * Math.PI * rand());
    const order = fit.map((f, i) => i).sort((a, b) => fit[b] - fit[a]);
    const pick = () => { let b = Math.floor(rand() * pop.length); for (let k = 0; k < 2; k++) { const i = Math.floor(rand() * pop.length); if (fit[i] > fit[b]) b = i; } return pop[b]; };
    const next = order.slice(0, elite).map((i) => pop[i].slice());
    while (next.length < pop.length) {
      const a = pick(), b = pick();
      next.push(a.map((x, i) => {
        let g = rand() < 0.5 ? x : b[i];
        if (rand() < mutation) g = Math.min(1, Math.max(0, g + 0.15 * normal()));
        if (rand() < mutation / 4) g = rand();
        return g;
      }));
    }
    return next;
  }

  const Growth = { RULES, BANDS, geneCount, grow, voxelize, project, slice, downsample, makeSim, absorptionLength, roomFor, evaluate, panelBaseline, nextGeneration, rng, bounds, envelope };
  if (typeof module !== "undefined" && module.exports) module.exports = Growth;
  else root.Growth = Growth;
})(typeof self !== "undefined" ? self : this);
