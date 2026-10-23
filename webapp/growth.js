/* Biothing Sound Lab growth: free forms grown by natural rules, tested with a wave simulator.
 *
 * A form is a list of parts in metres, z up: capsules {a, b, r} (a sphere when a == b) and
 * flat sheets {c, u, v, t} (centre, two half-edge vectors, absolute thickness in metres).
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
  // Felt sheet thickness is a manufacturing value: 15-80 mm, kept as is when a form is scaled.
  const sheetT = (u) => lerp(0.015, 0.08, u);
  const mmT = (u) => Math.round(sheetT(u) * 1000) + " mm";
  const isSheet = (p) => p.c !== undefined;

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
    // ---- hybrid rules: flat felt sheets (efficient, like panels) arranged by natural growth.
    leaves: {
      group: "hybrid", name: "Leaf canopy", imitates: "a tree canopy, ginkgo, lily pads",
      idea: "Thin branches carry flat felt leaves. Every leaf works on both faces and the canopy stays open, so sound reaches each one.",
      genes: ["branches", "spread", "length", "levels", "leaf", "tilt", "twist", "sheet"],
      build(g, down) {
        const n = 2 + Math.floor(g[0] * 1.999), spread = lerp(20, 70, g[1]) * Math.PI / 180, lr = lerp(0.55, 0.8, g[2]);
        const levels = 2 + Math.floor(g[3] * 2.999), leaf = lerp(0.15, 0.45, g[4]), tilt = lerp(0, Math.PI / 2, g[5]), twist = lerp(0, Math.PI, g[6]);
        const t = sheetT(g[7]), axis = down ? [0, 0, -1] : [0, 0, 1], parts = [];
        const grow = (pt, d, L, lvl, az) => {
          const q = add(pt, mul(d, L));
          parts.push({ a: pt, b: q, r: 0.015 });
          if (lvl === levels) {
            const side = perp(d), up = norm(add(mul(d, Math.cos(tilt)), mul(cross(d, side), Math.sin(tilt))));
            parts.push({ c: add(q, mul(up, leaf)), u: mul(up, leaf), v: mul(side, leaf * 0.75), t });
            return;
          }
          const u = perp(d), v = cross(d, u);
          for (let k = 0; k < n; k++) {
            const ang = az + (k * 2 * Math.PI) / n, sd = add(mul(u, Math.cos(ang)), mul(v, Math.sin(ang)));
            grow(q, norm(add(mul(d, Math.cos(spread)), mul(sd, Math.sin(spread)))), L * lr, lvl + 1, az + twist);
          }
        };
        grow([0, 0, 0], axis, 1, 0, 0);
        return parts;
      },
      describe(g) { return `${Math.pow(2 + Math.floor(g[0] * 1.999), 2 + Math.floor(g[3] * 2.999))} felt leaves ${mmT(g[7])} thick on a branching frame.`; },
    },
    gills: {
      group: "hybrid", name: "Mushroom gills", imitates: "gills under a mushroom cap, fish gills, turbine fins",
      idea: "Radial fins packed around a centre, like the gills under a mushroom. Lots of sheet surface in a small volume, open between the fins.",
      genes: ["fins", "tiers", "inner", "twist", "height", "cap", "sheet"],
      build(g, down) {
        const fins = 8 + Math.floor(g[0] * 32), tiers = 1 + Math.floor(g[1] * 2.999), rin = lerp(0, 0.4, g[2]), twist = lerp(0, Math.PI / 3, g[3]);
        const h = lerp(0.3, 0.9, g[4]), t = sheetT(g[6]), sgn = down ? -1 : 1, parts = [];
        for (let j = 0; j < tiers; j++) {
          const zc = sgn * (j + 0.5) * h * 1.08;
          for (let k = 0; k < fins; k++) {
            const th = (k * 2 * Math.PI) / fins + j * Math.PI / fins, rad = [Math.cos(th), Math.sin(th), 0], tang = [-Math.sin(th), Math.cos(th), 0];
            const up = norm(add(mul([0, 0, 1], Math.cos(twist)), mul(tang, Math.sin(twist))));
            parts.push({ c: add(mul(rad, (rin + 1) / 2), [0, 0, zc]), u: mul(rad, (1 - rin) / 2), v: mul(up, h / 2), t });
          }
        }
        if (g[5] > 0.5) parts.push({ c: [0, 0, sgn * -0.02], u: [0.75, 0, 0], v: [0, 0.75, 0], t });
        return parts;
      },
      describe(g) { return `${8 + Math.floor(g[0] * 32)} fins in ${1 + Math.floor(g[1] * 2.999)} tier(s), ${mmT(g[6])} felt${g[3] > 0.3 ? ", spiralling" : ""}${g[5] > 0.5 ? ", under a cap" : ""}.`; },
    },
    cone: {
      group: "hybrid", name: "Pine cone", imitates: "pine cones, artichokes, seed heads",
      idea: "Scales on a spiral (the golden angle) that open outward. Open scales catch sound on both faces and shade each other only a little.",
      genes: ["count", "aspect", "scale", "opening", "core", "sheet"],
      build(g) {
        const n = 20 + Math.floor(g[0] * 100), aspect = lerp(0.6, 2.2, g[1]), sc = lerp(0.12, 0.35, g[2]), open = lerp(0, 1.4, g[3]);
        const core = lerp(0.2, 0.5, g[4]), t = sheetT(g[5]), golden = Math.PI * (3 - Math.sqrt(5)), parts = [];
        parts.push({ a: [0, 0, -aspect * core], b: [0, 0, aspect * core], r: core * 0.7 });
        for (let i = 0; i < n; i++) {
          const z = 1 - (2 * (i + 0.5)) / n, r = Math.sqrt(1 - z * z), th = i * golden;
          const P = [r * Math.cos(th) * core, r * Math.sin(th) * core, z * aspect * core];
          const N = norm([P[0], P[1], P[2] / (aspect * aspect)]), side = norm(cross([0, 0, 1], N)), merid = norm(cross(N, side));
          const lift = norm(add(mul(merid, Math.cos(open)), mul(N, Math.sin(open))));
          parts.push({ c: add(P, mul(lift, sc)), u: mul(lift, sc), v: mul(side, sc * 0.8), t });
        }
        return parts;
      },
      describe(g) { return `${20 + Math.floor(g[0] * 100)} scales ${mmT(g[5])} thick, opened ${Math.round(lerp(0, 80, g[3]))}°, ${lerp(0.6, 2.2, g[1]) > 1.3 ? "tall" : "round"} cone.`; },
    },
    lamellae: {
      group: "hybrid", name: "Twisted lamellae", imitates: "bone lamellae, beetle shells (helicoidal layers), book gills",
      idea: "Stacked sheets, each turned a little from the last, like the layers in a beetle's shell. Gaps between sheets let sound in; the twist catches it from every direction.",
      genes: ["layers", "gradient", "twist", "aspect", "taper", "upright", "sheet"],
      build(g, down) {
        const layers = 3 + Math.floor(g[0] * 13), grad = lerp(-0.5, 0.5, g[1]), twist = lerp(0, Math.PI / 3, g[2]), aspect = lerp(1, 3, g[3]);
        const taper = lerp(0.4, 1, g[4]), upright = g[5] > 0.5, t = sheetT(g[6]), parts = [];
        let pos = 0;
        for (let j = 0; j < layers; j++) {
          const f = layers > 1 ? j / (layers - 1) : 0, size = lerp(1, taper, f), ang = j * twist;
          const a = [Math.cos(ang) * aspect * 0.5 * size, Math.sin(ang) * aspect * 0.5 * size, 0], b = [-Math.sin(ang) * 0.5 * size, Math.cos(ang) * 0.5 * size, 0];
          const c = upright ? [pos, 0, 0] : [0, 0, (down ? -1 : 1) * pos];
          const rot = (q) => (upright ? [q[2], q[0], q[1]] : q);
          parts.push({ c, u: rot(a), v: rot(b), t });
          pos += 0.15 * (1 + grad * (f * 2 - 1));
        }
        return parts;
      },
      describe(g) { return `${3 + Math.floor(g[0] * 13)} ${g[5] > 0.5 ? "upright" : "stacked"} sheets ${mmT(g[6])} thick, each turned ${Math.round(lerp(0, 60, g[2]))}° from the last.`; },
    },
  };

  // Every rule also gets four "stance" genes (used when the form stands on its own: root flare,
  // holdfast type, base size, number of roots) and a last gene choosing the soft material.
  const STANCE = 4;
  const geneCount = (rule) => RULES[rule].genes.length + STANCE + 1;
  const MIN_TIP_ANGLE = 10; // degrees a freestanding piece must tilt before it topples

  const corners = (p) => [add(add(p.c, p.u), p.v), add(sub(p.c, p.u), p.v), sub(sub(p.c, p.u), p.v), sub(add(p.c, p.u), p.v)];
  function bounds(prims) {
    const min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
    for (const p of prims) {
      const pts = isSheet(p) ? corners(p) : [p.a, p.b], pad = isSheet(p) ? p.t / 2 : p.r;
      for (const q of pts) for (let i = 0; i < 3; i++) { min[i] = Math.min(min[i], q[i] - pad); max[i] = Math.max(max[i], q[i] + pad); }
    }
    return { min, max };
  }

  // Scale and place a form inside its envelope: hanging from z = 0, or standing on z = 0.
  // Sheets keep their real thickness; everything else scales.
  function fit(prims, { width, height, mount }) {
    const { min, max } = bounds(prims.map((p) => (isSheet(p) ? { ...p, t: 0 } : p)));
    const s = Math.min(width / (max[0] - min[0] || 1e-9), width / (max[1] - min[1] || 1e-9), height / (max[2] - min[2] || 1e-9));
    const cx = (min[0] + max[0]) / 2, cy = (min[1] + max[1]) / 2, cz = mount === "hanging" ? max[2] : min[2];
    const tr = (q) => [(q[0] - cx) * s, (q[1] - cy) * s, (q[2] - cz) * s];
    return prims.map((p) => (isSheet(p) ? { c: tr(p.c), u: mul(p.u, s), v: mul(p.v, s), t: p.t } : { a: tr(p.a), b: tr(p.b), r: p.r * s }));
  }

  function envelope(size, mount) {
    return { width: size, height: mount === "hanging" ? size : Math.min(2.2, size * 1.5), mount };
  }

  function grow(rule, genes, opts) {
    const env = envelope(opts.size, opts.mount), raw = RULES[rule].build(genes, opts.mount === "hanging");
    const parts = fit(raw, env);
    return opts.mount === "standing" ? standUp(parts, genes.slice(RULES[rule].genes.length, RULES[rule].genes.length + STANCE), env) : parts;
  }

  // Freestanding adaptation, the way trees and corals do it: a flare that widens the form towards
  // the ground, then a holdfast (splayed roots, or a weighted plinth) with a central stem.
  // Holdfast parts are structure (rigid, not absorbing); everything else stays felt.
  function standUp(parts, st, env) {
    const flare = lerp(0, 0.6, st[0]), H = env.height;
    const k = (z) => 1 + flare * (1 - (2 * Math.max(0, Math.min(H, z))) / H);
    const fl = (q) => [q[0] * k(q[2]), q[1] * k(q[2]), q[2]];
    let out = parts.map((p) => (isSheet(p) ? { ...p, c: fl(p.c) } : { ...p, a: fl(p.a), b: fl(p.b) }));
    const b = bounds(out), wide = Math.max(b.max[0] - b.min[0], b.max[1] - b.min[1]);
    if (wide > env.width) {
      const s = env.width / wide, sx = (q) => [q[0] * s, q[1] * s, q[2]];
      out = out.map((p) => (isSheet(p) ? { ...p, c: sx(p.c), u: [p.u[0] * s, p.u[1] * s, p.u[2]], v: [p.v[0] * s, p.v[1] * s, p.v[2]] } : { ...p, a: sx(p.a), b: sx(p.b) }));
    }
    const R = lerp(0.25, 0.5, st[2]) * env.width, roots = st[1] < 0.5, legs = 3 + Math.floor(st[3] * 3.999), top = 0.45 * H;
    if (roots) {
      const hr = 0.12 * H;
      out.push({ a: [0, 0, hr], b: [0, 0, top], r: 0.015, structure: true });
      for (let i = 0; i < legs; i++) {
        const ang = (i * 2 * Math.PI) / legs;
        out.push({ a: [0, 0, hr], b: [R * Math.cos(ang), R * Math.sin(ang), 0.01], r: 0.01, structure: true });
      }
    } else {
      out.push({ c: [0, 0, 0.01], u: [R, 0, 0], v: [0, R, 0], t: 0.02, structure: true });
      out.push({ a: [0, 0, 0.02], b: [0, 0, top], r: 0.015, structure: true });
    }
    return out;
  }

  // Tipping check: the angle the piece can lean before its centre of mass leaves the footprint.
  function stability(parts, softCom, softMass, structureDensity) {
    let m = softMass, cx = softCom[0] * softMass, cy = softCom[1] * softMass, cz = softCom[2] * softMass, structureMass = 0, support = 0;
    for (const p of parts) {
      if (!p.structure) continue;
      let vol, c;
      if (isSheet(p)) {
        vol = 4 * Math.hypot(...cross(p.u, p.v)) * p.t; c = p.c;
        support = Math.max(support, Math.min(Math.hypot(...p.u), Math.hypot(...p.v)));
      } else {
        vol = Math.PI * p.r * p.r * Math.hypot(...sub(p.b, p.a)); c = mul(add(p.a, p.b), 0.5);
        if (p.b[2] < 0.02) support = Math.max(support, Math.hypot(p.b[0], p.b[1]));
      }
      const pm = vol * structureDensity;
      structureMass += pm; m += pm; cx += c[0] * pm; cy += c[1] * pm; cz += c[2] * pm;
    }
    const legs = parts.filter((p) => p.structure && !isSheet(p) && p.b[2] < 0.02).length;
    if (legs >= 3) support *= Math.cos(Math.PI / legs); // inscribed radius of the root polygon
    const off = Math.hypot(cx / m, cy / m), h = Math.max(cz / m, 1e-3);
    const angle = (Math.atan2(support - off, h) * 180) / Math.PI;
    return { angle, ok: angle >= MIN_TIP_ANGLE, structureMass, support, comHeight: h };
  }

  // Compact transfer format (11 numbers per part) for workers and storage.
  // Type code: 0 capsule, 1 sheet, +2 when the part is structure.
  const STRIDE = 11;
  function pack(prims) {
    const f = new Float32Array(prims.length * STRIDE);
    prims.forEach((p, i) => f.set(isSheet(p) ? [p.structure ? 3 : 1, ...p.c, ...p.u, ...p.v, p.t] : [p.structure ? 2 : 0, ...p.a, ...p.b, p.r, 0, 0, 0], i * STRIDE));
    return f;
  }
  function unpack(f) {
    const out = [];
    for (let i = 0; i < f.length; i += STRIDE) {
      const x = Array.from(f.subarray(i + 1, i + STRIDE));
      const part = f[i] % 2 === 1 ? { c: x.slice(0, 3), u: x.slice(3, 6), v: x.slice(6, 9), t: x[9] } : { a: x.slice(0, 3), b: x.slice(3, 6), r: x[6] };
      if (f[i] >= 2) part.structure = true;
      out.push(part);
    }
    return out;
  }
  const count = (f) => f.length / STRIDE;

  // ------------------------------------------------------------- voxels ----
  function segDist(c, a, b) {
    const ab = sub(b, a), t = Math.max(0, Math.min(1, dot(sub(c, a), ab) / (dot(ab, ab) || 1e-12)));
    return Math.hypot(...sub(c, add(a, mul(ab, t))));
  }
  // Solid fraction per voxel (thin parts spread over neighbouring voxels so their volume is kept).
  function voxelize(all, dx) {
    const prims = all.filter((p) => !p.structure);
    const { min, max } = bounds(prims);
    const o = min.map((m) => m - dx), n = max.map((m, i) => Math.max(1, Math.ceil((m + dx - o[i]) / dx)));
    const frac = new Float32Array(n[0] * n[1] * n[2]), R = 0.75 * dx, NORM = 0.589;
    const range = (lo3, hi3) => [[0, 1, 2].map((i) => Math.max(0, Math.floor((lo3[i] - o[i]) / dx))), [0, 1, 2].map((i) => Math.min(n[i] - 1, Math.ceil((hi3[i] - o[i]) / dx)))];
    const soft = (x) => Math.max(0, Math.min(1, x / dx + 0.5));
    for (const p of prims) {
      if (isSheet(p)) {
        const b = bounds([p]), [lo, hi] = range(b.min.map((m) => m - dx), b.max.map((m) => m + dx));
        const Lu = Math.hypot(...p.u), Lv = Math.hypot(...p.v), uh = mul(p.u, 1 / (Lu || 1)), vh = mul(p.v, 1 / (Lv || 1)), nh = norm(cross(p.u, p.v));
        for (let k = lo[2]; k <= hi[2]; k++) for (let j = lo[1]; j <= hi[1]; j++) for (let i = lo[0]; i <= hi[0]; i++) {
          const d = sub([o[0] + (i + 0.5) * dx, o[1] + (j + 0.5) * dx, o[2] + (k + 0.5) * dx], p.c);
          const inside = soft(Lu - Math.abs(dot(d, uh))) * soft(Lv - Math.abs(dot(d, vh)));
          if (!inside) continue;
          const dn = Math.abs(dot(d, nh));
          const f = p.t >= dx ? soft(p.t / 2 - dn) : (p.t / dx) * Math.max(0, 1 - dn / dx);
          if (f > 0) frac[i + n[0] * (j + n[1] * k)] += f * inside;
        }
        continue;
      }
      const [lo, hi] = range([0, 1, 2].map((i) => Math.min(p.a[i], p.b[i]) - p.r - R), [0, 1, 2].map((i) => Math.max(p.a[i], p.b[i]) + p.r + R));
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
    const com = [0, 0, 0];
    for (let k = 0, idx = 0; k < n[2]; k++) for (let j = 0; j < n[1]; j++) for (let i = 0; i < n[0]; i++, idx++) {
      if (frac[idx] > 1) frac[idx] = 1;
      const f = frac[idx];
      if (!f) continue;
      vol += f; com[0] += f * (o[0] + (i + 0.5) * dx); com[1] += f * (o[1] + (j + 0.5) * dx); com[2] += f * (o[2] + (k + 0.5) * dx);
    }
    return { frac, n, o, dx, volume: vol * dx ** 3, com: com.map((c) => c / Math.max(vol, 1e-9)) };
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
  function evaluate({ rule, genes, mount, size, material, budget, weights, structureDensity = 680 }) {
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
    let fitness = score - (over ? score * (mass / budget - 1) + 0.05 : 0);
    const stance = mount === "standing" ? stability(prims, vox.com, mass, structureDensity) : null;
    if (stance && !stance.ok) fitness -= score * 0.8 * Math.min(1, (MIN_TIP_ANGLE - stance.angle) / MIN_TIP_ANGLE) + 0.05;
    return { rule, genes, mount, size, prims, mass, volume: vox.volume, area, score, fitness, perKg: score / Math.max(mass, 1e-6), stance };
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

  const Growth = { RULES, BANDS, STRIDE, STANCE, MIN_TIP_ANGLE, stability, pack, unpack, count, isSheet, geneCount, grow, voxelize, project, slice, downsample, makeSim, absorptionLength, roomFor, evaluate, panelBaseline, nextGeneration, rng, bounds, envelope };
  if (typeof module !== "undefined" && module.exports) module.exports = Growth;
  else root.Growth = Growth;
})(typeof self !== "undefined" ? self : this);
