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

  // ------------------------------------------------- shared by sound rules
  // Plan of a sheet-like piece along its length s in [-L/2, L/2]: straight, arc, S-curve or zigzag
  // (curved plans stand more easily). Returns the point, unit tangent and unit normal in plan.
  function planCurve(plan, curv, L) {
    const A = curv * L * 0.25, bend = Math.max(1e-3, curv * Math.PI * 0.9), Rc = L / bend;
    return (s2) => {
      let x = s2, y = 0, dy = 0;
      if (plan === 1) { x = Rc * Math.sin(s2 / Rc); y = Rc * (1 - Math.cos(s2 / Rc)) - Rc * (1 - Math.cos(L / 2 / Rc)) / 2; dy = Math.tan(s2 / Rc); }
      else if (plan === 2) { y = A * Math.sin((2 * Math.PI * s2) / L); dy = A * (2 * Math.PI / L) * Math.cos((2 * Math.PI * s2) / L); }
      else if (plan === 3) { const ph = ((s2 / L + 0.5) * 4) % 2; y = A * (ph < 1 ? ph : 2 - ph) - A / 2; dy = (ph < 1 ? 1 : -1) * A * 4 / L; }
      const tl = Math.hypot(1, dy), t = plan === 1 ? [Math.cos(s2 / Rc), Math.sin(s2 / Rc)] : [1 / tl, dy / tl];
      return { p: [x, y], t, n: [-t[1], t[0]] };
    };
  }
  // Bilinear sample of a field (nx x nz, values 0..1) at fractions u (along) and v (up), both 0..1.
  function sampleField(F, nx, nz, u0, v0) {
    if (nx < 2 || nz < 2) return F[0];
    const u = Math.max(0, Math.min(nx - 1, u0 * (nx - 1))), v = Math.max(0, Math.min(nz - 1, v0 * (nz - 1)));
    const i = Math.min(nx - 2, Math.floor(u)), k = Math.min(nz - 2, Math.floor(v)), fu = u - i, fv = v - k, q = (a, b) => F[a + nx * b];
    return (q(i, k) * (1 - fu) + q(i + 1, k) * fu) * (1 - fv) + (q(i, k + 1) * (1 - fu) + q(i + 1, k + 1) * fu) * fv;
  }
  // Structure: feet across the piece with posts (freestanding), or hangers from the ceiling.
  function footings(parts, curve, L, H, feet, footLen, down) {
    for (let k = 0; k < feet; k++) {
      const c = curve(-L / 2 + (L * (k + 0.5)) / feet), P = c.p;
      if (down) { parts.push({ a: [P[0], P[1], 0], b: [P[0], P[1], -H * 0.85], r: 0.006, structure: true }); continue; }
      parts.push({ a: [P[0] - c.n[0] * footLen, P[1] - c.n[1] * footLen, 0.012], b: [P[0] + c.n[0] * footLen, P[1] + c.n[1] * footLen, 0.012], r: 0.012, structure: true });
      parts.push({ a: [P[0], P[1], 0.012], b: [P[0], P[1], H * 0.85], r: 0.01, structure: true });
    }
    return down ? `${feet} hangers` : `${feet} feet ${Math.round(footLen * 200)} cm long`;
  }
  const flatField = { data: new Float32Array(4).fill(0.5), nx: 2, nz: 2 };

  // ------------------------------------------------------------ growth rules
  // Each rule turns genes (numbers in [0, 1]) into capsules at an arbitrary scale.
  // `down` is true for hanging forms (growth points downwards).
  // Felt sheet thickness is a manufacturing value: 15-80 mm, kept as is when a form is scaled.
  const sheetT = (u) => lerp(0.015, 0.08, u);
  const mmT = (u) => Math.round(sheetT(u) * 1000) + " mm";
  const isSheet = (p) => p.c !== undefined;
  const isTri = (p) => p.tri !== undefined;
  // Places: "standing" a single freestanding module, "hanging" a ceiling piece, "wall" a double-sided
  // partition standing on the ground, "wallhung" the same partition hung from the ceiling.
  const onGround = (m) => m === "standing" || m === "wall", hung = (m) => m === "hanging" || m === "wallhung", isWall = (m) => m === "wall" || m === "wallhung";
  const simMount = (m) => (m === "wall" ? "standing" : m === "wallhung" ? "hanging" : m);
  // Apply an affine point map to every part (sheets: centre and edge vectors; triangles: corners).
  function mapParts(parts, fn) {
    const out = parts.map((p) => {
      if (isTri(p)) return { ...p, tri: p.tri.map(fn) };
      if (isSheet(p)) { const c = fn(p.c); return { ...p, c, u: sub(fn(add(p.c, p.u)), c), v: sub(fn(add(p.c, p.v)), c) }; }
      return { ...p, a: fn(p.a), b: fn(p.b) };
    });
    for (const k of ["label", "field", "cutPlan", "stripPlan"]) if (parts[k] !== undefined) out[k] = parts[k];
    return out;
  }

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
    // ---- sound-grown (after Biothing's Genware): a cellular felt membrane grown from a sound field.
    membrane: {
      group: "sound", direct: true, name: "Sound-grown membrane", imitates: "Biothing's Genware: sound patterns as structural fields, cellular membranes",
      idea: "Sound draws the structure. Where the air moves most, cells grow smaller, thicker and closer; where it is calm, they open up. Each round the field is lowered where the membrane already absorbs, so growth spreads to where sound is still strong.",
      genes: ["plan", "curve", "cellMin", "cellMax", "contrast", "open", "thickMin", "thickMax", "fold", "louvre", "coupling", "rhythm", "feet", "foot"],
      build(g, down, ctx) {
        const env = ctx.env, L = env.width, H = env.height, field = ctx.field || flatField;
        const plan = Math.floor(g[0] * 3.999), curv = g[1], rMin = lerp(0.05, 0.12, g[2]), rMax = lerp(0.14, 0.3, g[3]), gamma = lerp(0.5, 2.5, g[4]);
        const open = g[5], tMin = lerp(0.01, 0.03, g[6]), tMax = lerp(0.03, 0.08, g[7]), fold = lerp(0, 0.15, g[8]), louvre = lerp(0, 1.2, g[9]);
        const coupling = lerp(0, 2, g[10]), rhythmGain = g[11], feet = 2 + Math.floor(g[12] * 2.999), footLen = lerp(0.15, 0.45, g[13]);
        const rand = rng(seedOf(g)), { nx, nz } = field, F0 = Float32Array.from(field.data);
        const bpm = (field.rhythm && field.rhythm.bpm) || 0, strength = (field.rhythm && field.rhythm.strength) || 0;
        const period = bpm ? H / Math.max(2, Math.min(6, bpm / 30)) : H;
        const curve = planCurve(plan, curv, L), at = (Fm, s2, z) => sampleField(Fm, nx, nz, (s2 / L + 0.5), z / H);
        const makeTiles = (Fm) => {
          const tiles = [];
          for (let a = 0; a < 2600 && tiles.length < 700; a++) {
            const m = rMax * 0.6, s2 = (rand() - 0.5) * (L - 2 * m), z = 0.06 + m + rand() * (H - 0.08 - 2 * m), f = Math.pow(Math.max(0, Math.min(1, at(Fm, s2, z))), gamma);
            const r = lerp(rMax, rMin, f) * (1 + rhythmGain * strength * 0.4 * Math.sin((2 * Math.PI * z) / period));
            if (tiles.some((o) => Math.hypot(o.s - s2, o.z - z) < (o.r + r) / 2)) continue;
            tiles.push({ s: s2, z, r, f });
          }
          return tiles;
        };
        const normalise = (arr) => { let lo = Infinity, hi = -Infinity; for (const x of arr) { lo = Math.min(lo, x); hi = Math.max(hi, x); } return arr.map((x) => (hi > lo ? (x - lo) / (hi - lo) : 0.5)); };
        let Fm = F0, tiles = makeTiles(Fm);
        for (let round = 1; round < 3 && coupling > 0; round++) {
          const M = new Float32Array(nx * nz);
          for (const t of tiles) {
            const i = Math.round((t.s / L + 0.5) * (nx - 1)), k = Math.round((t.z / H) * (nz - 1));
            for (let dk = -1; dk <= 1; dk++) for (let di = -1; di <= 1; di++) { const ii = i + di, kk = k + dk; if (ii >= 0 && ii < nx && kk >= 0 && kk < nz) M[ii + nx * kk] += (1 - open * (1 - t.f)) * (di || dk ? 0.5 : 1); }
          }
          const mmax = Math.max(1e-6, ...M);
          Fm = normalise(F0.map((x, i) => x * Math.exp((-coupling * M[i]) / mmax)));
          tiles = makeTiles(Fm);
        }
        const parts = [], sgn = down ? -1 : 1, zOf = (z) => (down ? z - H : z);
        const centres = tiles.map((t) => {
          const c = curve(t.s), solid = 1 - open * (1 - t.f), half = Math.max(0.015, t.r * 0.5 * solid);
          const dfds = (at(Fm, t.s + 0.02, t.z) - at(Fm, t.s - 0.02, t.z)) / 0.04, ang = louvre * Math.max(-1, Math.min(1, dfds * L * 0.25));
          const off = fold * (t.f - 0.5) * 2, ud = [c.t[0] * Math.cos(ang) + c.n[0] * Math.sin(ang), c.t[1] * Math.cos(ang) + c.n[1] * Math.sin(ang), 0];
          const ctr = [c.p[0] + c.n[0] * off, c.p[1] + c.n[1] * off, zOf(t.z)];
          parts.push({ c: ctr, u: mul(ud, half), v: [0, 0, half], t: lerp(tMin, tMax, t.f) });
          return ctr;
        });
        // Felt stitches between neighbouring cells hold the membrane together.
        centres.forEach((a, i) => {
          const near = centres.map((b, j) => [j, Math.hypot(...sub(a, b))]).filter(([j, d]) => j > i && d < (tiles[i].r + tiles[j].r) * 0.75).sort((x, y) => x[1] - y[1]).slice(0, 2);
          for (const [j] of near) parts.push({ a, b: centres[j], r: 0.004 });
        });
        parts.label = footings(parts, curve, L, H, feet, footLen, down);
        parts.field = Fm;
        return parts;
      },
      describe(g) {
        const plans = ["straight", "arched", "S-curved", "zigzag"];
        return `A ${plans[Math.floor(g[0] * 3.999)]} cellular membrane, cells ${Math.round(lerp(5, 12, g[2]))}-${Math.round(lerp(14, 30, g[3]))} cm, felt ${Math.round(lerp(10, 30, g[6]))}-${Math.round(lerp(30, 80, g[7]))} mm${g[10] > 0.2 ? ", grown over three rounds with the sound" : ""}.`;
      },
    },
    // ---- after Biothing's New Museum ceiling: one flat laser-cut felt sheet whose cut pattern is drawn
    // by an interference field. Flaps cut on three sides open out of the sheet, most where the air
    // moves most. The felt's own stiffness decides how far a flap really stays open (material computation).
    cutsheet: {
      group: "sound", direct: true, name: "Interference cut-sheet", imitates: "Biothing's laser-cut ceiling at the New Museum: wave interference cut into flat sheets",
      idea: "One flat felt sheet, cut only in 2D. Where sound is strong the flaps are larger and open further, so sound reaches both faces and more felt edge; where it is calm the sheet stays whole. Gravity and the felt's stiffness finish the shape.",
      genes: ["plan", "curve", "cell", "cut", "open", "contrast", "sheet", "alternate", "hinge", "feet", "foot"],
      build(g, down, ctx) {
        const env = ctx.env, L = env.width, H = env.height, field = ctx.field || flatField, mat = ctx.material || { density: 120 };
        const plan = Math.floor(g[0] * 3.999), a = lerp(0.08, 0.25, g[2]), cutMax = lerp(0.5, 0.92, g[3]), open = lerp(0.2, 1.4, g[4]), gamma = lerp(0.5, 2.5, g[5]);
        const t = lerp(0.006, 0.02, g[6]), alternate = g[7] > 0.5, side = g[8] > 0.5, feet = 2 + Math.floor(g[9] * 2.999), footLen = lerp(0.15, 0.45, g[10]);
        const curve = planCurve(plan, g[1], L), base = 0.06, nxc = Math.max(1, Math.round(L / a)), nzc = Math.max(1, Math.round((H - base) / a));
        const cw = L / nxc, ch = (H - base) / nzc, zOf = (z) => (down ? z - H : z), parts = [], cells = [];
        // Felt as a cantilever under its own weight: end slope 2·rho·g·l³ / (E·t²), with E rising with density
        // (about 2 MPa for 200 kg/m³ PET felt). Top-hinged flaps sag closed; side hinges carry the load in plane.
        const E = 50 * mat.density * mat.density, sag = (l) => (2 * mat.density * 9.81 * l * l * l) / (E * t * t);
        const sheet = (s0, s1, z0, z1, extra) => {
          if (s1 - s0 < 1e-3 || z1 - z0 < 1e-3) return;
          const c = curve((s0 + s1) / 2);
          parts.push({ c: [c.p[0], c.p[1], zOf((z0 + z1) / 2)], u: [c.t[0] * (s1 - s0) / 2, c.t[1] * (s1 - s0) / 2, 0], v: [0, 0, (z1 - z0) / 2], t, ...extra });
        };
        sheet(-L / 2, L / 2, 0, base);
        let cutLength = 0, droopSum = 0, flaps = 0;
        for (let k = 0; k < nzc; k++) for (let i = 0; i < nxc; i++) {
          const s0 = -L / 2 + i * cw, z0 = base + k * ch, sc = s0 + cw / 2, zc = z0 + ch / 2;
          const f = Math.pow(Math.max(0, Math.min(1, sampleField(field.data, field.nx, field.nz, sc / L + 0.5, zc / H))), gamma);
          const frac = cutMax * f;
          if (frac < 0.2) { sheet(s0, s0 + cw, z0, z0 + ch); cells.push({ s: sc, z: zc, w: 0, h: 0 }); continue; }
          const w = cw * frac, h = ch * frac, wanted = open * f, lever = side ? 0 : h;
          const ang = wanted / (1 + sag(lever) / Math.max(wanted, 0.05)), dir = alternate && (i + k) % 2 ? -1 : 1;
          // the frame around the hole
          sheet(s0, s0 + cw, z0, zc - h / 2); sheet(s0, s0 + cw, zc + h / 2, z0 + ch);
          sheet(s0, sc - w / 2, zc - h / 2, zc + h / 2); sheet(sc + w / 2, s0 + cw, zc - h / 2, zc + h / 2);
          // the flap, turned out of the sheet about its hinge
          const c = curve(sc), n = [c.n[0] * dir, c.n[1] * dir];
          if (side) {
            const hc = curve(sc - w / 2), d = [c.t[0] * Math.cos(ang) + n[0] * Math.sin(ang), c.t[1] * Math.cos(ang) + n[1] * Math.sin(ang), 0];
            parts.push({ c: [hc.p[0] + d[0] * w / 2, hc.p[1] + d[1] * w / 2, zOf(zc)], u: mul(d, w / 2), v: [0, 0, h / 2], t });
          } else {
            const d = [n[0] * Math.sin(ang), n[1] * Math.sin(ang), -Math.cos(ang)], top = zOf(zc + h / 2);
            parts.push({ c: [c.p[0] + d[0] * h / 2, c.p[1] + d[1] * h / 2, top + d[2] * h / 2], u: [c.t[0] * w / 2, c.t[1] * w / 2, 0], v: mul(d, h / 2), t });
          }
          cutLength += side ? 2 * w + h : w + 2 * h; droopSum += 1 - ang / Math.max(wanted, 1e-6); flaps++;
          cells.push({ s: sc, z: zc, w, h, open: ang, side, dir });
        }
        parts.label = footings(parts, curve, L, H, feet, footLen, down);
        parts.cutPlan = { L, H, t, base, side, cells, cutLength, flaps, droop: flaps ? droopSum / flaps : 0 };
        return parts;
      },
      describe(g) {
        const plans = ["flat", "arched", "S-curved", "zigzag"];
        return `A ${plans[Math.floor(g[0] * 3.999)]} ${Math.round(lerp(6, 20, g[6]))} mm felt sheet cut into ${Math.round(lerp(8, 25, g[2]))} cm cells, flaps ${g[8] > 0.5 ? "hinged at the side" : "hinged at the top"}${g[7] > 0.5 ? ", opening to both faces" : ""}.`;
      },
    },
  };

  // ======================================================= pattern lab ======
  // Generators for modules that can go on a wall, hang from the ceiling or stand on the floor, with a
  // parameter schema (`params`) so the designer can drive every value by hand, lock it, or give
  // evolution a range to explore (a graph like Grasshopper or nTop). All are graded by a field.
  // Geometry is built in a "wall frame": the face is x (across, centred) by z (up, 0..Hm), the depth
  // grows out along +y from y = 0. mountPattern then puts it on the wall, the ceiling or the floor.
  const opt = (label, options, node) => ({ label, options, node });
  const num = (label, min, max, unit, node, round) => ({ label, min, max, unit, node, round });
  function paramValue(rule, g, i) {
    const p = RULES[rule].params[i], x = Math.max(0, Math.min(1, g[i]));
    if (p.options) return Math.min(p.options.length - 1, Math.floor(x * p.options.length));
    const v = lerp(p.min, p.max, x);
    return p.round ? Math.round(v) : v;
  }
  const clamp01 = (x) => Math.max(0, Math.min(1, x));
  function frameOf(env, mount) {
    const W = env.width, Hm = mount === "hanging" ? env.width : env.height - 0.1;
    return { W, Hm };
  }
  function fieldFn(field, W, Hm, gamma) {
    const F = field || flatField;
    return (x, z) => Math.pow(clamp01(sampleField(F.data, F.nx, F.nz, x / W + 0.5, z / Hm)), gamma);
  }
  function mountPattern(parts, mount, { W, Hm, depth, mirror = true, feet = 2, footLen = 0.3 }) {
    if (mount === "hanging") {
      const drop = 0.25, out = mapParts(parts, (q) => [q[0], q[2] - Hm / 2, -drop - q[1]]);
      for (const sx of [-1, 1]) for (const sy of [-1, 1]) out.push({ a: [(sx * W) / 3, (sy * Hm) / 3, 0], b: [(sx * W) / 3, (sy * Hm) / 3, -drop], r: 0.003, structure: true });
      out.label = "4 hanging cables";
      return out;
    }
    // freestanding module or partition: both faces work; on feet, or hung from the ceiling on cables
    let out = mirror ? mapParts(parts, (q) => [q[0], q[1] + 0.0045, q[2]]) : mapParts(parts, (q) => [q[0], q[1] - depth / 2, q[2]]);
    if (mirror) out = out.concat(mapParts(out, (q) => [q[0], -q[1], q[2]]));
    const z0 = mount === "wallhung" ? -Hm - 0.25 : 0.1;
    out = mapParts(out, (q) => [q[0], q[1], q[2] + z0]);
    if (mirror) out.push({ c: [0, 0, z0 + Hm / 2], u: [W / 2, 0, 0], v: [0, 0, Hm / 2], t: 0.009, structure: true });
    if (mount === "wallhung") {
      for (const sx of [-1, 1]) out.push({ a: [(sx * W) / 3, 0, 0], b: [(sx * W) / 3, 0, -0.25], r: 0.003, structure: true });
      out.label = "2 hanging cables";
    } else out.label = footings(out, planCurve(0, 0, W), W, Hm + 0.1, feet, footLen, false);
    return out;
  }
  // Cell boundaries along a length, each cell as large as sizeAt(fraction) says, stretched to fit exactly.
  function warp1D(len, sizeAt) {
    const b = [0];
    let x = 0;
    while (x < len - 1e-9 && b.length < 400) { x += Math.max(0.01, sizeAt(Math.min(1, x / len))); b.push(x); }
    if (b.length < 2) b.push(len);
    const k = len / b[b.length - 1];
    return b.map((v) => v * k);
  }
  const atWarp = (b, u) => { const i = Math.max(0, Math.min(b.length - 2, Math.floor(u))); return b[i] + (b[i + 1] - b[i]) * (u - i); };
  const invWarp = (b, x) => { let i = 0; while (i < b.length - 2 && b[i + 1] < x) i++; return i + (x - b[i]) / (b[i + 1] - b[i] || 1); };

  // 2D Voronoi cells in a rectangle (half-plane clipping, nearest seeds first).
  function clipPoly(poly, nx, nz, c) {
    const out = [];
    for (let i = 0; i < poly.length; i++) {
      const a = poly[i], b = poly[(i + 1) % poly.length], da = nx * a[0] + nz * a[1] - c, db = nx * b[0] + nz * b[1] - c;
      if (da <= 0) out.push(a);
      if ((da < 0) !== (db < 0) && da !== db) { const t = da / (da - db); out.push([a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t]); }
    }
    return out;
  }
  function voronoiCells(pts, x0, x1, z0, z1) {
    return pts.map((p, i) => {
      let poly = [[x0, z0], [x1, z0], [x1, z1], [x0, z1]];
      const others = pts.map((q, j) => [j, (q[0] - p[0]) ** 2 + (q[1] - p[1]) ** 2]).filter(([j]) => j !== i).sort((a, b) => a[1] - b[1]);
      for (const [j, d2] of others) {
        let R2 = 0;
        for (const v of poly) R2 = Math.max(R2, (v[0] - p[0]) ** 2 + (v[1] - p[1]) ** 2);
        if (d2 > 4 * R2) break;
        const q = pts[j];
        poly = clipPoly(poly, q[0] - p[0], q[1] - p[1], (q[0] * q[0] + q[1] * q[1] - p[0] * p[0] - p[1] * p[1]) / 2);
        if (poly.length < 3) break;
      }
      return poly;
    });
  }
  // Seeds over the face: relaxed (spacing follows the field), honeycomb (graded hex rows) or a sunflower spiral.
  function faceSeeds(kind, W, Hm, sizeAt, rand) {
    const pts = [], inRect = (x, z) => x > -W / 2 && x < W / 2 && z > 0 && z < Hm;
    if (kind === 1) {
      const rows = warp1D(Hm, (v) => sizeAt(0, v * Hm) * 0.87);
      rows.forEach((z, j) => {
        if (j === rows.length - 1) return;
        const zc = (z + rows[j + 1]) / 2, cols = warp1D(W, (u) => sizeAt(-W / 2 + u * W, zc));
        for (let i = 0; i < cols.length - 1; i++) { const x = -W / 2 + (cols[i] + cols[i + 1]) / 2 + (j % 2 ? (cols[i + 1] - cols[i]) / 2 : 0); if (inRect(x, zc)) pts.push([x, zc]); }
      });
      return pts;
    }
    if (kind === 2) {
      const golden = Math.PI * (3 - Math.sqrt(5)), cx = 0, cz = Hm / 2;
      let A = 0;
      for (let k = 0; k < 600; k++) {
        const r = Math.sqrt(A / Math.PI), th = k * golden, x = cx + r * Math.cos(th), z = cz + r * Math.sin(th);
        if (r > Math.hypot(W, Hm) / 2) break;
        const s = sizeAt(Math.max(-W / 2, Math.min(W / 2, x)), Math.max(0, Math.min(Hm, z)));
        A += s * s * 0.9;
        if (inRect(x, z)) pts.push([x, z]);
      }
      return pts;
    }
    for (let t = 0; t < 5000 && pts.length < 400; t++) {
      const x = (rand() - 0.5) * W, z = rand() * Hm, s = sizeAt(x, z);
      if (pts.every((q) => Math.hypot(q[0] - x, q[1] - z) > s * 0.5 + sizeAt(q[0], q[1]) * 0.5)) pts.push([x, z]);
    }
    return pts;
  }
  // A face quad (parallelogram) through four points given in (x, y, z), as a sheet of thickness t.
  const quad = (A, B, C, D, t, extra) => ({ c: mul(add(add(A, B), add(C, D)), 0.25), u: mul(add(sub(B, A), sub(C, D)), 0.25), v: mul(add(sub(D, A), sub(C, B)), 0.25), t, ...extra });

  // Unit cells of strut lattices as bonds in the unit cube.
  const UNIT = (() => {
    const corners = [];
    for (const a of [0, 1]) for (const b of [0, 1]) for (const c of [0, 1]) corners.push([a, b, c]);
    const faces = [[0.5, 0.5, 0], [0.5, 0.5, 1], [0.5, 0, 0.5], [0.5, 1, 0.5], [0, 0.5, 0.5], [1, 0.5, 0.5]];
    const make = (pts, L) => { const out = []; for (let i = 0; i < pts.length; i++) for (let j = i + 1; j < pts.length; j++) if (Math.abs(Math.hypot(...sub(pts[i], pts[j])) - L) < 1e-6) out.push([pts[i], pts[j]]); return out; };
    return [make(corners, 1), make([...corners, [0.5, 0.5, 0.5]], Math.sqrt(3) / 2), make([...corners, ...faces], Math.SQRT1_2),
      make([...corners, ...faces, [0.25, 0.25, 0.25], [0.25, 0.75, 0.75], [0.75, 0.25, 0.75], [0.75, 0.75, 0.25]], Math.sqrt(3) / 4)];
  })();
  // Struts of a lattice on warped cells. map(u, v, w) gives a point in metres from cell coordinates;
  // keep(p) can drop struts (an envelope); radius(p) grades them. Type 4 is a honeycomb with rungs.
  function latticeStruts(type, nu, nv, nw, map, keep, radius) {
    const seen = new Set(), out = [], key = (p) => p.map((x) => Math.round(x * 4)).join(",");
    const strut = (P, Q) => {
      const k1 = key(P), k2 = key(Q), k = k1 < k2 ? k1 + "|" + k2 : k2 + "|" + k1;
      if (seen.has(k)) return;
      seen.add(k);
      const a = map(...P), b = map(...Q);
      if (keep && (!keep(a) || !keep(b))) return;
      out.push({ a, b, r: radius(mul(add(a, b), 0.5)), structure: true });
    };
    if (type === 4) {
      const s3 = Math.sqrt(3), verts = new Map();
      for (let j = -1; j <= Math.ceil(nw / 1.5) + 1; j++) for (let i = -1; i <= Math.ceil(nu / s3) + 1; i++) {
        const cx = s3 * (i + (j & 1) / 2), cz = 1.5 * j;
        for (let e = 0; e < 6; e++) {
          const a1 = Math.PI / 6 + (e * Math.PI) / 3, a2 = a1 + Math.PI / 3;
          const p = [(cx + Math.cos(a1)) / s3 * 1, (cz + Math.sin(a1)) / 1.5], q = [(cx + Math.cos(a2)) / s3, (cz + Math.sin(a2)) / 1.5];
          if ([p, q].some(([u, w]) => u < -1e-6 || u > nu + 1e-6 || w < -1e-6 || w > nw + 1e-6)) continue;
          for (let v = 0; v <= nv; v++) strut([p[0], v, p[1]], [q[0], v, q[1]]);
          verts.set(key([p[0], 0, p[1]]), p);
        }
      }
      for (const p of verts.values()) for (let v = 0; v < nv; v++) strut([p[0], v, p[1]], [p[0], v + 1, p[1]]);
      return out;
    }
    const bonds = UNIT[type];
    for (let k = 0; k < nw; k++) for (let j = 0; j < nv; j++) for (let i = 0; i < nu; i++) for (const [A, B] of bonds) strut([i + A[0], j + A[1], k + A[2]], [i + B[0], j + B[1], k + B[2]]);
    return out;
  }
  // Triply periodic minimal surfaces (value 0 is the surface) and their largest value, in cell coordinates.
  const TPMS = [
    { f: (X, Y, Z) => Math.sin(X) * Math.cos(Y) + Math.sin(Y) * Math.cos(Z) + Math.sin(Z) * Math.cos(X), max: 1.5 },
    { f: (X, Y, Z) => Math.cos(X) + Math.cos(Y) + Math.cos(Z), max: 3 },
    { f: (X, Y, Z) => Math.sin(X) * Math.sin(Y) * Math.sin(Z) + Math.sin(X) * Math.cos(Y) * Math.cos(Z) + Math.cos(X) * Math.sin(Y) * Math.cos(Z) + Math.cos(X) * Math.cos(Y) * Math.sin(Z), max: 1.42 },
  ];
  // Marching tetrahedra over a grid: triangles of the surface value(x, y, z) = 0.
  function isoTriangles(value, x0, y0, z0, hx, hy, hz, nx, ny, nz, keepTri) {
    const V = new Float32Array((nx + 1) * (ny + 1) * (nz + 1)), id = (i, j, k) => i + (nx + 1) * (j + (ny + 1) * k), P = (i, j, k) => [x0 + i * hx, y0 + j * hy, z0 + k * hz];
    for (let k = 0; k <= nz; k++) for (let j = 0; j <= ny; j++) for (let i = 0; i <= nx; i++) V[id(i, j, k)] = value(...P(i, j, k));
    const C = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]];
    const TETS = [[0, 5, 1, 6], [0, 1, 2, 6], [0, 2, 3, 6], [0, 3, 7, 6], [0, 7, 4, 6], [0, 4, 5, 6]], out = [];
    const cut = (pa, va, pb, vb) => { const t = va / (va - vb); return [pa[0] + (pb[0] - pa[0]) * t, pa[1] + (pb[1] - pa[1]) * t, pa[2] + (pb[2] - pa[2]) * t]; };
    for (let k = 0; k < nz; k++) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
      const cv = C.map((c) => V[id(i + c[0], j + c[1], k + c[2])]);
      if (cv.every((v) => v < 0) || cv.every((v) => v >= 0)) continue;
      const cp = C.map((c) => P(i + c[0], j + c[1], k + c[2]));
      for (const T of TETS) {
        const ins = T.filter((q) => cv[q] < 0), outs = T.filter((q) => cv[q] >= 0);
        if (!ins.length || !outs.length) continue;
        if (ins.length === 1 || outs.length === 1) {
          const lone = ins.length === 1 ? ins[0] : outs[0], rest = ins.length === 1 ? outs : ins;
          out.push(rest.map((q) => cut(cp[lone], cv[lone], cp[q], cv[q])));
        } else {
          const [a, b] = ins, [c, d] = outs, p1 = cut(cp[a], cv[a], cp[c], cv[c]), p2 = cut(cp[a], cv[a], cp[d], cv[d]), p3 = cut(cp[b], cv[b], cp[d], cv[d]), p4 = cut(cp[b], cv[b], cp[c], cv[c]);
          out.push([p1, p2, p3], [p1, p3, p4]);
        }
      }
    }
    return keepTri ? out.filter(keepTri) : out;
  }

  const CELL = "Cell", RESP = "Response to the field", HOLD = "Holdfast (freestanding)";
  const holdParams = [opt("Feet", ["2 feet", "3 feet", "4 feet"], HOLD), num("Foot length", 0.15, 0.45, "cm", HOLD)];
  const LAB = {
    voronoi: {
      group: "pattern", lab: true, direct: true, mounts: true, name: "Voronoi cells", imitates: "cracked mud, giraffe skin, dragonfly wings, foam, the cells of plant tissue",
      idea: "The face is divided into cells, each the area closest to one seed. Where sound is strong the cells grow smaller and deeper, so there are more felt walls and edges there.",
      params: [opt("Seeds", ["relaxed", "honeycomb", "sunflower spiral"], CELL), opt("Style", ["open cells", "funnels", "stepped web"], CELL),
        num("Cell size where sound is strong", 0.05, 0.2, "cm", RESP), num("Cell size where calm", 0.12, 0.4, "cm", RESP),
        num("Depth where calm", 0.02, 0.15, "cm", RESP), num("Depth where sound is strong", 0.05, 0.3, "cm", RESP),
        num("Opening (funnel or web hole)", 0.1, 0.8, "%", CELL), num("Felt thickness", 0.005, 0.02, "mm", CELL), opt("Back", ["open back", "felt back"], CELL),
        num("Contrast", 0.5, 3, "", RESP), ...holdParams],
      build(g, down, ctx) {
        const v = (i) => paramValue("voronoi", g, i), mount = ctx.mount, { W, Hm } = frameOf(ctx.env, mount), rand = rng(seedOf(g));
        const style = v(1), cMin = v(2), cMax = Math.max(v(2), v(3)), dLo = v(4), dHi = Math.min(ctx.env.depth || 0.3, Math.max(v(4), v(5))), open = v(6), t = v(7), f = fieldFn(ctx.field, W, Hm, v(9));
        const sizeAt = (x, z) => lerp(cMax, cMin, f(x, z));
        const seeds = faceSeeds(v(0), W, Hm, sizeAt, rand), cells = voronoiCells(seeds, -W / 2, W / 2, 0, Hm), parts = [];
        const P3 = (q, y) => [q[0], y, q[1]];
        const edges = new Map();
        cells.forEach((poly, i) => {
          if (poly.length < 3) return;
          const s = seeds[i], fi = f(s[0], s[1]), d = lerp(dLo, dHi, fi), n = poly.length;
          const hole = poly.map((q) => [s[0] + (q[0] - s[0]) * open, s[1] + (q[1] - s[1]) * open]);
          for (let e = 0; e < n; e++) {
            const a = poly[e], b = poly[(e + 1) % n];
            if (style === 0) {
              const k = [a, b].map((q) => q.map((x) => Math.round(x * 2000)).join(",")).sort().join("|");
              if (!edges.has(k) || edges.get(k).d < d) edges.set(k, { a, b, d });
            } else if (style === 1) parts.push(quad(P3(hole[e], 0), P3(hole[(e + 1) % n], 0), P3(b, d), P3(a, d), t));
            else for (let L = 0; L < 3; L++) {
              const sc = 1 - ((1 - open) * (L + 1)) / 3, hA = [s[0] + (a[0] - s[0]) * sc, s[1] + (a[1] - s[1]) * sc], hB = [s[0] + (b[0] - s[0]) * sc, s[1] + (b[1] - s[1]) * sc];
              const y = d * (1 - (L + 0.5) / 3);
              parts.push(quad(P3(a, y), P3(b, y), P3(hB, y), P3(hA, y), d / 3));
            }
          }
        });
        for (const { a, b, d } of edges.values()) parts.push(quad(P3(a, 0), P3(b, 0), P3(b, d), P3(a, d), t));
        if (v(8) === 1 || style === 1) parts.push({ c: [0, -t / 2, Hm / 2], u: [W / 2, 0, 0], v: [0, 0, Hm / 2], t });
        return mountPattern(parts, mount, { W, Hm, depth: dHi, feet: 2 + v(10), footLen: v(11) });
      },
      describe(g) { const v = (i) => paramValue("voronoi", g, i); return `${["Relaxed", "Honeycomb", "Sunflower"][v(0)]} Voronoi ${["open cells", "funnels", "stepped web"][v(1)]}, cells ${Math.round(v(2) * 100)}-${Math.round(Math.max(v(2), v(3)) * 100)} cm, ${Math.round(v(4) * 100)}-${Math.round(Math.max(v(4), v(5)) * 100)} cm deep, ${Math.round(v(7) * 1000)} mm felt.`; },
    },
    petals: {
      group: "pattern", lab: true, direct: true, mounts: true, name: "Felt petals", imitates: "sunflower heads, artichokes, succulents, scales of a cone",
      idea: "Felt petals on sunflower spirals, larger towards the edge. Where sound is strong they open further out of the face, so both faces of every petal catch sound.",
      params: [opt("Spiral centres", ["1 centre", "2 centres", "3 centres"], CELL), num("Petals", 60, 400, "", CELL, true),
        num("Petal length at the centre", 0.04, 0.15, "cm", CELL), num("Petal length at the edge", 0.08, 0.3, "cm", CELL), num("Width / length", 0.3, 0.9, "%", CELL),
        num("Opening where calm", 5, 45, "°", RESP), num("Opening where sound is strong", 30, 85, "°", RESP), num("Curl", 0, 60, "°", CELL),
        num("Felt thickness", 0.003, 0.012, "mm", CELL), opt("Back", ["open back", "felt back"], CELL), num("Contrast", 0.5, 3, "", RESP), ...holdParams],
      build(g, down, ctx) {
        const v = (i) => paramValue("petals", g, i), mount = ctx.mount, { W, Hm } = frameOf(ctx.env, mount), f = fieldFn(ctx.field, W, Hm, v(10));
        const nC = 1 + v(0), count = v(1), Lc = v(2), Le = Math.max(v(2), v(3)), ratio = v(4), oLo = v(5), oHi = Math.max(v(5), v(6)), curl = (v(7) * Math.PI) / 180, t = v(8);
        // spiral centres at the strongest places of the field, kept apart
        const cand = [];
        for (let i = 0; i < 9; i++) for (let k = 0; k < 9; k++) { const x = -W / 2 + ((i + 0.5) / 9) * W, z = ((k + 0.5) / 9) * Hm; cand.push([x, z, f(x, z)]); }
        cand.sort((a, b) => b[2] - a[2]);
        const centres = [];
        for (const c of cand) if (centres.length < nC && centres.every((q) => Math.hypot(q[0] - c[0], q[1] - c[1]) > Math.min(W, Hm) * 0.4)) centres.push(c);
        const golden = Math.PI * (3 - Math.sqrt(5)), R = Math.hypot(W, Hm) / (2 * Math.sqrt(nC)), per = Math.ceil(count / nC), parts = [], depthMax = [0];
        for (const C of centres) for (let k = 1; k <= per; k++) {
          const fr = Math.sqrt(k / per), x = C[0] + R * fr * Math.cos(k * golden), z = C[1] + R * fr * Math.sin(k * golden);
          if (x < -W / 2 || x > W / 2 || z < 0 || z > Hm) continue;
          if (centres.some((o) => o !== C && Math.hypot(o[0] - x, o[1] - z) < Math.hypot(C[0] - x, C[1] - z))) continue;
          const L = lerp(Lc, Le, fr), w = L * ratio, ang = ((lerp(oLo, oHi, f(x, z)) * Math.PI) / 180), dl = Math.hypot(x - C[0], z - C[1]) || 1;
          const dir = [(x - C[0]) / dl, 0, (z - C[1]) / dl], side = [dir[2], 0, -dir[0]], Y = [0, 1, 0];
          const e1 = add(mul(dir, Math.cos(ang)), mul(Y, Math.sin(ang))), a2 = Math.min(Math.PI * 0.95, ang + curl), e2 = add(mul(dir, Math.cos(a2)), mul(Y, Math.sin(a2)));
          const p0 = [x, t, z], p1 = add(p0, mul(e1, L / 2));
          parts.push({ c: add(p0, mul(e1, L / 4)), u: mul(e1, L / 4), v: mul(side, w / 2), t });
          parts.push({ c: add(p1, mul(e2, L / 4)), u: mul(e2, L / 4), v: mul(side, w * 0.3), t });
          depthMax[0] = Math.max(depthMax[0], L * Math.sin(ang));
        }
        if (v(9) === 1) parts.push({ c: [0, t / 2, Hm / 2], u: [W / 2, 0, 0], v: [0, 0, Hm / 2], t });
        return mountPattern(parts, mount, { W, Hm, depth: depthMax[0], feet: 2 + v(11), footLen: v(12) });
      },
      describe(g) { const v = (i) => paramValue("petals", g, i); return `${v(1)} felt petals ${Math.round(v(8) * 1000)} mm thick on ${1 + v(0)} spiral${v(0) ? "s" : ""}, ${Math.round(v(2) * 100)}-${Math.round(Math.max(v(2), v(3)) * 100)} cm long, opening ${Math.round(v(5))}-${Math.round(Math.max(v(5), v(6)))}°.`; },
    },
    strips: {
      group: "pattern", lab: true, direct: true, mounts: true, name: "Slotted felt strips", imitates: "the gills of a bracket fungus, honeycomb walls, an egg-crate",
      idea: "Flat felt strips slotted into each other like an egg-crate: a 3D lattice made only from 2D cuts. Strips are closer and deeper where sound is strong. Flat-packs for touring.",
      params: [opt("Pattern", ["square", "diagonal"], CELL), num("Spacing where sound is strong", 0.05, 0.2, "cm", RESP), num("Spacing where calm", 0.1, 0.35, "cm", RESP),
        num("Depth where calm", 0.03, 0.15, "cm", RESP), num("Depth where sound is strong", 0.06, 0.3, "cm", RESP), num("Felt thickness", 0.004, 0.015, "mm", CELL),
        num("Contrast", 0.5, 3, "", RESP), ...holdParams],
      build(g, down, ctx) {
        const v = (i) => paramValue("strips", g, i), mount = ctx.mount, { W, Hm } = frameOf(ctx.env, mount), f = fieldFn(ctx.field, W, Hm, v(6));
        const sLo = v(1), sHi = Math.max(v(1), v(2)), dLo = v(3), dHi = Math.min(ctx.env.depth || 0.3, Math.max(v(3), v(4))), t = v(5), diag = v(0) === 1;
        const col = (u) => { let s = 0; for (let k = 0; k < 8; k++) s += f(-W / 2 + u * W, ((k + 0.5) / 8) * Hm); return s / 8; };
        const row = (w) => { let s = 0; for (let k = 0; k < 8; k++) s += f(-W / 2 + ((k + 0.5) / 8) * W, w * Hm); return s / 8; };
        const depth = (x, z) => lerp(dLo, dHi, f(x, z)), parts = [], strips = [];
        const xs = warp1D(W, (u) => lerp(sHi, sLo, col(u))).map((x) => x - W / 2), zs = warp1D(Hm, (w) => lerp(sHi, sLo, row(w)));
        const seg = (A, B, name) => {
          const m = mul(add(A, B), 0.5), d = depth(m[0], m[1]);
          parts.push({ c: [m[0], d / 2, m[1]], u: [(B[0] - A[0]) / 2, 0, (B[1] - A[1]) / 2], v: [0, d / 2, 0], t });
          return d;
        };
        if (!diag) {
          xs.forEach((x, i) => { if (i === 0 || i === xs.length - 1) return; const hs = []; for (let k = 0; k < zs.length - 1; k++) hs.push(seg([x, zs[k]], [x, zs[k + 1]])); strips.push({ name: "V" + i, along: zs, heights: hs, slotsFrom: "back" }); });
          zs.forEach((z, k) => { const hs = []; for (let i = 0; i < xs.length - 1; i++) hs.push(seg([xs[i], z], [xs[i + 1], z])); strips.push({ name: "H" + (k + 1), along: xs.map((x) => x + W / 2), heights: hs, slotsFrom: "front" }); });
        } else {
          // two families of lines at +-45 degrees, cut where they leave the face
          const step = (sLo + sHi) / 2, inside = (q) => q[0] >= -W / 2 - 1e-9 && q[0] <= W / 2 + 1e-9 && q[1] >= -1e-9 && q[1] <= Hm + 1e-9;
          for (const sg of [1, -1]) for (let c = -W - Hm; c <= W + Hm; c += step) {
            const pts = [];
            for (let s2 = -W - Hm; s2 <= W + Hm; s2 += step) { const q = [(c + s2) / 2, sg * (s2 - c) / 2 + Hm / 2]; pts.push(q); }
            for (let k = 0; k < pts.length - 1; k++) if (inside(pts[k]) && inside(pts[k + 1])) seg(pts[k], pts[k + 1]);
          }
        }
        const out = mountPattern(parts, mount, { W, Hm, depth: dHi, feet: 2 + v(7), footLen: v(8) });
        if (!diag) out.stripPlan = { t, strips, W, Hm };
        return out;
      },
      describe(g) { const v = (i) => paramValue("strips", g, i); return `${v(0) ? "Diagonal" : "Square"} egg-crate of ${Math.round(v(5) * 1000)} mm felt strips, ${Math.round(v(1) * 100)}-${Math.round(Math.max(v(1), v(2)) * 100)} cm apart, ${Math.round(v(3) * 100)}-${Math.round(Math.max(v(3), v(4)) * 100)} cm deep.`; },
    },
    strutlattice: {
      group: "lattice", lab: true, direct: true, mounts: true, name: "Strut lattice", imitates: "spongy bone, the Venus flower basket sponge, radiolarians",
      idea: "A printed lattice of thin struts, cells smaller and struts thicker where sound is strong, its surface rising and falling with the field. The struts alone hardly absorb (their gaps are far too big to rub the air): felt balls in the cells do the absorbing.",
      params: [opt("Cell type", ["cubic", "body-centred", "octet", "diamond", "honeycomb"], CELL), num("Cell size where sound is strong", 0.04, 0.12, "cm", RESP), num("Cell size where calm", 0.06, 0.25, "cm", RESP),
        num("Strut radius where calm", 0.001, 0.004, "mm", RESP), num("Strut radius where sound is strong", 0.002, 0.008, "mm", RESP), num("Lattice depth", 0.05, 0.3, "cm", CELL),
        num("Surface wave", 0, 0.8, "%", RESP), num("Cells filled with felt", 0, 1, "%", RESP), num("Felt ball size", 0.3, 0.9, "%", CELL), num("Contrast", 0.5, 3, "", RESP), ...holdParams],
      build(g, down, ctx) { return buildLattice("strutlattice", g, ctx, false); },
      describe(g) { const v = (i) => paramValue("strutlattice", g, i); return `A ${["cubic", "body-centred", "octet", "diamond", "honeycomb"][v(0)]} lattice, cells ${Math.round(v(1) * 100)}-${Math.round(Math.max(v(1), v(2)) * 100)} cm, ${Math.round(v(5) * 100)} cm deep, ${Math.round(v(7) * 100)}% of cells filled with felt.`; },
    },
    tpms: {
      group: "lattice", lab: true, direct: true, mounts: true, name: "Minimal surface (TPMS)", imitates: "butterfly wing scales (gyroids), sea urchin skeletons, cell membranes",
      idea: "A continuous curved printed shell that splits space into two interwoven labyrinths, graded by the field. The shell itself barely absorbs; felt balls in one labyrinth do, densest where sound is strong.",
      params: [opt("Surface", ["gyroid", "Schwarz P", "diamond"], CELL), num("Cell size where sound is strong", 0.05, 0.12, "cm", RESP), num("Cell size where calm", 0.08, 0.25, "cm", RESP),
        num("Shell thickness where calm", 0.0008, 0.003, "mm", RESP), num("Shell thickness where sound is strong", 0.001, 0.005, "mm", RESP), num("Lattice depth", 0.06, 0.3, "cm", CELL),
        num("Surface wave", 0, 0.8, "%", RESP), num("Cells filled with felt", 0, 1, "%", RESP), num("Felt ball size", 0.3, 0.9, "%", CELL), num("Contrast", 0.5, 3, "", RESP), ...holdParams],
      build(g, down, ctx) { return buildLattice("tpms", g, ctx, true); },
      describe(g) { const v = (i) => paramValue("tpms", g, i); return `A ${["gyroid", "Schwarz P", "diamond"][v(0)]} shell, cells ${Math.round(v(1) * 100)}-${Math.round(Math.max(v(1), v(2)) * 100)} cm, ${Math.round(v(5) * 100)} cm deep, ${Math.round(v(7) * 100)}% of its labyrinth filled with felt.`; },
    },
    spatial: {
      group: "lattice", lab: true, direct: true, mounts: true, name: "Spatial lattice", imitates: "Biothing's robotically printed clouds and pergolas, coral heads, a tree flaring into its canopy",
      idea: "A printed strut lattice that fills a grown body: a column that flares into a canopy, or a cloud. The body grows widest where sound is strong; felt pockets inside do the absorbing.",
      params: [opt("Body", ["column into canopy", "cloud"], CELL), opt("Cell type", ["cubic", "body-centred", "octet", "diamond"], CELL), num("Cell size", 0.05, 0.16, "cm", CELL),
        num("Strut radius", 0.0015, 0.006, "mm", CELL), num("Base or waist", 0.15, 0.5, "%", CELL), num("Flare or bulge", 0, 1, "%", RESP), num("Twist", 0, 90, "°", CELL),
        num("Cells filled with felt", 0, 1, "%", RESP), num("Felt ball size", 0.3, 0.9, "%", CELL), num("Contrast", 0.5, 3, "", RESP), ...holdParams],
      build(g, down, ctx) {
        const v = (i) => paramValue("spatial", g, i), mount = ctx.mount, env = ctx.env, W = env.width, H = mount === "hanging" ? env.width : env.height;
        const f = fieldFn(ctx.field, W, H, v(9)), body = v(0), a = v(2), waist = v(4), flare = v(5), twist = (v(6) * Math.PI) / 180;
        const fz = (z) => { let s = 0; for (let k = 0; k < 6; k++) s += f(-W / 2 + ((k + 0.5) / 6) * W, z); return s / 6; };
        // radius of the body at height z (0 at the floor or ceiling attachment, H at the far end)
        const radius = (z) => {
          const u = z / H;
          if (body === 0) return (W / 2) * (waist + (1 - waist) * Math.pow(u, 2.2) * (0.5 + 0.5 * flare + 0.4 * flare * fz(z)));
          return (W / 2) * Math.max(0.05, waist + (1 - waist) * Math.sin(Math.PI * u) * (0.6 + 0.4 * flare * fz(z)));
        };
        const inside = (p) => { const z = p[2], u = z / H, rot = twist * u, x = p[0] * Math.cos(rot) + p[1] * Math.sin(rot), y = -p[0] * Math.sin(rot) + p[1] * Math.cos(rot); return z >= -1e-6 && z <= H + 1e-6 && Math.max(Math.abs(x), Math.abs(y)) <= radius(z) + 1e-6; };
        let cell = a, n = 0;
        for (let it = 0; it < 4; it++) { n = Math.ceil(W / cell) ** 2 * Math.ceil(H / cell); if (n * UNIT[Math.min(3, v(1))].length * 0.35 < 14000) break; cell *= 1.25; }
        const nu = Math.ceil(W / cell), nw = Math.ceil(H / cell), map = (u, vv, w) => [-W / 2 + (u / nu) * W, -W / 2 + (vv / nu) * W, (w / nw) * H];
        const parts = latticeStruts(Math.min(3, v(1)), nu, nu, nw, map, inside, () => v(3));
        const balls = [];
        for (let k = 0; k < nw; k++) for (let j = 0; j < nu; j++) for (let i = 0; i < nu; i++) { const c = map(i + 0.5, j + 0.5, k + 0.5); if (inside(c)) balls.push([c, f(c[0], c[2])]); }
        balls.sort((p, q) => q[1] - p[1]);
        for (const [c] of balls.slice(0, Math.round(balls.length * v(7)))) parts.push({ a: c, b: c, r: (cell / 2) * v(8) * 0.95 });
        let out = parts;
        if (isWall(mount)) { const k = Math.min(1, (env.depth || 0.3) / W); out = mapParts(out, (q) => [q[0], q[1] * k, q[2]]); }
        if (hung(mount)) out = mapParts(out, (q) => [q[0], q[1], -q[2]]);
        const plate = Math.max(0.25, radius(0) * 1.6), py = isWall(mount) ? Math.max(0.12, env.depth || 0.3) : plate;
        if (onGround(mount)) { out.push({ c: [0, 0, 0.01], u: [plate, 0, 0], v: [0, py, 0], t: 0.02, structure: true }); out.label = "a weighted base plate"; }
        else { out.push({ c: [0, 0, -0.01], u: [radius(0), 0, 0], v: [0, isWall(mount) ? py / 2 : radius(0), 0], t: 0.02, structure: true }); out.label = "a ceiling plate"; }
        return out;
      },
      describe(g) { const v = (i) => paramValue("spatial", g, i); return `A ${["column flaring into a canopy", "cloud"][v(0)]} filled with a ${["cubic", "body-centred", "octet", "diamond"][v(1)]} lattice of ${Math.round(v(2) * 100)} cm cells, ${Math.round(v(7) * 100)}% felt pockets.`; },
    },
  };
  // Module lattices (struts or a minimal surface) on cells warped by the field, conforming to a wavy surface.
  function buildLattice(rule, g, ctx, surface) {
    const v = (i) => paramValue(rule, g, i), mount = ctx.mount, { W, Hm } = frameOf(ctx.env, mount), f = fieldFn(ctx.field, W, Hm, v(9));
    const cLo = v(1), cHi = Math.max(v(1), v(2)), D = Math.min(ctx.env.depth || 0.3, v(5)), wave = v(6), infill = v(7), plug = v(8);
    const col = (u) => { let s = 0; for (let k = 0; k < 8; k++) s += f(-W / 2 + u * W, ((k + 0.5) / 8) * Hm); return s / 8; };
    const row = (w) => { let s = 0; for (let k = 0; k < 8; k++) s += f(-W / 2 + ((k + 0.5) / 8) * W, w * Hm); return s / 8; };
    let grow = 1, xs, zs, nv;
    for (let it = 0; it < 5; it++) {
      xs = warp1D(W, (u) => lerp(cHi, cLo, col(u)) * grow); zs = warp1D(Hm, (w) => lerp(cHi, cLo, row(w)) * grow);
      nv = Math.max(1, Math.round(D / ((cLo + cHi) / 2 / grow)));
      const cells = (xs.length - 1) * (zs.length - 1) * nv;
      if (surface ? cells < 2500 : cells * (v(0) === 4 ? 6 : UNIT[v(0)].length) * 0.4 < 12000) break;
      grow *= 1.2;
    }
    const nu = xs.length - 1, nw = zs.length - 1, Dsurf = (x, z) => D * (1 - wave * (1 - f(x, z)));
    const map = (u, vv, w) => { const x = atWarp(xs, u) - W / 2, z = atWarp(zs, w); return [x, (vv / nv) * Dsurf(x, z), z]; };
    let parts;
    if (!surface) parts = latticeStruts(v(0), nu, nv, nw, map, null, (p) => lerp(v(3), Math.max(v(3), v(4)), f(p[0], p[2])));
    else {
      const T = TPMS[v(0)], tLo = v(3), tHi = Math.max(v(3), v(4)), cell = (cLo + cHi) / 2 / grow, h = Math.max(cLo / grow / 4, Math.cbrt((W * Hm * D) / 7000));
      const gx = Math.ceil(W / h), gy = Math.max(2, Math.ceil(D / h)), gz = Math.ceil(Hm / h);
      const val = (x, y, z) => { const ds = Dsurf(x, z); return T.f(2 * Math.PI * invWarp(xs, x + W / 2), (2 * Math.PI * nv * y) / Math.max(ds, 1e-3), 2 * Math.PI * invWarp(zs, z)); };
      const tris = isoTriangles(val, -W / 2, 0, 0, W / gx, D / gy, Hm / gz, gx, gy, gz, (tr) => { const c = mul(add(add(tr[0], tr[1]), tr[2]), 1 / 3); return c[1] <= Dsurf(c[0], c[2]) + 1e-6; });
      parts = tris.map((tr) => { const c = mul(add(add(tr[0], tr[1]), tr[2]), 1 / 3); return { tri: tr, t: lerp(tLo, tHi, f(c[0], c[2])), structure: true }; });
      void cell;
    }
    // felt balls in the cells (strut lattices) or in one labyrinth (minimal surfaces), strongest field first
    const balls = [];
    for (let k = 0; k < nw; k++) for (let j = 0; j < nv; j++) for (let i = 0; i < nu; i++) {
      const c = map(i + 0.5, j + 0.5, k + 0.5), sz = Math.min(xs[i + 1] - xs[i], zs[k + 1] - zs[k], Dsurf(c[0], c[2]) / nv);
      if (surface && TPMS[v(0)].f(2 * Math.PI * (i + 0.5) + Math.PI / 2, 2 * Math.PI * (j + 0.5), 2 * Math.PI * (k + 0.5)) > 0 && v(0) === 1) continue;
      balls.push([c, f(c[0], c[2]), sz]);
    }
    balls.sort((p, q) => q[1] - p[1]);
    for (const [c, , sz] of balls.slice(0, Math.round(balls.length * infill))) parts.push({ a: c, b: c, r: (sz / 2) * plug * (surface ? 0.55 : 0.9) });
    return mountPattern(parts, mount, { W, Hm, depth: D, mirror: false, feet: 2 + v(10), footLen: v(11) });
  }
  // Cellular automata (Andrasek's first example, Wolfram's "new kind of science"): a 3D grid of cells,
  // each born or dying by how many of its neighbours are alive, generation after generation. Birth is
  // likelier where the field is strong. Parts come out in birth order, so replaying them replays the growth.
  LAB.cellular = {
    group: "pattern", lab: true, direct: true, mounts: true, name: "Cellular automata", imitates: "coral and lichen growing by accretion, frost dendrites, sponges",
    idea: "Each cell only looks at its 26 neighbours: with the right number alive it is born, with too many it dies. Run over many generations, these dumb local rules grow intricate porous bodies, densest where sound is strong.",
    params: [opt("Rule", ["coral (accretes, never dies)", "sponge (born 1-4, dies when crowded)", "crystal (grows from faces only)"], CELL), num("Cell size", 0.03, 0.1, "cm", CELL),
      num("Generations", 4, 40, "", CELL, true), num("Seeds", 1, 30, "", CELL, true), num("Birth chance", 0.05, 0.8, "%", RESP), num("Contrast", 0.5, 3, "", RESP),
      num("Depth", 0.05, 0.3, "cm", CELL), opt("Made as", ["felt blocks", "printed bonds", "printed bonds with felt"], CELL), num("Block size", 0.4, 1, "%", CELL), ...holdParams],
    build(g, down, ctx) {
      const v = (i) => paramValue("cellular", g, i), mount = ctx.mount, { W, Hm } = frameOf(ctx.env, mount), f = fieldFn(ctx.field, W, Hm, v(5)), rand = rng(seedOf(g) + 5);
      const rule = v(0), D = Math.min(ctx.env.depth || 0.3, v(6)), style = v(7);
      let a = v(1);
      while ((W / a) * (Hm / a) * (D / a) > 40000) a *= 1.15;
      const nx = Math.max(2, Math.round(W / a)), ny = Math.max(2, Math.round(D / a)), nz = Math.max(2, Math.round(Hm / a)), N = nx * ny * nz;
      const id = (i, j, k) => i + nx * (j + ny * k), born = new Int16Array(N).fill(-1), pos = (i, j, k) => [-W / 2 + (i + 0.5) * (W / nx), (j + 0.5) * (D / ny), (k + 0.5) * (Hm / nz)];
      const fCell = new Float32Array(N);
      for (let k = 0; k < nz; k++) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) { const p = pos(i, j, k); fCell[id(i, j, k)] = f(p[0], p[2]); }
      // seeds on the back layer at the strongest places of the field
      const back = [];
      for (let k = 0; k < nz; k++) for (let i = 0; i < nx; i++) back.push([i, k, fCell[id(i, 0, k)] + rand() * 0.15]);
      back.sort((p, q) => q[2] - p[2]);
      for (const [i, k] of back.slice(0, v(3))) born[id(i, 0, k)] = 0;
      const gens = v(2), chance = v(4);
      for (let gen = 1; gen <= gens; gen++) {
        const alive = Uint8Array.from(born, (b) => (b >= 0 ? 1 : 0)), next = [], dead = [];
        for (let k = 0; k < nz; k++) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
          let n26 = 0, n6 = 0;
          for (let dk = -1; dk <= 1; dk++) for (let dj = -1; dj <= 1; dj++) for (let di = -1; di <= 1; di++) {
            if (!di && !dj && !dk) continue;
            const ii = i + di, jj = j + dj, kk = k + dk;
            if (ii < 0 || jj < 0 || kk < 0 || ii >= nx || jj >= ny || kk >= nz || !alive[id(ii, jj, kk)]) continue;
            n26++; if (Math.abs(di) + Math.abs(dj) + Math.abs(dk) === 1) n6++;
          }
          const c = id(i, j, k);
          if (alive[c]) { if (rule === 1 && n26 > 12) dead.push(c); continue; }
          const ok = rule === 0 ? n26 >= 1 && n26 <= 2 : rule === 1 ? n26 >= 1 && n26 <= 4 : n6 === 1 && n26 <= 3;
          if (ok && rand() < chance * (0.15 + 0.85 * fCell[c])) next.push(c);
        }
        for (const c of next) born[c] = gen;
        for (const c of dead) born[c] = -2;
      }
      const cells = [];
      for (let c = 0; c < N; c++) if (born[c] >= 0) cells.push(c);
      cells.sort((p, q) => born[p] - born[q]);
      const coord = (c) => [c % nx, Math.floor(c / nx) % ny, Math.floor(c / (nx * ny))], parts = [], cw = W / nx, cd = D / ny, ch = Hm / nz, bs = v(8);
      for (const c of cells) {
        const [i, j, k] = coord(c), p = pos(i, j, k);
        if (style !== 1) parts.push({ c: p, u: [(cw / 2) * bs, 0, 0], v: [0, 0, (ch / 2) * bs], t: cd * bs });
        if (style !== 0) for (const [di, dj, dk] of [[1, 0, 0], [0, 1, 0], [0, 0, 1]]) {
          const ii = i + di, jj = j + dj, kk = k + dk;
          if (ii < nx && jj < ny && kk < nz && born[id(ii, jj, kk)] >= 0) parts.push({ a: p, b: pos(ii, jj, kk), r: Math.min(cw, ch) * 0.06, structure: true });
        }
      }
      if (!parts.length) parts.push({ c: pos(0, 0, 0), u: [cw / 2, 0, 0], v: [0, 0, ch / 2], t: cd });
      return mountPattern(parts, mount, { W, Hm, depth: D, mirror: style === 0, feet: 2 + v(9), footLen: v(10) });
    },
    describe(g) { const v = (i) => paramValue("cellular", g, i); return `A ${["coral", "sponge", "crystal"][v(0)]} cellular automaton, ${Math.round(v(1) * 100)} cm cells grown over ${v(2)} generations from ${v(3)} seeds, made as ${["felt blocks", "printed bonds", "printed bonds with felt"][v(7)]}.`; },
  };
  // Liquid bioform: blobs that flow into one smooth body (an implicit surface of blended spheres),
  // made as a thin shell of heat-formed felt and perforated where sound is strong, so sound gets inside.
  LAB.bioform = {
    group: "pattern", lab: true, direct: true, mounts: true, name: "Liquid bioform", imitates: "drops and splashes, cell clusters, seed pods, organic design's fluid bodies",
    idea: "Blobs placed along a flowing path melt into one smooth body. It is made as a thin shell of heat-formed felt; where sound is strong the shell opens into holes, so sound reaches the felt inside as well as outside.",
    params: [num("Blobs", 3, 24, "", CELL, true), num("Blob size", 0.12, 0.4, "cm", CELL), num("Flow (how far the path wanders)", 0, 1, "%", CELL), num("Melt (smoothness)", 0.3, 1, "%", CELL),
      num("Holes", 0, 120, "", RESP, true), num("Hole size where sound is strong", 0.03, 0.15, "cm", RESP), num("Felt shell thickness", 0.003, 0.012, "mm", CELL), num("Contrast", 0.5, 3, "", RESP),
      opt("Body", ["flowing blobs", "pulled from a slab"], CELL), ...holdParams],
    build(g, down, ctx) {
      if (paramValue("bioform", g, 8) === 1) return pulledBody(g, ctx);
      const v = (i) => paramValue("bioform", g, i), mount = ctx.mount, env = ctx.env, W = env.width, H = mount === "hanging" ? env.width : env.height - 0.05;
      const D = isWall(mount) ? env.depth : W, f = fieldFn(ctx.field, W, H, v(7)), rand = rng(seedOf(g) + 11), nb = v(0), R0 = v(1), flow = v(2), melt = v(3), t = v(6);
      // a wandering path up the body, pulled sideways towards where the field is strong
      const blobs = [];
      let x = 0, y = 0;
      for (let k = 0; k < nb; k++) {
        const z = H * (0.12 + (0.76 * k) / Math.max(1, nb - 1));
        let bx = x, bs = -1;
        for (let s2 = 0; s2 < 5; s2++) { const cx = Math.max(-W / 3, Math.min(W / 3, x + (rand() - 0.5) * W * 0.5 * flow)), sc = f(cx, z) + rand() * 0.2; if (sc > bs) { bs = sc; bx = cx; } }
        x = bx; y = Math.max(-D / 3, Math.min(D / 3, y + (rand() - 0.5) * D * 0.4 * flow));
        blobs.push([x, isWall(mount) ? Math.max(-D / 4, Math.min(D / 4, y)) : y, z, R0 * (0.6 + 0.6 * rand())]);
      }
      const k2 = 1 / (melt * melt), h0 = Math.max(0.03, Math.cbrt((W * (isWall(mount) ? env.depth : W) * H) / 22000));
      for (const b of blobs) b[3] = Math.max(b[3], (2 * h0) / (0.83 * melt)); // never smaller than the mesh can see
      const val = (px, py, pz) => { let s2 = 0; for (const [bx, by, bz, r] of blobs) s2 += Math.exp((-k2 * ((px - bx) ** 2 + (py - by) ** 2 + (pz - bz) ** 2)) / (r * r)); return 0.5 - s2; };
      const y0 = -D / 2, h = Math.max(0.03, Math.cbrt((W * D * H) / 22000)), gx = Math.ceil(W / h), gy = Math.max(2, Math.ceil(D / h)), gz = Math.ceil(H / h);
      let tris = isoTriangles(val, -W / 2, y0, 0, W / gx, D / gy, H / gz, gx, gy, gz);
      // holes: random points on the surface, open widest where the field is strong
      const holes = [];
      for (let k = 0; k < v(4) && tris.length; k++) { const tr = tris[Math.floor(rand() * tris.length)], c = mul(add(add(tr[0], tr[1]), tr[2]), 1 / 3); holes.push([c, v(5) * (0.25 + 0.75 * f(c[0], c[2]))]); }
      if (holes.length) tris = tris.filter((tr) => { const c = mul(add(add(tr[0], tr[1]), tr[2]), 1 / 3); return holes.every(([q, r]) => Math.hypot(...sub(c, q)) > r); });
      let parts = tris.map((tr) => ({ tri: tr, t }));
      if (hung(mount)) { parts = mapParts(parts, (q) => [q[0], q[1], q[2] - H]); parts.push({ a: [0, 0, 0], b: [0, 0, -H * 0.15], r: 0.004, structure: true }); parts.label = "a hanging cable"; return parts; }
      parts.push({ c: [0, 0, 0.01], u: [W * 0.3, 0, 0], v: [0, isWall(mount) ? Math.max(0.12, D) : W * 0.3, 0], t: 0.02, structure: true }, { a: [0, 0, 0.02], b: [0, 0, H * 0.5], r: 0.012, structure: true });
      parts.label = "a weighted base with a stem";
      return parts;
    },
    describe(g) { const v = (i) => paramValue("bioform", g, i); if (v(8) === 1) return `A slab with ${v(0)} pulls stretched out of it (leaving liquid voids behind), a ${Math.round(v(6) * 1000)} mm heat-formed felt shell with ${v(4)} holes.`; return `${v(0)} blobs melted into one body, a ${Math.round(v(6) * 1000)} mm heat-formed felt shell with ${v(4)} holes up to ${Math.round(v(5) * 100)} cm.`; },
  };
  // Agent fibres: a swarm of agents crawls over the face, steered by a flow direction, by their
  // neighbours (they line up, like a flock) and towards where sound is strong, weaving in depth.
  // Their trails become cords; where sound is strong they sprout short hooked tufts.
  LAB.swarm = {
    group: "pattern", lab: true, direct: true, mounts: true, name: "Agent fibres", imitates: "mycelium, root mats, coir and sisal, a flock's paths, fibre printed by a robot",
    idea: "Many simple agents crawl over the face, following a flow, lining up with their neighbours and drifting towards strong sound. Their trails become a mat of soft cords, with hooked tufts where sound is strongest. A mat of fibres is how felt absorbs, here at a larger scale.",
    params: [num("Agents", 20, 300, "", CELL, true), num("Steps per agent", 20, 200, "", CELL, true), num("Flow angle", 0, 180, "°", CELL), num("Follow the flow", 0, 1, "%", CELL),
      num("Line up with neighbours", 0, 1, "%", CELL), num("Drift towards sound", 0, 1, "%", RESP), num("Wander", 0, 1, "%", CELL), num("Weave depth", 0.02, 0.2, "cm", CELL),
      num("Cord thickness", 0.002, 0.012, "mm", CELL), num("Tufts where sound is strong", 0, 1, "%", RESP), num("Tuft length", 0.02, 0.12, "cm", RESP), opt("Made as", ["soft cord (absorbs)", "printed"], CELL),
      num("Contrast", 0.5, 3, "", RESP), ...holdParams],
    build(g, down, ctx) {
      const v = (i) => paramValue("swarm", g, i), mount = ctx.mount, { W, Hm } = frameOf(ctx.env, mount), f = fieldFn(ctx.field, W, Hm, v(12)), rand = rng(seedOf(g) + 31);
      const D = Math.min(ctx.env.depth || 0.3, v(7)), r = v(8), printed = v(11) === 1, a0 = (v(2) * Math.PI) / 180, base = [Math.cos(a0), Math.sin(a0)];
      let nA = v(0), steps = v(1);
      while (nA * steps > 14000) { nA = Math.round(nA * 0.9); steps = Math.round(steps * 0.95); }
      const stepLen = Math.max(0.012, Math.min(W, Hm) / 60), R = 0.08, parts = [], extra = printed ? { structure: true } : {};
      const spawn = () => { const side = rand(); return side < 0.5 ? { x: -W / 2 + rand() * W, z: base[1] > 0 ? 0.001 : Hm - 0.001 } : { x: base[0] > 0 ? -W / 2 + 0.001 : W / 2 - 0.001, z: rand() * Hm }; };
      const agents = Array.from({ length: nA }, () => ({ ...{ x: -W / 2 + rand() * W, z: rand() * Hm }, d: [...base], ph: rand() * 6.283, y: 0 }));
      const order = [];
      for (let s2 = 0; s2 < steps; s2++) {
        const grid = new Map(), cellOf = (a) => Math.floor((a.x + W) / R) + "," + Math.floor(a.z / R);
        for (const a of agents) { const k = cellOf(a); if (!grid.has(k)) grid.set(k, []); grid.get(k).push(a); }
        for (const a of agents) {
          // line up with neighbours
          let ax = 0, az = 0;
          const cx = Math.floor((a.x + W) / R), cz = Math.floor(a.z / R);
          for (let i = -1; i <= 1; i++) for (let k = -1; k <= 1; k++) for (const b of grid.get(cx + i + "," + (cz + k)) || []) if (b !== a) { ax += b.d[0]; az += b.d[1]; }
          const al = Math.hypot(ax, az), e = 0.02, gx = (f(a.x + e, a.z) - f(a.x - e, a.z)) / (2 * e), gz = (f(a.x, a.z + e) - f(a.x, a.z - e)) / (2 * e), gl = Math.hypot(gx, gz) || 1;
          const turn = (rand() - 0.5) * 2 * v(6) * 0.9;
          let dx = a.d[0] + v(3) * 0.3 * (base[0] - a.d[0]) + (al ? v(4) * 0.5 * (ax / al - a.d[0]) : 0) + v(5) * 0.4 * (gx / gl);
          let dz = a.d[1] + v(3) * 0.3 * (base[1] - a.d[1]) + (al ? v(4) * 0.5 * (az / al - a.d[1]) : 0) + v(5) * 0.4 * (gz / gl);
          const c = Math.cos(turn), sn = Math.sin(turn);
          [dx, dz] = [dx * c - dz * sn, dx * sn + dz * c];
          const dl = Math.hypot(dx, dz) || 1;
          a.d = [dx / dl, dz / dl];
          const P = [a.x, D / 2 + (D / 2 - r) * Math.sin(a.ph), a.z];
          a.x += a.d[0] * stepLen; a.z += a.d[1] * stepLen; a.ph += stepLen / 0.06;
          if (a.x < -W / 2 || a.x > W / 2 || a.z < 0 || a.z > Hm) { Object.assign(a, spawn(), { d: [...base] }); continue; }
          const Q = [a.x, D / 2 + (D / 2 - r) * Math.sin(a.ph), a.z], fq = f(a.x, a.z);
          parts.push({ a: P, b: Q, r, ...extra });
          if (rand() < v(9) * fq * 0.15) {
            const L = v(10) * (0.5 + fq), hook = [a.d[0] * 0.5, 0, a.d[1] * 0.5], t1 = add(Q, [0, L * 0.7, 0]);
            parts.push({ a: Q, b: [t1[0], Math.min(D + L, t1[1]), t1[2]], r: r * 0.8, ...extra }, { a: [t1[0], Math.min(D + L, t1[1]), t1[2]], b: add([t1[0], Math.min(D + L, t1[1]), t1[2]], mul(hook, L * 0.6)), r: r * 0.7, ...extra });
          }
        }
      }
      return mountPattern(parts, mount, { W, Hm, depth: D, mirror: false, feet: 2 + v(13), footLen: v(14) });
    },
    describe(g) { const v = (i) => paramValue("swarm", g, i); return `${v(0)} agents over ${v(1)} steps, flowing at ${Math.round(v(2))}°, ${Math.round(v(8) * 1000)} mm ${v(11) ? "printed fibres" : "soft cords"} woven ${Math.round(v(7) * 100)} cm deep${v(9) > 0.1 ? ", hooked tufts where sound is strong" : ""}.`; },
  };
  // "Pulled from a slab" (after Ross Lovegrove's Liquid Table Bioform, whose legs seem stretched out of the
  // metal, leaving a liquid void in the top): a slab with pulls stretched out of it towards where the field
  // is strong, and a dimple left above each pull. On a wall the pulls come out of the wall; under a ceiling
  // they hang like stalactites; standing, the slab rests on its pulls.
  function pulledBody(g, ctx) {
    const v = (i) => paramValue("bioform", g, i), mount = ctx.mount, env = ctx.env, W = env.width, rand = rng(seedOf(g) + 23);
    const span = isWall(mount) ? env.height - 0.1 : W, L = isWall(mount) ? env.depth : mount === "hanging" ? env.width * 0.8 : env.height - 0.05;
    const f = fieldFn(ctx.field, W, span, v(7)), np = v(0), r0 = Math.max(0.03, v(1) * 0.35), melt = v(3), th = 0.04, t = v(6), flow = v(2);
    // pulls at the strongest places of the field, kept apart; local frame: a across, b along the slab, n out of it
    const cand = [];
    for (let i = 0; i < 10; i++) for (let k = 0; k < 10; k++) { const a = -W / 2 + ((i + 0.5) / 10) * W, b = ((k + 0.5) / 10) * span; cand.push([a, b, f(a, b) + rand() * 0.25]); }
    cand.sort((p, q) => q[2] - p[2]);
    const pulls = [];
    for (const c of cand) if (pulls.length < np && pulls.every((q) => Math.hypot(q[0] - c[0], q[1] - c[1]) > r0 * 2.2)) pulls.push([c[0], c[1], mount === "standing" ? L : L * (0.45 + 0.55 * (1 - flow * rand()))]);
    const k2 = 1 / (melt * melt);
    const seg = (a, b, n, P) => { const tt = Math.max(0, Math.min(1, n / P[2])), d2 = (a - P[0]) ** 2 + (b - P[1]) ** 2 + (n - tt * P[2]) ** 2, r = r0 * (1 - 0.55 * tt); return Math.exp((-k2 * d2) / (r * r)); };
    const val = (a, b, n) => {
      let s2 = Math.exp(-((n / th) ** 2) * 2);
      for (const P of pulls) { s2 += seg(a, b, n, P); s2 -= 0.9 * Math.exp((-((a - P[0]) ** 2 + (b - P[1]) ** 2 + (n + th * 0.6) ** 2)) / (r0 * r0 * 0.5)); }
      return 0.5 - s2;
    };
    const h = Math.max(0.025, Math.cbrt((W * span * (L + 0.1)) / 22000)), ga = Math.ceil(W / h), gb = Math.ceil(span / h), gn = Math.ceil((L + 0.12) / h);
    let tris = isoTriangles(val, -W / 2, 0, -0.08, W / ga, span / gb, (L + 0.12) / gn, ga, gb, gn);
    const holes = [];
    for (let k = 0; k < v(4) && tris.length; k++) { const tr = tris[Math.floor(rand() * tris.length)], c = mul(add(add(tr[0], tr[1]), tr[2]), 1 / 3); holes.push([c, v(5) * (0.25 + 0.75 * f(c[0], c[1]))]); }
    if (holes.length) tris = tris.filter((tr) => { const c = mul(add(add(tr[0], tr[1]), tr[2]), 1 / 3); return holes.every(([q, r]) => Math.hypot(...sub(c, q)) > r); });
    const local = tris.map((tr) => ({ tri: tr, t }));
    // local (a, b, n) to the room
    let parts;
    if (isWall(mount)) {
      // a partition with pulls stretched out of one face
      const z0 = mount === "wallhung" ? -span - 0.25 : 0.1;
      parts = mapParts(local, (q) => [q[0], q[2] - L / 2, q[1] + z0]);
      if (mount === "wallhung") { for (const sx of [-1, 1]) parts.push({ a: [(sx * W) / 3, 0, 0], b: [(sx * W) / 3, 0, -0.25], r: 0.003, structure: true }); parts.label = "2 hanging cables"; }
      else parts.label = footings(parts, planCurve(0, 0, W), W, span + 0.1, 2, 0.3, false);
    }
    else if (mount === "hanging") { parts = mapParts(local, (q) => [q[0], q[1] - span / 2, -q[2] - 0.1]); for (const sx of [-1, 1]) parts.push({ a: [(sx * W) / 3, 0, 0], b: [(sx * W) / 3, 0, -0.1], r: 0.003, structure: true }); parts.label = "2 hanging cables"; }
    else {
      parts = mapParts(local, (q) => [q[0], q[1] - span / 2, Math.max(0, L - q[2])]);
      for (const P of pulls) parts.push({ c: [P[0], P[1] - span / 2, 0.006], u: [r0, 0, 0], v: [0, r0, 0], t: 0.012, structure: true });
      parts.label = `its ${pulls.length} pulls as legs, on weighted pads`;
    }
    return parts;
  }
  for (const [k, r] of Object.entries(LAB)) { r.genes = r.params.map((p) => p.label); RULES[k] = r; }

  // Geometric fields to grow from when you want to draw by hand rather than by sound: an attractor
  // point, a gradient, rings, the interference of two waves (Biothing's ceiling) or soft noise.
  function patternField(spec, nx = 24, nz = 36) {
    const data = new Float32Array(nx * nz), { kind = "waves", x = 0.5, z = 0.5, angle = 0, f1 = 3, f2 = 5, seed = 1 } = spec;
    const a = (angle * Math.PI) / 180, rand = rng(seed * 977 + 3), knots = Array.from({ length: 49 }, () => rand());
    const noise = (u, w) => { const U = u * 6, Wv = w * 6, i = Math.min(5, Math.floor(U)), k = Math.min(5, Math.floor(Wv)), fu = U - i, fw = Wv - k, q = (p, r) => knots[p + 7 * r], s = (t) => t * t * (3 - 2 * t);
      return (q(i, k) * (1 - s(fu)) + q(i + 1, k) * s(fu)) * (1 - s(fw)) + (q(i, k + 1) * (1 - s(fu)) + q(i + 1, k + 1) * s(fu)) * s(fw); };
    for (let k = 0; k < nz; k++) for (let i = 0; i < nx; i++) {
      const u = i / (nx - 1), w = k / (nz - 1);
      let val;
      if (kind === "attractor") val = Math.exp(-((u - x) ** 2 + (w - z) ** 2) / (2 * (0.08 + 0.3 / f1) ** 2));
      else if (kind === "gradient") val = (u - 0.5) * Math.cos(a) + (w - 0.5) * Math.sin(a);
      else if (kind === "radial") val = 0.5 + 0.5 * Math.cos(2 * Math.PI * f1 * Math.hypot(u - x, w - z));
      else if (kind === "noise") val = noise(u, w);
      else val = Math.sin(2 * Math.PI * f1 * u) + Math.sin(2 * Math.PI * f2 * (u * Math.cos(a) + w * Math.sin(a)));
      data[i + nx * k] = val;
    }
    return { data: normalise01(data), nx, nz };
  }

  // Quick facts for a live preview: felt mass, print mass and (standing) stability, no wave simulation.
  function inspect(prims, softDensity, structureDensity, mount) {
    const vox = voxelize(prims, 0.05), mass = vox.volume * softDensity;
    let print = 0;
    for (const p of prims) {
      if (!p.structure) continue;
      if (isTri(p)) print += 0.5 * Math.hypot(...cross(sub(p.tri[1], p.tri[0]), sub(p.tri[2], p.tri[0]))) * p.t;
      else if (isSheet(p)) print += 4 * Math.hypot(...cross(p.u, p.v)) * p.t;
      else print += Math.PI * p.r * p.r * Math.hypot(...sub(p.b, p.a));
    }
    return { mass, printMass: print * structureDensity, stance: onGround(mount) ? stability(prims, vox.com, mass, structureDensity) : null };
  }

  // Every rule also gets four "stance" genes (used when the form stands on its own: root flare,
  // holdfast type, base size, number of roots) and a last gene choosing the soft material.
  const STANCE = 4;
  const geneCount = (rule) => RULES[rule].genes.length + STANCE + 1;
  const MIN_TIP_ANGLE = 10; // degrees a freestanding piece must tilt before it topples

  const corners = (p) => [add(add(p.c, p.u), p.v), add(sub(p.c, p.u), p.v), sub(sub(p.c, p.u), p.v), sub(add(p.c, p.u), p.v)];
  function bounds(prims) {
    const min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
    for (const p of prims) {
      const pts = isTri(p) ? p.tri : isSheet(p) ? corners(p) : [p.a, p.b], pad = isTri(p) || isSheet(p) ? p.t / 2 : p.r;
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
    if (isWall(mount)) return { width: size, height: Math.min(2.2, size * 1.5), depth: Math.min(0.3, Math.max(0.1, 0.25 * size)), mount };
    return { width: size, height: mount === "hanging" ? size : Math.min(2.2, size * 1.5), depth: 0.3, mount };
  }
  // A form grown for the ceiling, turned on its side and pressed into a double-sided relief no thicker
  // than a partition (like a tree trained flat on a frame); on feet, or hung from the ceiling.
  function toPartition(parts, env, mount) {
    const r = mapParts(parts, (q) => [q[0], -q[2], q[1]]), b = bounds(r.map((p) => (isSheet(p) ? { ...p, t: 0 } : p))), Hm = env.height - 0.1;
    const sf = Math.min(env.width / (b.max[0] - b.min[0] || 1e-9), Hm / (b.max[2] - b.min[2] || 1e-9)), sd = Math.min(sf, env.depth / (b.max[1] - b.min[1] || 1e-9));
    const cx = (b.min[0] + b.max[0]) / 2, cy = (b.min[1] + b.max[1]) / 2, z0 = mount === "wallhung" ? -Hm - 0.25 : 0.1;
    const tr = (q) => [(q[0] - cx) * sf, (q[1] - cy) * sd, (q[2] - b.min[2]) * sf + z0];
    const out = r.map((p) => (isSheet(p) ? { ...p, c: tr(p.c), u: [p.u[0] * sf, p.u[1] * sd, p.u[2] * sf], v: [p.v[0] * sf, p.v[1] * sd, p.v[2] * sf] } : isTri(p) ? { ...p, tri: p.tri.map(tr) } : { ...p, a: tr(p.a), b: tr(p.b), r: p.r * sd }));
    if (mount === "wallhung") {
      for (const sx of [-1, 1]) out.push({ a: [(sx * env.width) / 3, 0, 0], b: [(sx * env.width) / 3, 0, -0.25], r: 0.003, structure: true }, { a: [(sx * env.width) / 3, 0, -0.25], b: [0, 0, z0 + Hm * 0.9], r: 0.004, structure: true });
      out.label = "2 hanging cables";
    } else out.label = footings(out, planCurve(0, 0, env.width), env.width, Hm + 0.1, 2, Math.max(0.2, env.depth), false);
    return out;
  }
  function grow(rule, genes, opts) {
    const env = envelope(opts.size, opts.mount), R = RULES[rule], ctx = { env, field: opts.field, material: opts.material, mount: opts.mount };
    if (R.direct) return R.build(genes, hung(opts.mount), ctx);
    if (isWall(opts.mount)) return toPartition(R.build(genes, true), env, opts.mount);
    const parts = fit(RULES[rule].build(genes, opts.mount === "hanging"), env);
    return opts.mount === "standing" ? standUp(parts, genes.slice(RULES[rule].genes.length, RULES[rule].genes.length + STANCE), env) : parts;
  }

  // Freestanding adaptation, the way trees and corals do it: a flare that widens the form towards
  // the ground, then a holdfast (splayed roots, or a weighted plinth) with a central stem.
  // Holdfast parts are structure (rigid, not absorbing); everything else stays felt.
  function standUp(parts, st, env) {
    const flare = lerp(0, 0.6, st[0]), H = env.height;
    const k = (z) => 1 + flare * (1 - (2 * Math.max(0, Math.min(H, z))) / H);
    const fl = (q) => [q[0] * k(q[2]), q[1] * k(q[2]), q[2]];
    let out = parts.map((p) => (isTri(p) ? { ...p, tri: p.tri.map(fl) } : isSheet(p) ? { ...p, c: fl(p.c) } : { ...p, a: fl(p.a), b: fl(p.b) }));
    const b = bounds(out), wide = Math.max(b.max[0] - b.min[0], b.max[1] - b.min[1]);
    if (wide > env.width) {
      const s = env.width / wide, sx = (q) => [q[0] * s, q[1] * s, q[2]];
      out = out.map((p) => (isTri(p) ? { ...p, tri: p.tri.map(sx) } : isSheet(p) ? { ...p, c: sx(p.c), u: [p.u[0] * s, p.u[1] * s, p.u[2]], v: [p.v[0] * s, p.v[1] * s, p.v[2]] } : { ...p, a: sx(p.a), b: sx(p.b) }));
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
    out.label = (roots ? `${legs} splayed roots` : "a weighted plinth") + (flare > 0.18 ? ", widening towards the floor" : "");
    return out;
  }

  // Tipping check: the angle the piece can lean before its centre of mass leaves the footprint
  // (the convex hull of the structure touching the floor).
  function hull(pts) {
    const P = pts.slice().sort((a, b) => a[0] - b[0] || a[1] - b[1]);
    if (P.length < 3) return P;
    const cr = (o, a, b) => (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
    const lo = [], hi = [];
    for (const p of P) { while (lo.length >= 2 && cr(lo[lo.length - 2], lo[lo.length - 1], p) <= 0) lo.pop(); lo.push(p); }
    for (const p of P.slice().reverse()) { while (hi.length >= 2 && cr(hi[hi.length - 2], hi[hi.length - 1], p) <= 0) hi.pop(); hi.push(p); }
    return lo.slice(0, -1).concat(hi.slice(0, -1));
  }
  function stability(parts, softCom, softMass, structureDensity) {
    let m = softMass, cx = softCom[0] * softMass, cy = softCom[1] * softMass, cz = softCom[2] * softMass, structureMass = 0;
    const ground = [];
    for (const p of parts) {
      if (!p.structure) continue;
      let vol, c;
      if (isTri(p)) {
        const [A, B, C] = p.tri;
        vol = 0.5 * Math.hypot(...cross(sub(B, A), sub(C, A))) * p.t; c = mul(add(add(A, B), C), 1 / 3);
        for (const q of p.tri) if (q[2] < 0.03) ground.push([q[0], q[1]]);
      } else if (isSheet(p)) {
        vol = 4 * Math.hypot(...cross(p.u, p.v)) * p.t; c = p.c;
        if (p.c[2] < 0.03) for (const q of corners(p)) ground.push([q[0], q[1]]);
      } else {
        vol = Math.PI * p.r * p.r * Math.hypot(...sub(p.b, p.a)); c = mul(add(p.a, p.b), 0.5);
        for (const q of [p.a, p.b]) if (q[2] < 0.03) ground.push([q[0], q[1]]);
      }
      const pm = vol * structureDensity;
      structureMass += pm; m += pm; cx += c[0] * pm; cy += c[1] * pm; cz += c[2] * pm;
    }
    const H = hull(ground), com = [cx / m, cy / m], h = Math.max(cz / m, 1e-3);
    let inside = -Infinity, across = 0;
    if (H.length >= 3) {
      inside = Infinity;
      for (let i = 0; i < H.length; i++) {
        const a = H[i], b = H[(i + 1) % H.length], ex = b[0] - a[0], ey = b[1] - a[1], len = Math.hypot(ex, ey) || 1e-9;
        inside = Math.min(inside, ((com[0] - a[0]) * ey - (com[1] - a[1]) * ex) / -len);
      }
      for (const a of H) for (const b of H) across = Math.max(across, Math.hypot(a[0] - b[0], a[1] - b[1]));
    }
    const angle = H.length >= 3 ? (Math.atan2(inside, h) * 180) / Math.PI : -90;
    return { angle, ok: angle >= MIN_TIP_ANGLE, structureMass, across, comHeight: h };
  }

  // Compact transfer format (11 numbers per part) for workers and storage.
  // Type code: 0 capsule, 1 sheet, +2 when the part is structure; 4 triangle, 5 structure triangle.
  const STRIDE = 11;
  function pack(prims) {
    const f = new Float32Array(prims.length * STRIDE);
    prims.forEach((p, i) => f.set(isTri(p) ? [p.structure ? 5 : 4, ...p.tri[0], ...p.tri[1], ...p.tri[2], p.t] : isSheet(p) ? [p.structure ? 3 : 1, ...p.c, ...p.u, ...p.v, p.t] : [p.structure ? 2 : 0, ...p.a, ...p.b, p.r, 0, 0, 0], i * STRIDE));
    return f;
  }
  function unpack(f) {
    const out = [];
    for (let i = 0; i < f.length; i += STRIDE) {
      const x = Array.from(f.subarray(i + 1, i + STRIDE));
      if (f[i] >= 4) { const tp = { tri: [x.slice(0, 3), x.slice(3, 6), x.slice(6, 9)], t: x[9] }; if (f[i] === 5) tp.structure = true; out.push(tp); continue; }
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
    if (!prims.length) return { frac: new Float32Array(1), n: [1, 1, 1], o: [0, 0, 0], dx, volume: 0, com: [0, 0, 0] };
    const { min, max } = bounds(prims);
    const o = min.map((m) => m - dx), n = max.map((m, i) => Math.max(1, Math.ceil((m + dx - o[i]) / dx)));
    const frac = new Float32Array(n[0] * n[1] * n[2]), R = 0.75 * dx, NORM = 0.589;
    const range = (lo3, hi3) => [[0, 1, 2].map((i) => Math.max(0, Math.floor((lo3[i] - o[i]) / dx))), [0, 1, 2].map((i) => Math.min(n[i] - 1, Math.ceil((hi3[i] - o[i]) / dx)))];
    const soft = (x) => Math.max(0, Math.min(1, x / dx + 0.5));
    for (const p of prims) {
      if (isTri(p)) {
        const c = mul(add(add(p.tri[0], p.tri[1]), p.tri[2]), 1 / 3), i = Math.floor((c[0] - o[0]) / dx), j = Math.floor((c[1] - o[1]) / dx), k = Math.floor((c[2] - o[2]) / dx);
        if (i >= 0 && j >= 0 && k >= 0 && i < n[0] && j < n[1] && k < n[2]) frac[i + n[0] * (j + n[1] * k)] += (0.5 * Math.hypot(...cross(sub(p.tri[1], p.tri[0]), sub(p.tri[2], p.tri[0]))) * p.t) / dx ** 3;
        continue;
      }
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
    if (lo === Infinity) return { grid: new Float32Array(1), nh: 1, nz: 1, dx, extent: 0, x0: vox.o[h], z0: vox.o[2] };
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
  function roomFor(mount) { return hung(mount) ? { w: 6, h: 4, cx: 3.4 } : { w: 6, h: 3.5, cx: 3.4 }; }

  function makeSim(pr, sigma, freq, room, mount, opts = {}) {
    const dx = pr.dx, nx = Math.round(room.w / dx), nz = Math.round(room.h / dx), dt = (0.68 * dx) / C0;
    const p = new Float32Array(nx * nz), vx = new Float32Array((nx + 1) * nz), vz = new Float32Array(nx * (nz + 1)), s = new Float32Array(nx * nz);
    const ox = Math.round((room.cx + pr.x0) / dx), oz = hung(mount) ? Math.round((room.h + pr.z0) / dx) : Math.round(pr.z0 / dx);
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
    const sx = Math.max(1, Math.min(nx - 2, Math.round((opts.src ? opts.src[0] : 0.4) / dx))), sz = Math.max(1, Math.min(nz - 2, Math.round((opts.src ? opts.src[1] : 0.4) / dx)));
    const src = sx + nx * sz, sig = opts.signal, cycles = 6, tb = sig ? sig.length * dt : cycles / freq, pc = RHO0 * C0 * C0 * dt / dx;
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
          if (sig) { if (step < sig.length) p[src] += sig[step]; }
          else if (t < tb) p[src] += Math.sin(2 * Math.PI * freq * t) * Math.sin((Math.PI * t) / tb) ** 2;
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
      // Adds the air-motion energy at every cell centre into `acc` (where porous material works hardest).
      addKinetic(acc) {
        for (let k = 0; k < nz; k++) for (let i = 0; i < nx; i++) {
          const ux = (vx[i + (nx + 1) * k] + vx[i + 1 + (nx + 1) * k]) / 2, uz = (vz[i + nx * k] + vz[i + nx * (k + 1)]) / 2;
          acc[i + nx * k] += ux * ux + uz * uz;
        }
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
  function evaluate({ rule, genes, mount, size, material, budget, weights, structureDensity = 680, field = null, skin = null }) {
    let prims = grow(rule, genes, { mount, size, field, material });
    if (skin) {
      const tufts = growSkin(prims, skin.field, skin.genes, envelope(size, mount), mount), merged = prims.concat(tufts);
      for (const k of ["label", "field", "cutPlan"]) merged[k] = prims[k];
      merged.skinCount = tufts.length;
      prims = merged;
    }
    const vox = voxelize(prims, C0 / BANDS[BANDS.length - 1] / CELLS_PER_WAVE);
    const mass = vox.volume * material.density;
    // Strip method: absorption per metre in two slices per direction, times the extent.
    const room = roomFor(mount), slices = [], axes = [0, 1];
    for (const axis of axes) for (const at of [0.3, 0.7]) slices.push({ axis, sl: slice(vox, axis, at) });
    const area = vox.volume === 0 ? BANDS.map(() => 0) : BANDS.map((f, bi) => {
      const factor = 2 ** (BANDS.length - 1 - bi);
      const perAxis = axes.map((axis) => {
        const mine = slices.filter((x) => x.axis === axis);
        const mean = mine.reduce((acc, { sl }) => acc + absorptionLength(downsample(sl, factor), material.sigma, f, room, mount).length, 0) / mine.length;
        return mean * mine[0].sl.extent;
      });
      return perAxis.reduce((a, x) => a + x, 0) / perAxis.length;
    });
    const w = weights || [1 / 3, 1 / 3, 1 / 3];
    const score = area.reduce((a, x, i) => a + x * w[i], 0);
    const over = mass > budget;
    let fitness = score - (over ? score * (mass / budget - 1) + 0.05 : 0);
    const stance = onGround(mount) ? stability(prims, vox.com, mass, structureDensity) : null;
    if (stance && !stance.ok) fitness -= score * 0.8 * Math.min(1, (MIN_TIP_ANGLE - stance.angle) / MIN_TIP_ANGLE) + 0.05;
    const printMass = inspect(prims, 0, structureDensity, "wall").printMass;
    if (RULES[rule].group === "lattice" && printMass > 3 * budget) fitness -= score * (printMass / (3 * budget) - 1);
    return { rule, genes, mount, size, prims, mass, volume: vox.volume, area, score, fitness, perKg: score / Math.max(mass, 1e-6), stance,
      holdfast: prims.label || null, fieldMap: prims.field || null, cutPlan: prims.cutPlan || null, stripPlan: prims.stripPlan || null, skinCount: prims.skinCount || 0, printMass };
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

  // ------------------------------------------------------------- fields ----
  // Structural fields from sound: where the air moves most in front of the piece (normalised 0..1),
  // sampled on an nx x nz grid over the piece (width along its length, height from the floor).
  function patch(acc, sim, x0, width, height, nx, nz, zBase = 0) {
    const data = new Float32Array(nx * nz);
    for (let k = 0; k < nz; k++) for (let i = 0; i < nx; i++) {
      const x = x0 - width / 2 + (i / (nx - 1)) * width, z = zBase + (k / (nz - 1)) * height;
      const ci = Math.max(0, Math.min(sim.nx - 1, Math.round(x / sim.dx))), ck = Math.max(0, Math.min(sim.nz - 1, Math.round(z / sim.dx)));
      data[i + nx * k] = acc[ci + sim.nx * ck];
    }
    let lo = Infinity, hi = -Infinity;
    for (const v of data) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
    return data.map((v) => (hi > lo ? (v - lo) / (hi - lo) : 0.5));
  }
  // Tone bursts (weighted bands) or a recorded signal (sampled at the solver step), from a talker
  // position in a room section; the field is the air motion accumulated over the direct sound and
  // early reflections.
  function fieldFromSound({ room, weights = [1 / 3, 1 / 3, 1 / 3], signal = null, x0, width, height, nx = 24, nz = 36, src, zBase = 0 }) {
    const dx = C0 / BANDS[BANDS.length - 1] / CELLS_PER_WAVE, empty = { grid: new Float32Array(1), nh: 1, nz: 1, dx, x0: 0, z0: 0, extent: 1 };
    const total = new Float32Array(Math.round(room.w / dx) * Math.round(room.h / dx)), runs = signal ? [[0, 1]] : BANDS.map((f, i) => [f, weights[i]]);
    let sim = null;
    for (const [f, w] of runs) {
      if (!w) continue;
      sim = makeSim(empty, 0, f || 500, room, "standing", { src, signal });
      const acc = new Float32Array(sim.nx * sim.nz), t1 = sim.burstEnd + 0.1;
      while (sim.t < t1) { sim.step(4); sim.addKinetic(acc); }
      let sum = 0;
      for (const v of acc) sum += v;
      for (let i = 0; i < acc.length; i++) total[i] += (w * acc[i]) / (sum || 1);
    }
    return { data: patch(total, sim, x0, width, height, nx, nz, zBase), nx, nz, width, height };
  }
  // Standing waves of a box room (its resonances): the Chladni-like pattern of air motion on a
  // vertical plane through the piece, summed over the room modes between fmin and fmax (the range of
  // voice fundamentals, whose waves are short enough to draw a pattern across a piece).
  function roomModes(L, W, H, fmin, fmax, maxModes) {
    const modes = [], N = (len) => Math.ceil((2 * len * fmax) / C0);
    for (let a = 0; a <= N(L); a++) for (let b = 0; b <= N(W); b++) for (let c = 0; c <= N(H); c++) {
      if (!a && !b && !c) continue;
      const f = (C0 / 2) * Math.hypot(a / L, b / W, c / H);
      if (f >= fmin && f <= fmax) modes.push({ n: [a, b, c], f });
    }
    modes.sort((p, q) => p.f - q.f);
    const stride = Math.max(1, modes.length / maxModes);
    return Array.from({ length: Math.min(maxModes, modes.length) }, (_, i) => modes[Math.floor(i * stride)]);
  }
  // Squared air velocity of one mode at (x, y, z), up to a constant.
  function modeVelocity2(n, L, W, H, x, y, z) {
    const kx = (n[0] * Math.PI) / L, ky = (n[1] * Math.PI) / W, kz = (n[2] * Math.PI) / H;
    const cx = Math.cos(kx * x), sx = Math.sin(kx * x), cy = Math.cos(ky * y), sy = Math.sin(ky * y), cz = Math.cos(kz * z), sz = Math.sin(kz * z);
    return (kx * sx * cy * cz) ** 2 + (ky * cx * sy * cz) ** 2 + (kz * cx * cy * sz) ** 2;
  }
  const normalise01 = (data) => {
    let lo = Infinity, hi = -Infinity;
    for (const v of data) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
    return data.map((v) => (hi > lo ? (v - lo) / (hi - lo) : 0.5));
  };
  // ceilingZ: sample a horizontal plane at that height instead (a piece hanging under the ceiling).
  function fieldFromModes({ L, W, H, x0, y0, width, height, nx = 24, nz = 36, fmin = 70, fmax = 200, maxModes = 60, zBase = 0, ceilingZ = null }) {
    const use = roomModes(L, W, H, fmin, fmax, maxModes), data = new Float32Array(nx * nz);
    for (const { n } of use) for (let k = 0; k < nz; k++) for (let i = 0; i < nx; i++) {
      const x = x0 - width / 2 + (i / (nx - 1)) * width, v = k / (nz - 1);
      data[i + nx * k] += ceilingZ === null ? modeVelocity2(n, L, W, H, x, y0, zBase + v * height) : modeVelocity2(n, L, W, H, x, Math.max(0, Math.min(W, y0 - height / 2 + v * height)), ceilingZ);
    }
    return { data: normalise01(data), nx, nz, width, height, modes: use.map((m) => ({ n: m.n, f: m.f })) };
  }

  // ------------------------------------------------------ room placement ----
  // "Probabilistic programme" (Biothing's tower): instead of one fixed layout, a map of the floor plan
  // saying where a piece is likely to help, and layouts drawn from it. Two parts, weighted by what the
  // room needs (need = shares of 125 ... 4000 Hz):
  //  - low pitches (125-250 Hz): where the room's standing waves move the air most at ear height,
  //    because porous felt only works where air moves (not against a wall, where it stands still);
  //  - speech pitches (500 Hz up): a band 0.4-1.6 m in from the walls, where pieces meet the first
  //    reflections and leave the middle free for visitors.
  function placementMap({ lx, ly, H, need, nx = 40, z = 1.2, piece = 1.2 }) {
    const ny = Math.max(4, Math.round((nx * ly) / lx)), low = new Float32Array(nx * ny), mid = new Float32Array(nx * ny);
    const modes = roomModes(lx, ly, H, 30, 250, 80), at = (i, j) => [((i + 0.5) / nx) * lx, ((j + 0.5) / ny) * ly];
    for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
      const [x, y] = at(i, j);
      let v = 0;
      for (const { n } of modes) v += modeVelocity2(n, lx, ly, H, x, y, z);
      low[i + nx * j] = v;
      const d = Math.min(x, lx - x, y, ly - y);
      mid[i + nx * j] = 0.35 + 0.65 * Math.exp(-(((d - 1.0) / 0.7) ** 2));
    }
    const L1 = normalise01(low), M1 = normalise01(mid), wl = (need[0] || 0) + (need[1] || 0), wm = Math.max(0, 1 - wl);
    const data = new Float32Array(nx * ny);
    for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
      const [x, y] = at(i, j), d = Math.min(x, lx - x, y, ly - y), idx = i + nx * j;
      data[idx] = d < piece / 2 + 0.1 ? 0 : wl * L1[idx] + wm * M1[idx];
    }
    return { data: normalise01(data), nx, ny, lx, ly, low: L1, mid: M1, modes: modes.length };
  }
  // Draw n spots from the map, each at least `gap` from the others (relaxed if the room is full).
  // Returns spots in metres from the room's corner and how much better they sit than an even spread.
  function sampleLayout(map, n, seed, gap, sharp = 3) {
    const rand = rng(seed), { nx, ny, lx, ly, data } = map, spots = [];
    const cells = [];
    for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) if (data[i + nx * j] > 0) cells.push({ x: ((i + 0.5) / nx) * lx, y: ((j + 0.5) / ny) * ly, w: Math.pow(data[i + nx * j], sharp), v: data[i + nx * j] });
    let g = gap;
    while (spots.length < n && cells.length) {
      const ok = cells.filter((c) => spots.every((p) => Math.hypot(p.x - c.x, p.y - c.y) >= g));
      if (!ok.length) { g *= 0.8; if (g < 0.2) break; continue; }
      let r = rand() * ok.reduce((a, c) => a + c.w, 0), pick = ok[ok.length - 1];
      for (const c of ok) { r -= c.w; if (r <= 0) { pick = c; break; } }
      spots.push(pick);
    }
    const mean = cells.reduce((a, c) => a + c.v, 0) / Math.max(1, cells.length) || 1;
    const got = spots.reduce((a, c) => a + c.v, 0) / Math.max(1, spots.length);
    return { spots: spots.map((c) => [c.x, c.y]), gain: got / mean, gap: g };
  }

  // ------------------------------------------------------------- skin ------
  // Growing one behaviour on another (Genware): a second, finer growth of felt tufts on the surfaces
  // of an existing form, densest where a second field (usually the speech-band sound) is strongest.
  // genes: density, length, thickness, contrast.
  function growSkin(parts, field, sg, env, mount) {
    const rand = rng(seedOf(sg) + 17), density = lerp(15, 160, sg[0]), len = lerp(0.03, 0.12, sg[1]), r = lerp(0.004, 0.012, sg[2]), gamma = lerp(0.5, 3, sg[3]);
    const H = env.height, f = field || flatField, out = [], MAX = 1200;
    const fAt = (q) => Math.pow(Math.max(0, Math.min(1, sampleField(f.data, f.nx, f.nz, q[0] / env.width + 0.5, (hung(mount) ? q[2] + H : q[2]) / H))), gamma);
    const tuft = (q, nrm) => { if (out.length < MAX) out.push({ a: q, b: add(q, mul(nrm, len)), r }); };
    const tries = (area) => { const x = area * density; return Math.floor(x) + (rand() < x - Math.floor(x) ? 1 : 0); };
    for (const p of parts) {
      if (p.structure || isTri(p)) continue;
      if (isSheet(p)) {
        const nh = norm(cross(p.u, p.v)), area = 8 * Math.hypot(...cross(p.u, p.v));
        for (let k = tries(area); k > 0; k--) {
          const sgn = rand() < 0.5 ? -1 : 1, q = add(add(p.c, mul(p.u, rand() * 2 - 1)), add(mul(p.v, rand() * 2 - 1), mul(nh, (sgn * p.t) / 2)));
          if (rand() < fAt(q)) tuft(q, mul(nh, sgn));
        }
        continue;
      }
      const d = sub(p.b, p.a), Ls = Math.hypot(...d), sphere = Ls < 1e-6, area = sphere ? 4 * Math.PI * p.r * p.r : 2 * Math.PI * p.r * Ls;
      const ax = sphere ? [0, 0, 1] : mul(d, 1 / Ls), u = perp(ax), v = cross(ax, u);
      for (let k = tries(area); k > 0; k--) {
        let nrm;
        if (sphere) { do nrm = [rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1]; while (dot(nrm, nrm) > 1 || dot(nrm, nrm) < 1e-4); nrm = norm(nrm); }
        else { const th = rand() * 2 * Math.PI; nrm = add(mul(u, Math.cos(th)), mul(v, Math.sin(th))); }
        const q = add(add(p.a, mul(d, sphere ? 0 : rand())), mul(nrm, p.r));
        if (rand() < fAt(q)) tuft(q, nrm);
      }
    }
    return out;
  }
  const SKIN_GENES = ["density", "length", "thickness", "contrast"];

  const solverStep = () => (0.68 * (C0 / BANDS[BANDS.length - 1] / CELLS_PER_WAVE)) / C0;

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

  const Growth = { onGround, hung, isWall, simMount, paramValue, patternField, inspect, mapParts, isTri, placementMap, sampleLayout, growSkin, SKIN_GENES, roomModes, sampleField, RULES, BANDS, STRIDE, STANCE, MIN_TIP_ANGLE, stability, hull, fieldFromSound, fieldFromModes, solverStep, pack, unpack, count, isSheet, geneCount, grow, voxelize, project, slice, downsample, makeSim, absorptionLength, roomFor, evaluate, panelBaseline, nextGeneration, rng, bounds, envelope };
  if (typeof module !== "undefined" && module.exports) module.exports = Growth;
  else root.Growth = Growth;
})(typeof self !== "undefined" ? self : this);
