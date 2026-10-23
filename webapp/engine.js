/* Biothing Sound Lab engine: browser port of biothing_acoustics/layers.py and forms.py.
 * Kept numerically identical to the Python version (checked by webapp/engine.test.js).
 * Convention exp(+j*omega*t): lossy wavenumbers have a negative imaginary part. */
(function (root) {
  "use strict";

  // ---------------------------------------------------------- complex ----
  const C = (re, im = 0) => ({ re, im });
  const add = (a, b) => C(a.re + b.re, a.im + b.im);
  const sub = (a, b) => C(a.re - b.re, a.im - b.im);
  const mul = (a, b) => C(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
  const scale = (a, k) => C(a.re * k, a.im * k);
  const div = (a, b) => {
    const d = b.re * b.re + b.im * b.im;
    return C((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d);
  };
  const inv = (a) => div(C(1), a);
  const abs2 = (a) => a.re * a.re + a.im * a.im;
  const csqrt = (a) => { // principal branch, like numpy
    const r = Math.hypot(a.re, a.im);
    const re = Math.sqrt((r + a.re) / 2);
    const im = Math.sqrt(Math.max(0, (r - a.re) / 2));
    return C(re, a.im < 0 ? -im : im);
  };
  const csin = (a) => C(Math.sin(a.re) * Math.cosh(a.im), Math.cos(a.re) * Math.sinh(a.im));
  const ccos = (a) => C(Math.cos(a.re) * Math.cosh(a.im), -Math.sin(a.re) * Math.sinh(a.im));
  const ctan = (a) => div(csin(a), ccos(a));
  const ctanh = (a) => {
    const d = Math.cosh(2 * a.re) + Math.cos(2 * a.im);
    return C(Math.sinh(2 * a.re) / d, Math.sin(2 * a.im) / d);
  };
  const J = C(0, 1);

  // --------------------------------------------------------- constants ----
  const RHO0 = 1.204, C0 = 343.2, Z0 = RHO0 * C0;
  const ETA = 1.81e-5, PRANDTL = 0.71, GAMMA = 1.4, K0 = RHO0 * C0 * C0;
  const FREQS = [100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000];
  const N_THETA = 46;
  const THETA = Array.from({ length: N_THETA }, (_, i) => (i / (N_THETA - 1)) * (89 * Math.PI / 180));
  const N_SUB = 14, WALL = 0.002, PARTITION_SPACING = 0.1;

  // Paris weights: trapezoid of sin(2θ), normalised.
  const PARIS = (() => {
    const w = THETA.map((t) => Math.sin(2 * t));
    const h = THETA[1] - THETA[0];
    const tw = w.map((v, i) => v * h * (i === 0 || i === N_THETA - 1 ? 0.5 : 1));
    const total = tw.reduce((a, b) => a + b, 0);
    return tw.map((v) => v / total);
  })();

  // ---------------------------------------------------------- materials ---
  // Indicative values; same table as biothing_acoustics/materials.py.
  const FELTS = {
    pet_felt_light: { name: "Recycled PET felt, light", sigma: 9000, density: 35, cost: 5.0, carbon: 1.6, group: "Recycled", formable: true, note: "Soft and very light. Made from plastic bottles; can be pressed into shapes." },
    pet_felt_medium: { name: "Recycled PET felt, medium", sigma: 20000, density: 60, cost: 5.0, carbon: 1.6, group: "Recycled", formable: true, note: "The everyday acoustic felt. Can be heat-pressed into reliefs." },
    pet_felt_dense: { name: "Recycled PET felt, dense", sigma: 45000, density: 120, cost: 5.5, carbon: 1.6, group: "Recycled", formable: true, note: "Tough surface, best at middle and high pitches." },
    recycled_cotton: { name: "Recycled cotton fibre", sigma: 12000, density: 45, cost: 3.5, carbon: 0.8, group: "Recycled", formable: false, note: "From old clothes. Low carbon; needs fire treatment." },
    wool_felt: { name: "Wool felt", sigma: 25000, density: 100, cost: 14.0, carbon: 10.0, group: "Animal", formable: false, note: "Dense natural felt, beautiful surface, expensive." },
    sheep_wool_batt: { name: "Sheep wool batt", sigma: 6000, density: 25, cost: 8.0, carbon: 1.0, group: "Animal", formable: false, note: "Fluffy wool padding. Excellent absorber; needs a cover fabric." },
    hemp_fibre: { name: "Hemp fibre mat", sigma: 8000, density: 40, cost: 4.0, carbon: 0.5, group: "Plant", formable: false, note: "Fast-growing plant, very low carbon, slightly rough look." },
    coconut_coir: { name: "Coconut coir mat", sigma: 3000, density: 60, cost: 3.0, carbon: 0.4, group: "Plant", formable: false, note: "Coarse husk fibres. Very open: good as a front layer." },
    cork_granulate: { name: "Bonded cork granulate", sigma: 15000, density: 110, cost: 6.0, carbon: 0.3, group: "Plant", formable: true, note: "Cork-oak bark, harvested without cutting the tree. Can be moulded." },
    cellulose_fibre: { name: "Recycled paper cellulose", sigma: 20000, density: 45, cost: 1.5, carbon: 0.3, group: "Plant", formable: true, note: "Very cheap; can be moulded like papier-mache." },
    wood_wool: { name: "Wood wool board", sigma: 1500, density: 400, cost: 1.5, carbon: 0.5, group: "Plant", formable: false, note: "Wood shavings bound with cement. Robust but heavy." },
    mycelium: { name: "Mycelium composite", sigma: 60000, density: 90, cost: 9.0, carbon: 0.6, group: "Fungus", formable: true, note: "Mushroom roots grown in a mould: a literally grown material." },
    mineral_wool: { name: "Stone wool", sigma: 15000, density: 50, cost: 2.5, carbon: 1.3, group: "Mineral & synthetic", formable: false, note: "Melted rock fibres. Does not burn; must be fully covered." },
    melamine_foam: { name: "Melamine foam", sigma: 10000, density: 9, cost: 15.0, carbon: 4.0, group: "Mineral & synthetic", formable: true, note: "Feather-light foam, cuts into any shape. Fossil-based." },
  };
  const RIGIDS = {
    rigid_pet: { name: "Rigid recycled PET", density: 900, cost: 6.0, carbon: 1.8, note: "Same plastic as PET felt: recycles together." },
    birch_plywood: { name: "Birch plywood", density: 680, cost: 3.5, carbon: 0.7, note: "Strong, easy to CNC-cut." },
    solid_cardboard: { name: "Solid recycled cardboard", density: 700, cost: 1.5, carbon: 0.9, note: "Cheap and light, fine for short exhibitions." },
    pla_print: { name: "3D-printed PLA", density: 1240, cost: 25.0, carbon: 3.0, note: "Plant-starch plastic, any shape; slow to print." },
    fired_clay: { name: "Fired clay", density: 1900, cost: 4.0, carbon: 0.4, note: "Like the sound pots built into ancient theatres and medieval churches." },
  };
  const FELT_KEYS = Object.keys(FELTS);
  const RIGID_KEYS = Object.keys(RIGIDS);

  // ------------------------------------------------------------ targets ---
  const TARGETS = {
    voices: { 125: 0.02, 250: 0.10, 500: 0.20, 1000: 0.28, 2000: 0.25, 4000: 0.15 },
    crowd: { 125: 0.15, 250: 0.25, 500: 0.25, 1000: 0.20, 2000: 0.10, 4000: 0.05 },
    exhibition: { 125: 0.085, 250: 0.175, 500: 0.225, 1000: 0.24, 2000: 0.175, 4000: 0.10 },
    low_noise: { 125: 0.35, 250: 0.30, 500: 0.20, 1000: 0.10, 2000: 0.05, 4000: 0.0 },
  };

  function bandWeights(target) {
    const octaves = Object.keys(TARGETS[target]).map(Number).sort((a, b) => a - b);
    const nearest = FREQS.map((f) => {
      let best = octaves[0], d = Infinity;
      for (const o of octaves) { const x = Math.abs(Math.log2(f / o)); if (x < d) { d = x; best = o; } }
      return best;
    });
    const count = {};
    nearest.forEach((o) => { count[o] = (count[o] || 0) + 1; });
    const w = nearest.map((o) => TARGETS[target][o] / count[o]);
    const s = w.reduce((a, b) => a + b, 0);
    return w.map((v) => v / s);
  }

  // ------------------------------------------------------------- media ----
  function miki(f, sigma) {
    const x = f / sigma, a = Math.pow(x, -0.632), b = Math.pow(x, -0.618);
    return { zc: C(Z0 * (1 + 0.070 * a), -Z0 * 0.107 * a), kc: C((2 * Math.PI * f / C0) * (1 + 0.109 * b), -(2 * Math.PI * f / C0) * 0.160 * b) };
  }
  function porousMedium(f, sigma, frameDensity) {
    const w = 2 * Math.PI * f, { zc, kc } = miki(f, sigma);
    let rho = scale(mul(zc, kc), 1 / w);
    const k = div(scale(zc, w), kc);
    if (frameDensity != null) {
      const rt = frameDensity + RHO0;
      rho = div(sub(scale(rho, rt), C(RHO0 * RHO0)), add(rho, C(rt - 2 * RHO0)));
    }
    return { rho, k };
  }
  function mixedMedium(f, sigma, s, frameDensity) {
    s = Math.min(1, Math.max(0, s));
    const p = porousMedium(f, sigma, frameDensity);
    const rho = inv(add(scale(inv(p.rho), s), C((1 - s) / RHO0)));
    const k = inv(add(scale(inv(p.k), s), C((1 - s) / K0)));
    return { rho, k };
  }
  const AIR = { rho: C(RHO0), k: C(K0) };
  function slitMedium(f, width) {
    const w = 2 * Math.PI * f, h = width / 2, rot = C(Math.SQRT1_2, Math.SQRT1_2); // sqrt(j)
    const a = scale(rot, h * Math.sqrt(w * RHO0 / ETA));
    const b = scale(rot, h * Math.sqrt(w * RHO0 * PRANDTL / ETA));
    const rho = div(C(RHO0), sub(C(1), div(ctanh(a), a)));
    const k = div(C(GAMMA * 101325), add(C(1), scale(div(ctanh(b), b), GAMMA - 1)));
    return { rho, k };
  }
  function perforatedPlate(f, d, t, eps) {
    const w = 2 * Math.PI * f, kp = d / 2 * Math.sqrt(w * RHO0 / ETA);
    const r = 32 * ETA * t / (eps * RHO0 * C0 * d * d) * (Math.sqrt(1 + kp * kp / 32) + Math.SQRT2 / 32 * kp * d / t);
    const x = w * t / (eps * C0) * (1 + 1 / Math.sqrt(9 + kp * kp / 2) + 0.85 * d / t);
    return C(Z0 * r, Z0 * x);
  }

  // ----------------------------------------------------- transfer matrices
  function layer(f, theta, m, d) {
    const w = 2 * Math.PI * f, kx = w / C0 * Math.sin(theta);
    let kz = csqrt(sub(scale(div(m.rho, m.k), w * w), C(kx * kx)));
    if (kz.im > 0) kz = C(-kz.re, -kz.im);
    const zz = div(scale(m.rho, w), kz);
    const arg = scale(kz, d), c = ccos(arg), s = csin(arg);
    return [c, mul(J, mul(zz, s)), mul(J, div(s, zz)), c];
  }
  const mmul = (a, b) => [add(mul(a[0], b[0]), mul(a[1], b[2])), add(mul(a[0], b[1]), mul(a[1], b[3])),
    add(mul(a[2], b[0]), mul(a[3], b[2])), add(mul(a[2], b[1]), mul(a[3], b[3]))];
  const stack = (ms) => ms.reduce((acc, m) => (acc ? mmul(acc, m) : m), null);
  const admittanceOn = (m, yb) => div(add(m[2], mul(m[3], yb)), add(m[0], mul(m[1], yb)));
  function alphaFromY(y, theta) {
    const c = Math.cos(theta), zy = scale(y, Z0);
    return Math.min(1, Math.max(0, 1 - abs2(div(sub(C(c), zy), add(C(c), zy)))));
  }
  function alphaFree(m, theta) {
    const z = Z0 / Math.cos(theta);
    const a = add(m[0], scale(m[1], 1 / z)), b = add(scale(m[2], z), m[3]);
    const den = add(a, b);
    const r = div(sub(a, b), den), t = div(C(2), den);
    return Math.min(1, Math.max(0, 1 - abs2(r) - abs2(t)));
  }

  // Random-incidence absorption: perFreq(f) prepares what does not depend on angle
  // and returns the absorption as a function of angle.
  function randomIncidence(perFreq) {
    return FREQS.map((f) => {
      const atAngle = perFreq(f);
      return THETA.reduce((acc, t, i) => acc + PARIS[i] * atAngle(t), 0);
    });
  }

  const mid = () => Array.from({ length: N_SUB }, (_, i) => (i + 0.5) / N_SUB);
  const gradedMedia = (f, sigma, fr, fd) => fr.map((s) => mixedMedium(f, sigma, s, fd));
  const graded = (f, t, media, depth) => stack(media.map((m) => layer(f, t, m, depth / media.length)));

  // --------------------------------------------------------------- models -
  function spikes(p) {
    const felt = FELTS[p.porous], d = p.depth / 1000, taper = 1 - p.base;
    const s = mid().map((z) => (taper > 0 ? (z >= taper ? 1 : Math.pow(z / taper, p.taper)) : 1));
    const alpha = randomIncidence((f) => {
      const media = gradedMedia(f, felt.sigma, s);
      return (t) => alphaFromY(admittanceOn(graded(f, t, media, d), C(0)), t);
    });
    return { alpha, mass: mean(s) * d * felt.density, fractions: s };
  }
  function sponge(p) {
    const felt = FELTS[p.porous], d = p.depth / 1000, ds = d * (1 - p.gap);
    const s = mid().map((z) => p.outer + (p.inner - p.outer) * z);
    const alpha = randomIncidence((f) => {
      const media = gradedMedia(f, felt.sigma, s);
      return (t) => {
        let m = graded(f, t, media, ds);
        if (p.gap > 0) m = mmul(m, layer(f, t, AIR, d * p.gap));
        return alphaFromY(admittanceOn(m, C(0)), t);
      };
    });
    return { alpha, mass: mean(s) * ds * felt.density, fractions: s };
  }
  function cloud(p) {
    const felt = FELTS[p.porous], d = p.depth / 1000;
    const half = mid().map((z) => p.outer + (p.inner - p.outer) * z);
    const s = half.concat(half.slice().reverse());
    const alpha = randomIncidence((f) => {
      const media = gradedMedia(f, felt.sigma, s, felt.density);
      return (t) => 2 * alphaFree(graded(f, t, media, d), t);
    });
    return { alpha, mass: mean(s) * d * felt.density, fractions: s };
  }
  function resonatorZ(f, p, depth, fill, sigma) {
    const plate = perforatedPlate(f, p.hole / 1000, p.plate / 1000, p.porosity);
    const m = mmul(layer(f, 0, porousMedium(f, sigma), depth * fill), layer(f, 0, AIR, depth * (1 - fill)));
    return add(plate, inv(admittanceOn(m, C(0))));
  }
  function bottles(p) {
    const felt = FELTS[p.porous], rigid = RIGIDS[p.rigid], d = p.depth / 1000, depths = [d, d * p.ratio];
    const yf = FREQS.map((f) => depths.reduce((acc, x) => add(acc, scale(inv(resonatorZ(f, p, x, p.fill, felt.sigma)), 0.5)), C(0)));
    const alpha = FREQS.map((f, i) => THETA.reduce((acc, t, j) => acc + PARIS[j] * alphaFromY(yf[i], t), 0));
    const md = mean(depths);
    const mass = p.plate / 1000 * (1 - p.porosity) * rigid.density + p.fill * md * felt.density + 2 / PARTITION_SPACING * md * WALL * rigid.density;
    return { alpha, mass, depths };
  }
  function coil(p) {
    const rigid = RIGIDS[p.rigid], d = p.depth / 1000, w = p.width / 1000, a = p.cell / 1000;
    const lMax = p.coil * a * d / (w + WALL);
    const lengths = [0, 1, 2, 3].map((i) => lMax * (1 - p.grading * i / 3));
    const eps = w / a;
    const yf = FREQS.map((f) => {
      const m = slitMedium(f, w), zt = csqrt(mul(m.rho, m.k)), kt = scale(csqrt(div(m.rho, m.k)), 2 * Math.PI * f);
      const end = C(0, 2 * Math.PI * f * RHO0 * 0.6 * w);
      return lengths.reduce((acc, l) => add(acc, scale(inv(add(scale(div(mul(J, zt), ctan(scale(kt, l))), -1), end)), 0.25 * eps)), C(0));
    });
    const alpha = FREQS.map((f, i) => THETA.reduce((acc, t, j) => acc + PARIS[j] * alphaFromY(yf[i], t), 0));
    const mass = (WALL / (w + WALL) * d * p.coil + 0.002 * (1 - eps)) * rigid.density;
    return { alpha, mass, lengths };
  }
  function hybrid(p) {
    const felt = FELTS[p.porous], rigid = RIGIDS[p.rigid], d = p.depth / 1000, front = d * p.front;
    const s = mid().map((z) => p.outer + (p.inner - p.outer) * z);
    const cavity = Math.max(d - front - p.plate / 1000, 0.001);
    const alpha = randomIncidence((f) => {
      const media = gradedMedia(f, felt.sigma, s), y = inv(resonatorZ(f, p, cavity, 0, felt.sigma));
      return (t) => alphaFromY(admittanceOn(graded(f, t, media, front), y), t);
    });
    const mass = mean(s) * front * felt.density + p.plate / 1000 * (1 - p.porosity) * rigid.density + 2 / PARTITION_SPACING * cavity * WALL * rigid.density;
    return { alpha, mass, fractions: s, cavity };
  }
  const mean = (a) => a.reduce((x, y) => x + y, 0) / a.length;

  const DEPTH = ["depth", "depth", [0.15, 1.0]];
  const FAMILIES = [
    { key: "spikes", name: "Spikes and wedges", imitates: "grass, fur, owl-feather velvet", mounting: "wall",
      idea: "Sound slides in between the points and never meets a hard surface to bounce off.",
      genes: [DEPTH, ["base", "float", [0, 0.7]], ["taper", "float", [0.4, 3]], ["porous", "porous"]], model: spikes },
    { key: "sponge", name: "Gradient sponge", imitates: "sea sponge, moss, bone", mounting: "wall",
      idea: "Open and airy at the surface so sound gets in, denser inside so it gets trapped.",
      genes: [DEPTH, ["outer", "float", [0.05, 1]], ["inner", "float", [0.05, 1]], ["gap", "float", [0, 0.6]], ["porous", "porous"]], model: sponge },
    { key: "cloud", name: "Hanging cloud", imitates: "hanging moss, a forest canopy", mounting: "hanging",
      idea: "Hangs free in the air, so it catches sound from both sides.",
      genes: [DEPTH, ["outer", "float", [0.05, 1]], ["inner", "float", [0.05, 1]], ["porous", "porous"]], model: cloud },
    { key: "bottles", name: "Bottle resonators", imitates: "seashells, gourds, hollow stems", mounting: "wall",
      idea: "Air in a small neck bounces on the air in the cavity behind it and swallows a low pitch.",
      genes: [DEPTH, ["hole", "float", [2, 20]], ["plate", "float", [2, 15]], ["porosity", "float", [0.005, 0.2]], ["fill", "float", [0, 1]],
        ["ratio", "float", [0.2, 1]], ["porous", "porous"], ["rigid", "rigid"]], model: bottles },
    { key: "coil", name: "Coiled channels", imitates: "snail shell, the cochlea in your ear", mounting: "wall",
      idea: "A very long, narrow tube folded into a thin object catches low pitches without needing depth.",
      genes: [DEPTH, ["width", "float", [2, 15]], ["cell", "float", [20, 80]], ["grading", "float", [0, 0.9]], ["coil", "float", [0.3, 1]], ["rigid", "rigid"]], model: coil },
    { key: "hybrid", name: "Moss on a shell", imitates: "moss growing on a shell or bark", mounting: "wall",
      idea: "A soft sponge catches the middle and high pitches; bottle cavities behind it catch the low ones.",
      genes: [DEPTH, ["front", "float", [0.2, 0.8]], ["outer", "float", [0.05, 1]], ["inner", "float", [0.05, 1]], ["hole", "float", [2, 20]],
        ["plate", "float", [2, 10]], ["porosity", "float", [0.005, 0.2]], ["porous", "porous"], ["rigid", "rigid"]], model: hybrid },
  ];
  const FAMILY = Object.fromEntries(FAMILIES.map((f) => [f.key, f]));

  // ---------------------------------------------------------- evolution --
  const round = (x, n) => Math.round(x * 10 ** n) / 10 ** n;
  function decode(fam, v, brief) {
    const porous = brief.porous && brief.porous.length ? brief.porous : FELT_KEYS;
    const rigid = brief.rigid && brief.rigid.length ? brief.rigid : RIGID_KEYS;
    const p = {};
    fam.genes.forEach(([name, kind, spec], i) => {
      const u = Math.min(1, Math.max(0, v[i]));
      if (kind === "porous") p[name] = porous[Math.min(Math.floor(u * porous.length), porous.length - 1)];
      else if (kind === "rigid") p[name] = rigid[Math.min(Math.floor(u * rigid.length), rigid.length - 1)];
      else if (kind === "depth") p[name] = round((spec[0] + u * (spec[1] - spec[0])) * brief.maxDepth, 1);
      else p[name] = round(spec[0] + u * (spec[1] - spec[0]), 4);
    });
    return p;
  }
  function assess(fam, params, brief, weights) {
    weights = weights || bandWeights(brief.target);
    const r = fam.model(params);
    const score = r.alpha.reduce((acc, a, i) => acc + a * weights[i], 0);
    const over = r.mass > brief.massBudget;
    const fitness = score - (Math.max(0, r.mass / brief.massBudget - 1) + (over ? 0.1 : 0));
    return { family: fam.key, params, alpha: r.alpha, mass: r.mass, score, fitness, extra: r };
  }

  function rng(seed) { // mulberry32
    let a = seed >>> 0;
    return () => {
      a = (a + 0x6d2b79f5) >>> 0;
      let t = a;
      t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  // A breeder advances one generation per step() call so a page can animate it.
  function breeder(fam, brief, { population = 24, generations = 16, seed = 1, mutation = 0.15, elite = 2 } = {}) {
    const rand = rng(seed + fam.key.length * 101);
    const normal = () => Math.sqrt(-2 * Math.log(1 - rand())) * Math.cos(2 * Math.PI * rand());
    const n = fam.genes.length, weights = bandWeights(brief.target), archive = new Map();
    let pop = Array.from({ length: population }, () => Array.from({ length: n }, rand));
    let gen = 0, best = null;
    const evaluate = (v) => {
      const p = decode(fam, v, brief), key = JSON.stringify(p);
      if (!archive.has(key)) {
        const r = assess(fam, p, brief, weights);
        archive.set(key, r);
        if (!best || r.fitness > best.fitness) best = r;
      }
      return archive.get(key).fitness;
    };
    const tournament = (fit) => {
      let bi = Math.floor(rand() * fit.length);
      for (let k = 1; k < 3; k++) { const i = Math.floor(rand() * fit.length); if (fit[i] > fit[bi]) bi = i; }
      return bi;
    };
    return {
      get best() { return best; },
      get generation() { return gen; },
      get done() { return gen > generations; },
      generations,
      step() {
        const fit = pop.map(evaluate);
        gen += 1;
        if (gen > generations) return best;
        const order = fit.map((f, i) => i).sort((a, b) => fit[b] - fit[a]);
        const next = order.slice(0, elite).map((i) => pop[i].slice());
        while (next.length < population) {
          const a = pop[tournament(fit)], b = pop[tournament(fit)];
          next.push(a.map((x, i) => {
            let g = rand() < 0.5 ? x : b[i];
            if (rand() < mutation) g = Math.min(1, Math.max(0, g + 0.15 * normal()));
            if (rand() < mutation / 5) g = rand();
            return g;
          }));
        }
        pop = next;
        return best;
      },
    };
  }

  // Plain flat layer on a wall: the fair way to compare materials on their own.
  function flatLayer(key, thicknessMm) {
    const felt = FELTS[key], d = thicknessMm / 1000;
    const alpha = randomIncidence((f) => {
      const m = porousMedium(f, felt.sigma);
      return (t) => alphaFromY(admittanceOn(layer(f, t, m, d), C(0)), t);
    });
    return { alpha, mass: d * felt.density };
  }

  const quarterWave = (depthMm) => C0 / (4 * depthMm / 1000);

  const Engine = { FREQS, TARGETS, FELTS, RIGIDS, FAMILIES, FAMILY, N_SUB, bandWeights, decode, assess, breeder, flatLayer, quarterWave };
  if (typeof module !== "undefined" && module.exports) module.exports = Engine;
  else root.Engine = Engine;
})(this);
