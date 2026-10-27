// Run: node webapp/growth.test.js
const assert = require("assert");
const G = require("./growth.js");

const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);

// Voxelised volumes of a sphere and a thin rod match their real volume.
const sphere = G.voxelize([{ a: [0, 0, 0], b: [0, 0, 0], r: 0.2 }], 0.02);
close(sphere.volume, (4 / 3) * Math.PI * 0.2 ** 3, 0.03 * (4 / 3) * Math.PI * 0.2 ** 3, "sphere volume");
const rod = G.voxelize([{ a: [0, 0, 0], b: [1, 0.3, 0.2], r: 0.005 }], 0.043);
const rodVol = Math.PI * 0.005 ** 2 * Math.hypot(1, 0.3, 0.2);
close(rod.volume, rodVol, 0.25 * rodVol, "thin rod volume");

// Sheets: a thick and a thin felt sheet keep their volume, and pack/unpack round-trips.
const thick = G.voxelize([{ c: [0, 0, 0], u: [0.4, 0.1, 0], v: [-0.05, 0.3, 0.1], t: 0.06 }], 0.043);
const sheetArea = (u, v) => Math.hypot(u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]) * 4;
const tv = sheetArea([0.4, 0.1, 0], [-0.05, 0.3, 0.1]) * 0.06;
close(thick.volume, tv, 0.15 * tv, "thick sheet volume");
const thin = G.voxelize([{ c: [0, 0, 0], u: [0.5, 0, 0], v: [0, 0.35, 0.2], t: 0.015 }], 0.057);
const thv = sheetArea([0.5, 0, 0], [0, 0.35, 0.2]) * 0.015;
close(thin.volume, thv, 0.15 * thv, "thin sheet volume");
const parts = [{ a: [0, 0, 0], b: [1, 2, 3], r: 0.1 }, { c: [1, 1, 1], u: [0.2, 0, 0], v: [0, 0.3, 0], t: 0.03 }];
const back = G.unpack(G.pack(parts));
close(back[1].t, 0.03, 1e-6, "pack keeps sheets"); close(back[0].b[2], 3, 1e-6, "pack keeps capsules");

// Every rule grows inside its envelope: on a wall, hanging from the ceiling and standing.
for (const rule of Object.keys(G.RULES)) for (const mount of ["hanging", "standing", "wall"]) {
  const genes = Array.from({ length: G.geneCount(rule) }, (_, i) => (i * 0.37) % 1);
  const prims = G.grow(rule, genes, { mount, size: 1.2 });
  const { min, max } = G.bounds(prims), env = G.envelope(1.2, mount), R = G.RULES[rule];
  assert.ok(prims.length >= 3 && prims.length < 80000, `${rule} ${mount} part count ${prims.length}`);
  const slack = R.group === "sound" ? 0.5 : R.lab ? 0.1 : R.group === "hybrid" || mount !== "hanging" ? 0.081 : 1e-6; // feet reach out; sheets keep their real thickness when scaled
  assert.ok(max[0] - min[0] <= env.width + slack && max[2] - min[2] <= env.height + slack + (R.lab && mount === "hanging" ? 0.3 : 0), `${rule} ${mount} fits`);
  if (mount === "hanging") close(max[2], 0, slack, `${rule} hangs from the ceiling`);
  else if (mount === "standing") close(min[2], 0, slack, `${rule} stands on the floor`);
  else assert.ok(min[1] >= -0.02 && max[1] <= env.depth + 0.12 && min[2] >= -0.05, `${rule} sits on the wall (${min[1].toFixed(2)}..${max[1].toFixed(2)} m deep)`);
}
// Wall pieces are tested against the wall; a wall evaluation runs.
const wall = G.evaluate({ rule: "voronoi", genes: Array.from({ length: G.geneCount("voronoi") }, () => 0.5), mount: "wall", size: 1.2, material: { sigma: 20000, density: 60 }, budget: 5 });
assert.ok(wall.score > 0 && !wall.stance && wall.holdfast, "wall module evaluates");
// Lattices: struts and shells are print (structure); only their felt infill absorbs.
const lat = (rule, fill) => { const g = Array.from({ length: G.geneCount(rule) }, () => 0.5); g[7] = fill; return G.evaluate({ rule, genes: g, mount: "wall", size: 1, material: { sigma: 20000, density: 60 }, budget: 5, structureDensity: 1240 }); };
const empty0 = lat("strutlattice", 0), filled = lat("strutlattice", 1);
assert.ok(empty0.score < 0.02 && filled.score > empty0.score && filled.printMass > 0, `printed struts alone barely absorb (${empty0.score.toFixed(3)} vs ${filled.score.toFixed(2)} m2 with felt)`);
const gy = G.grow("tpms", Array.from({ length: G.geneCount("tpms") }, () => 0.2), { mount: "wall", size: 1 });
assert.ok(gy.filter(G.isTri).length > 1000, "gyroid shell is meshed");
const gyBack = G.unpack(G.pack(gy));
assert.ok(gyBack.filter(G.isTri).length === gy.filter(G.isTri).length && gyBack.find(G.isTri).structure, "triangles survive packing");
// Parameters map genes to real values; geometric fields are normalised.
close(G.paramValue("voronoi", [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], 2), 0.2, 1e-9, "parameter value");
for (const kind of ["attractor", "gradient", "radial", "waves", "noise"]) { const pf = G.patternField({ kind }); assert.ok(Math.max(...pf.data) === 1 && Math.min(...pf.data) === 0, `${kind} field`); }
console.log(`lab ok: wall voronoi ${wall.area.map((a) => a.toFixed(2)).join("/")} m2, gyroid ${gy.length} parts`);

// The wave solver: an empty room keeps its energy; a porous block absorbs, more material absorbs more.
const room = G.roomFor("standing");
const empty = { grid: new Float32Array(1), nh: 1, nz: 1, dx: 0.086, x0: 0, z0: 0, extent: 1 };
const none = G.absorptionLength(empty, 20000, 500, room, "standing");
assert.ok(none.length < 0.02, `empty room barely decays (${none.length})`);
const block = (fill) => ({ grid: new Float32Array(10 * 12).fill(fill), nh: 10, nz: 12, dx: 0.086, x0: -0.43, z0: 0, extent: 1 });
const light = G.absorptionLength(block(0.1), 20000, 500, room, "standing");
const dense = G.absorptionLength(block(0.4), 20000, 500, room, "standing");
assert.ok(light.length > 0.1, `porous block absorbs (${light.length.toFixed(2)} m)`);
assert.ok(dense.length > light.length, "more material absorbs more here");
// Order of magnitude: a 0.86 m x 1 m porous block cannot absorb more than its perimeter.
assert.ok(dense.length < 2 * (0.86 + 1.03) * 1.2, `plausible absorption length (${dense.length.toFixed(2)} m)`);

const t0 = Date.now();
const r = G.evaluate({ rule: "tree", genes: [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0], mount: "hanging", size: 1.2,
  material: { sigma: 20000, density: 60 }, budget: 5 });
assert.ok(r.score > 0 && r.mass > 0, "evaluation runs");
console.log(`growth ok: tree ${r.prims.length} parts, ${r.mass.toFixed(2)} kg, ${r.area.map((a) => a.toFixed(2)).join("/")} m2 at 250/500/1000 Hz, ${Date.now() - t0} ms`);

// Freestanding: standing forms get a holdfast of structure parts, and a stability verdict.
for (const rule of ["snow", "leaves"]) {
  const genes = Array.from({ length: G.geneCount(rule) }, () => 0.5);
  const r = G.evaluate({ rule, genes, mount: "standing", size: 1.2, material: { sigma: 20000, density: 60 }, budget: 5 });
  assert.ok(r.stance && r.prims.some((p) => p.structure), `${rule} has a holdfast`);
  assert.ok(Number.isFinite(r.stance.angle) && r.stance.structureMass > 0, `${rule} stability computed`);
  const back = G.unpack(G.pack(r.prims));
  assert.ok(back.filter((p) => p.structure).length === r.prims.filter((p) => p.structure).length, "structure survives packing");
}
// A wide heavy plinth under a light form is stable; a tiny base under a tall heavy form is not.
const tall = [{ a: [0, 0, 0.2], b: [0, 0, 2.0], r: 0.2 }];
const plinth = (R) => [...tall, { c: [0, 0, 0.01], u: [R, 0, 0], v: [0, R, 0], t: 0.02, structure: true }];
const vt = G.voxelize(tall, 0.05);
assert.ok(G.stability(plinth(0.6), vt.com, vt.volume * 60, 680).ok, "wide plinth stands");
assert.ok(!G.stability(plinth(0.08), vt.com, vt.volume * 60, 680).ok, "narrow base topples");
console.log("freestanding checks ok");

// Sound fields: tone bursts, a recorded signal and room modes all give a normalised field.
const roomSec = { w: 6, h: 3.5, cx: 3.4 };
const fa = G.fieldFromSound({ room: roomSec, x0: 3.4, width: 1.2, height: 1.8, src: [1.4, 1.5] });
assert.ok(fa.data.length === 24 * 36 && Math.max(...fa.data) === 1 && Math.min(...fa.data) === 0, "tone-burst field normalised");
const dtS = G.solverStep(), sig = Float32Array.from({ length: Math.round(0.05 / dtS) }, (_, i) => Math.sin(2 * Math.PI * 300 * i * dtS) * Math.exp(-i * dtS * 40));
const fb = G.fieldFromSound({ room: roomSec, signal: sig, x0: 3.4, width: 1.2, height: 1.8, src: [1.4, 1.5] });
assert.ok(Math.max(...fb.data) === 1, "signal-driven field");
const fc = G.fieldFromModes({ L: 12, W: 8, H: 4.5, x0: 7, y0: 3, width: 1.2, height: 1.8 });
assert.ok(fc.modes.length > 20 && fc.modes.every((m) => m.f >= 70 && m.f <= 200), `room modes in the voice range (${fc.modes.length})`);
const m0 = fc.modes[0];
close(m0.f, (343.2 / 2) * Math.hypot(m0.n[0] / 12, m0.n[1] / 8, m0.n[2] / 4.5), 1e-9, "mode frequency formula");
const sd = (d) => { const mu = d.reduce((x, y) => x + y, 0) / d.length; return Math.sqrt(d.reduce((x, y) => x + (y - mu) ** 2, 0) / d.length); };
assert.ok(sd(fc.data) > 0.15, "the resonance field has a pattern across the piece");
// The membrane follows its field: more cells where the field is strong.
const left = { data: Float32Array.from({ length: 24 * 36 }, (_, i) => ((i % 24) < 12 ? 1 : 0)), nx: 24, nz: 36 };
const mem = G.grow("membrane", Array.from({ length: G.geneCount("membrane") }, () => 0.5), { mount: "standing", size: 1.2, field: left });
const tiles = mem.filter((q) => G.isSheet(q) && !q.structure);
const onLeft = tiles.filter((q) => q.c[0] < 0).length;
assert.ok(onLeft > tiles.length * 0.6, `membrane is denser where sound is strong (${onLeft}/${tiles.length} cells on the loud side)`);
assert.ok(mem.label && mem.field, "membrane reports its holdfast and final field");
const st = G.evaluate({ rule: "membrane", genes: Array.from({ length: G.geneCount("membrane") }, () => 0.5), mount: "standing", size: 1.2, material: { sigma: 20000, density: 60 }, budget: 6, field: fa });
assert.ok(st.stance && Number.isFinite(st.stance.angle) && st.stance.across > 0.2, `membrane stability (${st.stance.angle.toFixed(1)}°)`);
console.log(`genware ok: ${tiles.length} cells, ${st.mass.toFixed(2)} kg, ${st.area.map((a) => a.toFixed(2)).join("/")} m2, tips at ${st.stance.angle.toFixed(0)}°`);

// Interference cut-sheet: flaps open where the field is strong; soft light felt sags, stiff felt holds.
const csGenes = Array.from({ length: G.geneCount("cutsheet") }, () => 0.5);
csGenes[8] = 0.2; // top hinges, so gravity matters
const cs = (density) => G.grow("cutsheet", csGenes, { mount: "standing", size: 1.2, field: left, material: { density, sigma: 20000 } });
const stiff = cs(250), soft = cs(40), plan = stiff.cutPlan;
const flapsLeft = plan.cells.filter((c) => c.w > 0 && c.s < 0).length, flapsRight = plan.cells.filter((c) => c.w > 0 && c.s > 0).length;
assert.ok(flapsLeft > flapsRight, `cut-sheet cuts where the sound is strong (${flapsLeft} vs ${flapsRight} flaps)`);
assert.ok(soft.cutPlan.droop > stiff.cutPlan.droop, `soft felt sags more (${soft.cutPlan.droop.toFixed(2)} vs ${stiff.cutPlan.droop.toFixed(2)})`);
assert.ok(plan.cutLength > 0 && stiff.label, "cut length and holdfast reported");
const ce = G.evaluate({ rule: "cutsheet", genes: csGenes, mount: "standing", size: 1.2, material: { sigma: 20000, density: 120 }, budget: 6, field: fa });
assert.ok(ce.cutPlan && ce.score > 0 && ce.stance, "cut-sheet evaluates with its cut plan");
console.log(`cut-sheet ok: ${plan.flaps} flaps, ${plan.cutLength.toFixed(1)} m of cuts, ${ce.mass.toFixed(2)} kg, ${ce.area.map((a) => a.toFixed(2)).join("/")} m2`);

// Skin: a second growth on a form, densest where its own field is strong, and counted in the evaluation.
const host = G.grow("leaves", Array.from({ length: G.geneCount("leaves") }, () => 0.5), { mount: "standing", size: 1.2 });
const tufts = G.growSkin(host, left, [0.8, 0.5, 0.5, 0.5], G.envelope(1.2, "standing"), "standing");
assert.ok(tufts.length > 20, `skin grows tufts (${tufts.length})`);
assert.ok(tufts.filter((q) => q.a[0] < 0).length > tufts.length * 0.6, "skin is densest where its field is strong");
const skinned = G.evaluate({ rule: "leaves", genes: Array.from({ length: G.geneCount("leaves") }, () => 0.5), mount: "standing", size: 1.2, material: { sigma: 20000, density: 60 }, budget: 8, skin: { field: fa, genes: [0.8, 0.5, 0.5, 0.5] } });
const bare = G.evaluate({ rule: "leaves", genes: Array.from({ length: G.geneCount("leaves") }, () => 0.5), mount: "standing", size: 1.2, material: { sigma: 20000, density: 60 }, budget: 8 });
assert.ok(skinned.skinCount > 0 && skinned.mass > bare.mass, "skin adds felt");

// Placement map: no spots inside the walls, layouts keep their spacing and favour the map.
const pm = G.placementMap({ lx: 12, ly: 8, H: 4.5, need: [0.3, 0.3, 0.2, 0.1, 0.05, 0.05] });
assert.ok(pm.data[0] === 0 && Math.max(...pm.data) === 1, "placement map normalised, walls excluded");
const lay = G.sampleLayout(pm, 10, 3, 1.5);
assert.ok(lay.spots.length === 10 && lay.gain > 1, `layout drawn from the map (${lay.gain.toFixed(2)}x an even spread)`);
for (let i = 0; i < lay.spots.length; i++) for (let j = i + 1; j < lay.spots.length; j++) assert.ok(Math.hypot(lay.spots[i][0] - lay.spots[j][0], lay.spots[i][1] - lay.spots[j][1]) >= lay.gap - 1e-9, "spacing kept");
console.log(`skin and placement ok: ${tufts.length} tufts, layout ${lay.gain.toFixed(2)}x`);
