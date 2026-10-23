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

// Every rule grows inside its envelope, hanging and standing.
for (const rule of Object.keys(G.RULES)) for (const mount of ["hanging", "standing"]) {
  const genes = Array.from({ length: G.geneCount(rule) }, (_, i) => (i * 0.37) % 1);
  const prims = G.grow(rule, genes, { mount, size: 1.2 });
  const { min, max } = G.bounds(prims), env = G.envelope(1.2, mount);
  assert.ok(prims.length >= 3 && prims.length < 6000, `${rule} primitive count ${prims.length}`);
  const slack = G.RULES[rule].group === "sound" ? 0.5 : G.RULES[rule].group === "hybrid" || mount === "standing" ? 0.081 : 1e-6; // membrane feet reach out sideways // sheets keep their real thickness when scaled
  assert.ok(max[0] - min[0] <= env.width + slack && max[2] - min[2] <= env.height + slack, `${rule} ${mount} fits`);
  if (mount === "hanging") close(max[2], 0, slack, `${rule} hangs from the ceiling`);
  else close(min[2], 0, slack, `${rule} stands on the floor`);
}

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
