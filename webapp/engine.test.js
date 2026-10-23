// Checks the browser engine against reference results from the Python engine.
// Run: python3 webapp/make_reference.py && node webapp/engine.test.js
const assert = require("assert");
const Engine = require("./engine.js");
const ref = require("./reference.json");

const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), `${what}: ${a} vs ${b}`);

for (const [target, w] of Object.entries(ref.weights)) {
  Engine.bandWeights(target).forEach((x, i) => close(x, w[i], 1e-9, `weights ${target}[${i}]`));
}
const brief = { target: "exhibition", maxDepth: 200, massBudget: 6 };
for (const [key, c] of Object.entries(ref.cases)) {
  const r = Engine.assess(Engine.FAMILY[key], c.params, brief);
  r.alpha.forEach((a, i) => close(a, c.alpha[i], 1e-6, `${key} alpha[${Engine.FREQS[i]} Hz]`));
  close(r.mass, c.mass, 1e-9, `${key} mass`);
  close(r.score, c.score, 1e-6, `${key} score`);
}
const flat = Engine.flatLayer("hemp_fibre", 50);
const asSponge = Engine.assess(Engine.FAMILY.sponge, { depth: 50, outer: 1, inner: 1, gap: 0, porous: "hemp_fibre" }, brief);
flat.alpha.forEach((a, i) => close(a, asSponge.alpha[i], 1e-9, `flat layer [${i}]`));

const t0 = Date.now();
const b = Engine.breeder(Engine.FAMILY.hybrid, brief, { population: 24, generations: 16 });
while (!b.done) b.step();
assert.ok(b.best.fitness > 0.3, "breeding finds a reasonable hybrid");
console.log(`engine matches Python; one family bred in ${Date.now() - t0} ms`);
