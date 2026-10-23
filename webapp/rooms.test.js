// Run: node webapp/rooms.test.js
const assert = require("assert");
const Rooms = require("./rooms.js");

const close = (a, b, tol, what) => assert.ok(Math.abs(a - b) <= tol, `${what}: ${a} vs ${b}`);

// A 10 x 8 x 4 m box: volume, areas.
const info = Rooms.analyze(Rooms.boxMesh(10, 8, 4));
close(info.volume, 320, 1e-9, "box volume");
close(info.areas.floor, 80, 1e-9, "floor");
close(info.areas.ceiling, 80, 1e-9, "ceiling");
close(info.areas.wall, 144, 1e-9, "walls");

// The same box written as an OBJ in millimetres with Y up, as many CAD tools export it.
const m = Rooms.boxMesh(10, 8, 4);
let obj = "", n = 0;
for (const t of m.tris) {
  for (const p of t) obj += `v ${p[0] * 1000} ${p[2] * 1000} ${-p[1] * 1000}\n`;
  obj += `f ${n + 1} ${n + 2} ${n + 3}\n`;
  n += 3;
}
const parsed = Rooms.parseOBJ(obj, { yUp: true });
const pinfo = Rooms.analyze(parsed);
close(parsed.unit, 0.001, 0, "millimetres detected");
close(pinfo.volume, 320, 1e-6, "OBJ volume");
close(pinfo.areas.floor, 80, 1e-6, "OBJ floor classified");
close(pinfo.areas.wall, 144, 1e-6, "OBJ walls classified");

// Sabine check against a hand calculation at 500 Hz.
const space = { volume: 320, areas: info.areas, finishes: { floor: "concrete", wall: "masonry", ceiling: "plaster" }, people: 0 };
const r = Rooms.reverb(space);
const a500 = 80 * 0.02 + 144 * 0.02 + 80 * 0.02 + 0.0026 * 320;
close(r.rtBefore[2], 0.161 * 320 / a500, 1e-9, "Sabine RT at 500 Hz");

// Adding the area the helper asks for reaches the target.
const alpha = [0.3, 0.6, 0.9, 0.95, 0.95, 0.95];
const area = Rooms.areaNeeded(space, alpha);
const after = Rooms.reverb(space, [{ area, alpha }]);
close(after.speechRtAfter, after.target, 1e-9, "area needed reaches target");
assert.ok(after.noiseDrop > 3, "absorbers lower the crowd noise");
close(Rooms.targetRT(1000), 0.32 * 3 - 0.17, 1e-12, "target RT");
console.log(`rooms ok: ${area.toFixed(1)} m2 needed in a 320 m3 concrete box, noise -${after.noiseDrop.toFixed(1)} dB`);
