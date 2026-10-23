/* Biothing Sound Lab rooms: exhibition spaces as 3D meshes, and their reverberation.
 * Statistical room acoustics (Sabine): good for comparing options in a space,
 * not a substitute for a ray-tracing study of a specific hall. */
(function (root) {
  "use strict";

  const OCTAVES = [125, 250, 500, 1000, 2000, 4000];
  // Air absorption 4m (1/m) at 20 C, 50 % RH (from ISO 9613-1 attenuation).
  const AIR_4M = [0.0004, 0.0012, 0.0026, 0.0045, 0.009, 0.027];
  // Equivalent absorption area per standing visitor, m2 (typical values).
  const PERSON = [0.15, 0.3, 0.5, 0.7, 0.8, 0.8];

  // Typical absorption coefficients per octave 125-4000 Hz (indicative textbook values).
  const SURFACES = {
    concrete: { name: "Concrete or screed", use: "floor", alpha: [0.01, 0.01, 0.02, 0.02, 0.02, 0.03] },
    stone: { name: "Stone or tiles", use: "floor", alpha: [0.01, 0.01, 0.01, 0.02, 0.02, 0.02] },
    wood_floor: { name: "Wooden floor", use: "floor", alpha: [0.15, 0.11, 0.1, 0.07, 0.06, 0.07] },
    carpet: { name: "Carpet", use: "floor", alpha: [0.02, 0.06, 0.14, 0.37, 0.6, 0.65] },
    masonry: { name: "Painted concrete or brick", use: "wall", alpha: [0.01, 0.01, 0.02, 0.02, 0.02, 0.03] },
    plasterboard: { name: "Plasterboard on studs", use: "wall", alpha: [0.29, 0.1, 0.05, 0.04, 0.07, 0.09] },
    glass: { name: "Large glass panes", use: "wall", alpha: [0.18, 0.06, 0.04, 0.03, 0.02, 0.02] },
    wood_panel: { name: "Wood panelling", use: "wall", alpha: [0.42, 0.21, 0.1, 0.08, 0.06, 0.06] },
    curtain: { name: "Heavy curtains", use: "wall", alpha: [0.07, 0.31, 0.49, 0.75, 0.7, 0.6] },
    plaster: { name: "Plaster ceiling", use: "ceiling", alpha: [0.01, 0.01, 0.02, 0.03, 0.04, 0.05] },
    metal_deck: { name: "Metal roof deck", use: "ceiling", alpha: [0.25, 0.15, 0.1, 0.08, 0.07, 0.07] },
    acoustic_tiles: { name: "Acoustic ceiling tiles", use: "ceiling", alpha: [0.5, 0.7, 0.6, 0.7, 0.7, 0.5] },
  };

  const PRESETS = [
    { key: "gallery", name: "Gallery room", size: [12, 8, 4.5], floor: "wood_floor", wall: "plasterboard", ceiling: "plaster", people: 30,
      note: "A classic white gallery with hard walls." },
    { key: "white_cube", name: "White cube", size: [20, 12, 5], floor: "concrete", wall: "masonry", ceiling: "plaster", people: 60,
      note: "Concrete floor, painted masonry: very echoey." },
    { key: "foyer", name: "Glass foyer", size: [30, 6, 4], floor: "stone", wall: "glass", ceiling: "plaster", people: 80,
      note: "A long entrance hall with a glass facade." },
    { key: "hall", name: "Exhibition hall", size: [40, 30, 10], floor: "concrete", wall: "plasterboard", ceiling: "metal_deck", people: 300,
      note: "A big fair hall under a metal roof." },
  ];

  // ---------------------------------------------------------------- meshes --
  // A mesh is {tris: [[x,y,z]x3, ...], tags: ["floor"|"wall"|"ceiling", ...]}, metres, z up.
  function boxMesh(w, d, h) {
    const p = (x, y, z) => [x, y, z];
    const quads = [
      ["floor", [p(0, 0, 0), p(0, d, 0), p(w, d, 0), p(w, 0, 0)]],
      ["ceiling", [p(0, 0, h), p(w, 0, h), p(w, d, h), p(0, d, h)]],
      ["wall", [p(0, 0, 0), p(w, 0, 0), p(w, 0, h), p(0, 0, h)]],
      ["wall", [p(w, 0, 0), p(w, d, 0), p(w, d, h), p(w, 0, h)]],
      ["wall", [p(w, d, 0), p(0, d, 0), p(0, d, h), p(w, d, h)]],
      ["wall", [p(0, d, 0), p(0, 0, 0), p(0, 0, h), p(0, d, h)]],
    ];
    const tris = [], tags = [];
    for (const [tag, q] of quads) {
      tris.push([q[0], q[1], q[2]], [q[0], q[2], q[3]]);
      tags.push(tag, tag);
    }
    return { tris, tags };
  }

  const subv = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
  const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const len = (a) => Math.hypot(a[0], a[1], a[2]);

  // Parse an OBJ file (vertices and faces only). Faces are fan-triangulated.
  // Units: if the model is larger than 500 units it is assumed to be in millimetres,
  // larger than 60 units in centimetres. `yUp` converts Y-up exports to z-up.
  function parseOBJ(text, { yUp = false } = {}) {
    const v = [], tris = [];
    for (const raw of text.split(/\r?\n/)) {
      const line = raw.trim();
      if (line.startsWith("v ")) {
        const [x, y, z] = line.slice(2).trim().split(/\s+/).map(Number);
        v.push(yUp ? [x, -z, y] : [x, y, z]);
      } else if (line.startsWith("f ")) {
        const idx = line.slice(2).trim().split(/\s+/).map((t) => {
          const i = parseInt(t.split("/")[0], 10);
          return i < 0 ? v.length + i : i - 1;
        });
        for (let k = 1; k + 1 < idx.length; k++) tris.push([v[idx[0]], v[idx[k]], v[idx[k + 1]]]);
      }
    }
    if (!tris.length || tris.some((t) => t.some((p) => !p))) throw new Error("No valid faces found in this OBJ file.");
    let min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
    for (const t of tris) for (const p of t) for (let i = 0; i < 3; i++) { min[i] = Math.min(min[i], p[i]); max[i] = Math.max(max[i], p[i]); }
    const span = Math.max(...max.map((m, i) => m - min[i]));
    const unit = span > 500 ? 0.001 : span > 60 ? 0.01 : 1;
    const moved = tris.map((t) => t.map((p) => [(p[0] - min[0]) * unit, (p[1] - min[1]) * unit, (p[2] - min[2]) * unit]));
    const h = (max[2] - min[2]) * unit;
    const tags = moved.map((t) => {
      const n = cross(subv(t[1], t[0]), subv(t[2], t[0])), l = len(n) || 1, nz = n[2] / l;
      const cz = (t[0][2] + t[1][2] + t[2][2]) / 3;
      if (Math.abs(nz) < 0.7) return "wall";
      return cz < h / 2 ? "floor" : "ceiling";
    });
    return { tris: moved, tags, unit };
  }

  function analyze(mesh) {
    let vol = 0;
    const areas = { floor: 0, wall: 0, ceiling: 0 };
    const min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
    mesh.tris.forEach((t, i) => {
      const [a, b, c] = t;
      vol += (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6;
      areas[mesh.tags[i]] += len(cross(subv(b, a), subv(c, a))) / 2;
      for (const p of t) for (let k = 0; k < 3; k++) { min[k] = Math.min(min[k], p[k]); max[k] = Math.max(max[k], p[k]); }
    });
    return { volume: Math.abs(vol), areas, total: areas.floor + areas.wall + areas.ceiling, min, max };
  }

  // ------------------------------------------------------------- acoustics --
  // One-third-octave absorption (100-5000 Hz, 18 bands) to octaves 125-4000.
  function thirdToOctave(alpha18) {
    return OCTAVES.map((_, i) => (alpha18[3 * i] + alpha18[3 * i + 1] + alpha18[3 * i + 2]) / 3);
  }

  // Recommended reverberation for rooms where people talk (after DIN 18041, group A3).
  const targetRT = (volume) => Math.max(0.3, 0.32 * Math.log10(volume) - 0.17);

  // space: {volume, areas, finishes: {floor, wall, ceiling}, people}
  // absorbers: [{area, alpha: [6]}] where area * alpha is the added absorption per octave.
  function reverb(space, absorbers = []) {
    const base = OCTAVES.map((_, i) =>
      space.areas.floor * SURFACES[space.finishes.floor].alpha[i]
      + space.areas.wall * SURFACES[space.finishes.wall].alpha[i]
      + space.areas.ceiling * SURFACES[space.finishes.ceiling].alpha[i]
      + space.people * PERSON[i]
      + AIR_4M[i] * space.volume);
    const added = OCTAVES.map((_, i) => absorbers.reduce((acc, a) => acc + a.area * a.alpha[i], 0));
    const rt = (A) => A.map((a) => 0.161 * space.volume / a);
    const after = base.map((a, i) => a + added[i]);
    const speech = (A) => (A[2] + A[3]) / 2; // 500-1000 Hz, where voices carry
    return {
      absorptionBefore: base, absorptionAfter: after,
      rtBefore: rt(base), rtAfter: rt(after),
      target: targetRT(space.volume),
      speechRtBefore: 0.161 * space.volume / speech(base), speechRtAfter: 0.161 * space.volume / speech(after),
      noiseDrop: 10 * Math.log10(speech(after) / speech(base)), // drop of the reverberant (crowd) level, dB
    };
  }

  // m2 of an absorber with octave absorption `alpha` needed to reach the target at 500-1000 Hz.
  function areaNeeded(space, alpha, existing = []) {
    const r = reverb(space, existing);
    const need = 0.161 * space.volume / r.target - (r.absorptionAfter[2] + r.absorptionAfter[3]) / 2;
    const per = (alpha[2] + alpha[3]) / 2;
    return need <= 0 ? 0 : need / per;
  }

  const Rooms = { OCTAVES, SURFACES, PRESETS, PERSON, boxMesh, parseOBJ, analyze, thirdToOctave, targetRT, reverb, areaNeeded };
  if (typeof module !== "undefined" && module.exports) module.exports = Rooms;
  else root.Rooms = Rooms;
})(this);
