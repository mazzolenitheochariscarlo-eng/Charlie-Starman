# Biothing acoustics: generative desk screens for temporary exhibitions

A Python pipeline that *grows* sound-absorbing desk screens instead of drawing
them, following the generative logic of Alisa Andrasek's Biothing: the designer
writes a brief and a genome, and an evolutionary algorithm breeds the form.

Each screen is a rigid core with a sculpted felt face on both sides. One felt
sheet is formed over a relief made of cells:

- **low cells** press the felt: denser, higher flow resistivity, absorbs mid and high frequencies;
- **high cells** let the felt stand off the core: an air cavity behind it pulls absorption down towards speech frequencies.

Natural growth patterns (phyllotaxis, relaxed Voronoi foam, honeycomb) and
relief fields (radial, vertical, gyroid, noise) spread these conditions over the
face. That gives every cell slightly different tuning, the way the multi-scale
scales of a moth wing absorb sound over a broad range.

The screens belong to a **modular exhibition system**: a small family of module
widths that tile common desk sizes, dry joints only (no glue or screws),
one-person module mass, flat-pack transport and reuse across many events.

## Quick start

```bash
pip install -r requirements.txt
python -m biothing_acoustics run                                     # default exhibition brief
python -m biothing_acoustics run --brief briefs/pet_mono_material.json --out output/pet
python -m biothing_acoustics evaluate output/design.json             # re-score an edited design
python -m biothing_acoustics catalogue                               # materials, joints, supports
python -m pytest
```

A run takes about 30 s and writes to `output/`:

| File | Content |
|---|---|
| `report.md` | Winning design, acoustic ratings, exhibition bill of materials, shortlist of alternatives |
| `module_<w>x<h>.svg` | 1:1 face pattern for each module width used; cell outlines are the forming lines for the mould or CNC relief, and darker cells stand further off the core |
| `module_<w>x<h>.obj` | 3D mesh of the module (core and both relief faces) in mm, for Rhino or Blender |
| `bill_of_materials.csv` | Screen runs per desk island, module counts, joints and supports |
| `design.json` | Genome and all metrics, editable and re-scorable with `evaluate` |
| `*.png` | Absorption curve, speech privacy, evolution and face relief plots |

`examples/pet_mono_material/` holds the output of the PET mono-material brief.

## How it works

```
brief.json ──► genome (12 genes) ──► geometry ──► acoustics ──► speech privacy ──► fitness
                    ▲                cells,        Miki model,    barrier            │
                    │                relief,       cavities,      diffraction,       │
                    └──── genetic algorithm ◄──────────────────────────────────────┘
                                         │
                                         ▼
                 modular layout ──► bill of materials ──► SVG / OBJ / report
```

| Module | Role |
|---|---|
| `geometry.py` | Seed patterns, Voronoi cells (half-plane clipping), relief fields |
| `acoustics.py` | Miki porous model, felt compression, air cavities, random-incidence absorption, ISO 11654 alpha_w and class, NRC, mass law, leaks through slits |
| `speech.py` | Diffraction over the top and around the run ends, transmission through the panel, ISO 3382-3 speech spectrum, simplified speech intelligibility index |
| `design.py` | Genome, decoding, evaluation, fitness |
| `optimize.py` | Genetic algorithm (tournament selection, uniform crossover, Gaussian mutation, elitism) |
| `modular.py` | Module tiling, desk-island layouts, joints, supports, bill of materials |
| `materials.py` | Felt and core library |
| `brief.py` | Scenario, constraints, weights, layout and material whitelists |

### Genome

Pattern, cell density, relief field and exponent, felt type, sheet thickness,
minimum compression, cavity depth, core, height above the desk, joint and support.

### Fitness (weights set in the brief)

- weighted absorption coefficient alpha_w, and mean absorption at 500-2000 Hz;
- speech privacy: the drop in speech intelligibility at the facing desk;
- lightness, cost and embodied carbon;
- penalties for modules heavier than one person should handle, or screens thicker than the limit.

## Writing your own brief

Copy `briefs/exhibition_desk_screen.json` and edit it. The most useful settings:

- `layout`: your desk islands, each either `bench` (facing desks: a spine plus dividers) or `row` (a counter facing visitors: a front run plus dividers);
- `constraints`: module mass, thickness, minimum and maximum height above the desk;
- `weights`: what matters most;
- `allowed_felts` / `allowed_cores` / `allowed_joints` / `allowed_supports`: restrict the catalogue (see `pet_mono_material.json`);
- `scenario.background_noise`: the hall's noise spectrum (a quieter hall makes privacy harder).

## What the first runs show

- Screen height drives privacy most. With seated talkers about 1 m from the screen, a 650 mm screen (top at 1.39 m) cuts speech by about 6 dB(A). Low screens mainly absorb, and don't block much.
- The gap between screen and desk matters. Desk clamps that leave 20 mm open lose sound reduction compared with feet or saddles.
- Deep cavities on both faces raise absorption but quickly reach the thickness limit. The optimiser trades relief depth against thickness, mass and cost.
- With the default weights, a cheap, low-carbon flat felt can beat a formed one. Raise `weights.absorption`, or restrict to formable felts, if the relief is part of the concept.

## Limits

- Absorption is predicted for an infinite surface. Real panels measured in a reverberation room (ISO 354) usually score higher because of edge effects. Use the numbers to compare designs, then measure.
- Speech privacy includes direct, diffracted and transmitted paths only. Ceiling and wall reflections are not modelled.
- Material, cost, carbon and hardware values are indicative placeholders. Replace them with supplier data sheets and EPDs.
- Validate with impedance-tube samples (ISO 10534-2), then a reverberation-room test of a full module.
- Variable-density felt panels are covered by existing patents (for example Caimi Brevetti's Snowsound technology). Run a freedom-to-operate check before production.
