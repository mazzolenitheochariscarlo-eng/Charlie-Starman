# Mono-material PET exhibition desk screen: generated design

Bred over 30 generations; best fitness 0.527.

## Winning design

| Parameter | Value |
|---|---|
| Cell pattern | phyllotaxis, 230 cells/m2 (180 per face) |
| Relief field | radial, exponent 0.45 |
| Felt | Recycled PET felt, light, sheet 25.0 mm |
| Forming | pressed to 68% of sheet thickness, cavities up to 28.5 mm |
| Core | Rigid recycled PET board 9 mm |
| Height above desk | 650 mm (top at 1390 mm) |
| Total thickness | 116 mm |
| Joint | Slot and key, CNC-cut |
| Support | Weighted foot |

## Acoustic performance

- Weighted absorption **alpha_w = 0.45(H)**, class **D**, NRC 0.50
- Equivalent absorption area per 1200 mm module (both faces): 0.78 m2 Sabine
- Speech at the facing desk: 51.4 dB(A) without screen, 45.0 dB(A) with screen (-6.4 dB)
- Speech intelligibility index: 0.43 -> 0.17 (privacy gain 61%)
- Open area from joints and the slit above the desk: 0.9%

| Octave (Hz) | 125 | 250 | 500 | 1000 | 2000 | 4000 |
|---|---|---|---|---|---|---|
| alpha_p | 0.05 | 0.20 | 0.40 | 0.65 | 0.85 | 0.90 |

## Exhibition system

| Item | Value |
|---|---|
| Modules | 3 x 700 mm, 8 x 800 mm, 6 x 1400 mm, 3 x 1600 mm |
| Module mass | 700 mm: 2.15 kg, 800 mm: 2.45 kg, 1400 mm: 4.3 kg, 1600 mm: 4.91 kg |
| Joints / supports | 17 / 34 |
| Screen area | 14.1 m2 |
| Total mass incl. fittings | 108 kg |
| Materials and fittings cost | 599 EUR |
| Embodied carbon | 185 kgCO2e (12 per event over 15 events) |
| Flat-pack volume | 1.88 m3 (~2 pallet(s)) |
| Assembly effort | 0.4 person-hours |
| Mono-material (felt and core recycle together) | yes |

## Shortlist (distinct alternatives)

| # | Pattern | Felt | Core | Height | alpha_w | Speech -dB(A) | Heaviest module | Fitness |
|---|---|---|---|---|---|---|---|---|
| 1 | phyllotaxis | pet_felt_light | pet_board_9mm | 650 | 0.45 | 6.4 | 4.9 kg | 0.527 |
| 2 | voronoi | pet_felt_light | pet_board_9mm | 650 | 0.45 | 6.4 | 4.9 kg | 0.527 |
| 3 | hex | pet_felt_light | pet_board_9mm | 650 | 0.45 | 6.4 | 4.9 kg | 0.526 |

## Files

- `module_1400x650.svg`
- `module_1400x650.obj`

- `bill_of_materials.csv`
- `design.json`
- `absorption.png`
- `speech_privacy.png`
- `evolution.png`
- `face_pattern.png`

## Limits of the models

- Absorption is predicted for an infinite surface (Miki model, Paris integration). Measured values on real panels (ISO 354) are usually higher because of edge effects.
- Speech privacy counts the direct, diffracted and transmitted paths only. Ceiling and wall reflections are not modelled.
- Material, cost and carbon values are indicative. Replace them with supplier data and EPDs.
- Confirm with impedance-tube samples (ISO 10534-2) and then a reverberation-room test.
