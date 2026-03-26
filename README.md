# Dravya

**Dravya** (द्रव्य — Sanskrit for "substance, matter") — material science engine for the [AGNOS](https://github.com/MacCracken/agnosticos) ecosystem.

## Features

- **Materials** — 9 presets (steel, aluminum, copper, titanium, glass, rubber, concrete, oak, carbon fiber) with E, v, sigma_y, sigma_uts, rho, alpha; derived elastic properties and thermal strain
- **Stress** — symmetric tensor (Voigt), von Mises, Tresca, principal stresses, invariants (I1/I2/I3/J2), deviatoric, octahedral shear, arithmetic ops
- **Strain** — engineering/true strain, volumetric, deviatoric, effective (equivalent) strain, arithmetic ops
- **Elasticity** — Hooke's law, bulk/shear modulus, Lame parameters, reverse conversions (G,v->E; K,G->E,v), plane stress/strain moduli, P-wave modulus
- **Yield** — von Mises and Tresca yield checks and safety factors, Drucker-Prager (concrete/rock/soil)
- **Beams** — deflection (point load, UDL, fixed-fixed), bending/shear/torsion stress, Euler buckling, moment of inertia (rect, circle, hollow), polar moment, section modulus
- **Fatigue** — Basquin's law (cycle and reversal forms), Miner's rule, Goodman/Gerber/Soderberg mean-stress corrections, endurance limit estimation, stress ratio utilities

## Quick Start

```rust
use dravya::{material::Material, stress::StressTensor, yield_criteria};

let steel = Material::steel();
let s = StressTensor::uniaxial(200e6);
let sf = yield_criteria::safety_factor(&s, steel.yield_strength);
println!("Safety factor: {sf:.2}"); // 1.25
```

## License

GPL-3.0
