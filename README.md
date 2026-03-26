# Dravya

**Dravya** (द्रव्य — Sanskrit for "substance, matter") — material science engine for the [AGNOS](https://github.com/MacCracken/agnosticos) ecosystem.

## Features

- **Materials** — 9 presets (steel, aluminum, copper, titanium, glass, rubber, concrete, oak, carbon fiber) with E, ν, σ_y, ρ, α
- **Stress** — symmetric tensor, von Mises, Tresca max shear, principal stresses, hydrostatic
- **Strain** — engineering strain, true strain, volumetric strain
- **Elasticity** — Hooke's law, bulk/shear modulus, Lamé parameters
- **Yield** — von Mises and Tresca yield checks, safety factor
- **Beams** — cantilever/simply-supported deflection, bending stress, moment of inertia (rect, circle)
- **Fatigue** — Basquin's law, Miner's rule, endurance limit estimation

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
