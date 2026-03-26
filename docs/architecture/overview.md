# Dravya Architecture

```
dravya
├── material.rs       — 9 engineering material presets, derived elastic properties, thermal strain
├── stress.rs         — Symmetric stress tensor, von Mises, principals, invariants (I1/I2/I3/J2), deviatoric, arithmetic
├── strain.rs         — Engineering/true strain, volumetric, deviatoric, effective strain, arithmetic
├── elastic.rs        — Hooke's law, bulk/shear modulus, Lame, reverse conversions, plane stress/strain
├── yield_criteria.rs — von Mises/Tresca yield, safety factors, Drucker-Prager
├── beam.rs           — Deflection (point/UDL/fixed), bending/shear/torsion stress, Euler buckling, section properties
├── fatigue.rs        — Basquin, Miner's rule, endurance limit, Goodman/Gerber/Soderberg corrections
├── error.rs          — DravyaError enum, Result alias
└── logging.rs        — Optional tracing subscriber init (feature: logging)
```

## Dependencies

- **hisab** — eigenvalue decomposition (principal stresses), EPSILON_F64 constant
- **serde** — serialization for Material, StressTensor, StrainTensor
- **thiserror** — error derive
- **tracing** — structured logging events

## Consumers

- **impetus** — collision materials, wave speeds
- **soorat** — PBR from physical properties
- **kiran/joshua** — destructible environments (Drucker-Prager for brittle materials)
- **ushma** — thermal expansion coupling
