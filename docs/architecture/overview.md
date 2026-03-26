# Dravya Architecture

```
dravya
├── material.rs       — 9 engineering material presets
├── stress.rs         — Symmetric stress tensor, von Mises, principals
├── strain.rs         — Engineering/true strain, volumetric
├── elastic.rs        — Hooke's law, bulk/shear modulus, Lamé
├── yield_criteria.rs — von Mises/Tresca yield, safety factor
├── beam.rs           — Deflection, bending stress, moment of inertia
└── fatigue.rs        — Basquin, Miner's rule, endurance limit
```

Consumers: impetus, soorat, kiran/joshua, ushma
