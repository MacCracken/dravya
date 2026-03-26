# Changelog

## [Unreleased]

### Added
- **material** — `ultimate_tensile_strength` field on `Material`, convenience methods `shear_modulus()`, `bulk_modulus()`, `lame_lambda()`, `thermal_strain()`, `thermal_stress()`; `Default` and `Display` trait implementations
- **stress** — deviatoric stress tensor, stress invariants I1/I2/I3, J2 deviatoric invariant, octahedral shear stress, hydrostatic/pure_shear constructors, `Add`/`Sub`/`Mul` arithmetic, `scale()`, `Default`/`Display`/`PartialEq` traits
- **strain** — deviatoric strain, effective (von Mises equivalent) strain, `Add`/`Sub`/`Mul` arithmetic, `scale()`, `Default`/`Display`/`PartialEq` traits; documented Voigt shear convention
- **elastic** — reverse conversions: `youngs_from_shear`, `youngs_from_bulk_shear`, `poisson_from_bulk_shear`; plane stress/strain moduli; P-wave modulus; `shear_modulus` division-by-zero guard
- **beam** — distributed-load deflections (cantilever UDL, simply-supported UDL), fixed-fixed beam deflection, transverse shear stress (`VQ/Ib`), torsional stress and angle of twist, Euler critical buckling load, hollow circular/rectangular moment of inertia, polar moment of inertia, section modulus for circular sections
- **yield_criteria** — Tresca safety factor, Drucker-Prager yield check with Mohr-Coulomb parameter conversion (for concrete, rock, soil)
- **fatigue** — reversal-based Basquin form (`basquin_cycles_reversals`), Goodman/Gerber/Soderberg mean-stress corrections, stress ratio and amplitude/mean decomposition helpers; endurance limit 700 MPa cap
- **error** — `DivisionByZero` and `InvalidParameter` variants; `Clone` derive
- **benchmarks** — 9 new benchmarks: max_shear, j2, deviatoric, euler_buckling, safety_factor, basquin_cycles, miners_rule_100, goodman_correction, effective_strain

### Changed
- **stress** — principal stress solver replaced with hisab's Jacobi eigenvalue decomposition (`eigen_symmetric`), improving numerical accuracy; von Mises refactored through J2 invariant
- **all modules** — epsilon guards now use `hisab::EPSILON_F64` (1e-12) instead of `f64::EPSILON`
- **lib.rs** — expanded re-exports: 25+ items now accessible from crate root

### Fixed
- **material** — copper yield strength corrected: 210 MPa -> 62 MPa (annealed C11000); concrete yield_strength corrected: 3 MPa -> 30 MPa (compressive strength); titanium density: 4507 -> 4430 kg/m^3; carbon fiber UTS: 3500 -> 1800 MPa (composite, not bare fiber); rubber CTE: 200e-6 -> 120e-6; aluminum CTE: 23e-6 -> 23.6e-6; added alloy/condition documentation to all preset names
- **fatigue** — `basquin_cycles` now guards against negative stress amplitude (was NaN); `endurance_limit_estimate` caps at 700 MPa per standard practice
- **elastic** — `shear_modulus` now guards against v = -1.0 (division by zero)
- **logging** — `init()` uses `try_init()` to avoid panic on double-call

## [0.1.0] - 2026-03-25

### Added
- **material** — 9 engineering material presets with full mechanical properties (steel, aluminum, copper, titanium, glass, rubber, concrete, oak, carbon fiber)
- **stress** — symmetric tensor (Voigt notation), von Mises equivalent stress, Tresca max shear, principal stresses, hydrostatic stress
- **strain** — engineering strain, true (logarithmic) strain, volumetric strain tensor
- **elastic** — Hooke's law, bulk/shear modulus, Lame parameters, inverse Hooke
- **yield_criteria** — von Mises/Tresca yield checks, safety factor
- **beam** — cantilever/simply-supported deflection, bending stress, moment of inertia (rectangular, circular), section modulus
- **fatigue** — Basquin's law (cycles to failure), Miner's rule (cumulative damage), endurance limit estimation
