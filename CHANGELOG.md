# Changelog

## [1.0.0] - 2026-03-26

### Performance
- **stress** — principal stresses: 860 ns -> 90 ns (9.6x faster) via closed-form Lode angle solver replacing iterative Jacobi eigenvalue decomposition with heap allocations
- **stress** — max_shear: 868 ns -> 90 ns (9.6x faster) — same improvement, propagated
- **fatigue** — rainflow: 1900 ns -> 1533 ns (19% faster) via O(1) pop/push replacing O(n) Vec::remove
- **material** — `Material::name` changed from `String` to `Cow<'static, str>`, eliminating heap allocation for presets and temperature interpolation

### Added
- **material** — `ultimate_tensile_strength` field on `Material`, convenience methods `shear_modulus()`, `bulk_modulus()`, `lame_lambda()`, `thermal_strain()`, `thermal_stress()`; `Default` and `Display` trait implementations
- **stress** — deviatoric stress tensor, stress invariants I1/I2/I3, J2 deviatoric invariant, octahedral shear stress, hydrostatic/pure_shear constructors, `Add`/`Sub`/`Mul` arithmetic, `scale()`, `Default`/`Display`/`PartialEq` traits
- **strain** — deviatoric strain, effective (von Mises equivalent) strain, `Add`/`Sub`/`Mul` arithmetic, `scale()`, `Default`/`Display`/`PartialEq` traits; documented Voigt shear convention
- **elastic** — reverse conversions: `youngs_from_shear`, `youngs_from_bulk_shear`, `poisson_from_bulk_shear`; plane stress/strain moduli; P-wave modulus; `shear_modulus` division-by-zero guard
- **beam** — distributed-load deflections (cantilever UDL, simply-supported UDL), fixed-fixed beam deflection, transverse shear stress (`VQ/Ib`), torsional stress and angle of twist, Euler critical buckling load, hollow circular/rectangular moment of inertia, polar moment of inertia, section modulus for circular sections
- **yield_criteria** — Tresca safety factor, Drucker-Prager yield check with Mohr-Coulomb parameter conversion (for concrete, rock, soil)
- **fatigue** — reversal-based Basquin form (`basquin_cycles_reversals`), Goodman/Gerber/Soderberg mean-stress corrections, stress ratio and amplitude/mean decomposition helpers; endurance limit 700 MPa cap
- **error** — `DivisionByZero` and `InvalidParameter` variants; `Clone` derive
- **constitutive** — new module: isotropic stiffness matrix C and compliance matrix S (6x6 Voigt), `stress_from_strain` / `strain_from_stress` (generalized 3D Hooke's law), elastic-perfectly-plastic uniaxial model, bilinear hardening with tangent modulus, Ramberg-Osgood nonlinear stress-strain (forward + Newton-Raphson inverse)
- **fracture** — new module: Mode I stress intensity factors (center crack infinite/finite, edge crack, penny-shaped, crack at hole/Bowie), fracture toughness check, critical crack length, fracture stress, KIc from energy release rate, Paris law crack growth rate and life prediction
- **fatigue** — Coffin-Manson low-cycle fatigue (strain-life equation, transition life), Marin endurance limit modification factors (surface finish, size, reliability, corrected endurance)
- **composite** — new module: `Lamina` type with orthotropic properties and presets (carbon/epoxy, glass/epoxy), reduced stiffness Q and transformed Q-bar matrices, Classical Laminate Theory ABD matrix and inverse, ply stress transformation, failure criteria (max stress, Tsai-Hill, Tsai-Wu with custom f* interaction parameter)
- **material** — 4 new presets: stainless steel 304, gray cast iron, brass C36000, HDPE (13 total)
- **thermal** — new module (feature-gated `thermal`): ushma coupling with `Material`-to-`ThermalMaterial` conversion, thermal strain/stress tensors, constrained thermal stress, mechanical strain extraction, 1D/2D thermal grid stress field computation, thermal yield detection, max safe temperature change
- **constitutive** — Johnson-Cook rate-dependent plasticity (with copper/4340/Ti-6Al-4V presets), Neo-Hookean hyperelasticity (strain energy, uniaxial stress/tangent), orthotropic 3D stiffness tensor (9 independent constants)
- **composite** — Hashin 2D failure criterion (4 modes: fiber/matrix tension/compression), maximum strain failure criterion, progressive failure analysis (ply discount with Hashin), ABD inverse, custom f* for Tsai-Wu
- **fatigue** — rainflow periodic variant for repeating signals, turning point extraction utility
- **benchmarks** — 12 new benchmarks: max_shear, j2, deviatoric, euler_buckling, safety_factor, basquin_cycles, miners_rule_100, goodman_correction, effective_strain, stress_from_strain, stiffness_matrix, elastic_perfectly_plastic

### Changed
- **stress** — principal stress solver: closed-form Lode angle method (zero-alloc, 9.6x faster than Jacobi); von Mises refactored through J2 invariant
- **material** — `Material::name` type changed from `String` to `Cow<'static, str>`
- **error** — `ComputationError(String)` replaced with `SolverNoConvergence { method, iterations }` (no heap allocation)
- **all modules** — epsilon guards now use `hisab::EPSILON_F64` (1e-12) instead of `f64::EPSILON`; silent 0.0 fallbacks now emit `tracing::warn!`
- **lib.rs** — expanded re-exports: 100+ items now accessible from crate root

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
