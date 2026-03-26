# Dravya Roadmap

## Status
**v0.1.0** — Initial scaffold with real material science, hardened via P(-1).

## Completed (P(-1) Hardening)
- Material property verification against engineering handbooks
- Stress invariants (I1, I2, I3, J2), deviatoric stress, octahedral shear
- Strain deviatoric and effective (von Mises equivalent) strain
- Elastic constant reverse conversions and plane stress/strain moduli
- Drucker-Prager yield criterion for pressure-dependent materials
- Beam distributed loads, torsion, Euler buckling, hollow sections
- Goodman/Gerber/Soderberg mean-stress corrections
- Tensor arithmetic (Add, Sub, Mul)
- Thermal strain/stress on Material
- hisab integration (eigen_symmetric for principal stresses, EPSILON_F64)

## Future Features (demand-gated)

### Plasticity
- Bilinear hardening model
- Power-law hardening (Ramberg-Osgood)
- Isotropic/kinematic hardening

### Fracture Mechanics
- Stress intensity factor (KI, KII, KIII)
- Fracture toughness (KIc)
- Crack growth (Paris law)
- J-integral

### Composites
- Classical laminate theory
- Ply failure criteria (Tsai-Hill, Tsai-Wu)
- Laminate stiffness matrices

### Thermal
- ushma coupling for transient thermal analysis

### Fatigue (advanced)
- Coffin-Manson low-cycle fatigue
- Rainflow cycle counting
- SN curve interpolation from tabulated data
- Marin endurance limit modification factors
- Neuber's rule for notch effects

### Constitutive
- Generalized 3D Hooke's law (stiffness/compliance matrices)
- Stress-strain tensor conversion via constitutive law
- Elastic-perfectly-plastic model

## v1.0.0 Criteria
- API frozen, zero unwrap/panic, 90%+ coverage, benchmark golden numbers
