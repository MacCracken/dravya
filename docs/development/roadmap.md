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
- ~~Stress intensity factor (KI)~~ done (center/edge/penny/finite-width/hole geometries)
- ~~Fracture toughness (KIc)~~ done (check, critical crack length, energy release)
- ~~Crack growth (Paris law)~~ done (rate + life prediction)
- Stress intensity factors KII, KIII (mode II/III)
- J-integral

### Composites
- Classical laminate theory
- Ply failure criteria (Tsai-Hill, Tsai-Wu)
- Laminate stiffness matrices

### Thermal
- ushma coupling for transient thermal analysis

### Fatigue (advanced)
- ~~Coffin-Manson low-cycle fatigue~~ done
- Rainflow cycle counting
- SN curve interpolation from tabulated data
- ~~Marin endurance limit modification factors~~ done
- Neuber's rule for notch effects

### Constitutive (advanced)
- ~~Generalized 3D Hooke's law (stiffness/compliance matrices)~~ done
- ~~Stress-strain tensor conversion via constitutive law~~ done
- ~~Elastic-perfectly-plastic model~~ done
- ~~Ramberg-Osgood nonlinear stress-strain~~ done
- ~~Bilinear hardening (elastic-plastic with tangent modulus)~~ done
- Temperature-dependent material properties

## v1.0.0 Criteria
- API frozen, zero unwrap/panic, 90%+ coverage, benchmark golden numbers
