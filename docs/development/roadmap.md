# Dravya Roadmap

## Status
**v1.0.0** — API stable. 13 modules, 278 tests, 26 benchmarks, zero panics. Hardened via P(-1) with 9.6x principal stress speedup, Johnson-Cook, Neo-Hookean, Hashin, progressive failure, orthotropic 3D.

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
- ~~Bilinear hardening model~~ done
- ~~Power-law hardening (Ramberg-Osgood)~~ done
- ~~Isotropic/kinematic hardening~~ done (isotropic, kinematic/Prager, combined)

### Fracture Mechanics
- ~~Stress intensity factor (KI)~~ done (center/edge/penny/finite-width/hole geometries)
- ~~Fracture toughness (KIc)~~ done (check, critical crack length, energy release)
- ~~Crack growth (Paris law)~~ done (rate + life prediction)
- ~~Stress intensity factors KII, KIII (mode II/III)~~ done
- ~~J-integral~~ done (from SIFs, mode I, mixed-mode, K-from-J roundtrip)

### Composites
- ~~Classical laminate theory~~ done (ABD matrix from ply layup)
- ~~Ply failure criteria (Tsai-Hill, Tsai-Wu)~~ done (+ max stress)
- ~~Laminate stiffness matrices~~ done (Q, Q-bar, ABD)

### Thermal
- ~~ushma coupling for transient thermal analysis~~ done (feature-gated `thermal` module)

### Fatigue (advanced)
- ~~Coffin-Manson low-cycle fatigue~~ done
- ~~Rainflow cycle counting~~ done
- ~~SN curve interpolation from tabulated data~~ done (log-log interpolation)
- ~~Marin endurance limit modification factors~~ done
- ~~Neuber's rule for notch effects~~ done (product + Ramberg-Osgood solver)

### Constitutive (advanced)
- ~~Generalized 3D Hooke's law (stiffness/compliance matrices)~~ done
- ~~Stress-strain tensor conversion via constitutive law~~ done
- ~~Elastic-perfectly-plastic model~~ done
- ~~Ramberg-Osgood nonlinear stress-strain~~ done
- ~~Bilinear hardening (elastic-plastic with tangent modulus)~~ done
- ~~Temperature-dependent material properties~~ done (TempDependentMaterial with linear interpolation)

## v1.0.0 Criteria — MET
- API frozen
- Zero unwrap/panic in library code
- 243 tests (233 unit + 10 integration)
- 26 golden benchmarks captured
- All cleanliness gates pass (fmt, clippy, audit, deny, doc)
