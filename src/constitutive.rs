//! Constitutive models: stress-strain relationships beyond simple Hooke's law.
//!
//! Provides isotropic stiffness/compliance matrices (Voigt notation),
//! tensor conversion, material models (EPP, bilinear, Ramberg-Osgood),
//! and hardening laws (isotropic, kinematic, combined).

use serde::{Deserialize, Serialize};

use crate::material::Material;
use crate::strain::StrainTensor;
use crate::stress::StressTensor;

/// 6x6 isotropic stiffness matrix C in Voigt notation.
///
/// Maps strain to stress: {σ} = \[C\]{ε}
///
/// Note: shear strain components must be **engineering** shear strain
/// (γ = 2ε_ij), which matches [`StrainTensor`]'s convention.
///
/// ```text
/// C = | λ+2μ   λ     λ    0   0   0  |
///     |  λ    λ+2μ    λ    0   0   0  |
///     |  λ      λ   λ+2μ   0   0   0  |
///     |  0      0     0    μ   0   0  |
///     |  0      0     0    0   μ   0  |
///     |  0      0     0    0   0   μ  |
/// ```
#[must_use]
pub fn stiffness_matrix(youngs_modulus: f64, poisson_ratio: f64) -> [[f64; 6]; 6] {
    let e = youngs_modulus;
    let v = poisson_ratio;
    let denom = (1.0 + v) * (1.0 - 2.0 * v);
    if denom.abs() < hisab::EPSILON_F64 {
        return [[0.0; 6]; 6];
    }
    let factor = e / denom;
    let c11 = factor * (1.0 - v);
    let c12 = factor * v;
    let c44 = factor * (1.0 - 2.0 * v) / 2.0; // = G = E / (2(1+v))

    let mut c = [[0.0; 6]; 6];
    c[0][0] = c11;
    c[1][1] = c11;
    c[2][2] = c11;
    c[0][1] = c12;
    c[0][2] = c12;
    c[1][0] = c12;
    c[1][2] = c12;
    c[2][0] = c12;
    c[2][1] = c12;
    c[3][3] = c44;
    c[4][4] = c44;
    c[5][5] = c44;
    c
}

/// 6x6 isotropic compliance matrix S in Voigt notation.
///
/// Maps stress to strain: {ε} = \[S\]{σ}
///
/// ```text
/// S = | 1/E   -v/E  -v/E   0      0      0     |
///     | -v/E   1/E  -v/E   0      0      0     |
///     | -v/E  -v/E   1/E   0      0      0     |
///     |  0      0     0    1/G    0      0     |
///     |  0      0     0     0    1/G    0     |
///     |  0      0     0     0     0    1/G   |
/// ```
#[must_use]
pub fn compliance_matrix(youngs_modulus: f64, poisson_ratio: f64) -> [[f64; 6]; 6] {
    if youngs_modulus.abs() < hisab::EPSILON_F64 {
        return [[0.0; 6]; 6];
    }
    let e = youngs_modulus;
    let v = poisson_ratio;
    let s11 = 1.0 / e;
    let s12 = -v / e;
    let g = e / (2.0 * (1.0 + v));
    let s44 = if g.abs() < hisab::EPSILON_F64 {
        0.0
    } else {
        1.0 / g
    };

    let mut s = [[0.0; 6]; 6];
    s[0][0] = s11;
    s[1][1] = s11;
    s[2][2] = s11;
    s[0][1] = s12;
    s[0][2] = s12;
    s[1][0] = s12;
    s[1][2] = s12;
    s[2][0] = s12;
    s[2][1] = s12;
    s[3][3] = s44;
    s[4][4] = s44;
    s[5][5] = s44;
    s
}

/// Convert strain to stress using the isotropic stiffness matrix.
///
/// σ = C ε (generalized 3D Hooke's law).
#[must_use]
pub fn stress_from_strain_3d(material: &Material, strain: &StrainTensor) -> StressTensor {
    let c = stiffness_matrix(material.youngs_modulus, material.poisson_ratio);
    let e = &strain.components;
    let mut s = [0.0; 6];
    for i in 0..6 {
        for j in 0..6 {
            s[i] += c[i][j] * e[j];
        }
    }
    StressTensor::new(s[0], s[1], s[2], s[3], s[4], s[5])
}

/// Convert stress to strain using the isotropic compliance matrix.
///
/// ε = S σ (inverse generalized Hooke's law).
#[must_use]
pub fn strain_from_stress_3d(material: &Material, stress: &StressTensor) -> StrainTensor {
    let s_mat = compliance_matrix(material.youngs_modulus, material.poisson_ratio);
    let sig = &stress.components;
    let mut e = [0.0; 6];
    for i in 0..6 {
        for j in 0..6 {
            e[i] += s_mat[i][j] * sig[j];
        }
    }
    StrainTensor::new(e[0], e[1], e[2], e[3], e[4], e[5])
}

/// Elastic-perfectly-plastic uniaxial response.
///
/// Given total strain, returns stress clamped at yield:
/// - |ε| <= ε_y: σ = E × ε  (elastic)
/// - |ε| > ε_y:  σ = ±σ_y   (plastic, constant stress)
///
/// Returns (stress, is_yielded).
#[must_use]
#[inline]
pub fn elastic_perfectly_plastic(
    youngs_modulus: f64,
    yield_strength: f64,
    strain: f64,
) -> (f64, bool) {
    let stress = youngs_modulus * strain;
    if stress.abs() <= yield_strength {
        (stress, false)
    } else {
        (stress.signum() * yield_strength, true)
    }
}

/// Elastic-perfectly-plastic response using material properties.
///
/// Convenience wrapper around [`elastic_perfectly_plastic`].
#[must_use]
#[inline]
pub fn elastic_perfectly_plastic_material(material: &Material, strain: f64) -> (f64, bool) {
    elastic_perfectly_plastic(material.youngs_modulus, material.yield_strength, strain)
}

/// Bilinear hardening uniaxial response.
///
/// Two-slope model:
/// - |ε| <= ε_y: σ = E × ε (elastic)
/// - |ε| > ε_y: σ = ±σ_y + E_t × (ε - ±ε_y) (plastic with tangent modulus)
///
/// `tangent_modulus` is the post-yield slope (0 = perfectly plastic,
/// E = fully elastic, typical metals: 0.01E-0.1E).
///
/// Returns (stress, plastic_strain, is_yielded).
#[must_use]
pub fn bilinear_hardening(
    youngs_modulus: f64,
    yield_strength: f64,
    tangent_modulus: f64,
    strain: f64,
) -> (f64, f64, bool) {
    let yield_strain = if youngs_modulus.abs() < hisab::EPSILON_F64 {
        return (0.0, 0.0, false);
    } else {
        yield_strength / youngs_modulus
    };

    let elastic_stress = youngs_modulus * strain;
    if elastic_stress.abs() <= yield_strength {
        (elastic_stress, 0.0, false)
    } else {
        let sign = strain.signum();
        let excess_strain = strain.abs() - yield_strain;
        let stress = sign * (yield_strength + tangent_modulus * excess_strain);
        let plastic_strain = sign * excess_strain * (1.0 - tangent_modulus / youngs_modulus);
        (stress, plastic_strain, true)
    }
}

/// Ramberg-Osgood nonlinear stress-strain model.
///
/// Total strain from stress:
/// ε = σ/E + (σ/K)^n
///
/// K = strength coefficient (Pa), n = hardening exponent (typically 5-50
/// for metals; higher n = more abrupt yield).
///
/// This is the forward direction (stress -> strain). For the inverse
/// (strain -> stress), use [`ramberg_osgood_stress`] which uses
/// Newton-Raphson iteration.
#[must_use]
#[inline]
pub fn ramberg_osgood_strain(
    youngs_modulus: f64,
    strength_coeff: f64,
    hardening_exp: f64,
    stress: f64,
) -> f64 {
    if youngs_modulus.abs() < hisab::EPSILON_F64 || strength_coeff.abs() < hisab::EPSILON_F64 {
        return 0.0;
    }
    let elastic = stress / youngs_modulus;
    let plastic = (stress.abs() / strength_coeff).powf(hardening_exp) * stress.signum();
    elastic + plastic
}

/// Inverse Ramberg-Osgood: find stress from total strain.
///
/// Uses [`hisab::num::newton_raphson`] to solve ε = σ/E + (σ/K)^n for σ.
///
/// Returns the converged stress, or `Err` if the solver fails.
pub fn ramberg_osgood_stress(
    youngs_modulus: f64,
    strength_coeff: f64,
    hardening_exp: f64,
    strain: f64,
    tol: f64,
    max_iter: usize,
) -> crate::Result<f64> {
    if youngs_modulus.abs() < hisab::EPSILON_F64 {
        return Ok(0.0);
    }

    let f = |sigma: f64| {
        ramberg_osgood_strain(youngs_modulus, strength_coeff, hardening_exp, sigma) - strain
    };
    let df = |sigma: f64| {
        let d_elastic = 1.0 / youngs_modulus;
        let d_plastic = if sigma.abs() < hisab::EPSILON_F64 {
            0.0
        } else {
            (hardening_exp / strength_coeff)
                * (sigma.abs() / strength_coeff).powf(hardening_exp - 1.0)
        };
        d_elastic + d_plastic
    };

    let x0 = youngs_modulus * strain;
    hisab::num::newton_raphson(f, df, x0, tol, max_iter).map_err(|_| {
        crate::DravyaError::SolverNoConvergence {
            method: "ramberg_osgood_stress",
            iterations: max_iter,
        }
    })
}

// ---------------------------------------------------------------------------
// Hardening models
// ---------------------------------------------------------------------------

/// Isotropic hardening state.
///
/// The yield surface expands uniformly as plastic strain accumulates.
/// Current yield stress: σ_y(ε_p) = σ_y0 + H * ε_p_accumulated
///
/// H = isotropic hardening modulus (Pa).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct IsotropicHardening {
    /// Initial yield strength (Pa).
    pub initial_yield: f64,
    /// Hardening modulus H (Pa). σ_y grows linearly with accumulated plastic strain.
    pub hardening_modulus: f64,
    /// Accumulated equivalent plastic strain (scalar, always >= 0).
    pub accumulated_plastic_strain: f64,
}

impl IsotropicHardening {
    /// Create a new isotropic hardening state.
    #[must_use]
    pub fn new(initial_yield: f64, hardening_modulus: f64) -> Self {
        Self {
            initial_yield,
            hardening_modulus,
            accumulated_plastic_strain: 0.0,
        }
    }

    /// Current yield stress.
    #[must_use]
    #[inline]
    pub fn current_yield(&self) -> f64 {
        self.initial_yield + self.hardening_modulus * self.accumulated_plastic_strain
    }

    /// Apply a uniaxial strain increment. Returns (stress, updated state).
    ///
    /// Uses a radial-return-like 1D algorithm:
    /// 1. Trial stress = E * (total strain)
    /// 2. If |trial| <= current yield: elastic
    /// 3. Else: plastic correction
    #[must_use]
    pub fn apply_uniaxial(&self, youngs_modulus: f64, total_strain: f64) -> (f64, Self) {
        let trial = youngs_modulus * total_strain;
        let current_y = self.current_yield();

        if trial.abs() <= current_y {
            (trial, *self)
        } else {
            // Plastic: solve for stress and plastic strain increment
            // σ = sign(trial) * (σ_y0 + H * (ε_p_old + Δε_p))
            // σ = sign(trial) * σ_y0 + E * (ε - ε_p_old - Δε_p)  ... from elastic strain
            // Combining: Δε_p = (|trial| - current_y) / (E + H)
            let denom = youngs_modulus + self.hardening_modulus;
            let d_ep = if denom.abs() < hisab::EPSILON_F64 {
                0.0
            } else {
                (trial.abs() - current_y) / denom
            };
            let new_ep = self.accumulated_plastic_strain + d_ep;
            let stress = trial.signum() * (self.initial_yield + self.hardening_modulus * new_ep);
            let new_state = Self {
                accumulated_plastic_strain: new_ep,
                ..*self
            };
            (stress, new_state)
        }
    }
}

/// Kinematic hardening state (Prager linear kinematic).
///
/// The yield surface translates in stress space without expanding.
/// Yield check: |σ - α| <= σ_y0
///
/// α = back-stress, evolves as dα = C * dε_p (Prager rule).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct KinematicHardening {
    /// Initial yield strength (Pa, constant — surface does not expand).
    pub yield_strength: f64,
    /// Kinematic hardening modulus C (Pa).
    pub hardening_modulus: f64,
    /// Current back-stress (Pa).
    pub back_stress: f64,
}

impl KinematicHardening {
    /// Create a new kinematic hardening state.
    #[must_use]
    pub fn new(yield_strength: f64, hardening_modulus: f64) -> Self {
        Self {
            yield_strength,
            hardening_modulus,
            back_stress: 0.0,
        }
    }

    /// Apply a uniaxial strain increment. Returns (stress, updated state).
    #[must_use]
    pub fn apply_uniaxial(&self, youngs_modulus: f64, total_strain: f64) -> (f64, Self) {
        let trial = youngs_modulus * total_strain;
        let shifted = trial - self.back_stress;

        if shifted.abs() <= self.yield_strength {
            (trial, *self)
        } else {
            let denom = youngs_modulus + self.hardening_modulus;
            let d_ep = if denom.abs() < hisab::EPSILON_F64 {
                0.0
            } else {
                (shifted.abs() - self.yield_strength) / denom
            };
            let sign = shifted.signum();
            let new_back = self.back_stress + sign * self.hardening_modulus * d_ep;
            let stress = new_back + sign * self.yield_strength;
            let new_state = Self {
                back_stress: new_back,
                ..*self
            };
            (stress, new_state)
        }
    }
}

/// Combined isotropic + kinematic hardening state.
///
/// Yield check: |σ - α| <= σ_y(ε_p)
///
/// Both yield surface expansion (isotropic) and translation (kinematic)
/// occur simultaneously.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CombinedHardening {
    /// Initial yield strength (Pa).
    pub initial_yield: f64,
    /// Isotropic hardening modulus H_iso (Pa).
    pub iso_modulus: f64,
    /// Kinematic hardening modulus H_kin (Pa).
    pub kin_modulus: f64,
    /// Accumulated equivalent plastic strain.
    pub accumulated_plastic_strain: f64,
    /// Current back-stress.
    pub back_stress: f64,
}

impl CombinedHardening {
    /// Create a new combined hardening state.
    #[must_use]
    pub fn new(initial_yield: f64, iso_modulus: f64, kin_modulus: f64) -> Self {
        Self {
            initial_yield,
            iso_modulus,
            kin_modulus,
            accumulated_plastic_strain: 0.0,
            back_stress: 0.0,
        }
    }

    /// Current yield radius (isotropic contribution).
    #[must_use]
    #[inline]
    pub fn current_yield(&self) -> f64 {
        self.initial_yield + self.iso_modulus * self.accumulated_plastic_strain
    }

    /// Apply a uniaxial strain increment. Returns (stress, updated state).
    #[must_use]
    pub fn apply_uniaxial(&self, youngs_modulus: f64, total_strain: f64) -> (f64, Self) {
        let trial = youngs_modulus * total_strain;
        let shifted = trial - self.back_stress;
        let current_y = self.current_yield();

        if shifted.abs() <= current_y {
            (trial, *self)
        } else {
            let h_total = self.iso_modulus + self.kin_modulus;
            let denom = youngs_modulus + h_total;
            let d_ep = if denom.abs() < hisab::EPSILON_F64 {
                0.0
            } else {
                (shifted.abs() - current_y) / denom
            };
            let sign = shifted.signum();
            let new_ep = self.accumulated_plastic_strain + d_ep;
            let new_back = self.back_stress + sign * self.kin_modulus * d_ep;
            let new_yield = self.initial_yield + self.iso_modulus * new_ep;
            let stress = new_back + sign * new_yield;
            let new_state = Self {
                accumulated_plastic_strain: new_ep,
                back_stress: new_back,
                ..*self
            };
            (stress, new_state)
        }
    }
}

// ---------------------------------------------------------------------------
// Johnson-Cook rate-dependent plasticity
// ---------------------------------------------------------------------------

/// Johnson-Cook flow stress model for rate-dependent plasticity.
///
/// σ = (A + B * ε_p^n) * (1 + C * ln(ε̇*)) * (1 - T*^m)
///
/// where ε̇* = ε̇ / ε̇_0 (normalized strain rate),
/// T* = (T - T_room) / (T_melt - T_room) (homologous temperature).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct JohnsonCook {
    /// Initial yield stress A (Pa).
    pub a: f64,
    /// Hardening coefficient B (Pa).
    pub b: f64,
    /// Strain hardening exponent n.
    pub n: f64,
    /// Strain rate sensitivity C.
    pub c: f64,
    /// Thermal softening exponent m.
    pub m_thermal: f64,
    /// Reference strain rate (1/s).
    pub strain_rate_ref: f64,
    /// Room temperature (K).
    pub t_room: f64,
    /// Melting temperature (K).
    pub t_melt: f64,
}

impl JohnsonCook {
    /// Compute flow stress for given plastic strain, strain rate, and temperature.
    #[must_use]
    pub fn flow_stress(&self, plastic_strain: f64, strain_rate: f64, temperature: f64) -> f64 {
        // Strain hardening
        let hardening = self.a + self.b * plastic_strain.max(0.0).powf(self.n);

        // Strain rate effect
        let rate_ratio = if self.strain_rate_ref.abs() < hisab::EPSILON_F64 {
            1.0
        } else {
            (strain_rate / self.strain_rate_ref).max(1.0)
        };
        let rate_factor = 1.0 + self.c * rate_ratio.ln();

        // Thermal softening
        let t_star = if (self.t_melt - self.t_room).abs() < hisab::EPSILON_F64 {
            0.0
        } else {
            ((temperature - self.t_room) / (self.t_melt - self.t_room)).clamp(0.0, 1.0)
        };
        let thermal_factor = 1.0 - t_star.powf(self.m_thermal);

        hardening * rate_factor * thermal_factor
    }

    /// Typical parameters for OFHC copper.
    #[must_use]
    pub fn copper() -> Self {
        Self {
            a: 90e6,
            b: 292e6,
            n: 0.31,
            c: 0.025,
            m_thermal: 1.09,
            strain_rate_ref: 1.0,
            t_room: 293.0,
            t_melt: 1356.0,
        }
    }

    /// Typical parameters for 4340 steel.
    #[must_use]
    pub fn steel_4340() -> Self {
        Self {
            a: 792e6,
            b: 510e6,
            n: 0.26,
            c: 0.014,
            m_thermal: 1.03,
            strain_rate_ref: 1.0,
            t_room: 293.0,
            t_melt: 1793.0,
        }
    }

    /// Typical parameters for Ti-6Al-4V.
    #[must_use]
    pub fn titanium_ti6al4v() -> Self {
        Self {
            a: 1098e6,
            b: 1092e6,
            n: 0.93,
            c: 0.014,
            m_thermal: 1.1,
            strain_rate_ref: 1.0,
            t_room: 293.0,
            t_melt: 1878.0,
        }
    }
}

// ---------------------------------------------------------------------------
// Neo-Hookean hyperelasticity
// ---------------------------------------------------------------------------

/// Neo-Hookean hyperelastic model.
///
/// Strain energy: W = C1 * (I1_bar - 3) + 1/D1 * (J - 1)^2
///
/// where C1 = mu/2, D1 = 2/kappa (mu = shear modulus, kappa = bulk modulus),
/// I1_bar = J^(-2/3) * I1 (isochoric first invariant),
/// J = det(F) (volume ratio).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct NeoHookean {
    /// Material parameter C1 = mu/2 (Pa).
    pub c1: f64,
    /// Material parameter D1 = 2/kappa (1/Pa). 0 = incompressible.
    pub d1: f64,
}

impl NeoHookean {
    /// Create from shear modulus and bulk modulus.
    #[must_use]
    pub fn from_moduli(shear_modulus: f64, bulk_modulus: f64) -> Self {
        Self {
            c1: shear_modulus / 2.0,
            d1: if bulk_modulus.abs() < hisab::EPSILON_F64 {
                0.0
            } else {
                2.0 / bulk_modulus
            },
        }
    }

    /// Create from Young's modulus and Poisson's ratio.
    #[must_use]
    pub fn from_elastic(youngs_modulus: f64, poisson_ratio: f64) -> Self {
        let g = crate::elastic::shear_modulus(youngs_modulus, poisson_ratio);
        let k = crate::elastic::bulk_modulus(youngs_modulus, poisson_ratio);
        Self::from_moduli(g, k)
    }

    /// Strain energy density for a given deformation state.
    ///
    /// `i1` = first invariant of left Cauchy-Green tensor (tr(B)),
    /// `j` = determinant of deformation gradient (volume ratio).
    #[must_use]
    #[inline]
    pub fn strain_energy(&self, i1: f64, j: f64) -> f64 {
        let i1_bar = j.powf(-2.0 / 3.0) * i1;
        let volumetric = if self.d1.abs() < hisab::EPSILON_F64 {
            0.0
        } else {
            (j - 1.0).powi(2) / self.d1
        };
        self.c1 * (i1_bar - 3.0) + volumetric
    }

    /// Cauchy stress for uniaxial stretch ratio λ (incompressible).
    ///
    /// σ = 2 * C1 * (λ^2 - 1/λ)
    #[must_use]
    #[inline]
    pub fn uniaxial_stress(&self, stretch: f64) -> f64 {
        if stretch.abs() < hisab::EPSILON_F64 {
            return 0.0;
        }
        2.0 * self.c1 * (stretch * stretch - 1.0 / stretch)
    }

    /// Tangent modulus for uniaxial stretch (incompressible).
    ///
    /// dσ/dλ = 2 * C1 * (2λ + 1/λ^2)
    #[must_use]
    #[inline]
    pub fn uniaxial_tangent(&self, stretch: f64) -> f64 {
        if stretch.abs() < hisab::EPSILON_F64 {
            return 0.0;
        }
        2.0 * self.c1 * (2.0 * stretch + 1.0 / (stretch * stretch))
    }
}

// ---------------------------------------------------------------------------
// Orthotropic 3D stiffness
// ---------------------------------------------------------------------------

/// 3D orthotropic elastic material properties.
///
/// 9 independent constants: E1, E2, E3, G12, G23, G13, v12, v23, v13.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Orthotropic3D {
    pub e1: f64,
    pub e2: f64,
    pub e3: f64,
    pub g12: f64,
    pub g23: f64,
    pub g13: f64,
    pub nu12: f64,
    pub nu23: f64,
    pub nu13: f64,
}

impl Orthotropic3D {
    /// Compute the 6x6 stiffness matrix in Voigt notation.
    ///
    /// Uses reciprocity: v_ij / E_i = v_ji / E_j.
    #[must_use]
    pub fn stiffness_matrix(&self) -> [[f64; 6]; 6] {
        let nu21 = self.nu12 * self.e2 / self.e1;
        let nu31 = self.nu13 * self.e3 / self.e1;
        let nu32 = self.nu23 * self.e3 / self.e2;

        let delta = 1.0
            - self.nu12 * nu21
            - self.nu23 * nu32
            - self.nu13 * nu31
            - 2.0 * self.nu12 * nu32 * self.nu13;

        if delta.abs() < hisab::EPSILON_F64 {
            return [[0.0; 6]; 6];
        }

        let c11 = self.e1 * (1.0 - self.nu23 * nu32) / delta;
        let c22 = self.e2 * (1.0 - self.nu13 * nu31) / delta;
        let c33 = self.e3 * (1.0 - self.nu12 * nu21) / delta;
        let c12 = self.e1 * (nu21 + nu31 * self.nu23) / delta;
        let c13 = self.e1 * (nu31 + nu21 * nu32) / delta;
        let c23 = self.e2 * (nu32 + self.nu12 * nu31) / delta;

        let mut c = [[0.0; 6]; 6];
        c[0][0] = c11;
        c[1][1] = c22;
        c[2][2] = c33;
        c[0][1] = c12;
        c[1][0] = c12;
        c[0][2] = c13;
        c[2][0] = c13;
        c[1][2] = c23;
        c[2][1] = c23;
        c[3][3] = self.g12;
        c[4][4] = self.g23;
        c[5][5] = self.g13;
        c
    }

    /// Create from isotropic properties (all directions equal).
    #[must_use]
    pub fn from_isotropic(youngs_modulus: f64, poisson_ratio: f64) -> Self {
        let g = crate::elastic::shear_modulus(youngs_modulus, poisson_ratio);
        Self {
            e1: youngs_modulus,
            e2: youngs_modulus,
            e3: youngs_modulus,
            g12: g,
            g23: g,
            g13: g,
            nu12: poisson_ratio,
            nu23: poisson_ratio,
            nu13: poisson_ratio,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stiffness_matrix_symmetry() {
        let c = stiffness_matrix(200e9, 0.30);
        for (i, row) in c.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    (val - c[j][i]).abs() < 1.0,
                    "C should be symmetric at [{i}][{j}]"
                );
            }
        }
    }

    #[test]
    fn compliance_matrix_symmetry() {
        let s = compliance_matrix(200e9, 0.30);
        for (i, row) in s.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    (val - s[j][i]).abs() < hisab::EPSILON_F64,
                    "S should be symmetric at [{i}][{j}]"
                );
            }
        }
    }

    #[test]
    fn stiffness_compliance_inverse() {
        // C * S should be identity (approximately)
        let c = stiffness_matrix(200e9, 0.30);
        let s = compliance_matrix(200e9, 0.30);
        for (i, c_row) in c.iter().enumerate() {
            for (j, _) in s.iter().enumerate() {
                let product: f64 = c_row
                    .iter()
                    .enumerate()
                    .map(|(k, &c_ik)| c_ik * s[k][j])
                    .sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (product - expected).abs() < 1e-6,
                    "C*S should be identity at [{i}][{j}], got {product}"
                );
            }
        }
    }

    #[test]
    fn uniaxial_stress_from_strain() {
        let steel = Material::steel();
        // Uniaxial strain εxx = 0.001, εyy = εzz = -v * εxx
        let v = steel.poisson_ratio;
        let eps = StrainTensor::new(0.001, -v * 0.001, -v * 0.001, 0.0, 0.0, 0.0);
        let sigma = stress_from_strain_3d(&steel, &eps);
        // Should produce uniaxial stress σxx = E * εxx = 200 MPa
        assert!(
            (sigma.components[0] - 200e6).abs() < 1e3,
            "σxx should be 200 MPa, got {}",
            sigma.components[0] / 1e6
        );
        assert!(
            sigma.components[1].abs() < 1e3,
            "σyy should be ~0, got {}",
            sigma.components[1]
        );
        assert!(
            sigma.components[2].abs() < 1e3,
            "σzz should be ~0, got {}",
            sigma.components[2]
        );
    }

    #[test]
    fn uniaxial_strain_from_stress() {
        let steel = Material::steel();
        let sigma = StressTensor::uniaxial(200e6);
        let eps = strain_from_stress_3d(&steel, &sigma);
        // εxx = σ/E = 0.001
        assert!(
            (eps.components[0] - 0.001).abs() < 1e-9,
            "εxx should be 0.001, got {}",
            eps.components[0]
        );
        // εyy = εzz = -v * σ/E = -0.0003
        let v = steel.poisson_ratio;
        let expected_lat = -v * 0.001;
        assert!(
            (eps.components[1] - expected_lat).abs() < 1e-9,
            "εyy should be {expected_lat}, got {}",
            eps.components[1]
        );
    }

    #[test]
    fn stress_strain_roundtrip() {
        let steel = Material::steel();
        let original = StressTensor::new(100e6, 50e6, 30e6, 20e6, 10e6, 5e6);
        let eps = strain_from_stress_3d(&steel, &original);
        let recovered = stress_from_strain_3d(&steel, &eps);
        for i in 0..6 {
            assert!(
                (original.components[i] - recovered.components[i]).abs() < 1.0,
                "component {i} mismatch: {} vs {}",
                original.components[i],
                recovered.components[i]
            );
        }
    }

    #[test]
    fn hydrostatic_stress_from_strain() {
        let steel = Material::steel();
        // Pure volumetric strain
        let eps = StrainTensor::new(0.001, 0.001, 0.001, 0.0, 0.0, 0.0);
        let sigma = stress_from_strain_3d(&steel, &eps);
        // All normal stresses should be equal (hydrostatic)
        assert!(
            (sigma.components[0] - sigma.components[1]).abs() < 1.0,
            "hydrostatic: σxx should equal σyy"
        );
        assert!(
            (sigma.components[1] - sigma.components[2]).abs() < 1.0,
            "hydrostatic: σyy should equal σzz"
        );
        // σ = K * 3ε_vol / 3 ... actually σ_ii = (λ + 2μ)ε + λε + λε = 3Kε
        let k = steel.bulk_modulus();
        let expected = 3.0 * k * 0.001;
        assert!(
            (sigma.components[0] - expected).abs() < 1e3,
            "hydrostatic stress should be 3Kε, got {}",
            sigma.components[0]
        );
    }

    #[test]
    fn pure_shear_stress_from_strain() {
        let steel = Material::steel();
        let g = steel.shear_modulus();
        // Engineering shear strain γxy = 0.001
        let eps = StrainTensor::new(0.0, 0.0, 0.0, 0.001, 0.0, 0.0);
        let sigma = stress_from_strain_3d(&steel, &eps);
        // τxy = G * γxy
        assert!(
            (sigma.components[3] - g * 0.001).abs() < 1.0,
            "τxy should be G*γ, got {}",
            sigma.components[3]
        );
    }

    #[test]
    fn epp_elastic_region() {
        let (stress, yielded) = elastic_perfectly_plastic(200e9, 250e6, 0.001);
        assert!(!yielded);
        assert!(
            (stress - 200e6).abs() < 1.0,
            "elastic region: σ = Eε = 200 MPa"
        );
    }

    #[test]
    fn epp_plastic_region() {
        let (stress, yielded) = elastic_perfectly_plastic(200e9, 250e6, 0.01);
        assert!(yielded);
        assert!(
            (stress - 250e6).abs() < 1.0,
            "plastic region: σ = σ_y = 250 MPa"
        );
    }

    #[test]
    fn epp_compression() {
        let (stress, yielded) = elastic_perfectly_plastic(200e9, 250e6, -0.01);
        assert!(yielded);
        assert!((stress - (-250e6)).abs() < 1.0, "compression: σ = -σ_y");
    }

    #[test]
    fn epp_at_yield() {
        let yield_strain = 250e6 / 200e9; // 0.00125
        let (stress, yielded) = elastic_perfectly_plastic(200e9, 250e6, yield_strain);
        assert!(!yielded, "exactly at yield should be elastic");
        assert!((stress - 250e6).abs() < 1.0);
    }

    #[test]
    fn epp_material_wrapper() {
        let steel = Material::steel();
        let (stress, yielded) = elastic_perfectly_plastic_material(&steel, 0.01);
        assert!(yielded);
        assert!((stress - steel.yield_strength).abs() < 1.0);
    }

    #[test]
    fn incompressible_guard() {
        let c = stiffness_matrix(200e9, 0.5);
        // Should return zeros (degenerate case)
        assert_eq!(c[0][0], 0.0);
    }

    #[test]
    fn zero_modulus_compliance_guard() {
        let s = compliance_matrix(0.0, 0.3);
        assert_eq!(s[0][0], 0.0);
    }

    #[test]
    fn stiffness_diagonal_positive() {
        let c = stiffness_matrix(200e9, 0.30);
        for (i, row) in c.iter().enumerate() {
            assert!(row[i] > 0.0, "C[{i}][{i}] should be positive");
        }
    }

    #[test]
    fn compliance_diagonal_positive() {
        let s = compliance_matrix(200e9, 0.30);
        for (i, row) in s.iter().enumerate() {
            assert!(row[i] > 0.0, "S[{i}][{i}] should be positive");
        }
    }

    // --- Bilinear hardening tests ---

    #[test]
    fn bilinear_elastic_region() {
        let (stress, plastic, yielded) = bilinear_hardening(200e9, 250e6, 20e9, 0.001);
        assert!(!yielded);
        assert!((stress - 200e6).abs() < 1.0);
        assert!(plastic.abs() < hisab::EPSILON_F64);
    }

    #[test]
    fn bilinear_plastic_region() {
        // E=200 GPa, σ_y=250 MPa, E_t=20 GPa, ε=0.005
        // ε_y = 250e6/200e9 = 0.00125
        // excess = 0.005 - 0.00125 = 0.00375
        // σ = 250e6 + 20e9 * 0.00375 = 250e6 + 75e6 = 325e6
        let (stress, plastic, yielded) = bilinear_hardening(200e9, 250e6, 20e9, 0.005);
        assert!(yielded);
        assert!(
            (stress - 325e6).abs() < 1e3,
            "stress should be 325 MPa, got {}",
            stress / 1e6
        );
        assert!(plastic > 0.0, "plastic strain should be positive");
    }

    #[test]
    fn bilinear_reduces_to_epp_at_zero_tangent() {
        let (stress_bl, _, yielded_bl) = bilinear_hardening(200e9, 250e6, 0.0, 0.01);
        let (stress_epp, yielded_epp) = elastic_perfectly_plastic(200e9, 250e6, 0.01);
        assert_eq!(yielded_bl, yielded_epp);
        assert!(
            (stress_bl - stress_epp).abs() < 1.0,
            "bilinear with E_t=0 should match EPP"
        );
    }

    #[test]
    fn bilinear_compression() {
        let (stress, plastic, yielded) = bilinear_hardening(200e9, 250e6, 20e9, -0.005);
        assert!(yielded);
        assert!(stress < 0.0, "stress should be negative");
        assert!(
            plastic < 0.0,
            "plastic strain should be negative in compression"
        );
    }

    // --- Ramberg-Osgood tests ---

    #[test]
    fn ramberg_osgood_elastic_dominates_low_stress() {
        // At low stress, plastic term is negligible
        let eps = ramberg_osgood_strain(200e9, 1000e6, 10.0, 50e6);
        let elastic = 50e6 / 200e9;
        assert!(
            (eps - elastic).abs() < elastic * 0.01,
            "at low stress, RO should be ~elastic: {eps} vs {elastic}"
        );
    }

    #[test]
    fn ramberg_osgood_higher_stress_more_strain() {
        let eps_low = ramberg_osgood_strain(200e9, 1000e6, 10.0, 200e6);
        let eps_high = ramberg_osgood_strain(200e9, 1000e6, 10.0, 400e6);
        assert!(eps_high > eps_low, "higher stress should give more strain");
    }

    #[test]
    fn ramberg_osgood_symmetry() {
        let eps_pos = ramberg_osgood_strain(200e9, 1000e6, 10.0, 300e6);
        let eps_neg = ramberg_osgood_strain(200e9, 1000e6, 10.0, -300e6);
        assert!(
            (eps_pos + eps_neg).abs() < hisab::EPSILON_F64,
            "RO should be antisymmetric"
        );
    }

    #[test]
    fn ramberg_osgood_inverse_roundtrip() {
        let original_stress = 300e6;
        let eps = ramberg_osgood_strain(200e9, 1000e6, 10.0, original_stress);
        let recovered = ramberg_osgood_stress(200e9, 1000e6, 10.0, eps, 1e-6, 50)
            .expect("Newton-Raphson should converge");
        assert!(
            (recovered - original_stress).abs() < 1e3,
            "roundtrip: expected {original_stress}, got {recovered}"
        );
    }

    #[test]
    fn ramberg_osgood_inverse_high_strain() {
        let original_stress = 800e6;
        let eps = ramberg_osgood_strain(200e9, 1000e6, 10.0, original_stress);
        let recovered =
            ramberg_osgood_stress(200e9, 1000e6, 10.0, eps, 1e-6, 50).expect("should converge");
        assert!(
            (recovered - original_stress).abs() < 1e3,
            "high strain roundtrip: expected {original_stress}, got {recovered}"
        );
    }

    #[test]
    fn ramberg_osgood_zero_stress() {
        let eps = ramberg_osgood_strain(200e9, 1000e6, 10.0, 0.0);
        assert!(eps.abs() < hisab::EPSILON_F64);
    }

    #[test]
    fn ramberg_osgood_guard_zero_coeff() {
        assert_eq!(ramberg_osgood_strain(200e9, 0.0, 10.0, 100e6), 0.0);
    }

    // --- Isotropic hardening tests ---

    #[test]
    fn iso_hardening_elastic() {
        let state = IsotropicHardening::new(250e6, 10e9);
        let (stress, new_state) = state.apply_uniaxial(200e9, 0.001);
        assert!((stress - 200e6).abs() < 1.0, "should be elastic");
        assert!(
            new_state.accumulated_plastic_strain.abs() < hisab::EPSILON_F64,
            "no plastic strain in elastic region"
        );
    }

    #[test]
    fn iso_hardening_plastic() {
        let state = IsotropicHardening::new(250e6, 10e9);
        let (stress, new_state) = state.apply_uniaxial(200e9, 0.005);
        assert!(stress > 250e6, "stress should exceed initial yield");
        assert!(
            new_state.accumulated_plastic_strain > 0.0,
            "should have plastic strain"
        );
        assert!(
            new_state.current_yield() > 250e6,
            "yield should have expanded"
        );
    }

    #[test]
    fn iso_hardening_yield_grows() {
        let state = IsotropicHardening::new(250e6, 10e9);
        let (_, state1) = state.apply_uniaxial(200e9, 0.005);
        let (_, state2) = state1.apply_uniaxial(200e9, 0.010);
        assert!(
            state2.current_yield() > state1.current_yield(),
            "yield should keep growing"
        );
    }

    // --- Kinematic hardening tests ---

    #[test]
    fn kin_hardening_elastic() {
        let state = KinematicHardening::new(250e6, 10e9);
        let (stress, new_state) = state.apply_uniaxial(200e9, 0.001);
        assert!((stress - 200e6).abs() < 1.0);
        assert!(
            new_state.back_stress.abs() < hisab::EPSILON_F64,
            "no back-stress shift in elastic"
        );
    }

    #[test]
    fn kin_hardening_plastic() {
        let state = KinematicHardening::new(250e6, 10e9);
        let (_, new_state) = state.apply_uniaxial(200e9, 0.005);
        assert!(
            new_state.back_stress > 0.0,
            "back-stress should shift in tension"
        );
    }

    #[test]
    fn kin_hardening_bauschinger() {
        // After tensile yielding, compressive yield should occur earlier
        let state = KinematicHardening::new(250e6, 10e9);
        let (_, state_after_tension) = state.apply_uniaxial(200e9, 0.005);
        // Compressive yield now at: -(σ_y - back_stress) relative to back_stress
        // i.e., yield in compression at α - σ_y (shifted down)
        let compressive_yield =
            state_after_tension.back_stress - state_after_tension.yield_strength;
        // Original compressive yield would have been -250e6
        assert!(
            compressive_yield > -250e6,
            "Bauschinger: compressive yield should be reduced, got {compressive_yield}"
        );
    }

    // --- Combined hardening tests ---

    #[test]
    fn combined_elastic() {
        let state = CombinedHardening::new(250e6, 5e9, 5e9);
        let (stress, new_state) = state.apply_uniaxial(200e9, 0.001);
        assert!((stress - 200e6).abs() < 1.0);
        assert!(new_state.accumulated_plastic_strain.abs() < hisab::EPSILON_F64);
        assert!(new_state.back_stress.abs() < hisab::EPSILON_F64);
    }

    #[test]
    fn combined_plastic() {
        let state = CombinedHardening::new(250e6, 5e9, 5e9);
        let (stress, new_state) = state.apply_uniaxial(200e9, 0.005);
        assert!(stress > 250e6);
        assert!(new_state.accumulated_plastic_strain > 0.0);
        assert!(new_state.back_stress > 0.0);
        assert!(new_state.current_yield() > 250e6);
    }

    #[test]
    fn combined_reduces_to_iso_when_kin_zero() {
        let iso = IsotropicHardening::new(250e6, 10e9);
        let combined = CombinedHardening::new(250e6, 10e9, 0.0);
        let (s_iso, _) = iso.apply_uniaxial(200e9, 0.005);
        let (s_comb, _) = combined.apply_uniaxial(200e9, 0.005);
        assert!(
            (s_iso - s_comb).abs() < 1.0,
            "combined with H_kin=0 should match pure isotropic"
        );
    }

    #[test]
    fn combined_reduces_to_kin_when_iso_zero() {
        let kin = KinematicHardening::new(250e6, 10e9);
        let combined = CombinedHardening::new(250e6, 0.0, 10e9);
        let (s_kin, _) = kin.apply_uniaxial(200e9, 0.005);
        let (s_comb, _) = combined.apply_uniaxial(200e9, 0.005);
        assert!(
            (s_kin - s_comb).abs() < 1.0,
            "combined with H_iso=0 should match pure kinematic"
        );
    }

    // --- Johnson-Cook tests ---

    #[test]
    fn jc_quasi_static_room_temp() {
        let jc = JohnsonCook::steel_4340();
        let sigma = jc.flow_stress(0.1, 1.0, 293.0);
        // At room temp, quasi-static: σ = A + B*0.1^n = 792e6 + 510e6*0.1^0.26
        assert!(sigma > jc.a, "should exceed initial yield");
        assert!(sigma < 2.0 * jc.a, "should be reasonable");
    }

    #[test]
    fn jc_rate_increases_stress() {
        let jc = JohnsonCook::copper();
        let s_low = jc.flow_stress(0.1, 1.0, 293.0);
        let s_high = jc.flow_stress(0.1, 1000.0, 293.0);
        assert!(s_high > s_low, "higher strain rate should increase stress");
    }

    #[test]
    fn jc_temperature_decreases_stress() {
        let jc = JohnsonCook::copper();
        let s_cold = jc.flow_stress(0.1, 1.0, 293.0);
        let s_hot = jc.flow_stress(0.1, 1.0, 800.0);
        assert!(s_hot < s_cold, "higher temperature should decrease stress");
    }

    #[test]
    fn jc_at_melt_zero_stress() {
        let jc = JohnsonCook::copper();
        let sigma = jc.flow_stress(0.1, 1.0, jc.t_melt);
        assert!(sigma.abs() < 1.0, "at melting point, stress should be ~0");
    }

    #[test]
    fn jc_zero_plastic_strain() {
        let jc = JohnsonCook::steel_4340();
        let sigma = jc.flow_stress(0.0, 1.0, 293.0);
        assert!(
            (sigma - jc.a).abs() < 1.0,
            "at zero plastic strain, stress should be A"
        );
    }

    // --- Neo-Hookean tests ---

    #[test]
    fn neo_hookean_from_moduli() {
        let nh = NeoHookean::from_moduli(76.9e9, 166.7e9);
        assert!((nh.c1 - 38.45e9).abs() < 0.1e9);
    }

    #[test]
    fn neo_hookean_zero_stretch_zero_stress() {
        let nh = NeoHookean::from_elastic(0.01e9, 0.49);
        let sigma = nh.uniaxial_stress(1.0); // λ = 1 means no deformation
        assert!(sigma.abs() < 1.0, "no deformation = no stress");
    }

    #[test]
    fn neo_hookean_tension() {
        let nh = NeoHookean::from_elastic(0.01e9, 0.49); // rubber-like
        let sigma = nh.uniaxial_stress(1.5); // 50% stretch
        assert!(sigma > 0.0, "tension should give positive stress");
    }

    #[test]
    fn neo_hookean_compression() {
        let nh = NeoHookean::from_elastic(0.01e9, 0.49);
        let sigma = nh.uniaxial_stress(0.8); // 20% compression
        assert!(sigma < 0.0, "compression should give negative stress");
    }

    #[test]
    fn neo_hookean_energy_zero_at_reference() {
        let nh = NeoHookean::from_elastic(0.01e9, 0.49);
        let w = nh.strain_energy(3.0, 1.0); // I1=3, J=1 = reference state
        assert!(w.abs() < 1.0, "energy at reference should be ~0");
    }

    #[test]
    fn neo_hookean_tangent_positive() {
        let nh = NeoHookean::from_elastic(0.01e9, 0.49);
        let t = nh.uniaxial_tangent(1.0);
        assert!(t > 0.0, "tangent at reference should be positive");
    }

    // --- Orthotropic 3D tests ---

    #[test]
    fn orthotropic_isotropic_matches_standard() {
        let ort = Orthotropic3D::from_isotropic(200e9, 0.30);
        let c_ort = ort.stiffness_matrix();
        let c_iso = stiffness_matrix(200e9, 0.30);
        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    (c_ort[i][j] - c_iso[i][j]).abs() < 1e3,
                    "orthotropic-as-isotropic should match isotropic at [{i}][{j}]: {} vs {}",
                    c_ort[i][j],
                    c_iso[i][j]
                );
            }
        }
    }

    #[test]
    fn orthotropic_symmetric() {
        let ort = Orthotropic3D {
            e1: 138e9,
            e2: 11e9,
            e3: 11e9,
            g12: 5.5e9,
            g23: 3.9e9,
            g13: 5.5e9,
            nu12: 0.28,
            nu23: 0.40,
            nu13: 0.28,
        };
        let c = ort.stiffness_matrix();
        for (i, row) in c.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    (val - c[j][i]).abs() < 1.0,
                    "should be symmetric at [{i}][{j}]"
                );
            }
        }
    }

    #[test]
    fn orthotropic_diagonal_positive() {
        let ort = Orthotropic3D::from_isotropic(200e9, 0.30);
        let c = ort.stiffness_matrix();
        for (i, row) in c.iter().enumerate() {
            assert!(row[i] > 0.0, "C[{i}][{i}] should be positive");
        }
    }
}
