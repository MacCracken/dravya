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
    hisab::num::newton_raphson(f, df, x0, tol, max_iter)
        .map_err(|e| crate::DravyaError::ComputationError(format!("ramberg_osgood_stress: {e}")))
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
}
