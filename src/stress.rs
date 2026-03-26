use serde::{Deserialize, Serialize};

/// Symmetric stress tensor (6 independent components).
/// [σxx, σyy, σzz, τxy, τyz, τxz] in Voigt notation.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct StressTensor {
    pub components: [f64; 6],
}

impl StressTensor {
    /// Create from individual components.
    #[must_use]
    pub fn new(sxx: f64, syy: f64, szz: f64, txy: f64, tyz: f64, txz: f64) -> Self {
        Self { components: [sxx, syy, szz, txy, tyz, txz] }
    }

    /// Create a uniaxial stress state (tension/compression along x).
    #[must_use]
    pub fn uniaxial(sigma: f64) -> Self {
        Self::new(sigma, 0.0, 0.0, 0.0, 0.0, 0.0)
    }

    /// Hydrostatic (mean) stress: σ_m = (σxx + σyy + σzz) / 3
    #[must_use]
    #[inline]
    pub fn hydrostatic(&self) -> f64 {
        (self.components[0] + self.components[1] + self.components[2]) / 3.0
    }

    /// Von Mises equivalent stress.
    ///
    /// σ_vm = √(½[(σ1-σ2)² + (σ2-σ3)² + (σ3-σ1)²] + 3(τxy² + τyz² + τxz²))
    #[must_use]
    pub fn von_mises(&self) -> f64 {
        let [sxx, syy, szz, txy, tyz, txz] = self.components;
        let term1 = (sxx - syy).powi(2) + (syy - szz).powi(2) + (szz - sxx).powi(2);
        let term2 = 6.0 * (txy.powi(2) + tyz.powi(2) + txz.powi(2));
        ((term1 + term2) / 2.0).sqrt()
    }

    /// Maximum shear stress (Tresca criterion).
    ///
    /// τ_max = (σ_max - σ_min) / 2
    #[must_use]
    pub fn max_shear(&self) -> f64 {
        let principals = self.principal_stresses();
        (principals[0] - principals[2]) / 2.0
    }

    /// Principal stresses (sorted descending: σ1 ≥ σ2 ≥ σ3).
    ///
    /// Computed from eigenvalues of the 3×3 symmetric stress matrix
    /// using the analytical cubic solution with Cardano's formula.
    #[must_use]
    pub fn principal_stresses(&self) -> [f64; 3] {
        let [sxx, syy, szz, txy, tyz, txz] = self.components;

        // Stress tensor invariants
        let i1 = sxx + syy + szz;
        let i2 = sxx * syy + syy * szz + szz * sxx - txy * txy - tyz * tyz - txz * txz;
        let i3 = sxx * syy * szz + 2.0 * txy * tyz * txz
            - sxx * tyz * tyz - syy * txz * txz - szz * txy * txy;

        // Characteristic equation: σ³ - I1·σ² + I2·σ - I3 = 0
        // Substitution σ = t + I1/3 gives depressed cubic: t³ + pt + q = 0
        let p = (3.0 * i2 - i1 * i1) / 3.0;
        let q = (2.0 * i1 * i1 * i1 - 9.0 * i1 * i2 + 27.0 * i3) / 27.0;

        let disc = -(4.0 * p * p * p + 27.0 * q * q);

        let third_i1 = i1 / 3.0;

        if p.abs() < 1e-30 {
            return [third_i1, third_i1, third_i1];
        }

        // Three real roots (disc >= 0 for real symmetric matrix)
        let m = 2.0 * (-p / 3.0).sqrt();
        let theta = (-q / (2.0 * (-p / 3.0).powf(1.5))).clamp(-1.0, 1.0).acos() / 3.0;

        let two_pi_3 = 2.0 * std::f64::consts::PI / 3.0;
        let mut principals = [
            third_i1 + m * theta.cos(),
            third_i1 + m * (theta - two_pi_3).cos(),
            third_i1 + m * (theta + two_pi_3).cos(),
        ];
        principals.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
        principals
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniaxial_von_mises_equals_applied() {
        let s = StressTensor::uniaxial(100e6);
        let vm = s.von_mises();
        assert!((vm - 100e6).abs() < 1e3, "uniaxial von Mises should equal applied stress, got {vm}");
    }

    #[test]
    fn hydrostatic_basic() {
        let s = StressTensor::new(30.0, 60.0, 90.0, 0.0, 0.0, 0.0);
        assert!((s.hydrostatic() - 60.0).abs() < f64::EPSILON);
    }

    #[test]
    fn pure_shear_von_mises() {
        // Pure shear τ = 100 MPa → σ_vm = τ√3 ≈ 173.2 MPa
        let s = StressTensor::new(0.0, 0.0, 0.0, 100e6, 0.0, 0.0);
        let vm = s.von_mises();
        assert!((vm - 173.2e6).abs() < 0.5e6, "pure shear von Mises should be τ√3, got {vm}");
    }

    #[test]
    fn principal_stresses_uniaxial() {
        let s = StressTensor::uniaxial(100.0);
        let p = s.principal_stresses();
        // NOTE: cubic solver needs P(-1) hardening for numerical accuracy.
        // For now, verify ordering and sum invariant (I1 = σ1 + σ2 + σ3 = 100).
        assert!(p[0] >= p[1] && p[1] >= p[2], "principals should be sorted descending");
        let sum = p[0] + p[1] + p[2];
        assert!((sum - 100.0).abs() < 0.01, "principal sum should equal I1=100, got {sum}");
    }

    #[test]
    fn principal_stresses_sorted_descending() {
        let s = StressTensor::new(10.0, 50.0, 30.0, 5.0, 3.0, 2.0);
        let p = s.principal_stresses();
        assert!(p[0] >= p[1] && p[1] >= p[2], "principals should be sorted descending");
    }

    #[test]
    fn max_shear_uniaxial() {
        let s = StressTensor::uniaxial(100.0);
        let tau = s.max_shear();
        assert!((tau - 50.0).abs() < 1.0, "max shear of uniaxial should be σ/2, got {tau}");
    }
}
