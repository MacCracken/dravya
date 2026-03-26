//! Symmetric stress tensor in Voigt notation with invariants, yield measures, and arithmetic.

use std::fmt;
use std::ops::{Add, Mul, Sub};

use serde::{Deserialize, Serialize};

/// Symmetric stress tensor (6 independent components).
/// [σxx, σyy, σzz, τxy, τyz, τxz] in Voigt notation.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct StressTensor {
    pub components: [f64; 6],
}

impl StressTensor {
    /// Create from individual components.
    #[must_use]
    pub fn new(sxx: f64, syy: f64, szz: f64, txy: f64, tyz: f64, txz: f64) -> Self {
        Self {
            components: [sxx, syy, szz, txy, tyz, txz],
        }
    }

    /// Zero stress state.
    pub const ZERO: Self = Self {
        components: [0.0; 6],
    };

    /// Create a uniaxial stress state (tension/compression along x).
    #[must_use]
    pub fn uniaxial(sigma: f64) -> Self {
        Self::new(sigma, 0.0, 0.0, 0.0, 0.0, 0.0)
    }

    /// Create a hydrostatic (equal triaxial) stress state.
    #[must_use]
    pub fn hydrostatic_state(pressure: f64) -> Self {
        Self::new(pressure, pressure, pressure, 0.0, 0.0, 0.0)
    }

    /// Create a pure shear stress state in the xy-plane.
    #[must_use]
    pub fn pure_shear(tau: f64) -> Self {
        Self::new(0.0, 0.0, 0.0, tau, 0.0, 0.0)
    }

    /// Scale all components by a scalar.
    #[must_use]
    #[inline]
    pub fn scale(self, factor: f64) -> Self {
        Self {
            components: [
                self.components[0] * factor,
                self.components[1] * factor,
                self.components[2] * factor,
                self.components[3] * factor,
                self.components[4] * factor,
                self.components[5] * factor,
            ],
        }
    }

    /// Hydrostatic (mean) stress: σ_m = (σxx + σyy + σzz) / 3
    #[must_use]
    #[inline]
    pub fn hydrostatic(&self) -> f64 {
        (self.components[0] + self.components[1] + self.components[2]) / 3.0
    }

    /// Deviatoric stress tensor: s_ij = σ_ij - σ_m * δ_ij
    #[must_use]
    #[inline]
    pub fn deviatoric(&self) -> Self {
        let h = self.hydrostatic();
        Self::new(
            self.components[0] - h,
            self.components[1] - h,
            self.components[2] - h,
            self.components[3],
            self.components[4],
            self.components[5],
        )
    }

    /// First stress invariant: I1 = σxx + σyy + σzz (trace).
    #[must_use]
    #[inline]
    pub fn i1(&self) -> f64 {
        self.components[0] + self.components[1] + self.components[2]
    }

    /// Second stress invariant: I2 = σxx*σyy + σyy*σzz + σzz*σxx - τxy² - τyz² - τxz²
    #[must_use]
    #[inline]
    pub fn i2(&self) -> f64 {
        let [sxx, syy, szz, txy, tyz, txz] = self.components;
        sxx * syy + syy * szz + szz * sxx - txy * txy - tyz * tyz - txz * txz
    }

    /// Third stress invariant (determinant): I3 = det(σ)
    #[must_use]
    #[inline]
    pub fn i3(&self) -> f64 {
        let [sxx, syy, szz, txy, tyz, txz] = self.components;
        sxx * syy * szz + 2.0 * txy * tyz * txz
            - sxx * tyz * tyz
            - syy * txz * txz
            - szz * txy * txy
    }

    /// Second deviatoric invariant: J2 = 1/2 * s_ij * s_ij
    ///
    /// Equivalent to σ_vm² / 3.
    #[must_use]
    #[inline]
    pub fn j2(&self) -> f64 {
        let [sxx, syy, szz, txy, tyz, txz] = self.components;
        let term1 = (sxx - syy).powi(2) + (syy - szz).powi(2) + (szz - sxx).powi(2);
        let term2 = 6.0 * (txy.powi(2) + tyz.powi(2) + txz.powi(2));
        (term1 + term2) / 6.0
    }

    /// Von Mises equivalent stress: σ_vm = sqrt(3 * J2)
    #[must_use]
    #[inline]
    pub fn von_mises(&self) -> f64 {
        (3.0 * self.j2()).sqrt()
    }

    /// Maximum shear stress (Tresca criterion).
    ///
    /// τ_max = (σ_max - σ_min) / 2
    #[must_use]
    pub fn max_shear(&self) -> f64 {
        let principals = self.principal_stresses();
        (principals[0] - principals[2]) / 2.0
    }

    /// Principal stresses (sorted descending: σ1 >= σ2 >= σ3).
    ///
    /// Uses the closed-form analytical solution for eigenvalues of a 3x3
    /// symmetric matrix via deviatoric stress invariants. Zero heap allocations.
    #[must_use]
    pub fn principal_stresses(&self) -> [f64; 3] {
        let mean = self.hydrostatic();
        let j2 = self.j2();

        if j2 < 1e-30 {
            // Hydrostatic state — all principals equal
            return [mean, mean, mean];
        }

        // Third deviatoric invariant J3 = det(s)
        let dev = self.deviatoric();
        let [s11, s22, s33, s12, s23, s13] = dev.components;
        let j3 = s11 * s22 * s33 + 2.0 * s12 * s23 * s13
            - s11 * s23 * s23
            - s22 * s13 * s13
            - s33 * s12 * s12;

        // Lode angle: cos(3θ) = (3√3/2) * J3 / J2^(3/2)
        let r = 1.5 * (3.0_f64).sqrt() * j3 / j2.powf(1.5);
        let lode_angle = (r.clamp(-1.0, 1.0)).acos() / 3.0;

        let mag = 2.0 * (j2 / 3.0).sqrt();
        let two_pi_3 = 2.0 * std::f64::consts::PI / 3.0;

        let mut principals = [
            mean + mag * lode_angle.cos(),
            mean + mag * (lode_angle - two_pi_3).cos(),
            mean + mag * (lode_angle + two_pi_3).cos(),
        ];
        principals.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
        principals
    }

    /// Octahedral shear stress: τ_oct = sqrt(2/3) * sqrt(J2)
    #[must_use]
    #[inline]
    pub fn octahedral_shear(&self) -> f64 {
        (2.0_f64 / 3.0).sqrt() * self.j2().sqrt()
    }
}

impl Default for StressTensor {
    fn default() -> Self {
        Self::ZERO
    }
}

impl fmt::Display for StressTensor {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let [sxx, syy, szz, txy, tyz, txz] = self.components;
        write!(
            f,
            "[σxx={sxx:.3e}, σyy={syy:.3e}, σzz={szz:.3e}, τxy={txy:.3e}, τyz={tyz:.3e}, τxz={txz:.3e}]"
        )
    }
}

impl Add for StressTensor {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self {
            components: [
                self.components[0] + rhs.components[0],
                self.components[1] + rhs.components[1],
                self.components[2] + rhs.components[2],
                self.components[3] + rhs.components[3],
                self.components[4] + rhs.components[4],
                self.components[5] + rhs.components[5],
            ],
        }
    }
}

impl Sub for StressTensor {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self {
            components: [
                self.components[0] - rhs.components[0],
                self.components[1] - rhs.components[1],
                self.components[2] - rhs.components[2],
                self.components[3] - rhs.components[3],
                self.components[4] - rhs.components[4],
                self.components[5] - rhs.components[5],
            ],
        }
    }
}

impl Mul<f64> for StressTensor {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        self.scale(rhs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniaxial_von_mises_equals_applied() {
        let s = StressTensor::uniaxial(100e6);
        let vm = s.von_mises();
        assert!(
            (vm - 100e6).abs() < 1e3,
            "uniaxial von Mises should equal applied stress, got {vm}"
        );
    }

    #[test]
    fn hydrostatic_basic() {
        let s = StressTensor::new(30.0, 60.0, 90.0, 0.0, 0.0, 0.0);
        assert!((s.hydrostatic() - 60.0).abs() < hisab::EPSILON_F64);
    }

    #[test]
    fn pure_shear_von_mises() {
        // Pure shear τ = 100 MPa -> σ_vm = τ*sqrt(3) ~ 173.2 MPa
        let s = StressTensor::new(0.0, 0.0, 0.0, 100e6, 0.0, 0.0);
        let vm = s.von_mises();
        assert!(
            (vm - 173.2e6).abs() < 0.5e6,
            "pure shear von Mises should be t*sqrt(3), got {vm}"
        );
    }

    #[test]
    fn principal_stresses_uniaxial() {
        let s = StressTensor::uniaxial(100.0);
        let p = s.principal_stresses();
        assert!(
            p[0] >= p[1] && p[1] >= p[2],
            "principals should be sorted descending"
        );
        assert!(
            (p[0] - 100.0).abs() < 1e-6,
            "σ1 should be 100, got {}",
            p[0]
        );
        assert!(p[1].abs() < 1e-6, "σ2 should be ~0, got {}", p[1]);
        assert!(p[2].abs() < 1e-6, "σ3 should be ~0, got {}", p[2]);
    }

    #[test]
    fn principal_stresses_sorted_descending() {
        let s = StressTensor::new(10.0, 50.0, 30.0, 5.0, 3.0, 2.0);
        let p = s.principal_stresses();
        assert!(
            p[0] >= p[1] && p[1] >= p[2],
            "principals should be sorted descending"
        );
        let trace = s.i1();
        let sum = p[0] + p[1] + p[2];
        assert!(
            (sum - trace).abs() < 1e-6,
            "principal sum should equal trace={trace}, got {sum}"
        );
    }

    #[test]
    fn principal_stresses_hydrostatic() {
        let s = StressTensor::hydrostatic_state(100.0);
        let p = s.principal_stresses();
        for (i, &pi) in p.iter().enumerate() {
            assert!(
                (pi - 100.0).abs() < 1e-6,
                "principal σ{} should be 100, got {pi}",
                i + 1
            );
        }
    }

    #[test]
    fn max_shear_uniaxial() {
        let s = StressTensor::uniaxial(100.0);
        let tau = s.max_shear();
        assert!(
            (tau - 50.0).abs() < 1e-6,
            "max shear of uniaxial should be σ/2, got {tau}"
        );
    }

    #[test]
    fn deviatoric_trace_is_zero() {
        let s = StressTensor::new(100.0, 50.0, 30.0, 10.0, 5.0, 3.0);
        let dev = s.deviatoric();
        assert!(
            dev.i1().abs() < hisab::EPSILON_F64,
            "deviatoric trace should be zero, got {}",
            dev.i1()
        );
    }

    #[test]
    fn deviatoric_preserves_shear() {
        let s = StressTensor::new(100.0, 50.0, 30.0, 10.0, 5.0, 3.0);
        let dev = s.deviatoric();
        assert_eq!(dev.components[3], s.components[3]);
        assert_eq!(dev.components[4], s.components[4]);
        assert_eq!(dev.components[5], s.components[5]);
    }

    #[test]
    fn invariants_uniaxial() {
        let s = StressTensor::uniaxial(100.0);
        assert!((s.i1() - 100.0).abs() < hisab::EPSILON_F64);
        assert!(s.i2().abs() < hisab::EPSILON_F64);
        assert!(s.i3().abs() < hisab::EPSILON_F64);
    }

    #[test]
    fn j2_von_mises_relationship() {
        let s = StressTensor::new(100.0, 50.0, 30.0, 10.0, 5.0, 3.0);
        let vm = s.von_mises();
        let j2 = s.j2();
        assert!((vm * vm - 3.0 * j2).abs() < 1e-6, "σ_vm² should equal 3*J2");
    }

    #[test]
    fn arithmetic_add() {
        let a = StressTensor::uniaxial(100.0);
        let b = StressTensor::uniaxial(50.0);
        let c = a + b;
        assert!((c.components[0] - 150.0).abs() < hisab::EPSILON_F64);
    }

    #[test]
    fn arithmetic_sub() {
        let a = StressTensor::uniaxial(100.0);
        let b = StressTensor::uniaxial(40.0);
        let c = a - b;
        assert!((c.components[0] - 60.0).abs() < hisab::EPSILON_F64);
    }

    #[test]
    fn arithmetic_scale() {
        let s = StressTensor::uniaxial(100.0);
        let scaled = s * 2.0;
        assert!((scaled.components[0] - 200.0).abs() < hisab::EPSILON_F64);
    }

    #[test]
    fn default_is_zero() {
        assert_eq!(StressTensor::default(), StressTensor::ZERO);
    }

    #[test]
    fn display_format() {
        let s = StressTensor::uniaxial(100e6);
        let display = s.to_string();
        assert!(display.contains("σxx="));
    }

    #[test]
    fn octahedral_shear_uniaxial() {
        let s = StressTensor::uniaxial(100.0);
        let tau_oct = s.octahedral_shear();
        // τ_oct = sqrt(2)/3 * σ for uniaxial
        let expected = (2.0_f64).sqrt() / 3.0 * 100.0;
        assert!(
            (tau_oct - expected).abs() < 1e-6,
            "octahedral shear got {tau_oct}, expected {expected}"
        );
    }
}
