//! Strain tensors (Voigt notation), engineering/true strain, and effective strain.

use std::fmt;
use std::ops::{Add, Mul, Sub};

use serde::{Deserialize, Serialize};

/// Symmetric strain tensor (6 independent components, Voigt notation).
///
/// Components: [εxx, εyy, εzz, γxy, γyz, γxz] where γ = engineering
/// shear strain (γ = 2ε_ij for shear components).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct StrainTensor {
    pub components: [f64; 6],
}

impl StrainTensor {
    /// Create from individual components (engineering shear strain convention).
    #[must_use]
    pub fn new(exx: f64, eyy: f64, ezz: f64, gxy: f64, gyz: f64, gxz: f64) -> Self {
        Self {
            components: [exx, eyy, ezz, gxy, gyz, gxz],
        }
    }

    /// Zero strain state.
    pub const ZERO: Self = Self {
        components: [0.0; 6],
    };

    /// Volumetric strain: ε_vol = εxx + εyy + εzz
    #[must_use]
    #[inline]
    pub fn volumetric(&self) -> f64 {
        self.components[0] + self.components[1] + self.components[2]
    }

    /// Deviatoric strain tensor.
    #[must_use]
    #[inline]
    pub fn deviatoric(&self) -> Self {
        let mean = self.volumetric() / 3.0;
        Self::new(
            self.components[0] - mean,
            self.components[1] - mean,
            self.components[2] - mean,
            self.components[3],
            self.components[4],
            self.components[5],
        )
    }

    /// Von Mises equivalent strain.
    ///
    /// ε_eq = sqrt(2/3 * (e_xx² + e_yy² + e_zz² + γxy²/2 + γyz²/2 + γxz²/2))
    /// where e_ij are deviatoric strain components.
    #[must_use]
    #[inline]
    pub fn effective_strain(&self) -> f64 {
        let dev = self.deviatoric();
        let [exx, eyy, ezz, gxy, gyz, gxz] = dev.components;
        let normal_sum = exx.powi(2) + eyy.powi(2) + ezz.powi(2);
        let shear_sum = (gxy.powi(2) + gyz.powi(2) + gxz.powi(2)) / 2.0;
        ((2.0 / 3.0) * (normal_sum + shear_sum)).sqrt()
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
}

impl Default for StrainTensor {
    fn default() -> Self {
        Self::ZERO
    }
}

impl fmt::Display for StrainTensor {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let [exx, eyy, ezz, gxy, gyz, gxz] = self.components;
        write!(
            f,
            "[εxx={exx:.3e}, εyy={eyy:.3e}, εzz={ezz:.3e}, γxy={gxy:.3e}, γyz={gyz:.3e}, γxz={gxz:.3e}]"
        )
    }
}

impl Add for StrainTensor {
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

impl Sub for StrainTensor {
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

impl Mul<f64> for StrainTensor {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        self.scale(rhs)
    }
}

/// Engineering strain: ε = (L - L0) / L0
#[must_use]
#[inline]
pub fn engineering_strain(original: f64, deformed: f64) -> f64 {
    if original.abs() < hisab::EPSILON_F64 {
        return 0.0;
    }
    (deformed - original) / original
}

/// True (logarithmic) strain: ε_true = ln(L / L0)
#[must_use]
#[inline]
pub fn true_strain(original: f64, deformed: f64) -> f64 {
    if original <= 0.0 || deformed <= 0.0 {
        return 0.0;
    }
    (deformed / original).ln()
}

/// Checked engineering strain. Returns `Err` for zero original length.
pub fn try_engineering_strain(original: f64, deformed: f64) -> crate::Result<f64> {
    if original.abs() < hisab::EPSILON_F64 {
        return Err(crate::DravyaError::DivisionByZero(
            "engineering_strain: zero original length",
        ));
    }
    Ok((deformed - original) / original)
}

/// Checked true strain. Returns `Err` for non-positive inputs.
pub fn try_true_strain(original: f64, deformed: f64) -> crate::Result<f64> {
    if original <= 0.0 {
        return Err(crate::DravyaError::InvalidParameter {
            name: "original",
            value: original,
            reason: "must be positive",
        });
    }
    if deformed <= 0.0 {
        return Err(crate::DravyaError::InvalidParameter {
            name: "deformed",
            value: deformed,
            reason: "must be positive",
        });
    }
    Ok((deformed / original).ln())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn engineering_strain_1_percent() {
        let e = engineering_strain(100.0, 101.0);
        assert!((e - 0.01).abs() < 1e-10, "1% extension, got {e}");
    }

    #[test]
    fn true_strain_1_percent() {
        let e = true_strain(100.0, 101.0);
        assert!(
            (e - 0.00995).abs() < 0.001,
            "true strain ~0.00995 for 1%, got {e}"
        );
    }

    #[test]
    fn true_strain_close_to_engineering_small() {
        let eng = engineering_strain(100.0, 100.5);
        let tru = true_strain(100.0, 100.5);
        assert!(
            (eng - tru).abs() < 0.001,
            "small strains should be nearly equal"
        );
    }

    #[test]
    fn volumetric_strain() {
        let s = StrainTensor::new(0.01, 0.02, 0.03, 0.0, 0.0, 0.0);
        assert!((s.volumetric() - 0.06).abs() < hisab::EPSILON_F64);
    }

    #[test]
    fn zero_length_safe() {
        assert_eq!(engineering_strain(0.0, 1.0), 0.0);
        assert_eq!(true_strain(0.0, 1.0), 0.0);
    }

    #[test]
    fn deviatoric_trace_zero() {
        let s = StrainTensor::new(0.01, 0.02, 0.03, 0.001, 0.002, 0.003);
        let dev = s.deviatoric();
        assert!(
            dev.volumetric().abs() < hisab::EPSILON_F64,
            "deviatoric volumetric strain should be zero"
        );
    }

    #[test]
    fn effective_strain_uniaxial() {
        // Uniaxial strain 0.01 with lateral contraction (v=0.3):
        // exx=0.01, eyy=ezz=-0.003
        let s = StrainTensor::new(0.01, -0.003, -0.003, 0.0, 0.0, 0.0);
        let eff = s.effective_strain();
        assert!(eff > 0.0, "effective strain should be positive");
    }

    #[test]
    fn arithmetic_add() {
        let a = StrainTensor::new(0.01, 0.0, 0.0, 0.0, 0.0, 0.0);
        let b = StrainTensor::new(0.02, 0.0, 0.0, 0.0, 0.0, 0.0);
        let c = a + b;
        assert!((c.components[0] - 0.03).abs() < hisab::EPSILON_F64);
    }

    #[test]
    fn arithmetic_scale() {
        let s = StrainTensor::new(0.01, 0.02, 0.03, 0.0, 0.0, 0.0);
        let scaled = s * 2.0;
        assert!((scaled.components[0] - 0.02).abs() < hisab::EPSILON_F64);
    }

    #[test]
    fn default_is_zero() {
        assert_eq!(StrainTensor::default(), StrainTensor::ZERO);
    }

    #[test]
    fn display_format() {
        let s = StrainTensor::new(0.001, 0.0, 0.0, 0.0, 0.0, 0.0);
        let display = s.to_string();
        assert!(display.contains("εxx="));
    }

    #[test]
    fn serde_roundtrip() {
        let s = StrainTensor::new(0.01, 0.02, 0.03, 0.001, 0.002, 0.003);
        let json = serde_json::to_string(&s).unwrap();
        let back: StrainTensor = serde_json::from_str(&json).unwrap();
        assert_eq!(s, back);
    }

    #[test]
    fn try_engineering_strain_zero_length() {
        let result = try_engineering_strain(0.0, 1.0);
        assert!(result.is_err());
    }

    #[test]
    fn try_engineering_strain_valid() {
        let result = try_engineering_strain(100.0, 101.0);
        assert!(result.is_ok());
        assert!((result.unwrap() - 0.01).abs() < 1e-10);
    }

    #[test]
    fn try_true_strain_negative_original() {
        let result = try_true_strain(-1.0, 1.0);
        assert!(result.is_err());
    }

    #[test]
    fn try_true_strain_negative_deformed() {
        let result = try_true_strain(1.0, -1.0);
        assert!(result.is_err());
    }

    #[test]
    fn try_true_strain_valid() {
        let result = try_true_strain(100.0, 101.0);
        assert!(result.is_ok());
    }
}
