use serde::{Deserialize, Serialize};

/// Symmetric strain tensor (6 independent components, Voigt notation).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct StrainTensor {
    pub components: [f64; 6],
}

impl StrainTensor {
    #[must_use]
    pub fn new(exx: f64, eyy: f64, ezz: f64, gxy: f64, gyz: f64, gxz: f64) -> Self {
        Self { components: [exx, eyy, ezz, gxy, gyz, gxz] }
    }

    /// Volumetric strain: ε_vol = εxx + εyy + εzz
    #[must_use]
    #[inline]
    pub fn volumetric(&self) -> f64 {
        self.components[0] + self.components[1] + self.components[2]
    }
}

/// Engineering strain: ε = (L - L₀) / L₀
#[must_use]
#[inline]
pub fn engineering_strain(original: f64, deformed: f64) -> f64 {
    if original.abs() < f64::EPSILON { return 0.0; }
    (deformed - original) / original
}

/// True (logarithmic) strain: ε_true = ln(L / L₀)
#[must_use]
#[inline]
pub fn true_strain(original: f64, deformed: f64) -> f64 {
    if original <= 0.0 || deformed <= 0.0 { return 0.0; }
    (deformed / original).ln()
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
        assert!((e - 0.00995).abs() < 0.001, "true strain ~0.00995 for 1%, got {e}");
    }

    #[test]
    fn true_strain_close_to_engineering_small() {
        let eng = engineering_strain(100.0, 100.5);
        let tru = true_strain(100.0, 100.5);
        assert!((eng - tru).abs() < 0.001, "small strains should be nearly equal");
    }

    #[test]
    fn volumetric_strain() {
        let s = StrainTensor::new(0.01, 0.02, 0.03, 0.0, 0.0, 0.0);
        assert!((s.volumetric() - 0.06).abs() < f64::EPSILON);
    }

    #[test]
    fn zero_length_safe() {
        assert_eq!(engineering_strain(0.0, 1.0), 0.0);
        assert_eq!(true_strain(0.0, 1.0), 0.0);
    }
}
