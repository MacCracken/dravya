use crate::stress::StressTensor;

/// Check if stress exceeds yield via von Mises criterion.
#[must_use]
#[inline]
pub fn von_mises_check(stress: &StressTensor, yield_strength: f64) -> bool {
    stress.von_mises() >= yield_strength
}

/// Check if stress exceeds yield via Tresca criterion (max shear).
#[must_use]
#[inline]
pub fn tresca_check(stress: &StressTensor, yield_strength: f64) -> bool {
    stress.max_shear() >= yield_strength / 2.0
}

/// Safety factor: σ_yield / σ_applied
#[must_use]
#[inline]
pub fn safety_factor(stress: &StressTensor, yield_strength: f64) -> f64 {
    let vm = stress.von_mises();
    if vm.abs() < f64::EPSILON { return f64::INFINITY; }
    yield_strength / vm
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn below_yield() {
        let s = StressTensor::uniaxial(100e6);
        assert!(!von_mises_check(&s, 250e6), "100 MPa should not yield at 250 MPa");
    }

    #[test]
    fn above_yield() {
        let s = StressTensor::uniaxial(300e6);
        assert!(von_mises_check(&s, 250e6), "300 MPa should yield at 250 MPa");
    }

    #[test]
    fn tresca_uniaxial() {
        let s = StressTensor::uniaxial(250e6);
        assert!(tresca_check(&s, 250e6));
    }

    #[test]
    fn safety_factor_basic() {
        let s = StressTensor::uniaxial(100e6);
        let sf = safety_factor(&s, 250e6);
        assert!((sf - 2.5).abs() < 0.01, "safety factor should be 2.5, got {sf}");
    }

    #[test]
    fn safety_factor_zero_stress() {
        let s = StressTensor::uniaxial(0.0);
        assert!(safety_factor(&s, 250e6).is_infinite());
    }
}
