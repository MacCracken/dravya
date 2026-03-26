//! Yield criteria: von Mises, Tresca, Drucker-Prager, and safety factors.

use crate::stress::StressTensor;

/// Check if stress exceeds yield via von Mises criterion.
///
/// Returns `true` when von Mises stress >= yield strength (at or beyond yield).
#[must_use]
#[inline]
pub fn von_mises_check(stress: &StressTensor, yield_strength: f64) -> bool {
    stress.von_mises() >= yield_strength
}

/// Check if stress exceeds yield via Tresca criterion (max shear).
///
/// Returns `true` when max shear stress >= yield_strength / 2.
#[must_use]
#[inline]
pub fn tresca_check(stress: &StressTensor, yield_strength: f64) -> bool {
    stress.max_shear() >= yield_strength / 2.0
}

/// Von Mises safety factor: σ_yield / σ_vm
#[must_use]
#[inline]
pub fn safety_factor(stress: &StressTensor, yield_strength: f64) -> f64 {
    let vm = stress.von_mises();
    if vm.abs() < hisab::EPSILON_F64 {
        return f64::INFINITY;
    }
    yield_strength / vm
}

/// Tresca safety factor: (σ_yield / 2) / τ_max
#[must_use]
#[inline]
pub fn safety_factor_tresca(stress: &StressTensor, yield_strength: f64) -> f64 {
    let tau = stress.max_shear();
    if tau.abs() < hisab::EPSILON_F64 {
        return f64::INFINITY;
    }
    (yield_strength / 2.0) / tau
}

/// Drucker-Prager yield check for pressure-dependent materials
/// (concrete, rock, soil).
///
/// f = sqrt(J2) + alpha * I1 - k
///
/// Yields when f >= 0. Parameters alpha and k relate to cohesion
/// and friction angle of the material:
/// - alpha = 2 sin(phi) / (sqrt(3) * (3 - sin(phi)))
/// - k = 6 c cos(phi) / (sqrt(3) * (3 - sin(phi)))
///
/// where phi = friction angle, c = cohesion.
#[must_use]
#[inline]
pub fn drucker_prager_check(stress: &StressTensor, alpha: f64, k: f64) -> bool {
    stress.j2().sqrt() + alpha * stress.i1() >= k
}

/// Compute Drucker-Prager parameters from Mohr-Coulomb friction angle
/// and cohesion (outer cone approximation matching triaxial compression).
///
/// Returns (alpha, k).
#[must_use]
pub fn drucker_prager_from_mohr_coulomb(friction_angle_rad: f64, cohesion: f64) -> (f64, f64) {
    let sin_phi = friction_angle_rad.sin();
    let cos_phi = friction_angle_rad.cos();
    let denom = 3.0_f64.sqrt() * (3.0 - sin_phi);
    let alpha = 2.0 * sin_phi / denom;
    let k = 6.0 * cohesion * cos_phi / denom;
    (alpha, k)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn below_yield() {
        let s = StressTensor::uniaxial(100e6);
        assert!(
            !von_mises_check(&s, 250e6),
            "100 MPa should not yield at 250 MPa"
        );
    }

    #[test]
    fn above_yield() {
        let s = StressTensor::uniaxial(300e6);
        assert!(
            von_mises_check(&s, 250e6),
            "300 MPa should yield at 250 MPa"
        );
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
        assert!(
            (sf - 2.5).abs() < 0.01,
            "safety factor should be 2.5, got {sf}"
        );
    }

    #[test]
    fn safety_factor_zero_stress() {
        let s = StressTensor::uniaxial(0.0);
        assert!(safety_factor(&s, 250e6).is_infinite());
    }

    #[test]
    fn tresca_safety_factor_uniaxial() {
        let s = StressTensor::uniaxial(100e6);
        // Tresca: τ_max = 50 MPa, allowable = 125 MPa, SF = 2.5
        let sf = safety_factor_tresca(&s, 250e6);
        assert!(
            (sf - 2.5).abs() < 0.01,
            "Tresca SF should be 2.5 for uniaxial, got {sf}"
        );
    }

    #[test]
    fn tresca_safety_factor_zero_stress() {
        let s = StressTensor::uniaxial(0.0);
        assert!(safety_factor_tresca(&s, 250e6).is_infinite());
    }

    #[test]
    fn von_mises_tresca_ordering() {
        // For any stress state, von Mises SF >= Tresca SF
        // (Tresca is more conservative)
        let s = StressTensor::new(100e6, 50e6, 0.0, 30e6, 0.0, 0.0);
        let sf_vm = safety_factor(&s, 500e6);
        let sf_tr = safety_factor_tresca(&s, 500e6);
        assert!(
            sf_tr <= sf_vm + 1e-6,
            "Tresca should be more conservative: VM={sf_vm}, Tresca={sf_tr}"
        );
    }

    #[test]
    fn drucker_prager_hydrostatic() {
        // Pure hydrostatic compression should be safe with reasonable parameters
        let s = StressTensor::hydrostatic_state(-10e6); // compression
        let phi = std::f64::consts::FRAC_PI_6; // 30 deg
        let (alpha, k) = drucker_prager_from_mohr_coulomb(phi, 5e6);
        // Hydrostatic compression: I1 = -30 MPa (negative), J2 ~ 0
        // f = 0 + alpha * (-30e6) - k < 0 (safe)
        assert!(
            !drucker_prager_check(&s, alpha, k),
            "hydrostatic compression should not yield"
        );
    }

    #[test]
    fn drucker_prager_parameters() {
        let phi = std::f64::consts::FRAC_PI_6; // 30 deg
        let (alpha, k) = drucker_prager_from_mohr_coulomb(phi, 5e6);
        assert!(alpha > 0.0, "alpha should be positive");
        assert!(k > 0.0, "k should be positive");
    }
}
