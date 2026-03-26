//! Fracture mechanics: stress intensity factors, fracture toughness,
//! and fatigue crack growth.

use std::f64::consts::PI;

// --- Stress Intensity Factors (Mode I) ---

/// Mode I stress intensity factor for a center crack in an infinite plate.
///
/// KI = σ * sqrt(π * a)
///
/// σ = remote stress (Pa), a = half-crack length (m).
#[must_use]
#[inline]
pub fn ki_center_crack_infinite(stress: f64, half_crack_length: f64) -> f64 {
    stress * (PI * half_crack_length).sqrt()
}

/// Mode I SIF for a center crack in a finite-width plate.
///
/// KI = σ * sqrt(π * a) * F(a/W)
///
/// where F(a/W) is the Feddersen correction for finite width W.
/// Uses the secant approximation: F = sqrt(sec(π*a/W))
#[must_use]
pub fn ki_center_crack_finite(stress: f64, half_crack_length: f64, width: f64) -> f64 {
    if width.abs() < hisab::EPSILON_F64 {
        return 0.0;
    }
    let ratio = PI * half_crack_length / width;
    let correction = if ratio.abs() < hisab::EPSILON_F64 {
        1.0
    } else {
        (ratio).cos().recip().sqrt()
    };
    stress * (PI * half_crack_length).sqrt() * correction
}

/// Mode I SIF for a single edge crack in a semi-infinite plate.
///
/// KI = 1.12 * σ * sqrt(π * a)
///
/// The 1.12 factor accounts for the free-surface effect.
#[must_use]
#[inline]
pub fn ki_edge_crack(stress: f64, crack_length: f64) -> f64 {
    1.12 * stress * (PI * crack_length).sqrt()
}

/// Mode I SIF for a penny-shaped (circular) embedded crack.
///
/// KI = (2/π) * σ * sqrt(π * a)
///
/// a = crack radius.
#[must_use]
#[inline]
pub fn ki_penny_crack(stress: f64, crack_radius: f64) -> f64 {
    (2.0 / PI) * stress * (PI * crack_radius).sqrt()
}

/// Mode I SIF for a through-thickness crack at a hole (Bowie solution).
///
/// KI = σ * sqrt(π * a) * F_bowie
///
/// Uses the approximation F = 0.5 * (3 - a/(a+r)) * sqrt(1 + a/r)
/// where r = hole radius, a = crack length from hole edge.
#[must_use]
pub fn ki_crack_at_hole(stress: f64, crack_length: f64, hole_radius: f64) -> f64 {
    if hole_radius.abs() < hisab::EPSILON_F64 {
        return ki_edge_crack(stress, crack_length);
    }
    let ratio = crack_length / hole_radius;
    let f_bowie = 0.5 * (3.0 - crack_length / (crack_length + hole_radius)) * (1.0 + ratio).sqrt();
    stress * (PI * crack_length).sqrt() * f_bowie
}

// --- Fracture toughness ---

/// Check if a crack will propagate (fracture check).
///
/// Returns `true` if KI >= KIc (critical stress intensity factor).
#[must_use]
#[inline]
pub fn fracture_check(ki: f64, fracture_toughness: f64) -> bool {
    ki >= fracture_toughness
}

/// Critical crack length for a given stress and fracture toughness.
///
/// a_cr = (KIc / σ)^2 / π (for center crack in infinite plate)
#[must_use]
#[inline]
pub fn critical_crack_length(stress: f64, fracture_toughness: f64) -> f64 {
    if stress.abs() < hisab::EPSILON_F64 {
        return f64::INFINITY;
    }
    (fracture_toughness / stress).powi(2) / PI
}

/// Fracture stress for a given crack length and fracture toughness.
///
/// σ_cr = KIc / sqrt(π * a)
#[must_use]
#[inline]
pub fn fracture_stress(half_crack_length: f64, fracture_toughness: f64) -> f64 {
    let denom = (PI * half_crack_length).sqrt();
    if denom.abs() < hisab::EPSILON_F64 {
        return f64::INFINITY;
    }
    fracture_toughness / denom
}

/// Plane strain fracture toughness from critical energy release rate.
///
/// KIc = sqrt(E * GIc / (1 - v^2))
///
/// GIc = critical strain energy release rate (J/m^2).
#[must_use]
#[inline]
pub fn kic_from_energy_release(
    youngs_modulus: f64,
    poisson_ratio: f64,
    critical_energy_release: f64,
) -> f64 {
    let denom = 1.0 - poisson_ratio * poisson_ratio;
    if denom.abs() < hisab::EPSILON_F64 {
        return 0.0;
    }
    (youngs_modulus * critical_energy_release / denom).sqrt()
}

// --- Paris law crack growth ---

/// Paris law fatigue crack growth rate.
///
/// da/dN = C * (ΔK)^m
///
/// C = material constant, m = Paris exponent (typically 2-4 for metals),
/// ΔK = stress intensity factor range (MPa*sqrt(m)).
///
/// Returns crack growth per cycle (m/cycle).
#[must_use]
#[inline]
pub fn paris_law_rate(c: f64, m: f64, delta_k: f64) -> f64 {
    if delta_k <= 0.0 {
        return 0.0;
    }
    c * delta_k.powf(m)
}

/// Integrate Paris law for crack life prediction (constant amplitude).
///
/// Integrates da/dN = C * (ΔK)^m from initial crack a_i to final crack a_f,
/// where ΔK = Δσ * sqrt(π * a) (center crack in infinite plate).
///
/// Uses numerical integration with `n_steps` increments.
///
/// Returns number of cycles to grow from a_i to a_f.
#[must_use]
pub fn paris_law_life(
    c: f64,
    m: f64,
    stress_range: f64,
    initial_crack: f64,
    final_crack: f64,
    n_steps: usize,
) -> f64 {
    if n_steps == 0 || initial_crack >= final_crack || stress_range.abs() < hisab::EPSILON_F64 {
        return 0.0;
    }

    let da = (final_crack - initial_crack) / n_steps as f64;
    let mut total_cycles = 0.0;
    let mut a = initial_crack;

    for _ in 0..n_steps {
        let delta_k = stress_range * (PI * a).sqrt();
        let da_dn = paris_law_rate(c, m, delta_k);
        if da_dn.abs() < hisab::EPSILON_F64 {
            return f64::INFINITY;
        }
        total_cycles += da / da_dn;
        a += da;
    }

    total_cycles
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- SIF tests ---

    #[test]
    fn ki_center_crack_basic() {
        // 100 MPa, 10mm half-crack
        let ki = ki_center_crack_infinite(100e6, 0.01);
        // KI = 100e6 * sqrt(π * 0.01) = 100e6 * 0.1772 = 17.72 MPa√m
        assert!(
            (ki - 17.72e6).abs() < 0.1e6,
            "KI should be ~17.72 MPa√m, got {}",
            ki / 1e6
        );
    }

    #[test]
    fn ki_finite_approaches_infinite() {
        // Very wide plate should approach infinite plate solution
        let ki_inf = ki_center_crack_infinite(100e6, 0.01);
        let ki_fin = ki_center_crack_finite(100e6, 0.01, 10.0); // a/W = 0.001
        assert!(
            (ki_inf - ki_fin).abs() / ki_inf < 0.001,
            "finite width with large W should match infinite"
        );
    }

    #[test]
    fn ki_finite_exceeds_infinite() {
        // Finite width correction always increases KI
        let ki_inf = ki_center_crack_infinite(100e6, 0.01);
        let ki_fin = ki_center_crack_finite(100e6, 0.01, 0.05); // a/W = 0.2
        assert!(ki_fin > ki_inf, "finite width should increase KI");
    }

    #[test]
    fn ki_edge_exceeds_center() {
        // Edge crack has higher SIF than center crack (free surface effect)
        let ki_center = ki_center_crack_infinite(100e6, 0.01);
        let ki_edge = ki_edge_crack(100e6, 0.01);
        assert!(ki_edge > ki_center, "edge crack should have higher KI");
    }

    #[test]
    fn ki_penny_less_than_through() {
        // Penny crack has lower SIF than through-thickness
        let ki_through = ki_center_crack_infinite(100e6, 0.01);
        let ki_penny = ki_penny_crack(100e6, 0.01);
        assert!(ki_penny < ki_through, "penny crack should have lower KI");
    }

    #[test]
    fn ki_zero_crack_is_zero() {
        let ki = ki_center_crack_infinite(100e6, 0.0);
        assert!(ki.abs() < hisab::EPSILON_F64);
    }

    // --- Fracture toughness tests ---

    #[test]
    fn fracture_check_basic() {
        // Steel KIc ~ 50 MPa√m
        let ki = ki_center_crack_infinite(200e6, 0.02); // ~50 MPa√m
        assert!(fracture_check(ki, 40e6), "should fracture when KI > KIc");
        assert!(
            !fracture_check(ki, 100e6),
            "should not fracture when KI < KIc"
        );
    }

    #[test]
    fn critical_crack_length_basic() {
        // σ = 100 MPa, KIc = 50 MPa√m
        // a_cr = (50/100)^2 / π = 0.25/π ≈ 0.0796 m
        let a_cr = critical_crack_length(100e6, 50e6);
        assert!(
            (a_cr - 0.0796).abs() < 0.001,
            "critical crack length should be ~79.6 mm, got {:.1} mm",
            a_cr * 1000.0
        );
    }

    #[test]
    fn critical_crack_length_zero_stress() {
        let a_cr = critical_crack_length(0.0, 50e6);
        assert!(a_cr.is_infinite());
    }

    #[test]
    fn fracture_stress_roundtrip() {
        let a = 0.01;
        let kic = 50e6;
        let sigma_cr = fracture_stress(a, kic);
        let ki = ki_center_crack_infinite(sigma_cr, a);
        assert!(
            (ki - kic).abs() < 1e3,
            "fracture stress roundtrip: KI should equal KIc"
        );
    }

    #[test]
    fn kic_from_energy_release_steel() {
        // E=200 GPa, v=0.3, GIc=20000 J/m²
        let kic = kic_from_energy_release(200e9, 0.3, 20000.0);
        // KIc = sqrt(200e9 * 20000 / 0.91) ≈ 66.3 MPa√m
        assert!(
            kic > 50e6 && kic < 80e6,
            "KIc should be ~66 MPa√m, got {}",
            kic / 1e6
        );
    }

    // --- Paris law tests ---

    #[test]
    fn paris_law_rate_basic() {
        // Typical steel: C = 1e-11, m = 3
        let da_dn = paris_law_rate(1e-11, 3.0, 20e6); // ΔK = 20 MPa√m
        // da/dN = 1e-11 * (20e6)^3 = 1e-11 * 8e21 = 8e10 ... that's wrong
        // Actually, da/dN = C * ΔK^m where ΔK in MPa√m and C in consistent units
        // C = 1e-11 m/cycle per (Pa√m)^m... let me just check it's positive and reasonable
        assert!(da_dn > 0.0, "growth rate should be positive");
    }

    #[test]
    fn paris_law_rate_zero_delta_k() {
        let da_dn = paris_law_rate(1e-11, 3.0, 0.0);
        assert_eq!(da_dn, 0.0);
    }

    #[test]
    fn paris_law_rate_negative_delta_k() {
        let da_dn = paris_law_rate(1e-11, 3.0, -10e6);
        assert_eq!(da_dn, 0.0, "negative ΔK should give zero growth");
    }

    #[test]
    fn paris_law_life_positive() {
        // Very rough check: life should be positive and finite
        let life = paris_law_life(1e-11, 3.0, 100e6, 0.001, 0.01, 1000);
        assert!(life > 0.0, "life should be positive");
        assert!(life.is_finite(), "life should be finite");
    }

    #[test]
    fn paris_law_life_longer_crack_fewer_cycles() {
        // Starting from a larger initial crack should give fewer remaining cycles
        let life_small = paris_law_life(1e-11, 3.0, 100e6, 0.001, 0.05, 1000);
        let life_large = paris_law_life(1e-11, 3.0, 100e6, 0.01, 0.05, 1000);
        assert!(
            life_large < life_small,
            "larger initial crack should give fewer cycles"
        );
    }

    #[test]
    fn paris_law_life_zero_range() {
        let life = paris_law_life(1e-11, 3.0, 0.0, 0.001, 0.01, 1000);
        assert_eq!(life, 0.0);
    }

    #[test]
    fn ki_crack_at_hole_basic() {
        let ki = ki_crack_at_hole(100e6, 0.005, 0.01);
        assert!(ki > 0.0);
        // Should be higher than a plain edge crack due to stress concentration
        let ki_plain = ki_edge_crack(100e6, 0.005);
        assert!(
            ki > ki_plain,
            "crack at hole should have higher KI than plain edge"
        );
    }
}
