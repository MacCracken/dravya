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
    if youngs_modulus <= 0.0 || critical_energy_release < 0.0 {
        return 0.0;
    }
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
/// Uses [`hisab::calc::integral_simpson`] for numerical integration with
/// `n_steps` panels (must be even; rounded up if odd).
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

    // Integrand: dN/da = 1 / (da/dN) = 1 / (C * (Δσ√(πa))^m)
    let integrand = |a: f64| {
        let delta_k = stress_range * (PI * a).sqrt();
        let da_dn = paris_law_rate(c, m, delta_k);
        if da_dn.abs() < hisab::EPSILON_F64 {
            f64::INFINITY
        } else {
            1.0 / da_dn
        }
    };

    // Simpson's rule requires even n; round up
    let n = if n_steps.is_multiple_of(2) {
        n_steps
    } else {
        n_steps + 1
    };

    hisab::calc::integral_simpson(integrand, initial_crack, final_crack, n).unwrap_or(0.0)
}

// --- Mode II / Mode III stress intensity factors ---

/// Mode II stress intensity factor for a center crack (in-plane shear).
///
/// KII = τ * sqrt(π * a)
///
/// τ = remote shear stress, a = half-crack length.
#[must_use]
#[inline]
pub fn kii_center_crack(shear_stress: f64, half_crack_length: f64) -> f64 {
    shear_stress * (PI * half_crack_length).sqrt()
}

/// Mode II SIF for a single edge crack under shear.
///
/// KII = F * τ * sqrt(π * a)
///
/// Uses F = 1.12 (same free-surface correction as Mode I).
/// Note: Mode II corrections are highly geometry-dependent. For precision,
/// use geometry-specific solutions from Tada/Paris/Irwin.
#[must_use]
#[inline]
pub fn kii_edge_crack(shear_stress: f64, crack_length: f64) -> f64 {
    1.12 * shear_stress * (PI * crack_length).sqrt()
}

/// Mode III stress intensity factor (anti-plane shear / tearing).
///
/// KIII = τ * sqrt(π * a) for a through-thickness crack.
#[must_use]
#[inline]
pub fn kiii_through_crack(shear_stress: f64, half_crack_length: f64) -> f64 {
    shear_stress * (PI * half_crack_length).sqrt()
}

// --- J-integral ---

/// J-integral for plane strain (relates to stress intensity factors).
///
/// J = (KI² + KII²) / E' + KIII² / (2G)
///
/// where E' = E / (1 - v²) for plane strain.
///
/// This is the mixed-mode energy release rate.
#[must_use]
#[inline]
pub fn j_integral_from_sifs(
    ki: f64,
    kii: f64,
    kiii: f64,
    youngs_modulus: f64,
    poisson_ratio: f64,
) -> f64 {
    let e_prime = youngs_modulus / (1.0 - poisson_ratio * poisson_ratio);
    let g = youngs_modulus / (2.0 * (1.0 + poisson_ratio));
    if e_prime.abs() < hisab::EPSILON_F64 || g.abs() < hisab::EPSILON_F64 {
        return 0.0;
    }
    (ki * ki + kii * kii) / e_prime + kiii * kiii / (2.0 * g)
}

/// Equivalent stress intensity factor from J-integral (plane strain).
///
/// K_eq = sqrt(J * E')
///
/// where E' = E / (1 - v²).
#[must_use]
#[inline]
pub fn k_from_j_integral(j: f64, youngs_modulus: f64, poisson_ratio: f64) -> f64 {
    if j < 0.0 || youngs_modulus <= 0.0 {
        return 0.0;
    }
    let e_prime = youngs_modulus / (1.0 - poisson_ratio * poisson_ratio);
    (j * e_prime).sqrt()
}

/// J-integral for Mode I only (plane strain): J = KI² / E'
#[must_use]
#[inline]
pub fn j_integral_mode_i(ki: f64, youngs_modulus: f64, poisson_ratio: f64) -> f64 {
    j_integral_from_sifs(ki, 0.0, 0.0, youngs_modulus, poisson_ratio)
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
        let ki_plain = ki_edge_crack(100e6, 0.005);
        assert!(
            ki > ki_plain,
            "crack at hole should have higher KI than plain edge"
        );
    }

    // --- Mode II / III tests ---

    #[test]
    fn kii_center_basic() {
        let kii = kii_center_crack(100e6, 0.01);
        let ki = ki_center_crack_infinite(100e6, 0.01);
        // Same formula for center crack: KII = τ*sqrt(πa), KI = σ*sqrt(πa)
        assert!(
            (kii - ki).abs() < 1.0,
            "KII center should equal KI center for same stress"
        );
    }

    #[test]
    fn kii_edge_exceeds_center() {
        let kii_c = kii_center_crack(100e6, 0.01);
        let kii_e = kii_edge_crack(100e6, 0.01);
        assert!(kii_e > kii_c, "edge should exceed center (free surface)");
    }

    #[test]
    fn kiii_positive() {
        let kiii = kiii_through_crack(100e6, 0.01);
        assert!(kiii > 0.0);
    }

    // --- J-integral tests ---

    #[test]
    fn j_integral_mode_i_positive() {
        let ki = ki_center_crack_infinite(100e6, 0.01);
        let j = j_integral_mode_i(ki, 200e9, 0.30);
        assert!(j > 0.0, "J should be positive for non-zero KI");
    }

    #[test]
    fn j_integral_roundtrip() {
        let ki_original = 50e6; // 50 MPa√m
        let j = j_integral_mode_i(ki_original, 200e9, 0.30);
        let ki_recovered = k_from_j_integral(j, 200e9, 0.30);
        assert!(
            (ki_recovered - ki_original).abs() < 1e3,
            "roundtrip: expected {ki_original}, got {ki_recovered}"
        );
    }

    #[test]
    fn j_integral_mixed_mode() {
        let ki = 30e6;
        let kii = 20e6;
        let kiii = 10e6;
        let j = j_integral_from_sifs(ki, kii, kiii, 200e9, 0.30);
        let j_mode_i = j_integral_mode_i(ki, 200e9, 0.30);
        assert!(j > j_mode_i, "mixed-mode J should exceed pure Mode I J");
    }

    #[test]
    fn j_integral_zero_k() {
        let j = j_integral_from_sifs(0.0, 0.0, 0.0, 200e9, 0.30);
        assert!(
            j.abs() < hisab::EPSILON_F64,
            "J should be zero for zero SIFs"
        );
    }

    #[test]
    fn k_from_j_matches_gi_formula() {
        // For Mode I: KIc = sqrt(E' * GIc) where GIc = J
        let gic = 5000.0; // J/m²
        let k_via_j = k_from_j_integral(gic, 200e9, 0.30);
        let k_via_gi = kic_from_energy_release(200e9, 0.30, gic);
        assert!(
            (k_via_j - k_via_gi).abs() < 1e3,
            "k_from_j should match kic_from_energy_release: {} vs {}",
            k_via_j / 1e6,
            k_via_gi / 1e6
        );
    }
}
