/// Basquin's law: cycles to failure (high-cycle fatigue).
///
/// N = (σ_a / A)^(1/b)
///
/// σ_a = stress amplitude, A = fatigue strength coefficient
/// (cycle-based), b = fatigue exponent (negative, typically -0.05 to -0.12).
///
/// Note: if using the standard reversal-based coefficient σ_f', convert
/// via A = σ_f' * 2^b, or use [`basquin_cycles_reversals`] directly.
#[must_use]
pub fn basquin_cycles(
    stress_amplitude: f64,
    fatigue_strength_coeff: f64,
    fatigue_exponent: f64,
) -> f64 {
    if fatigue_strength_coeff.abs() < hisab::EPSILON_F64
        || fatigue_exponent.abs() < hisab::EPSILON_F64
        || stress_amplitude <= 0.0
    {
        return 0.0;
    }
    (stress_amplitude / fatigue_strength_coeff).powf(1.0 / fatigue_exponent)
}

/// Basquin's law using the standard reversal-based form:
///
/// σ_a = σ_f' * (2N_f)^b
///
/// Solving for N_f: N_f = 0.5 * (σ_a / σ_f')^(1/b)
///
/// σ_f' = fatigue strength coefficient (from handbooks, reversal-based),
/// b = fatigue exponent (negative).
#[must_use]
pub fn basquin_cycles_reversals(
    stress_amplitude: f64,
    fatigue_strength_coeff: f64,
    fatigue_exponent: f64,
) -> f64 {
    if fatigue_strength_coeff.abs() < hisab::EPSILON_F64
        || fatigue_exponent.abs() < hisab::EPSILON_F64
        || stress_amplitude <= 0.0
    {
        return 0.0;
    }
    0.5 * (stress_amplitude / fatigue_strength_coeff).powf(1.0 / fatigue_exponent)
}

/// Miner's rule cumulative damage.
///
/// D = Σ(n_i / N_i)
///
/// Failure when D >= 1.0. Input: slice of (cycles_applied, cycles_to_failure) pairs.
#[must_use]
pub fn miners_rule(load_cycles: &[(f64, f64)]) -> f64 {
    load_cycles
        .iter()
        .map(
            |&(ni, n_total)| {
                if n_total <= 0.0 { 0.0 } else { ni / n_total }
            },
        )
        .sum()
}

/// Estimate endurance limit from ultimate tensile strength.
///
/// For wrought steels under rotating bending:
/// - S_e' = 0.5 * σ_uts  when σ_uts <= 1400 MPa
/// - S_e' = 700 MPa       when σ_uts > 1400 MPa
///
/// This is the unmodified endurance limit (laboratory specimen).
/// Apply Marin factors for real components.
#[must_use]
#[inline]
pub fn endurance_limit_estimate(ultimate_strength: f64) -> f64 {
    let raw = 0.5 * ultimate_strength;
    if raw > 700e6 { 700e6 } else { raw }
}

/// Check if cumulative damage exceeds threshold (failure predicted).
#[must_use]
#[inline]
pub fn is_fatigue_failure(damage: f64) -> bool {
    damage >= 1.0
}

/// Goodman mean-stress correction.
///
/// Returns the equivalent fully-reversed stress amplitude:
/// σ_ar = σ_a / (1 - σ_m / σ_uts)
///
/// σ_a = stress amplitude, σ_m = mean stress, σ_uts = ultimate tensile strength.
/// Returns 0.0 if σ_m >= σ_uts (static failure region).
#[must_use]
#[inline]
pub fn goodman_correction(stress_amplitude: f64, mean_stress: f64, uts: f64) -> f64 {
    let denom = 1.0 - mean_stress / uts;
    if denom <= hisab::EPSILON_F64 {
        return 0.0;
    }
    stress_amplitude / denom
}

/// Gerber mean-stress correction (parabolic).
///
/// σ_ar = σ_a / (1 - (σ_m / σ_uts)^2)
#[must_use]
#[inline]
pub fn gerber_correction(stress_amplitude: f64, mean_stress: f64, uts: f64) -> f64 {
    let ratio = mean_stress / uts;
    let denom = 1.0 - ratio * ratio;
    if denom <= hisab::EPSILON_F64 {
        return 0.0;
    }
    stress_amplitude / denom
}

/// Soderberg mean-stress correction (conservative, uses yield).
///
/// σ_ar = σ_a / (1 - σ_m / σ_y)
#[must_use]
#[inline]
pub fn soderberg_correction(stress_amplitude: f64, mean_stress: f64, yield_strength: f64) -> f64 {
    let denom = 1.0 - mean_stress / yield_strength;
    if denom <= hisab::EPSILON_F64 {
        return 0.0;
    }
    stress_amplitude / denom
}

/// Stress ratio R = σ_min / σ_max.
///
/// R = -1 is fully reversed, R = 0 is zero-to-tension, R = 1 is static.
#[must_use]
#[inline]
pub fn stress_ratio(sigma_min: f64, sigma_max: f64) -> f64 {
    if sigma_max.abs() < hisab::EPSILON_F64 {
        return 0.0;
    }
    sigma_min / sigma_max
}

/// Decompose stress range into amplitude and mean.
///
/// σ_a = (σ_max - σ_min) / 2, σ_m = (σ_max + σ_min) / 2
///
/// Returns (amplitude, mean).
#[must_use]
#[inline]
pub fn stress_amplitude_mean(sigma_min: f64, sigma_max: f64) -> (f64, f64) {
    let amplitude = (sigma_max - sigma_min) / 2.0;
    let mean = (sigma_max + sigma_min) / 2.0;
    (amplitude, mean)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn miners_rule_basic() {
        let d = miners_rule(&[(500.0, 1000.0), (400.0, 1000.0)]);
        assert!((d - 0.9).abs() < hisab::EPSILON_F64);
    }

    #[test]
    fn miners_rule_failure() {
        let d = miners_rule(&[(600.0, 1000.0), (500.0, 1000.0)]);
        assert!(is_fatigue_failure(d), "1.1 damage should indicate failure");
    }

    #[test]
    fn miners_rule_empty() {
        assert_eq!(miners_rule(&[]), 0.0);
    }

    #[test]
    fn endurance_limit_steel() {
        let se = endurance_limit_estimate(500e6);
        assert!((se - 250e6).abs() < hisab::EPSILON_F64);
    }

    #[test]
    fn endurance_limit_cap() {
        // High-strength steel: UTS = 2000 MPa -> capped at 700 MPa
        let se = endurance_limit_estimate(2000e6);
        assert!(
            (se - 700e6).abs() < hisab::EPSILON_F64,
            "endurance limit should cap at 700 MPa, got {}",
            se / 1e6
        );
    }

    #[test]
    fn basquin_higher_stress_fewer_cycles() {
        let n_low = basquin_cycles(200e6, 1000e6, -0.1);
        let n_high = basquin_cycles(400e6, 1000e6, -0.1);
        assert!(n_high < n_low, "higher stress should give fewer cycles");
    }

    #[test]
    fn basquin_negative_amplitude_safe() {
        assert_eq!(basquin_cycles(-100.0, 1000e6, -0.1), 0.0);
    }

    #[test]
    fn zero_coefficient_safe() {
        assert_eq!(basquin_cycles(100.0, 0.0, -0.1), 0.0);
    }

    #[test]
    fn basquin_reversals_half_of_cycles() {
        let n_cycles = basquin_cycles(200e6, 1000e6, -0.1);
        let n_rev = basquin_cycles_reversals(200e6, 1000e6, -0.1);
        assert!(
            (n_rev - 0.5 * n_cycles).abs() < 1.0,
            "reversal form should give half the cycles"
        );
    }

    #[test]
    fn goodman_correction_basic() {
        // σ_a = 100 MPa, σ_m = 200 MPa, σ_uts = 400 MPa
        // σ_ar = 100 / (1 - 200/400) = 100 / 0.5 = 200 MPa
        let sar = goodman_correction(100e6, 200e6, 400e6);
        assert!(
            (sar - 200e6).abs() < 1.0,
            "Goodman correction should be 200 MPa, got {}",
            sar / 1e6
        );
    }

    #[test]
    fn goodman_at_uts_returns_zero() {
        let sar = goodman_correction(100e6, 400e6, 400e6);
        assert_eq!(sar, 0.0, "at UTS, should return 0 (static failure)");
    }

    #[test]
    fn gerber_less_conservative_than_goodman() {
        let g = goodman_correction(100e6, 100e6, 400e6);
        let ge = gerber_correction(100e6, 100e6, 400e6);
        assert!(ge < g, "Gerber should be less conservative than Goodman");
    }

    #[test]
    fn soderberg_more_conservative_than_goodman() {
        let good = goodman_correction(100e6, 100e6, 400e6);
        let sod = soderberg_correction(100e6, 100e6, 250e6);
        assert!(
            sod > good,
            "Soderberg (using yield) should be more conservative"
        );
    }

    #[test]
    fn stress_ratio_fully_reversed() {
        let r = stress_ratio(-100e6, 100e6);
        assert!((r - (-1.0)).abs() < hisab::EPSILON_F64, "R should be -1");
    }

    #[test]
    fn stress_amplitude_mean_decomposition() {
        let (amp, mean) = stress_amplitude_mean(-50e6, 150e6);
        assert!((amp - 100e6).abs() < 1.0);
        assert!((mean - 50e6).abs() < 1.0);
    }
}
