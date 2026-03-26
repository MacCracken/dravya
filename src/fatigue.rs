/// Basquin's law: cycles to failure.
///
/// N = (σ_a / σ_f')^(1/b)
///
/// σ_a = stress amplitude, σ_f' = fatigue strength coefficient, b = fatigue exponent (negative).
#[must_use]
pub fn basquin_cycles(stress_amplitude: f64, fatigue_strength_coeff: f64, fatigue_exponent: f64) -> f64 {
    if fatigue_strength_coeff.abs() < f64::EPSILON || fatigue_exponent.abs() < f64::EPSILON {
        return 0.0;
    }
    (stress_amplitude / fatigue_strength_coeff).powf(1.0 / fatigue_exponent)
}

/// Miner's rule cumulative damage.
///
/// D = Σ(n_i / N_i)
///
/// Failure when D ≥ 1.0. Input: slice of (cycles_applied, cycles_to_failure) pairs.
#[must_use]
pub fn miners_rule(load_cycles: &[(f64, f64)]) -> f64 {
    load_cycles.iter()
        .map(|&(ni, n_total)| {
            if n_total <= 0.0 { 0.0 } else { ni / n_total }
        })
        .sum()
}

/// Estimate endurance limit from ultimate tensile strength.
///
/// σ_e ≈ 0.5 × σ_uts (for steels, σ_uts < 1400 MPa)
#[must_use]
#[inline]
pub fn endurance_limit_estimate(ultimate_strength: f64) -> f64 {
    0.5 * ultimate_strength
}

/// Check if cumulative damage exceeds threshold (failure predicted).
#[must_use]
#[inline]
pub fn is_fatigue_failure(damage: f64) -> bool {
    damage >= 1.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn miners_rule_basic() {
        // 50% of life at one load + 40% at another = 90% damage
        let d = miners_rule(&[(500.0, 1000.0), (400.0, 1000.0)]);
        assert!((d - 0.9).abs() < f64::EPSILON);
    }

    #[test]
    fn miners_rule_failure() {
        let d = miners_rule(&[(600.0, 1000.0), (500.0, 1000.0)]);
        assert!(is_fatigue_failure(d), "1.1 damage should indicate failure");
    }

    #[test]
    fn endurance_limit_steel() {
        let se = endurance_limit_estimate(500e6);
        assert!((se - 250e6).abs() < f64::EPSILON);
    }

    #[test]
    fn basquin_higher_stress_fewer_cycles() {
        let n_low = basquin_cycles(200e6, 1000e6, -0.1);
        let n_high = basquin_cycles(400e6, 1000e6, -0.1);
        assert!(n_high < n_low, "higher stress should give fewer cycles");
    }

    #[test]
    fn zero_coefficient_safe() {
        assert_eq!(basquin_cycles(100.0, 0.0, -0.1), 0.0);
    }
}
