//! Fatigue analysis: Basquin, Coffin-Manson, Miner's rule, mean-stress corrections, and Marin factors.

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

// --- Coffin-Manson (low-cycle fatigue) ---

/// Coffin-Manson total strain-life equation.
///
/// ε_a = (σ_f' / E) * (2N_f)^b + ε_f' * (2N_f)^c
///
/// Returns total strain amplitude for a given number of cycles.
/// σ_f' = fatigue strength coefficient, b = fatigue strength exponent (negative),
/// ε_f' = fatigue ductility coefficient, c = fatigue ductility exponent (negative).
#[must_use]
pub fn coffin_manson_strain(
    youngs_modulus: f64,
    fatigue_strength_coeff: f64,
    fatigue_strength_exp: f64,
    fatigue_ductility_coeff: f64,
    fatigue_ductility_exp: f64,
    cycles: f64,
) -> f64 {
    if youngs_modulus.abs() < hisab::EPSILON_F64 || cycles <= 0.0 {
        return 0.0;
    }
    let reversals = 2.0 * cycles;
    let elastic = (fatigue_strength_coeff / youngs_modulus) * reversals.powf(fatigue_strength_exp);
    let plastic = fatigue_ductility_coeff * reversals.powf(fatigue_ductility_exp);
    elastic + plastic
}

/// Transition life: number of cycles where elastic and plastic strain
/// amplitudes are equal.
///
/// 2N_t = (ε_f' * E / σ_f')^(1/(b-c))
#[must_use]
pub fn coffin_manson_transition_life(
    youngs_modulus: f64,
    fatigue_strength_coeff: f64,
    fatigue_strength_exp: f64,
    fatigue_ductility_coeff: f64,
    fatigue_ductility_exp: f64,
) -> f64 {
    let exp_diff = fatigue_strength_exp - fatigue_ductility_exp;
    if exp_diff.abs() < hisab::EPSILON_F64 || fatigue_strength_coeff.abs() < hisab::EPSILON_F64 {
        return 0.0;
    }
    let ratio = fatigue_ductility_coeff * youngs_modulus / fatigue_strength_coeff;
    let reversals = ratio.powf(1.0 / exp_diff);
    reversals / 2.0 // convert reversals to cycles
}

// --- Marin endurance limit modification factors ---

/// Surface finish factor k_a (Marin).
///
/// k_a = a * σ_uts^b
///
/// Common coefficients (a in MPa units, b dimensionless):
/// - Ground: a = 1.58, b = -0.085
/// - Machined/cold-drawn: a = 4.51, b = -0.265
/// - Hot-rolled: a = 57.7, b = -0.718
/// - As-forged: a = 272, b = -0.995
///
/// `uts` in Pa, `a_coeff`/`b_exp` from the table above (a in Pa-compatible units).
#[must_use]
#[inline]
pub fn marin_surface_factor(a_coeff: f64, b_exp: f64, uts: f64) -> f64 {
    if uts.abs() < hisab::EPSILON_F64 {
        return 1.0;
    }
    // Convert UTS to MPa for standard Marin coefficients
    let uts_mpa = uts / 1e6;
    a_coeff * uts_mpa.powf(b_exp)
}

/// Size factor k_b (Marin) for rotating round specimens.
///
/// - d <= 8 mm: k_b = 1.0
/// - 8 < d <= 250 mm: k_b = 1.189 * d^(-0.097)
///
/// `diameter` in meters.
#[must_use]
#[inline]
pub fn marin_size_factor(diameter: f64) -> f64 {
    let d_mm = diameter * 1000.0;
    if d_mm <= 8.0 {
        1.0
    } else {
        1.189 * d_mm.powf(-0.097)
    }
}

/// Reliability factor k_c (Marin).
///
/// Standard values:
/// - 50% reliability: 1.000
/// - 90%: 0.897
/// - 95%: 0.868
/// - 99%: 0.814
/// - 99.9%: 0.753
/// - 99.99%: 0.702
///
/// Uses the z-score approximation: k_c = 1 - 0.08 * z
/// where z is the standard normal quantile.
#[must_use]
#[inline]
pub fn marin_reliability_factor(reliability_z: f64) -> f64 {
    (1.0 - 0.08 * reliability_z).clamp(0.5, 1.0)
}

/// Apply Marin factors to unmodified endurance limit.
///
/// S_e = k_a * k_b * k_c * k_d * k_e * S_e'
///
/// `factors` is a slice of all applicable modification factors.
#[must_use]
pub fn marin_corrected_endurance(base_endurance: f64, factors: &[f64]) -> f64 {
    factors.iter().product::<f64>() * base_endurance
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

    // --- Coffin-Manson tests ---

    #[test]
    fn coffin_manson_strain_positive() {
        // Typical steel: σ_f'=900 MPa, b=-0.1, ε_f'=0.25, c=-0.6
        let eps = coffin_manson_strain(200e9, 900e6, -0.1, 0.25, -0.6, 1000.0);
        assert!(eps > 0.0, "strain amplitude should be positive");
    }

    #[test]
    fn coffin_manson_more_cycles_less_strain() {
        let eps_low = coffin_manson_strain(200e9, 900e6, -0.1, 0.25, -0.6, 100.0);
        let eps_high = coffin_manson_strain(200e9, 900e6, -0.1, 0.25, -0.6, 100_000.0);
        assert!(
            eps_high < eps_low,
            "more cycles should correspond to less strain amplitude"
        );
    }

    #[test]
    fn coffin_manson_transition_life_positive() {
        let nt = coffin_manson_transition_life(200e9, 900e6, -0.1, 0.25, -0.6);
        assert!(nt > 0.0, "transition life should be positive");
    }

    #[test]
    fn coffin_manson_at_transition_elastic_equals_plastic() {
        let nt = coffin_manson_transition_life(200e9, 900e6, -0.1, 0.25, -0.6);
        let reversals = 2.0 * nt;
        let elastic = (900e6 / 200e9) * reversals.powf(-0.1);
        let plastic = 0.25 * reversals.powf(-0.6);
        assert!(
            (elastic - plastic).abs() / elastic < 0.01,
            "at transition life, elastic and plastic should be equal: {elastic} vs {plastic}"
        );
    }

    // --- Marin factor tests ---

    #[test]
    fn marin_surface_factor_ground() {
        // Ground: a=1.58, b=-0.085, UTS=400 MPa
        let ka = marin_surface_factor(1.58, -0.085, 400e6);
        assert!(
            ka > 0.9 && ka < 1.0,
            "ground finish ka should be ~0.95, got {ka}"
        );
    }

    #[test]
    fn marin_surface_factor_hot_rolled() {
        // Hot-rolled: a=57.7, b=-0.718, UTS=400 MPa
        let ka = marin_surface_factor(57.7, -0.718, 400e6);
        assert!(
            ka > 0.4 && ka < 0.8,
            "hot-rolled ka should be ~0.5-0.7, got {ka}"
        );
    }

    #[test]
    fn marin_size_factor_small() {
        // 5mm diameter -> k_b = 1.0
        let kb = marin_size_factor(0.005);
        assert_eq!(kb, 1.0);
    }

    #[test]
    fn marin_size_factor_large() {
        // 50mm diameter -> k_b < 1.0
        let kb = marin_size_factor(0.050);
        assert!(kb < 1.0 && kb > 0.7, "50mm kb should be ~0.8, got {kb}");
    }

    #[test]
    fn marin_reliability_50_percent() {
        // z = 0 for 50% reliability
        let kc = marin_reliability_factor(0.0);
        assert!((kc - 1.0).abs() < hisab::EPSILON_F64);
    }

    #[test]
    fn marin_reliability_99_percent() {
        // z ≈ 2.326 for 99%
        let kc = marin_reliability_factor(2.326);
        assert!(
            (kc - 0.814).abs() < 0.01,
            "99% reliability kc should be ~0.814, got {kc}"
        );
    }

    #[test]
    fn marin_corrected_endurance_basic() {
        let se_prime = 250e6; // unmodified
        let ka = 0.9;
        let kb = 0.85;
        let kc = 0.897;
        let se = marin_corrected_endurance(se_prime, &[ka, kb, kc]);
        let expected = se_prime * ka * kb * kc;
        assert!((se - expected).abs() < 1.0);
    }

    #[test]
    fn marin_corrected_no_factors() {
        let se = marin_corrected_endurance(250e6, &[]);
        assert!((se - 250e6).abs() < hisab::EPSILON_F64);
    }
}
