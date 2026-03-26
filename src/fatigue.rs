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

// --- Rainflow cycle counting ---

/// Extract fatigue cycles from a load history using the rainflow algorithm
/// (ASTM E1049-85 three-point method).
///
/// Input: slice of stress peaks/valleys (turning points only — monotonic
/// segments should be pre-filtered to extrema).
///
/// Returns a list of (stress_range, mean_stress, count) tuples.
/// Half-cycles are counted as 0.5. For periodic signals, use
/// [`rainflow_count_periodic`] which eliminates residue.
#[must_use]
pub fn rainflow_count(peaks: &[f64]) -> Vec<(f64, f64, f64)> {
    if peaks.len() < 2 {
        return Vec::new();
    }

    let mut cycles = Vec::new();
    let mut stack: Vec<f64> = Vec::with_capacity(peaks.len());

    for &peak in peaks {
        stack.push(peak);
        while stack.len() >= 3 {
            let n = stack.len();
            let range_last = (stack[n - 1] - stack[n - 2]).abs();
            let range_prev = (stack[n - 2] - stack[n - 3]).abs();

            if range_last >= range_prev {
                let s_max = stack[n - 2].max(stack[n - 3]);
                let s_min = stack[n - 2].min(stack[n - 3]);
                cycles.push((s_max - s_min, (s_max + s_min) / 2.0, 1.0));
                // Pop the last, remove the two cycle points, re-push last
                let last = stack.pop().unwrap_or(0.0);
                stack.pop();
                stack.pop();
                stack.push(last);
            } else {
                break;
            }
        }
    }

    // Remaining points form half-cycles (residue)
    for window in stack.windows(2) {
        let range = (window[1] - window[0]).abs();
        let mean = (window[1] + window[0]) / 2.0;
        cycles.push((range, mean, 0.5));
    }

    cycles
}

/// Rainflow counting for periodic signals.
///
/// Rearranges the signal to start at the overall maximum, wraps the
/// residue, and re-counts to eliminate half-cycles. This is the standard
/// approach per ASTM E1049-85 for repeating load histories.
#[must_use]
pub fn rainflow_count_periodic(peaks: &[f64]) -> Vec<(f64, f64, f64)> {
    if peaks.len() < 2 {
        return Vec::new();
    }

    // Find the index of the overall maximum
    let max_idx = peaks
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i)
        .unwrap_or(0);

    // Rearrange: start from max, wrap around
    let mut rearranged = Vec::with_capacity(peaks.len() + 1);
    rearranged.extend_from_slice(&peaks[max_idx..]);
    rearranged.extend_from_slice(&peaks[..max_idx]);
    rearranged.push(peaks[max_idx]); // close the loop

    rainflow_count(&rearranged)
}

/// Extract turning points (peaks and valleys) from a raw load history.
///
/// Filters out intermediate points on monotonic segments, keeping only
/// local extrema. The first and last points are always included.
#[must_use]
pub fn extract_turning_points(signal: &[f64]) -> Vec<f64> {
    if signal.len() <= 2 {
        return signal.to_vec();
    }

    let mut points = Vec::with_capacity(signal.len() / 2);
    points.push(signal[0]);

    for i in 1..signal.len() - 1 {
        let prev = signal[i - 1];
        let curr = signal[i];
        let next = signal[i + 1];
        if (curr - prev) * (next - curr) < 0.0 {
            points.push(curr);
        }
    }

    points.push(signal[signal.len() - 1]);
    points
}

// --- SN curve interpolation ---

/// Interpolate cycles to failure from tabulated SN data.
///
/// `sn_data`: slice of (stress_amplitude, cycles_to_failure) sorted by
/// descending stress. Uses log-log linear interpolation between data points.
///
/// Returns cycles to failure for the given stress amplitude, or `None`
/// if outside the data range.
#[must_use]
pub fn sn_interpolate(sn_data: &[(f64, f64)], stress_amplitude: f64) -> Option<f64> {
    if sn_data.len() < 2 || stress_amplitude <= 0.0 {
        return None;
    }

    let log_s = stress_amplitude.ln();

    // Find bracketing interval (data sorted descending by stress)
    for window in sn_data.windows(2) {
        let (s_high, n_high) = window[0];
        let (s_low, n_low) = window[1];

        if stress_amplitude <= s_high && stress_amplitude >= s_low {
            if s_high <= 0.0 || s_low <= 0.0 || n_high <= 0.0 || n_low <= 0.0 {
                return None;
            }
            let log_s_high = s_high.ln();
            let log_s_low = s_low.ln();
            let log_n_high = n_high.ln();
            let log_n_low = n_low.ln();

            let denom = log_s_high - log_s_low;
            if denom.abs() < hisab::EPSILON_F64 {
                return Some(n_high);
            }
            let t = (log_s - log_s_low) / denom;
            let log_n = log_n_low + t * (log_n_high - log_n_low);
            return Some(log_n.exp());
        }
    }

    None
}

// --- Neuber's rule ---

/// Neuber's rule for elastic-plastic notch correction.
///
/// Given the elastic stress concentration factor Kt and the nominal
/// stress/strain, Neuber's rule gives:
///
/// σ_local * ε_local = Kt² * σ_nominal * ε_nominal
///
/// For the elastic case (σ = Eε):
/// σ_local * ε_local = Kt² * σ_nominal² / E
///
/// Returns the Neuber product (σ * ε at the notch root).
#[must_use]
#[inline]
pub fn neuber_product(kt: f64, nominal_stress: f64, youngs_modulus: f64) -> f64 {
    if youngs_modulus.abs() < hisab::EPSILON_F64 {
        return 0.0;
    }
    kt * kt * nominal_stress * nominal_stress / youngs_modulus
}

/// Solve Neuber's rule with Ramberg-Osgood material for notch-root stress.
///
/// Finds local stress σ satisfying:
/// σ * (σ/E + (σ/K)^n) = Kt² * S² / E
///
/// Uses [`hisab::num::newton_raphson`]. Returns (notch_stress, notch_strain),
/// or `Err` if the solver fails to converge.
pub fn neuber_ramberg_osgood(
    kt: f64,
    nominal_stress: f64,
    youngs_modulus: f64,
    strength_coeff: f64,
    hardening_exp: f64,
    tol: f64,
    max_iter: usize,
) -> crate::Result<(f64, f64)> {
    let target = neuber_product(kt, nominal_stress, youngs_modulus);
    if target.abs() < hisab::EPSILON_F64 {
        return Ok((0.0, 0.0));
    }

    let ro_strain = |sigma: f64| {
        crate::constitutive::ramberg_osgood_strain(
            youngs_modulus,
            strength_coeff,
            hardening_exp,
            sigma,
        )
    };

    let f = |sigma: f64| sigma * ro_strain(sigma) - target;
    let df = |sigma: f64| {
        let eps = ro_strain(sigma);
        let d_elastic = 1.0 / youngs_modulus;
        let d_plastic = if sigma.abs() < hisab::EPSILON_F64 {
            0.0
        } else {
            (hardening_exp / strength_coeff)
                * (sigma.abs() / strength_coeff).powf(hardening_exp - 1.0)
        };
        eps + sigma * (d_elastic + d_plastic)
    };

    let x0 = kt * nominal_stress;
    let sigma = hisab::num::newton_raphson(f, df, x0, tol, max_iter).map_err(|_| {
        crate::DravyaError::SolverNoConvergence {
            method: "neuber_ramberg_osgood",
            iterations: max_iter,
        }
    })?;
    let eps = ro_strain(sigma);
    Ok((sigma, eps))
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

    // --- Rainflow tests ---

    #[test]
    fn rainflow_simple_cycle() {
        // One full cycle: 0 -> 100 -> 0 -> 100 -> 0
        let peaks = vec![0.0, 100.0, 0.0, 100.0, 0.0];
        let cycles = rainflow_count(&peaks);
        assert!(!cycles.is_empty(), "should extract cycles");
        // Should find at least one full cycle with range=100
        let full_cycles: Vec<_> = cycles.iter().filter(|c| c.2 >= 1.0).collect();
        assert!(
            !full_cycles.is_empty(),
            "should have at least one full cycle"
        );
    }

    #[test]
    fn rainflow_empty() {
        assert!(rainflow_count(&[]).is_empty());
        assert!(rainflow_count(&[100.0]).is_empty());
    }

    #[test]
    fn rainflow_two_points() {
        let cycles = rainflow_count(&[0.0, 100.0]);
        assert_eq!(cycles.len(), 1);
        assert!(
            (cycles[0].0 - 100.0).abs() < hisab::EPSILON_F64,
            "range=100"
        );
        assert!((cycles[0].2 - 0.5).abs() < hisab::EPSILON_F64, "half-cycle");
    }

    #[test]
    fn rainflow_preserves_damage() {
        // Total damage should be consistent
        let peaks = vec![0.0, 200.0, 50.0, 180.0, 20.0, 150.0, 0.0];
        let cycles = rainflow_count(&peaks);
        let total_count: f64 = cycles.iter().map(|c| c.2).sum();
        assert!(total_count > 0.0, "should have positive cycle count");
    }

    // --- SN interpolation tests ---

    #[test]
    fn sn_interpolate_basic() {
        // Simple SN data (descending stress)
        let sn_data = vec![
            (300e6, 1e4), // 300 MPa -> 10^4 cycles
            (200e6, 1e5), // 200 MPa -> 10^5 cycles
            (100e6, 1e6), // 100 MPa -> 10^6 cycles
        ];
        let n = sn_interpolate(&sn_data, 200e6);
        assert!(n.is_some());
        assert!((n.unwrap() - 1e5).abs() < 1e3, "exact point should match");
    }

    #[test]
    fn sn_interpolate_between_points() {
        let sn_data = vec![(300e6, 1e4), (100e6, 1e6)];
        let n = sn_interpolate(&sn_data, 200e6);
        assert!(n.is_some());
        let cycles = n.unwrap();
        assert!(cycles > 1e4 && cycles < 1e6, "should be between bounds");
    }

    #[test]
    fn sn_interpolate_out_of_range() {
        let sn_data = vec![(300e6, 1e4), (100e6, 1e6)];
        assert!(sn_interpolate(&sn_data, 400e6).is_none());
        assert!(sn_interpolate(&sn_data, 50e6).is_none());
    }

    #[test]
    fn sn_interpolate_higher_stress_fewer_cycles() {
        let sn_data = vec![(300e6, 1e4), (200e6, 1e5), (100e6, 1e6)];
        let n_high = sn_interpolate(&sn_data, 250e6).unwrap();
        let n_low = sn_interpolate(&sn_data, 150e6).unwrap();
        assert!(n_high < n_low, "higher stress = fewer cycles");
    }

    // --- Neuber tests ---

    #[test]
    fn neuber_product_basic() {
        // Kt=2, σ_nom=100 MPa, E=200 GPa
        // Product = 4 * (100e6)^2 / 200e9 = 4e16 / 2e11 = 200000 Pa
        let np = neuber_product(2.0, 100e6, 200e9);
        assert!(
            (np - 200_000.0).abs() < 1.0,
            "Neuber product should be 200000, got {np}"
        );
    }

    #[test]
    fn neuber_product_kt1() {
        // Kt=1 (no notch): σ*ε = σ²/E (purely elastic)
        let np = neuber_product(1.0, 100e6, 200e9);
        let elastic = 100e6 * 100e6 / 200e9;
        assert!((np - elastic).abs() < hisab::EPSILON_F64);
    }

    #[test]
    fn neuber_ramberg_osgood_converges() {
        let (sigma, eps) = neuber_ramberg_osgood(2.0, 100e6, 200e9, 1000e6, 10.0, 1e-6, 50)
            .expect("should converge");
        assert!(sigma > 0.0);
        assert!(eps > 0.0);
        // Verify Neuber's rule: σ*ε ≈ Kt²*S²/E
        let target = neuber_product(2.0, 100e6, 200e9);
        let actual = sigma * eps;
        assert!(
            (actual - target).abs() / target < 0.001,
            "Neuber product mismatch: {actual} vs {target}"
        );
    }

    #[test]
    fn neuber_ramberg_osgood_notch_exceeds_nominal() {
        let (sigma, _) = neuber_ramberg_osgood(3.0, 100e6, 200e9, 1000e6, 10.0, 1e-6, 50)
            .expect("should converge");
        assert!(sigma > 100e6, "notch stress should exceed nominal");
    }
}
