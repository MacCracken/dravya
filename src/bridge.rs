//! Cross-crate bridges — convert primitive values from other AGNOS science crates
//! into dravya material science parameters and vice versa.
//!
//! Always available — takes primitive values (f64), no science crate deps.
//!
//! # Architecture
//!
//! ```text
//! sharira (physiology)          ──┐
//! impetus (physics)               ┼──> bridge ──> dravya material parameters
//! ushma   (thermodynamics)       ┘
//! ```

// ── Sharira bridges (physiology/biomechanics) ──────────────────────────────

/// Convert bone density (kg/m³) to estimated Young's modulus (Pa).
///
/// Carter & Hayes (1977) power-law: E = 3790 × ρ³ (for cancellous bone)
/// where ρ is apparent density in g/cm³. We convert from kg/m³.
#[must_use]
#[inline]
pub fn bone_density_to_youngs_modulus(density_kg_m3: f64) -> f64 {
    if density_kg_m3 <= 0.0 {
        return 0.0;
    }
    let density_g_cm3 = density_kg_m3 / 1000.0;
    3790.0e6 * density_g_cm3.powi(3) // Pa
}

/// Convert bone density (kg/m³) to estimated yield strength (Pa).
///
/// Empirical: σ_y ≈ 137 × ρ^1.88 (MPa), ρ in g/cm³.
/// Keller (1994) for trabecular bone compression.
#[must_use]
#[inline]
pub fn bone_density_to_yield_strength(density_kg_m3: f64) -> f64 {
    if density_kg_m3 <= 0.0 {
        return 0.0;
    }
    let density_g_cm3 = density_kg_m3 / 1000.0;
    137.0e6 * density_g_cm3.powf(1.88) // Pa
}

/// Convert muscle force (N) and tendon cross-section area (m²)
/// to tendon stress (Pa).
///
/// σ = F / A
#[must_use]
#[inline]
pub fn muscle_force_to_tendon_stress(force_n: f64, cross_section_m2: f64) -> f64 {
    if cross_section_m2 <= 0.0 {
        return 0.0;
    }
    force_n.abs() / cross_section_m2
}

/// Convert tendon strain and stiffness to tendon force (N).
///
/// F = E × A × ε (Hooke's law for uniaxial tension)
#[must_use]
#[inline]
pub fn tendon_strain_to_force(strain: f64, youngs_modulus_pa: f64, cross_section_m2: f64) -> f64 {
    youngs_modulus_pa * cross_section_m2 * strain
}

/// Estimate bone safety factor from applied stress (Pa) and bone density (kg/m³).
///
/// SF = σ_yield / σ_applied, where σ_yield is estimated from bone density.
#[must_use]
#[inline]
pub fn bone_safety_factor(applied_stress_pa: f64, density_kg_m3: f64) -> f64 {
    if applied_stress_pa.abs() < 1e-15 {
        return f64::INFINITY;
    }
    let yield_strength = bone_density_to_yield_strength(density_kg_m3);
    yield_strength / applied_stress_pa.abs()
}

/// Convert joint impact energy (J) and bone toughness (J/m³) to
/// estimated fracture volume (m³).
///
/// V_fracture = E_impact / toughness
#[must_use]
#[inline]
pub fn impact_to_fracture_volume(impact_energy_j: f64, toughness_j_per_m3: f64) -> f64 {
    if toughness_j_per_m3 <= 0.0 {
        return 0.0;
    }
    impact_energy_j.max(0.0) / toughness_j_per_m3
}

// ── Impetus bridges (physics) ──────────────────────────────────────────────

/// Convert collision force (N) and contact area (m²) to contact stress (Pa).
///
/// σ = F / A
#[must_use]
#[inline]
pub fn collision_to_contact_stress(force_n: f64, contact_area_m2: f64) -> f64 {
    if contact_area_m2 <= 0.0 {
        return 0.0;
    }
    force_n.abs() / contact_area_m2
}

/// Convert deformation velocity (m/s) and characteristic length (m)
/// to strain rate (1/s).
///
/// ε̇ = v / L
#[must_use]
#[inline]
pub fn velocity_to_strain_rate(velocity_ms: f64, length_m: f64) -> f64 {
    if length_m <= 0.0 {
        return 0.0;
    }
    velocity_ms.abs() / length_m
}

// ── Ushma bridges (thermodynamics) ─────────────────────────────────────────

/// Convert temperature difference (K) and thermal expansion coefficient (1/K)
/// to thermal strain (dimensionless).
///
/// ε = α × ΔT
#[must_use]
#[inline]
pub fn temperature_to_thermal_strain(delta_t_k: f64, expansion_coeff_per_k: f64) -> f64 {
    expansion_coeff_per_k * delta_t_k
}

/// Convert thermal gradient (K/m) to thermal stress (Pa) for a constrained body.
///
/// σ = E × α × ΔT
#[must_use]
#[inline]
pub fn thermal_gradient_to_stress(
    gradient_k_per_m: f64,
    youngs_modulus_pa: f64,
    expansion_coeff_per_k: f64,
    length_m: f64,
) -> f64 {
    youngs_modulus_pa * expansion_coeff_per_k * gradient_k_per_m * length_m
}

// ── Bijli bridges (electromagnetism) ───────────────────────────────────────

/// Convert electric field strength (V/m) to piezoelectric stress (Pa).
///
/// σ = e × E, where e is the piezoelectric stress constant (C/m²).
/// Typical PZT: e ≈ 5–25 C/m².
#[must_use]
#[inline]
pub fn e_field_to_piezo_stress(e_field_v_m: f64, piezo_coeff_c_m2: f64) -> f64 {
    e_field_v_m * piezo_coeff_c_m2
}

/// Convert magnetic flux density (T) to magnetostrictive strain.
///
/// ε = (3/2) × λ_s × (B/B_sat)² for cubic symmetry.
#[must_use]
#[inline]
pub fn magnetic_to_magnetostrictive_strain(
    flux_density_t: f64,
    saturation_magnetostriction: f64,
    saturation_flux_t: f64,
) -> f64 {
    if saturation_flux_t <= 0.0 {
        return 0.0;
    }
    let ratio = flux_density_t / saturation_flux_t;
    1.5 * saturation_magnetostriction * ratio * ratio
}

// ── Khanij bridges (geology) ──────────────────────────────────────────────

/// Convert mineral composition to bulk material density (kg/m³).
///
/// Weighted average: ρ = Σ(f_i × ρ_i) where f_i is volume fraction.
#[must_use]
pub fn mineral_fractions_to_density(fractions: &[f64], densities_kg_m3: &[f64]) -> f64 {
    let n = fractions.len().min(densities_kg_m3.len());
    let mut sum = 0.0;
    let mut frac_sum = 0.0;
    for i in 0..n {
        sum += fractions[i] * densities_kg_m3[i];
        frac_sum += fractions[i];
    }
    if frac_sum > 0.0 { sum / frac_sum } else { 0.0 }
}

/// Convert grain size (mm) to estimated fracture toughness scaling factor.
///
/// Hall-Petch-like: K_Ic ∝ d^(-1/2). Returns a dimensionless scale
/// relative to 1mm grain size.
#[must_use]
#[inline]
pub fn grain_size_to_toughness_scale(grain_size_mm: f64) -> f64 {
    if grain_size_mm <= 0.0 {
        return 0.0;
    }
    (1.0 / grain_size_mm).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Sharira ────────────────────────────────────────────────────────

    #[test]
    fn bone_youngs_cortical() {
        // Cortical bone: ~1800 kg/m³ → E ≈ 17-20 GPa
        let e = bone_density_to_youngs_modulus(1800.0);
        assert!(e > 10e9 && e < 30e9, "cortical bone E = {e} Pa");
    }

    #[test]
    fn bone_youngs_zero_density() {
        assert_eq!(bone_density_to_youngs_modulus(0.0), 0.0);
    }

    #[test]
    fn bone_yield_cortical() {
        // Cortical bone: ~1800 kg/m³ → σ_y ≈ 100-500 MPa (Keller power-law)
        let sy = bone_density_to_yield_strength(1800.0);
        assert!(sy > 50e6 && sy < 500e6, "cortical bone σ_y = {sy} Pa");
    }

    #[test]
    fn tendon_stress_basic() {
        // 1000N through 50mm² → 20 MPa
        let s = muscle_force_to_tendon_stress(1000.0, 50e-6);
        assert!((s - 20e6).abs() < 1e6);
    }

    #[test]
    fn tendon_stress_zero_area() {
        assert_eq!(muscle_force_to_tendon_stress(1000.0, 0.0), 0.0);
    }

    #[test]
    fn tendon_force_from_strain() {
        // E = 1.2 GPa (tendon), A = 50mm², ε = 0.03 → F = 1800 N
        let f = tendon_strain_to_force(0.03, 1.2e9, 50e-6);
        assert!((f - 1800.0).abs() < 1.0);
    }

    #[test]
    fn bone_safety_factor_safe() {
        // 10 MPa applied, cortical bone → SF >> 1
        let sf = bone_safety_factor(10e6, 1800.0);
        assert!(sf > 5.0, "should be safe, SF = {sf}");
    }

    #[test]
    fn bone_safety_factor_zero_stress() {
        let sf = bone_safety_factor(0.0, 1800.0);
        assert!(sf.is_infinite());
    }

    #[test]
    fn fracture_volume_basic() {
        let v = impact_to_fracture_volume(100.0, 1e6);
        assert!((v - 1e-4).abs() < 1e-10);
    }

    // ── Impetus ────────────────────────────────────────────────────────

    #[test]
    fn contact_stress_basic() {
        let s = collision_to_contact_stress(1000.0, 0.01);
        assert!((s - 100_000.0).abs() < 0.1);
    }

    #[test]
    fn strain_rate_basic() {
        let sr = velocity_to_strain_rate(10.0, 1.0);
        assert!((sr - 10.0).abs() < 0.001);
    }

    // ── Ushma ──────────────────────────────────────────────────────────

    #[test]
    fn thermal_strain_steel() {
        // Steel α ≈ 12e-6, ΔT = 100K → ε = 0.0012
        let e = temperature_to_thermal_strain(100.0, 12e-6);
        assert!((e - 0.0012).abs() < 1e-6);
    }

    #[test]
    fn thermal_stress_basic() {
        let s = thermal_gradient_to_stress(10.0, 200e9, 12e-6, 1.0);
        assert!((s - 24_000_000.0).abs() < 1.0);
    }

    // ── Bijli ──────────────────────────────────────────────────────────

    #[test]
    fn piezo_stress_basic() {
        let s = e_field_to_piezo_stress(1000.0, 15.0);
        assert!((s - 15000.0).abs() < 0.1);
    }

    #[test]
    fn magnetostriction_iron() {
        let e = magnetic_to_magnetostrictive_strain(1.0, 35e-6, 2.15);
        let expected = 1.5 * 35e-6 * (1.0 / 2.15_f64).powi(2);
        assert!((e - expected).abs() < 1e-10);
    }

    #[test]
    fn magnetostriction_zero_saturation() {
        assert_eq!(magnetic_to_magnetostrictive_strain(1.0, 35e-6, 0.0), 0.0);
    }

    // ── Khanij ─────────────────────────────────────────────────────────

    #[test]
    fn mineral_density_quartz_feldspar() {
        // 60% quartz (2650), 40% feldspar (2600) → ~2630
        let rho = mineral_fractions_to_density(&[0.6, 0.4], &[2650.0, 2600.0]);
        assert!((rho - 2630.0).abs() < 1.0);
    }

    #[test]
    fn mineral_density_empty() {
        assert_eq!(mineral_fractions_to_density(&[], &[]), 0.0);
    }

    #[test]
    fn grain_size_toughness_1mm() {
        let s = grain_size_to_toughness_scale(1.0);
        assert!((s - 1.0).abs() < 0.01);
    }

    #[test]
    fn grain_size_toughness_finer_is_tougher() {
        let fine = grain_size_to_toughness_scale(0.1);
        let coarse = grain_size_to_toughness_scale(1.0);
        assert!(fine > coarse);
    }

    #[test]
    fn grain_size_zero() {
        assert_eq!(grain_size_to_toughness_scale(0.0), 0.0);
    }
}
