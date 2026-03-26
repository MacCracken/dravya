use std::f64::consts::PI;

// --- Deflection formulas ---

/// Cantilever beam deflection at free end (point load): δ = FL^3 / (3EI)
#[must_use]
#[inline]
pub fn cantilever_deflection(
    force: f64,
    length: f64,
    youngs_modulus: f64,
    moment_of_inertia: f64,
) -> f64 {
    let denom = 3.0 * youngs_modulus * moment_of_inertia;
    if denom.abs() < hisab::EPSILON_F64 {
        return 0.0;
    }
    force * length.powi(3) / denom
}

/// Cantilever beam deflection at free end (uniform distributed load):
/// δ = wL^4 / (8EI)
#[must_use]
#[inline]
pub fn cantilever_deflection_udl(
    load_per_length: f64,
    length: f64,
    youngs_modulus: f64,
    moment_of_inertia: f64,
) -> f64 {
    let denom = 8.0 * youngs_modulus * moment_of_inertia;
    if denom.abs() < hisab::EPSILON_F64 {
        return 0.0;
    }
    load_per_length * length.powi(4) / denom
}

/// Simply supported beam deflection at midspan (point load):
/// δ = FL^3 / (48EI)
#[must_use]
#[inline]
pub fn simply_supported_deflection(
    force: f64,
    length: f64,
    youngs_modulus: f64,
    moment_of_inertia: f64,
) -> f64 {
    let denom = 48.0 * youngs_modulus * moment_of_inertia;
    if denom.abs() < hisab::EPSILON_F64 {
        return 0.0;
    }
    force * length.powi(3) / denom
}

/// Simply supported beam deflection at midspan (uniform distributed load):
/// δ = 5wL^4 / (384EI)
#[must_use]
#[inline]
pub fn simply_supported_deflection_udl(
    load_per_length: f64,
    length: f64,
    youngs_modulus: f64,
    moment_of_inertia: f64,
) -> f64 {
    let denom = 384.0 * youngs_modulus * moment_of_inertia;
    if denom.abs() < hisab::EPSILON_F64 {
        return 0.0;
    }
    5.0 * load_per_length * length.powi(4) / denom
}

/// Fixed-fixed beam deflection at midspan (point load):
/// δ = FL^3 / (192EI)
#[must_use]
#[inline]
pub fn fixed_fixed_deflection(
    force: f64,
    length: f64,
    youngs_modulus: f64,
    moment_of_inertia: f64,
) -> f64 {
    let denom = 192.0 * youngs_modulus * moment_of_inertia;
    if denom.abs() < hisab::EPSILON_F64 {
        return 0.0;
    }
    force * length.powi(3) / denom
}

// --- Stress ---

/// Bending stress: σ = My / I
#[must_use]
#[inline]
pub fn bending_stress(moment: f64, distance_from_neutral: f64, moment_of_inertia: f64) -> f64 {
    if moment_of_inertia.abs() < hisab::EPSILON_F64 {
        return 0.0;
    }
    moment * distance_from_neutral / moment_of_inertia
}

/// Transverse shear stress in a beam: τ = VQ / (Ib)
///
/// V = shear force, Q = first moment of area above the point,
/// I = moment of inertia, b = width at the point.
#[must_use]
#[inline]
pub fn shear_stress_beam(
    shear_force: f64,
    first_moment: f64,
    moment_of_inertia: f64,
    width: f64,
) -> f64 {
    let denom = moment_of_inertia * width;
    if denom.abs() < hisab::EPSILON_F64 {
        return 0.0;
    }
    shear_force * first_moment / denom
}

/// Torsional shear stress in a circular shaft: τ = Tr / J
#[must_use]
#[inline]
pub fn torsional_stress(torque: f64, radius: f64, polar_moment: f64) -> f64 {
    if polar_moment.abs() < hisab::EPSILON_F64 {
        return 0.0;
    }
    torque * radius / polar_moment
}

/// Angle of twist: θ = TL / (GJ) (radians)
#[must_use]
#[inline]
pub fn angle_of_twist(torque: f64, length: f64, shear_modulus: f64, polar_moment: f64) -> f64 {
    let denom = shear_modulus * polar_moment;
    if denom.abs() < hisab::EPSILON_F64 {
        return 0.0;
    }
    torque * length / denom
}

// --- Euler buckling ---

/// Euler critical buckling load: P_cr = π²EI / L_eff²
///
/// `effective_length` accounts for end conditions:
/// - Pinned-pinned: L_eff = L
/// - Fixed-free: L_eff = 2L
/// - Fixed-pinned: L_eff = 0.7L
/// - Fixed-fixed: L_eff = 0.5L
#[must_use]
#[inline]
pub fn euler_buckling_load(
    youngs_modulus: f64,
    moment_of_inertia: f64,
    effective_length: f64,
) -> f64 {
    let denom = effective_length * effective_length;
    if denom.abs() < hisab::EPSILON_F64 {
        return 0.0;
    }
    PI * PI * youngs_modulus * moment_of_inertia / denom
}

// --- Cross-section properties ---

/// Moment of inertia for a rectangular cross-section: I = bh^3/12
#[must_use]
#[inline]
pub fn moment_of_inertia_rect(width: f64, height: f64) -> f64 {
    width * height.powi(3) / 12.0
}

/// Moment of inertia for a circular cross-section: I = πr^4/4
#[must_use]
#[inline]
pub fn moment_of_inertia_circle(radius: f64) -> f64 {
    PI * radius.powi(4) / 4.0
}

/// Moment of inertia for a hollow circular section:
/// I = π/4 * (R_outer^4 - R_inner^4)
#[must_use]
#[inline]
pub fn moment_of_inertia_hollow_circle(outer_radius: f64, inner_radius: f64) -> f64 {
    PI / 4.0 * (outer_radius.powi(4) - inner_radius.powi(4))
}

/// Moment of inertia for a hollow rectangular section:
/// I = (b_o * h_o^3 - b_i * h_i^3) / 12
#[must_use]
#[inline]
pub fn moment_of_inertia_hollow_rect(
    outer_width: f64,
    outer_height: f64,
    inner_width: f64,
    inner_height: f64,
) -> f64 {
    (outer_width * outer_height.powi(3) - inner_width * inner_height.powi(3)) / 12.0
}

/// Polar moment of inertia for a solid circular section: J = πr^4/2
///
/// Used for torsion analysis.
#[must_use]
#[inline]
pub fn polar_moment_circle(radius: f64) -> f64 {
    PI * radius.powi(4) / 2.0
}

/// Polar moment of inertia for a hollow circular section:
/// J = π/2 * (R_outer^4 - R_inner^4)
#[must_use]
#[inline]
pub fn polar_moment_hollow_circle(outer_radius: f64, inner_radius: f64) -> f64 {
    PI / 2.0 * (outer_radius.powi(4) - inner_radius.powi(4))
}

/// Section modulus for rectangular section: S = bh^2/6
#[must_use]
#[inline]
pub fn section_modulus_rect(width: f64, height: f64) -> f64 {
    width * height.powi(2) / 6.0
}

/// Section modulus for circular section: S = πr^3/4
#[must_use]
#[inline]
pub fn section_modulus_circle(radius: f64) -> f64 {
    PI * radius.powi(3) / 4.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cantilever_basic() {
        let i = moment_of_inertia_rect(0.1, 0.01);
        let delta = cantilever_deflection(1000.0, 1.0, 200e9, i);
        assert!(delta > 0.0, "deflection should be positive");
        assert!(delta < 1.0, "deflection should be finite, got {delta}");
    }

    #[test]
    fn simply_supported_less_than_cantilever() {
        let i = moment_of_inertia_rect(0.1, 0.01);
        let cant = cantilever_deflection(1000.0, 1.0, 200e9, i);
        let ss = simply_supported_deflection(1000.0, 1.0, 200e9, i);
        assert!(
            ss < cant,
            "simply supported should deflect less than cantilever"
        );
    }

    #[test]
    fn fixed_fixed_less_than_simply_supported() {
        let i = moment_of_inertia_rect(0.1, 0.01);
        let ss = simply_supported_deflection(1000.0, 1.0, 200e9, i);
        let ff = fixed_fixed_deflection(1000.0, 1.0, 200e9, i);
        assert!(
            ff < ss,
            "fixed-fixed should deflect less than simply supported"
        );
    }

    #[test]
    fn udl_deflections() {
        let i = moment_of_inertia_rect(0.1, 0.01);
        let cant_udl = cantilever_deflection_udl(1000.0, 1.0, 200e9, i);
        let ss_udl = simply_supported_deflection_udl(1000.0, 1.0, 200e9, i);
        assert!(cant_udl > 0.0);
        assert!(ss_udl > 0.0);
        assert!(
            ss_udl < cant_udl,
            "SS UDL should deflect less than cantilever UDL"
        );
    }

    #[test]
    fn moment_of_inertia_rect_100x10mm() {
        let i = moment_of_inertia_rect(0.1, 0.01);
        assert!(
            (i - 8.333e-9).abs() < 1e-11,
            "I should be ~8.33e-9, got {i}"
        );
    }

    #[test]
    fn moment_of_inertia_circle_10mm() {
        let i = moment_of_inertia_circle(0.005);
        assert!(i > 0.0 && i < 1e-6);
    }

    #[test]
    fn hollow_circle_less_solid() {
        let solid = moment_of_inertia_circle(0.05);
        let hollow = moment_of_inertia_hollow_circle(0.05, 0.04);
        assert!(hollow < solid, "hollow should have less I than solid");
        assert!(hollow > 0.0);
    }

    #[test]
    fn hollow_rect_less_solid() {
        let solid = moment_of_inertia_rect(0.1, 0.1);
        let hollow = moment_of_inertia_hollow_rect(0.1, 0.1, 0.08, 0.08);
        assert!(hollow < solid);
        assert!(hollow > 0.0);
    }

    #[test]
    fn polar_moment_is_2x_area_moment() {
        let r = 0.05;
        let j = polar_moment_circle(r);
        let i = moment_of_inertia_circle(r);
        assert!(
            (j - 2.0 * i).abs() < hisab::EPSILON_F64,
            "J = 2I for circular section"
        );
    }

    #[test]
    fn bending_stress_basic() {
        let i = moment_of_inertia_rect(0.1, 0.01);
        let sigma = bending_stress(100.0, 0.005, i);
        assert!(sigma > 0.0);
    }

    #[test]
    fn torsional_stress_basic() {
        let j = polar_moment_circle(0.025);
        let tau = torsional_stress(500.0, 0.025, j);
        assert!(tau > 0.0);
    }

    #[test]
    fn euler_buckling_basic() {
        let i = moment_of_inertia_rect(0.05, 0.05);
        let p_cr = euler_buckling_load(200e9, i, 2.0);
        assert!(p_cr > 0.0, "critical load should be positive");
        // Shorter column -> higher buckling load
        let p_cr_short = euler_buckling_load(200e9, i, 1.0);
        assert!(p_cr_short > p_cr, "shorter column should have higher P_cr");
    }

    #[test]
    fn zero_inertia_safe() {
        assert_eq!(cantilever_deflection(100.0, 1.0, 200e9, 0.0), 0.0);
        assert_eq!(bending_stress(100.0, 0.005, 0.0), 0.0);
        assert_eq!(torsional_stress(100.0, 0.025, 0.0), 0.0);
    }

    #[test]
    fn section_modulus_circle_basic() {
        let s = section_modulus_circle(0.025);
        assert!(s > 0.0);
    }
}
