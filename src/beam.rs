use std::f64::consts::PI;

/// Cantilever beam deflection at free end: δ = FL³ / (3EI)
#[must_use]
#[inline]
pub fn cantilever_deflection(force: f64, length: f64, youngs_modulus: f64, moment_of_inertia: f64) -> f64 {
    let denom = 3.0 * youngs_modulus * moment_of_inertia;
    if denom.abs() < f64::EPSILON { return 0.0; }
    force * length.powi(3) / denom
}

/// Simply supported beam deflection at midspan: δ = FL³ / (48EI)
#[must_use]
#[inline]
pub fn simply_supported_deflection(force: f64, length: f64, youngs_modulus: f64, moment_of_inertia: f64) -> f64 {
    let denom = 48.0 * youngs_modulus * moment_of_inertia;
    if denom.abs() < f64::EPSILON { return 0.0; }
    force * length.powi(3) / denom
}

/// Bending stress: σ = My / I
#[must_use]
#[inline]
pub fn bending_stress(moment: f64, distance_from_neutral: f64, moment_of_inertia: f64) -> f64 {
    if moment_of_inertia.abs() < f64::EPSILON { return 0.0; }
    moment * distance_from_neutral / moment_of_inertia
}

/// Moment of inertia for a rectangular cross-section: I = bh³/12
#[must_use]
#[inline]
pub fn moment_of_inertia_rect(width: f64, height: f64) -> f64 {
    width * height.powi(3) / 12.0
}

/// Moment of inertia for a circular cross-section: I = πr⁴/4
#[must_use]
#[inline]
pub fn moment_of_inertia_circle(radius: f64) -> f64 {
    PI * radius.powi(4) / 4.0
}

/// Section modulus for rectangular section: S = bh²/6
#[must_use]
#[inline]
pub fn section_modulus_rect(width: f64, height: f64) -> f64 {
    width * height.powi(2) / 6.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cantilever_basic() {
        // F=1000N, L=1m, E=200GPa, I=8.33e-6 (100mm×10mm rect)
        let i = moment_of_inertia_rect(0.1, 0.01);
        let delta = cantilever_deflection(1000.0, 1.0, 200e9, i);
        assert!(delta > 0.0, "deflection should be positive");
        // 1kN on 100×10mm steel cantilever, 1m long: δ = FL³/3EI ≈ 0.2m (thin beam, large deflection)
        assert!(delta < 1.0, "deflection should be finite, got {delta}");
    }

    #[test]
    fn simply_supported_less_than_cantilever() {
        let i = moment_of_inertia_rect(0.1, 0.01);
        let cant = cantilever_deflection(1000.0, 1.0, 200e9, i);
        let ss = simply_supported_deflection(1000.0, 1.0, 200e9, i);
        assert!(ss < cant, "simply supported should deflect less than cantilever");
    }

    #[test]
    fn moment_of_inertia_rect_100x10mm() {
        let i = moment_of_inertia_rect(0.1, 0.01); // 100mm × 10mm
        // I = 0.1 × 0.01³ / 12 = 8.333e-9 m⁴
        assert!((i - 8.333e-9).abs() < 1e-11, "I should be ~8.33e-9, got {i}");
    }

    #[test]
    fn moment_of_inertia_circle_10mm() {
        let i = moment_of_inertia_circle(0.005); // 10mm diameter
        assert!(i > 0.0 && i < 1e-6);
    }

    #[test]
    fn bending_stress_basic() {
        let i = moment_of_inertia_rect(0.1, 0.01);
        let sigma = bending_stress(100.0, 0.005, i); // 100 N·m, 5mm from neutral axis
        assert!(sigma > 0.0);
    }

    #[test]
    fn zero_inertia_safe() {
        assert_eq!(cantilever_deflection(100.0, 1.0, 200e9, 0.0), 0.0);
        assert_eq!(bending_stress(100.0, 0.005, 0.0), 0.0);
    }
}
