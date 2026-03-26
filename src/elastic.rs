/// Hooke's law: σ = E × ε
#[must_use]
#[inline]
pub fn hookes_law(youngs_modulus: f64, strain: f64) -> f64 {
    youngs_modulus * strain
}

/// Bulk modulus: K = E / (3(1 - 2ν))
#[must_use]
#[inline]
pub fn bulk_modulus(youngs_modulus: f64, poisson_ratio: f64) -> f64 {
    let denom = 3.0 * (1.0 - 2.0 * poisson_ratio);
    if denom.abs() < f64::EPSILON { return 0.0; }
    youngs_modulus / denom
}

/// Shear modulus: G = E / (2(1 + ν))
#[must_use]
#[inline]
pub fn shear_modulus(youngs_modulus: f64, poisson_ratio: f64) -> f64 {
    youngs_modulus / (2.0 * (1.0 + poisson_ratio))
}

/// First Lamé parameter: λ = Eν / ((1+ν)(1-2ν))
#[must_use]
#[inline]
pub fn lame_lambda(youngs_modulus: f64, poisson_ratio: f64) -> f64 {
    let denom = (1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio);
    if denom.abs() < f64::EPSILON { return 0.0; }
    youngs_modulus * poisson_ratio / denom
}

/// Strain from stress (inverse Hooke): ε = σ / E
#[must_use]
#[inline]
pub fn strain_from_stress(youngs_modulus: f64, stress: f64) -> f64 {
    if youngs_modulus.abs() < f64::EPSILON { return 0.0; }
    stress / youngs_modulus
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hookes_law_steel() {
        // 200 GPa × 0.001 = 200 MPa
        let stress = hookes_law(200e9, 0.001);
        assert!((stress - 200e6).abs() < 1.0, "stress should be 200 MPa, got {stress}");
    }

    #[test]
    fn bulk_modulus_steel() {
        // K = 200e9 / (3 × 0.4) ≈ 166.7 GPa
        let k = bulk_modulus(200e9, 0.30);
        assert!((k - 166.7e9).abs() < 1e9, "steel bulk modulus should be ~167 GPa, got {k}");
    }

    #[test]
    fn shear_modulus_steel() {
        // G = 200e9 / (2 × 1.3) ≈ 76.9 GPa
        let g = shear_modulus(200e9, 0.30);
        assert!((g - 76.9e9).abs() < 0.5e9, "steel shear modulus should be ~77 GPa, got {g}");
    }

    #[test]
    fn lame_lambda_steel() {
        let lambda = lame_lambda(200e9, 0.30);
        assert!(lambda > 0.0);
    }

    #[test]
    fn strain_roundtrip() {
        let e = 200e9;
        let stress = hookes_law(e, 0.002);
        let strain_back = strain_from_stress(e, stress);
        assert!((strain_back - 0.002).abs() < 1e-12);
    }

    #[test]
    fn incompressible_bulk_modulus() {
        // ν = 0.5 → denominator = 0 → should return 0 (guard)
        let k = bulk_modulus(200e9, 0.5);
        assert_eq!(k, 0.0);
    }
}
