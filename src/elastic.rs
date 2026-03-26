//! Elastic constants and relationships: Hooke's law, moduli conversions, plane stress/strain.

/// Hooke's law: σ = E × ε
#[must_use]
#[inline]
pub fn hookes_law(youngs_modulus: f64, strain: f64) -> f64 {
    youngs_modulus * strain
}

/// Bulk modulus: K = E / (3(1 - 2v))
#[must_use]
#[inline]
pub fn bulk_modulus(youngs_modulus: f64, poisson_ratio: f64) -> f64 {
    let denom = 3.0 * (1.0 - 2.0 * poisson_ratio);
    if denom.abs() < hisab::EPSILON_F64 {
        tracing::warn!("bulk_modulus: degenerate v={poisson_ratio}, returning 0.0");
        return 0.0;
    }
    youngs_modulus / denom
}

/// Shear modulus: G = E / (2(1 + v))
#[must_use]
#[inline]
pub fn shear_modulus(youngs_modulus: f64, poisson_ratio: f64) -> f64 {
    let denom = 2.0 * (1.0 + poisson_ratio);
    if denom.abs() < hisab::EPSILON_F64 {
        tracing::warn!("shear_modulus: degenerate v={poisson_ratio}, returning 0.0");
        return 0.0;
    }
    youngs_modulus / denom
}

/// First Lame parameter: λ = Ev / ((1+v)(1-2v))
#[must_use]
#[inline]
pub fn lame_lambda(youngs_modulus: f64, poisson_ratio: f64) -> f64 {
    let denom = (1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio);
    if denom.abs() < hisab::EPSILON_F64 {
        tracing::warn!("lame_lambda: degenerate v={poisson_ratio}, returning 0.0");
        return 0.0;
    }
    youngs_modulus * poisson_ratio / denom
}

/// Strain from stress (inverse Hooke): ε = σ / E
#[must_use]
#[inline]
pub fn strain_from_stress(youngs_modulus: f64, stress: f64) -> f64 {
    if youngs_modulus.abs() < hisab::EPSILON_F64 {
        tracing::warn!("strain_from_stress: zero modulus, returning 0.0");
        return 0.0;
    }
    stress / youngs_modulus
}

/// Young's modulus from shear modulus and Poisson's ratio: E = 2G(1 + v)
#[must_use]
#[inline]
pub fn youngs_from_shear(shear_mod: f64, poisson_ratio: f64) -> f64 {
    2.0 * shear_mod * (1.0 + poisson_ratio)
}

/// Young's modulus from bulk and shear modulus: E = 9KG / (3K + G)
#[must_use]
#[inline]
pub fn youngs_from_bulk_shear(bulk_mod: f64, shear_mod: f64) -> f64 {
    let denom = 3.0 * bulk_mod + shear_mod;
    if denom.abs() < hisab::EPSILON_F64 {
        return 0.0;
    }
    9.0 * bulk_mod * shear_mod / denom
}

/// Poisson's ratio from bulk and shear modulus: v = (3K - 2G) / (2(3K + G))
#[must_use]
#[inline]
pub fn poisson_from_bulk_shear(bulk_mod: f64, shear_mod: f64) -> f64 {
    let denom = 2.0 * (3.0 * bulk_mod + shear_mod);
    if denom.abs() < hisab::EPSILON_F64 {
        return 0.0;
    }
    (3.0 * bulk_mod - 2.0 * shear_mod) / denom
}

/// Plane stress reduced modulus: E* = E / (1 - v^2)
#[must_use]
#[inline]
pub fn plane_stress_modulus(youngs_modulus: f64, poisson_ratio: f64) -> f64 {
    let denom = 1.0 - poisson_ratio * poisson_ratio;
    if denom.abs() < hisab::EPSILON_F64 {
        return 0.0;
    }
    youngs_modulus / denom
}

/// Plane strain reduced modulus: E* = E(1 - v) / ((1 + v)(1 - 2v))
#[must_use]
#[inline]
pub fn plane_strain_modulus(youngs_modulus: f64, poisson_ratio: f64) -> f64 {
    let denom = (1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio);
    if denom.abs() < hisab::EPSILON_F64 {
        return 0.0;
    }
    youngs_modulus * (1.0 - poisson_ratio) / denom
}

/// P-wave modulus: M = E(1-v) / ((1+v)(1-2v)) = K + 4G/3
///
/// Speed of longitudinal waves: c_p = sqrt(M / rho)
#[must_use]
#[inline]
pub fn p_wave_modulus(youngs_modulus: f64, poisson_ratio: f64) -> f64 {
    plane_strain_modulus(youngs_modulus, poisson_ratio)
}

/// Poisson's ratio from Young's modulus and shear modulus: v = E/(2G) - 1
#[must_use]
#[inline]
pub fn poisson_from_youngs_shear(youngs_modulus: f64, shear_mod: f64) -> f64 {
    if shear_mod.abs() < hisab::EPSILON_F64 {
        return 0.0;
    }
    youngs_modulus / (2.0 * shear_mod) - 1.0
}

// --- Checked (Result-returning) variants ---

/// Checked bulk modulus. Returns `Err` for incompressible (v=0.5) or degenerate inputs.
pub fn try_bulk_modulus(youngs_modulus: f64, poisson_ratio: f64) -> crate::Result<f64> {
    let denom = 3.0 * (1.0 - 2.0 * poisson_ratio);
    if denom.abs() < hisab::EPSILON_F64 {
        return Err(crate::DravyaError::DivisionByZero(
            "bulk_modulus: v=0.5 (incompressible)",
        ));
    }
    Ok(youngs_modulus / denom)
}

/// Checked shear modulus. Returns `Err` for v=-1.0 or degenerate inputs.
pub fn try_shear_modulus(youngs_modulus: f64, poisson_ratio: f64) -> crate::Result<f64> {
    let denom = 2.0 * (1.0 + poisson_ratio);
    if denom.abs() < hisab::EPSILON_F64 {
        return Err(crate::DravyaError::DivisionByZero("shear_modulus: v=-1.0"));
    }
    Ok(youngs_modulus / denom)
}

/// Checked strain from stress. Returns `Err` for zero modulus.
pub fn try_strain_from_stress(youngs_modulus: f64, stress: f64) -> crate::Result<f64> {
    if youngs_modulus.abs() < hisab::EPSILON_F64 {
        return Err(crate::DravyaError::DivisionByZero(
            "strain_from_stress: zero modulus",
        ));
    }
    Ok(stress / youngs_modulus)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hookes_law_steel() {
        // 200 GPa × 0.001 = 200 MPa
        let stress = hookes_law(200e9, 0.001);
        assert!(
            (stress - 200e6).abs() < 1.0,
            "stress should be 200 MPa, got {stress}"
        );
    }

    #[test]
    fn bulk_modulus_steel() {
        let k = bulk_modulus(200e9, 0.30);
        assert!(
            (k - 166.7e9).abs() < 1e9,
            "steel bulk modulus should be ~167 GPa, got {k}"
        );
    }

    #[test]
    fn shear_modulus_steel() {
        let g = shear_modulus(200e9, 0.30);
        assert!(
            (g - 76.9e9).abs() < 0.5e9,
            "steel shear modulus should be ~77 GPa, got {g}"
        );
    }

    #[test]
    fn shear_modulus_guard() {
        // v = -1.0 -> denominator = 0
        let g = shear_modulus(200e9, -1.0);
        assert_eq!(g, 0.0, "should guard against v = -1.0");
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
        let k = bulk_modulus(200e9, 0.5);
        assert_eq!(k, 0.0);
    }

    #[test]
    fn youngs_from_shear_steel() {
        let e = youngs_from_shear(76.9e9, 0.30);
        assert!(
            (e - 200e9).abs() < 1e9,
            "should recover E ~ 200 GPa, got {e}"
        );
    }

    #[test]
    fn youngs_from_bulk_shear_steel() {
        let k = bulk_modulus(200e9, 0.30);
        let g = shear_modulus(200e9, 0.30);
        let e = youngs_from_bulk_shear(k, g);
        assert!(
            (e - 200e9).abs() < 1e6,
            "should recover E = 200 GPa, got {e}"
        );
    }

    #[test]
    fn poisson_from_bulk_shear_steel() {
        let k = bulk_modulus(200e9, 0.30);
        let g = shear_modulus(200e9, 0.30);
        let v = poisson_from_bulk_shear(k, g);
        assert!((v - 0.30).abs() < 1e-6, "should recover v = 0.30, got {v}");
    }

    #[test]
    fn plane_stress_modulus_steel() {
        let es = plane_stress_modulus(200e9, 0.30);
        // E / (1 - 0.09) = 200 / 0.91 ~ 219.8 GPa
        assert!(
            (es - 219.8e9).abs() < 1e9,
            "plane stress modulus ~219.8 GPa, got {es}"
        );
    }

    #[test]
    fn poisson_from_youngs_shear_steel() {
        let v = poisson_from_youngs_shear(200e9, 76.923e9);
        assert!((v - 0.30).abs() < 0.01, "should recover v ~0.30, got {v}");
    }

    #[test]
    fn try_bulk_modulus_incompressible() {
        let result = try_bulk_modulus(200e9, 0.5);
        assert!(result.is_err(), "v=0.5 should return Err");
    }

    #[test]
    fn try_bulk_modulus_valid() {
        let result = try_bulk_modulus(200e9, 0.30);
        assert!(result.is_ok());
        assert!((result.unwrap() - 166.7e9).abs() < 1e9);
    }

    #[test]
    fn try_shear_modulus_degenerate() {
        let result = try_shear_modulus(200e9, -1.0);
        assert!(result.is_err());
    }

    #[test]
    fn try_strain_from_stress_zero_modulus() {
        let result = try_strain_from_stress(0.0, 100e6);
        assert!(result.is_err());
    }

    #[test]
    fn plane_strain_modulus_steel() {
        let ep = plane_strain_modulus(200e9, 0.30);
        // E(1-v) / ((1+v)(1-2v)) = 200*0.7 / (1.3*0.4) = 140/0.52 ~ 269.2 GPa
        assert!(
            (ep - 269.2e9).abs() < 1e9,
            "plane strain modulus ~269.2 GPa, got {ep}"
        );
    }
}
