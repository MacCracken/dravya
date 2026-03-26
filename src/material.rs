use serde::{Deserialize, Serialize};

/// Engineering material with mechanical properties.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Material {
    pub name: String,
    /// Young's modulus (Pa).
    pub youngs_modulus: f64,
    /// Poisson's ratio (dimensionless, typically 0.2–0.5).
    pub poisson_ratio: f64,
    /// Yield strength (Pa).
    pub yield_strength: f64,
    /// Density (kg/m³).
    pub density: f64,
    /// Coefficient of thermal expansion (1/K).
    pub thermal_expansion: f64,
}

impl Material {
    #[must_use] pub fn steel() -> Self {
        Self { name: "Steel".into(), youngs_modulus: 200e9, poisson_ratio: 0.30, yield_strength: 250e6, density: 7850.0, thermal_expansion: 12e-6 }
    }
    #[must_use] pub fn aluminum() -> Self {
        Self { name: "Aluminum".into(), youngs_modulus: 69e9, poisson_ratio: 0.33, yield_strength: 276e6, density: 2700.0, thermal_expansion: 23e-6 }
    }
    #[must_use] pub fn copper() -> Self {
        Self { name: "Copper".into(), youngs_modulus: 117e9, poisson_ratio: 0.34, yield_strength: 210e6, density: 8960.0, thermal_expansion: 17e-6 }
    }
    #[must_use] pub fn titanium() -> Self {
        Self { name: "Titanium".into(), youngs_modulus: 116e9, poisson_ratio: 0.32, yield_strength: 880e6, density: 4507.0, thermal_expansion: 8.6e-6 }
    }
    #[must_use] pub fn glass() -> Self {
        Self { name: "Glass".into(), youngs_modulus: 70e9, poisson_ratio: 0.22, yield_strength: 33e6, density: 2500.0, thermal_expansion: 9e-6 }
    }
    #[must_use] pub fn rubber() -> Self {
        Self { name: "Rubber".into(), youngs_modulus: 0.01e9, poisson_ratio: 0.49, yield_strength: 15e6, density: 1100.0, thermal_expansion: 200e-6 }
    }
    #[must_use] pub fn concrete() -> Self {
        Self { name: "Concrete".into(), youngs_modulus: 30e9, poisson_ratio: 0.20, yield_strength: 3e6, density: 2400.0, thermal_expansion: 12e-6 }
    }
    #[must_use] pub fn wood_oak() -> Self {
        Self { name: "Oak".into(), youngs_modulus: 12e9, poisson_ratio: 0.35, yield_strength: 40e6, density: 600.0, thermal_expansion: 5e-6 }
    }
    #[must_use] pub fn carbon_fiber() -> Self {
        Self { name: "Carbon Fiber".into(), youngs_modulus: 181e9, poisson_ratio: 0.27, yield_strength: 3500e6, density: 1600.0, thermal_expansion: -0.5e-6 }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn steel_properties() {
        let s = Material::steel();
        assert!((s.youngs_modulus - 200e9).abs() < 1e6);
        assert!((s.poisson_ratio - 0.30).abs() < 0.01);
        assert!((s.density - 7850.0).abs() < 1.0);
    }

    #[test]
    fn carbon_fiber_negative_thermal_expansion() {
        let cf = Material::carbon_fiber();
        assert!(cf.thermal_expansion < 0.0, "carbon fiber has negative CTE");
    }

    #[test]
    fn rubber_nearly_incompressible() {
        let r = Material::rubber();
        assert!(r.poisson_ratio > 0.48, "rubber should be nearly incompressible");
    }

    #[test]
    fn all_presets_positive_modulus() {
        let mats = [Material::steel(), Material::aluminum(), Material::copper(), Material::titanium(),
                     Material::glass(), Material::rubber(), Material::concrete(), Material::wood_oak(), Material::carbon_fiber()];
        for m in &mats {
            assert!(m.youngs_modulus > 0.0, "{} should have positive E", m.name);
            assert!(m.density > 0.0, "{} should have positive density", m.name);
        }
    }

    #[test]
    fn serde_roundtrip() {
        let m = Material::steel();
        let json = serde_json::to_string(&m).unwrap();
        let back: Material = serde_json::from_str(&json).unwrap();
        assert!((back.youngs_modulus - m.youngs_modulus).abs() < 1.0);
    }
}
