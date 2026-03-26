//! Engineering material definitions with mechanical, thermal, and elastic properties.

use std::borrow::Cow;
use std::fmt;

use serde::{Deserialize, Serialize};

use crate::elastic;

/// Engineering material with mechanical properties.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Material {
    /// Material name / identifier.
    pub name: Cow<'static, str>,
    /// Young's modulus (Pa).
    pub youngs_modulus: f64,
    /// Poisson's ratio (dimensionless, typically 0.2-0.5).
    pub poisson_ratio: f64,
    /// Yield strength (Pa).
    pub yield_strength: f64,
    /// Ultimate tensile strength (Pa).
    pub ultimate_tensile_strength: f64,
    /// Density (kg/m^3).
    pub density: f64,
    /// Coefficient of thermal expansion (1/K).
    pub thermal_expansion: f64,
}

impl Material {
    /// Create a validated material.
    ///
    /// Checks: E > 0, -1 < v < 0.5, sigma_y >= 0, UTS >= sigma_y, rho > 0.
    pub fn new(
        name: impl Into<Cow<'static, str>>,
        youngs_modulus: f64,
        poisson_ratio: f64,
        yield_strength: f64,
        ultimate_tensile_strength: f64,
        density: f64,
        thermal_expansion: f64,
    ) -> crate::Result<Self> {
        if youngs_modulus <= 0.0 {
            return Err(crate::DravyaError::InvalidParameter {
                name: "youngs_modulus",
                value: youngs_modulus,
                reason: "must be positive",
            });
        }
        if poisson_ratio <= -1.0 || poisson_ratio >= 0.5 {
            return Err(crate::DravyaError::InvalidParameter {
                name: "poisson_ratio",
                value: poisson_ratio,
                reason: "must be in (-1, 0.5)",
            });
        }
        if yield_strength < 0.0 {
            return Err(crate::DravyaError::InvalidParameter {
                name: "yield_strength",
                value: yield_strength,
                reason: "must be non-negative",
            });
        }
        if ultimate_tensile_strength < yield_strength {
            return Err(crate::DravyaError::InvalidParameter {
                name: "ultimate_tensile_strength",
                value: ultimate_tensile_strength,
                reason: "must be >= yield_strength",
            });
        }
        if density <= 0.0 {
            return Err(crate::DravyaError::InvalidParameter {
                name: "density",
                value: density,
                reason: "must be positive",
            });
        }
        Ok(Self {
            name: name.into(),
            youngs_modulus,
            poisson_ratio,
            yield_strength,
            ultimate_tensile_strength,
            density,
            thermal_expansion,
        })
    }

    /// Shear modulus G = E / (2(1 + v)).
    #[must_use]
    #[inline]
    pub fn shear_modulus(&self) -> f64 {
        elastic::shear_modulus(self.youngs_modulus, self.poisson_ratio)
    }

    /// Bulk modulus K = E / (3(1 - 2v)).
    #[must_use]
    #[inline]
    pub fn bulk_modulus(&self) -> f64 {
        elastic::bulk_modulus(self.youngs_modulus, self.poisson_ratio)
    }

    /// First Lame parameter.
    #[must_use]
    #[inline]
    pub fn lame_lambda(&self) -> f64 {
        elastic::lame_lambda(self.youngs_modulus, self.poisson_ratio)
    }

    /// Thermal strain for a given temperature change: e_th = alpha * dT.
    #[must_use]
    #[inline]
    pub fn thermal_strain(&self, delta_t: f64) -> f64 {
        self.thermal_expansion * delta_t
    }

    /// Thermal stress under full constraint: sigma = E * alpha * dT.
    #[must_use]
    #[inline]
    pub fn thermal_stress(&self, delta_t: f64) -> f64 {
        self.youngs_modulus * self.thermal_expansion * delta_t
    }

    /// Structural steel (ASTM A36 / S275).
    #[must_use]
    pub fn steel() -> Self {
        Self {
            name: Cow::Borrowed("Steel"),
            youngs_modulus: 200e9,
            poisson_ratio: 0.30,
            yield_strength: 250e6,
            ultimate_tensile_strength: 400e6,
            density: 7850.0,
            thermal_expansion: 12e-6,
        }
    }

    /// Aluminum 6061-T6.
    #[must_use]
    pub fn aluminum() -> Self {
        Self {
            name: Cow::Borrowed("Aluminum 6061-T6"),
            youngs_modulus: 69e9,
            poisson_ratio: 0.33,
            yield_strength: 276e6,
            ultimate_tensile_strength: 310e6,
            density: 2700.0,
            thermal_expansion: 23.6e-6,
        }
    }

    /// Copper (annealed, C11000).
    #[must_use]
    pub fn copper() -> Self {
        Self {
            name: Cow::Borrowed("Copper"),
            youngs_modulus: 117e9,
            poisson_ratio: 0.34,
            yield_strength: 62e6,
            ultimate_tensile_strength: 210e6,
            density: 8960.0,
            thermal_expansion: 17e-6,
        }
    }

    /// Titanium Ti-6Al-4V (annealed).
    #[must_use]
    pub fn titanium() -> Self {
        Self {
            name: Cow::Borrowed("Titanium Ti-6Al-4V"),
            youngs_modulus: 114e9,
            poisson_ratio: 0.33,
            yield_strength: 880e6,
            ultimate_tensile_strength: 950e6,
            density: 4430.0,
            thermal_expansion: 8.6e-6,
        }
    }

    /// Soda-lime glass.
    ///
    /// Note: glass is brittle — `yield_strength` represents tensile fracture
    /// strength, not a ductile yield point.
    #[must_use]
    pub fn glass() -> Self {
        Self {
            name: Cow::Borrowed("Glass"),
            youngs_modulus: 70e9,
            poisson_ratio: 0.22,
            yield_strength: 33e6,
            ultimate_tensile_strength: 33e6,
            density: 2500.0,
            thermal_expansion: 9e-6,
        }
    }

    /// Natural rubber (vulcanized).
    ///
    /// Note: rubber properties are highly compound-dependent.
    /// `yield_strength` represents approximate tensile strength.
    #[must_use]
    pub fn rubber() -> Self {
        Self {
            name: Cow::Borrowed("Rubber"),
            youngs_modulus: 0.01e9,
            poisson_ratio: 0.49,
            yield_strength: 15e6,
            ultimate_tensile_strength: 15e6,
            density: 1100.0,
            thermal_expansion: 120e-6,
        }
    }

    /// Normal-strength concrete (C30).
    ///
    /// Note: `yield_strength` is the compressive strength. Concrete has
    /// negligible tensile strength (~3 MPa) and is primarily a
    /// compressive material.
    #[must_use]
    pub fn concrete() -> Self {
        Self {
            name: Cow::Borrowed("Concrete"),
            youngs_modulus: 30e9,
            poisson_ratio: 0.20,
            yield_strength: 30e6,
            ultimate_tensile_strength: 30e6,
            density: 2400.0,
            thermal_expansion: 12e-6,
        }
    }

    /// Oak wood (along grain).
    ///
    /// Note: wood is anisotropic — all values are longitudinal (along grain).
    /// Cross-grain properties differ significantly.
    #[must_use]
    pub fn wood_oak() -> Self {
        Self {
            name: Cow::Borrowed("Oak"),
            youngs_modulus: 12e9,
            poisson_ratio: 0.35,
            yield_strength: 40e6,
            ultimate_tensile_strength: 60e6,
            density: 600.0,
            thermal_expansion: 5e-6,
        }
    }

    /// Carbon fiber composite (unidirectional, ~60% Vf, intermediate modulus).
    ///
    /// Note: values are longitudinal. Transverse properties are much lower.
    /// `yield_strength` represents ultimate tensile strength (composites
    /// are brittle — no ductile yield).
    #[must_use]
    pub fn carbon_fiber() -> Self {
        Self {
            name: Cow::Borrowed("Carbon Fiber"),
            youngs_modulus: 181e9,
            poisson_ratio: 0.27,
            yield_strength: 1800e6,
            ultimate_tensile_strength: 1800e6,
            density: 1600.0,
            thermal_expansion: -0.5e-6,
        }
    }

    /// Stainless steel 304 (austenitic).
    #[must_use]
    pub fn stainless_steel_304() -> Self {
        Self {
            name: Cow::Borrowed("Stainless Steel 304"),
            youngs_modulus: 193e9,
            poisson_ratio: 0.29,
            yield_strength: 215e6,
            ultimate_tensile_strength: 505e6,
            density: 8000.0,
            thermal_expansion: 17.3e-6,
        }
    }

    /// Gray cast iron (ASTM Class 40).
    ///
    /// Note: cast iron is brittle in tension. `yield_strength` is the
    /// 0.2% proof stress; tensile fracture occurs at ~293 MPa with
    /// minimal plastic deformation.
    #[must_use]
    pub fn cast_iron() -> Self {
        Self {
            name: Cow::Borrowed("Gray Cast Iron"),
            youngs_modulus: 130e9,
            poisson_ratio: 0.26,
            yield_strength: 276e6,
            ultimate_tensile_strength: 293e6,
            density: 7200.0,
            thermal_expansion: 10.8e-6,
        }
    }

    /// Brass (C36000, free-cutting).
    #[must_use]
    pub fn brass() -> Self {
        Self {
            name: Cow::Borrowed("Brass C36000"),
            youngs_modulus: 100e9,
            poisson_ratio: 0.31,
            yield_strength: 140e6,
            ultimate_tensile_strength: 340e6,
            density: 8500.0,
            thermal_expansion: 20.5e-6,
        }
    }

    /// High-density polyethylene (HDPE).
    #[must_use]
    pub fn hdpe() -> Self {
        Self {
            name: Cow::Borrowed("HDPE"),
            youngs_modulus: 1.1e9,
            poisson_ratio: 0.42,
            yield_strength: 26e6,
            ultimate_tensile_strength: 33e6,
            density: 950.0,
            thermal_expansion: 120e-6,
        }
    }
}

impl Default for Material {
    /// Defaults to structural steel.
    fn default() -> Self {
        Self::steel()
    }
}

/// Temperature-dependent material properties.
///
/// Stores material properties at discrete temperatures and linearly
/// interpolates for arbitrary temperature queries.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TempDependentMaterial {
    /// Material name / identifier.
    pub name: Cow<'static, str>,
    /// (temperature_K, material_at_that_temperature) pairs,
    /// sorted by ascending temperature.
    points: Vec<(f64, Material)>,
}

impl TempDependentMaterial {
    /// Create from a set of (temperature, material) data points.
    ///
    /// Points will be sorted by temperature. At least one point is required.
    pub fn new(
        name: impl Into<Cow<'static, str>>,
        mut points: Vec<(f64, Material)>,
    ) -> crate::Result<Self> {
        if points.is_empty() {
            return Err(crate::DravyaError::InvalidMaterial(
                "at least one temperature point required".into(),
            ));
        }
        points.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        Ok(Self {
            name: name.into(),
            points,
        })
    }

    /// Interpolate material properties at the given temperature (K).
    ///
    /// Clamps to the nearest data point if outside the defined range.
    #[must_use]
    pub fn at_temperature(&self, temp: f64) -> Material {
        if self.points.len() == 1 {
            return self.points[0].1.clone();
        }

        // Clamp to range
        let first = &self.points[0];
        let last = &self.points[self.points.len() - 1];
        if temp <= first.0 {
            return first.1.clone();
        }
        if temp >= last.0 {
            return last.1.clone();
        }

        // Find bracketing interval
        for window in self.points.windows(2) {
            let (t_lo, m_lo) = &window[0];
            let (t_hi, m_hi) = &window[1];
            if temp >= *t_lo && temp <= *t_hi {
                let dt = t_hi - t_lo;
                if dt.abs() < hisab::EPSILON_F64 {
                    return m_lo.clone();
                }
                let t = (temp - t_lo) / dt;
                return Self::lerp_material(m_lo, m_hi, t);
            }
        }

        last.1.clone()
    }

    /// Temperature range covered by the data.
    #[must_use]
    pub fn temperature_range(&self) -> (f64, f64) {
        (self.points[0].0, self.points[self.points.len() - 1].0)
    }

    fn lerp_material(a: &Material, b: &Material, t: f64) -> Material {
        let lerp = |x: f64, y: f64| x + t * (y - x);
        Material {
            name: a.name.clone(),
            youngs_modulus: lerp(a.youngs_modulus, b.youngs_modulus),
            poisson_ratio: lerp(a.poisson_ratio, b.poisson_ratio),
            yield_strength: lerp(a.yield_strength, b.yield_strength),
            ultimate_tensile_strength: lerp(
                a.ultimate_tensile_strength,
                b.ultimate_tensile_strength,
            ),
            density: lerp(a.density, b.density),
            thermal_expansion: lerp(a.thermal_expansion, b.thermal_expansion),
        }
    }
}

impl fmt::Display for Material {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}: E={:.1} GPa, v={:.2}, σ_y={:.0} MPa, ρ={:.0} kg/m³",
            self.name,
            self.youngs_modulus / 1e9,
            self.poisson_ratio,
            self.yield_strength / 1e6,
            self.density
        )
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
        assert!(s.ultimate_tensile_strength > s.yield_strength);
    }

    #[test]
    fn carbon_fiber_negative_thermal_expansion() {
        let cf = Material::carbon_fiber();
        assert!(cf.thermal_expansion < 0.0, "carbon fiber has negative CTE");
    }

    #[test]
    fn rubber_nearly_incompressible() {
        let r = Material::rubber();
        assert!(
            r.poisson_ratio > 0.48,
            "rubber should be nearly incompressible"
        );
    }

    #[test]
    fn all_presets_positive_modulus() {
        let mats = [
            Material::steel(),
            Material::aluminum(),
            Material::copper(),
            Material::titanium(),
            Material::glass(),
            Material::rubber(),
            Material::concrete(),
            Material::wood_oak(),
            Material::carbon_fiber(),
            Material::stainless_steel_304(),
            Material::cast_iron(),
            Material::brass(),
            Material::hdpe(),
        ];
        for m in &mats {
            assert!(m.youngs_modulus > 0.0, "{} should have positive E", m.name);
            assert!(m.density > 0.0, "{} should have positive density", m.name);
            assert!(
                m.ultimate_tensile_strength >= m.yield_strength,
                "{} UTS should be >= yield",
                m.name
            );
        }
    }

    #[test]
    fn serde_roundtrip() {
        let m = Material::steel();
        let json = serde_json::to_string(&m).unwrap();
        let back: Material = serde_json::from_str(&json).unwrap();
        assert_eq!(m, back);
    }

    #[test]
    fn derived_elastic_properties() {
        let s = Material::steel();
        let g = s.shear_modulus();
        assert!(
            (g - 76.9e9).abs() < 0.5e9,
            "steel G should be ~77 GPa, got {g}"
        );
        let k = s.bulk_modulus();
        assert!(
            (k - 166.7e9).abs() < 1e9,
            "steel K should be ~167 GPa, got {k}"
        );
    }

    #[test]
    fn thermal_strain_steel() {
        let s = Material::steel();
        let strain = s.thermal_strain(100.0); // 100 K rise
        // 12e-6 * 100 = 1.2e-3
        assert!(
            (strain - 1.2e-3).abs() < hisab::EPSILON_F64,
            "thermal strain should be 1.2e-3, got {strain}"
        );
    }

    #[test]
    fn thermal_stress_steel() {
        let s = Material::steel();
        let stress = s.thermal_stress(100.0);
        // 200e9 * 12e-6 * 100 = 240 MPa
        assert!(
            (stress - 240e6).abs() < 1.0,
            "thermal stress should be 240 MPa, got {stress}"
        );
    }

    #[test]
    fn default_is_steel() {
        assert_eq!(Material::default(), Material::steel());
    }

    #[test]
    fn display_format() {
        let s = Material::steel();
        let display = s.to_string();
        assert!(display.contains("Steel"));
        assert!(display.contains("200.0 GPa"));
    }

    #[test]
    fn copper_yield_reasonable() {
        // Annealed copper yield ~33-70 MPa
        let cu = Material::copper();
        assert!(
            cu.yield_strength < 100e6,
            "annealed copper yield should be < 100 MPa"
        );
    }

    #[test]
    fn concrete_compressive_strength() {
        let c = Material::concrete();
        assert!(
            c.yield_strength >= 20e6,
            "concrete compressive strength should be >= 20 MPa"
        );
    }

    #[test]
    fn new_valid() {
        let m = Material::new("Test", 200e9, 0.30, 250e6, 400e6, 7850.0, 12e-6);
        assert!(m.is_ok());
    }

    #[test]
    fn new_negative_modulus() {
        let m = Material::new("Bad", -1.0, 0.30, 250e6, 400e6, 7850.0, 12e-6);
        assert!(m.is_err());
    }

    #[test]
    fn new_invalid_poisson() {
        let m = Material::new("Bad", 200e9, 0.5, 250e6, 400e6, 7850.0, 12e-6);
        assert!(m.is_err());
        let m = Material::new("Bad", 200e9, -1.0, 250e6, 400e6, 7850.0, 12e-6);
        assert!(m.is_err());
    }

    #[test]
    fn new_uts_less_than_yield() {
        let m = Material::new("Bad", 200e9, 0.30, 400e6, 200e6, 7850.0, 12e-6);
        assert!(m.is_err());
    }

    #[test]
    fn new_zero_density() {
        let m = Material::new("Bad", 200e9, 0.30, 250e6, 400e6, 0.0, 12e-6);
        assert!(m.is_err());
    }

    #[test]
    fn new_negative_yield_ok_zero() {
        // Zero yield is valid (e.g., fluids)
        let m = Material::new("Fluid", 200e9, 0.30, 0.0, 0.0, 1000.0, 0.0);
        assert!(m.is_ok());
    }

    // --- Temperature-dependent material tests ---

    #[test]
    fn temp_dependent_single_point() {
        let tdm = TempDependentMaterial::new("Steel", vec![(293.0, Material::steel())]).unwrap();
        let m = tdm.at_temperature(500.0);
        assert_eq!(m.youngs_modulus, Material::steel().youngs_modulus);
    }

    #[test]
    fn temp_dependent_interpolation() {
        let m_low = Material::steel(); // E = 200 GPa at 293 K
        let mut m_high = Material::steel();
        m_high.youngs_modulus = 150e9; // E = 150 GPa at 800 K

        let tdm =
            TempDependentMaterial::new("Steel", vec![(293.0, m_low), (800.0, m_high)]).unwrap();

        // Midpoint: (293+800)/2 = 546.5 K -> E should be ~175 GPa
        let m = tdm.at_temperature(546.5);
        assert!(
            (m.youngs_modulus - 175e9).abs() < 1e9,
            "E at midpoint should be ~175 GPa, got {}",
            m.youngs_modulus / 1e9
        );
    }

    #[test]
    fn temp_dependent_clamps_low() {
        let tdm = TempDependentMaterial::new(
            "Steel",
            vec![(293.0, Material::steel()), (800.0, Material::steel())],
        )
        .unwrap();
        let m = tdm.at_temperature(100.0); // below range
        assert_eq!(m.youngs_modulus, Material::steel().youngs_modulus);
    }

    #[test]
    fn temp_dependent_clamps_high() {
        let mut m_high = Material::steel();
        m_high.youngs_modulus = 150e9;
        let tdm = TempDependentMaterial::new(
            "Steel",
            vec![(293.0, Material::steel()), (800.0, m_high.clone())],
        )
        .unwrap();
        let m = tdm.at_temperature(1000.0); // above range
        assert_eq!(m.youngs_modulus, m_high.youngs_modulus);
    }

    #[test]
    fn temp_dependent_empty_errors() {
        let result = TempDependentMaterial::new("Empty", vec![]);
        assert!(result.is_err());
    }

    #[test]
    fn temp_dependent_range() {
        let tdm = TempDependentMaterial::new(
            "Steel",
            vec![(293.0, Material::steel()), (800.0, Material::steel())],
        )
        .unwrap();
        assert_eq!(tdm.temperature_range(), (293.0, 800.0));
    }

    #[test]
    fn temp_dependent_modulus_decreases() {
        let mut m_high = Material::steel();
        m_high.youngs_modulus = 150e9;
        let tdm =
            TempDependentMaterial::new("Steel", vec![(293.0, Material::steel()), (800.0, m_high)])
                .unwrap();
        let m_400 = tdm.at_temperature(400.0);
        let m_600 = tdm.at_temperature(600.0);
        assert!(
            m_600.youngs_modulus < m_400.youngs_modulus,
            "E should decrease with temperature"
        );
    }
}
