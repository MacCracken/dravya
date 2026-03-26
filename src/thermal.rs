//! Thermal-mechanical coupling with [`ushma`].
//!
//! Requires the `thermal` feature: `dravya = { features = ["thermal"] }`
//!
//! Provides conversion between dravya [`Material`] and ushma `ThermalMaterial`,
//! thermal strain/stress tensors from temperature changes, and utilities for
//! computing mechanical response from thermal analysis results.

use crate::material::Material;
use crate::strain::StrainTensor;
use crate::stress::StressTensor;

/// Convert a dravya [`Material`] to an ushma [`ushma::material::ThermalMaterial`].
///
/// Maps shared properties (density) and uses the material name.
/// Thermal conductivity and specific heat are not stored in dravya's
/// `Material`, so they must be provided separately.
#[must_use]
pub fn to_thermal_material(
    material: &Material,
    conductivity: f64,
    specific_heat: f64,
) -> ushma::material::ThermalMaterial {
    ushma::material::ThermalMaterial {
        name: material.name.clone().into(),
        conductivity,
        specific_heat,
        density: material.density,
        melting_point: 0.0,
        boiling_point: 0.0,
    }
}

/// Isotropic thermal strain tensor from a uniform temperature change.
///
/// ε_th = α * ΔT on all normal components, zero shear.
#[must_use]
#[inline]
pub fn thermal_strain_tensor(material: &Material, delta_t: f64) -> StrainTensor {
    let e = material.thermal_expansion * delta_t;
    StrainTensor::new(e, e, e, 0.0, 0.0, 0.0)
}

/// Thermal stress tensor for a fully constrained body.
///
/// When all deformation is prevented (ε_total = 0), the thermal stress is:
/// σ = -\[C\] * ε_thermal
///
/// Uses the isotropic stiffness matrix to compute the full 3D thermal stress.
#[must_use]
pub fn constrained_thermal_stress(material: &Material, delta_t: f64) -> StressTensor {
    // For isotropic material under uniform thermal strain:
    // σ_ii = -(λ + 2μ)*α*ΔT - λ*α*ΔT - λ*α*ΔT = -(3λ + 2μ)*α*ΔT
    // Which equals -E*α*ΔT / (1 - 2v) for each normal component
    let v = material.poisson_ratio;
    let denom = 1.0 - 2.0 * v;
    let sigma = if denom.abs() < hisab::EPSILON_F64 {
        0.0
    } else {
        -material.youngs_modulus * material.thermal_expansion * delta_t / denom
    };
    StressTensor::new(sigma, sigma, sigma, 0.0, 0.0, 0.0)
}

/// Mechanical strain from total strain minus thermal strain.
///
/// ε_mech = ε_total - ε_thermal
///
/// The mechanical strain can then be used with constitutive models
/// to compute stress.
#[must_use]
#[inline]
pub fn mechanical_strain(
    total_strain: &StrainTensor,
    material: &Material,
    delta_t: f64,
) -> StrainTensor {
    let eth = thermal_strain_tensor(material, delta_t);
    *total_strain - eth
}

/// Compute thermal stress at each node of an ushma 1D thermal grid.
///
/// Given a reference temperature and a solved temperature field,
/// returns the constrained thermal stress at each node.
///
/// This assumes full constraint (no free expansion). For partial
/// constraint, use [`mechanical_strain`] with actual displacements.
#[must_use]
pub fn stress_from_thermal_grid_1d(
    material: &Material,
    grid: &ushma::numerical::ThermalGrid1D,
    reference_temp: f64,
) -> Vec<StressTensor> {
    grid.nodes
        .iter()
        .map(|&t| constrained_thermal_stress(material, t - reference_temp))
        .collect()
}

/// Compute thermal stress at each node of an ushma 2D thermal grid.
///
/// Returns a 2D vector of stress tensors matching the grid dimensions.
#[must_use]
pub fn stress_from_thermal_grid_2d(
    material: &Material,
    grid: &ushma::numerical::ThermalGrid2D,
    reference_temp: f64,
) -> Vec<Vec<StressTensor>> {
    grid.nodes
        .iter()
        .map(|row: &Vec<f64>| {
            row.iter()
                .map(|&t| constrained_thermal_stress(material, t - reference_temp))
                .collect()
        })
        .collect()
}

/// Check if thermal stress exceeds yield at any node in a 1D grid.
///
/// Returns the index of the first node that yields, or `None` if safe.
#[must_use]
pub fn find_thermal_yield_1d(
    material: &Material,
    grid: &ushma::numerical::ThermalGrid1D,
    reference_temp: f64,
) -> Option<usize> {
    for (i, &t) in grid.nodes.iter().enumerate() {
        let stress = constrained_thermal_stress(material, t - reference_temp);
        if stress.von_mises() >= material.yield_strength {
            return Some(i);
        }
    }
    None
}

/// Maximum temperature change before yielding under full constraint.
///
/// ΔT_max = σ_y * (1 - 2v) / (E * α)
///
/// Returns `f64::INFINITY` if thermal expansion is zero.
#[must_use]
#[inline]
pub fn max_thermal_delta_t(material: &Material) -> f64 {
    let denom = material.youngs_modulus * material.thermal_expansion;
    if denom.abs() < hisab::EPSILON_F64 {
        return f64::INFINITY;
    }
    let v = material.poisson_ratio;
    material.yield_strength * (1.0 - 2.0 * v) / denom
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thermal_strain_tensor_isotropic() {
        let steel = Material::steel();
        let ts = thermal_strain_tensor(&steel, 100.0);
        let expected = steel.thermal_expansion * 100.0; // 12e-6 * 100 = 1.2e-3
        assert!((ts.components[0] - expected).abs() < hisab::EPSILON_F64);
        assert!((ts.components[1] - expected).abs() < hisab::EPSILON_F64);
        assert!((ts.components[2] - expected).abs() < hisab::EPSILON_F64);
        // Shear components should be zero
        assert!(ts.components[3].abs() < hisab::EPSILON_F64);
        assert!(ts.components[4].abs() < hisab::EPSILON_F64);
        assert!(ts.components[5].abs() < hisab::EPSILON_F64);
    }

    #[test]
    fn constrained_thermal_stress_hydrostatic() {
        let steel = Material::steel();
        let stress = constrained_thermal_stress(&steel, 100.0);
        // Should be hydrostatic (equal normal stresses, zero shear)
        assert!((stress.components[0] - stress.components[1]).abs() < 1.0);
        assert!((stress.components[1] - stress.components[2]).abs() < 1.0);
        assert!(stress.components[3].abs() < 1.0);
        // σ = -E*α*ΔT / (1-2v) = -200e9 * 12e-6 * 100 / 0.4 = -600 MPa
        assert!(
            (stress.components[0] - (-600e6)).abs() < 1e3,
            "constrained thermal stress should be -600 MPa, got {}",
            stress.components[0] / 1e6
        );
    }

    #[test]
    fn constrained_thermal_stress_heating_is_compressive() {
        let steel = Material::steel();
        let stress = constrained_thermal_stress(&steel, 50.0);
        assert!(
            stress.components[0] < 0.0,
            "heating under constraint should produce compressive stress"
        );
    }

    #[test]
    fn constrained_thermal_stress_cooling_is_tensile() {
        let steel = Material::steel();
        let stress = constrained_thermal_stress(&steel, -50.0);
        assert!(
            stress.components[0] > 0.0,
            "cooling under constraint should produce tensile stress"
        );
    }

    #[test]
    fn mechanical_strain_subtracts_thermal() {
        let steel = Material::steel();
        let total = StrainTensor::new(0.002, 0.002, 0.002, 0.0, 0.0, 0.0);
        let mech = mechanical_strain(&total, &steel, 100.0);
        let expected = 0.002 - steel.thermal_expansion * 100.0;
        assert!((mech.components[0] - expected).abs() < hisab::EPSILON_F64);
    }

    #[test]
    fn to_thermal_material_preserves_density() {
        let steel = Material::steel();
        let tm = to_thermal_material(&steel, 43.0, 490.0);
        assert!((tm.density - steel.density).abs() < hisab::EPSILON_F64);
        assert!((tm.conductivity - 43.0).abs() < hisab::EPSILON_F64);
        assert!((tm.specific_heat - 490.0).abs() < hisab::EPSILON_F64);
    }

    #[test]
    fn max_thermal_delta_t_steel() {
        let steel = Material::steel();
        let dt_max = max_thermal_delta_t(&steel);
        // ΔT_max = 250e6 * 0.4 / (200e9 * 12e-6) = 100e6 / 2.4e6 = 41.67 K
        assert!(
            (dt_max - 41.67).abs() < 0.1,
            "max ΔT for steel should be ~41.7 K, got {dt_max}"
        );
    }

    #[test]
    fn max_thermal_delta_t_zero_expansion() {
        let mut m = Material::steel();
        m.thermal_expansion = 0.0;
        assert!(max_thermal_delta_t(&m).is_infinite());
    }

    #[test]
    fn stress_from_grid_1d_basic() {
        let steel = Material::steel();
        let grid = ushma::numerical::ThermalGrid1D::new(
            5,
            0.1,
            steel.thermal_expansion,
            350.0, // 350 K initial
            ushma::numerical::BoundaryCondition::Fixed(350.0),
            ushma::numerical::BoundaryCondition::Fixed(350.0),
        )
        .expect("grid creation should succeed");
        let stresses = stress_from_thermal_grid_1d(&steel, &grid, 293.0);
        assert_eq!(stresses.len(), 5);
        // All nodes at 350 K, ref 293 K -> ΔT = 57 K -> compressive
        for s in &stresses {
            assert!(s.components[0] < 0.0, "heating should be compressive");
        }
    }

    #[test]
    fn find_thermal_yield_safe() {
        let steel = Material::steel();
        let grid = ushma::numerical::ThermalGrid1D::new(
            5,
            0.1,
            steel.thermal_expansion,
            300.0, // small ΔT from ref
            ushma::numerical::BoundaryCondition::Fixed(300.0),
            ushma::numerical::BoundaryCondition::Fixed(300.0),
        )
        .expect("grid creation should succeed");
        assert!(
            find_thermal_yield_1d(&steel, &grid, 293.0).is_none(),
            "small ΔT should not yield"
        );
    }
}
