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
        name: material.name.clone(),
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

// ---------------------------------------------------------------------------
// Transient thermal-mechanical analysis
// ---------------------------------------------------------------------------

/// Result of a single transient thermal-mechanical time step.
#[derive(Debug, Clone)]
pub struct ThermalStepResult {
    /// Temperature at each node (K).
    pub temperatures: Vec<f64>,
    /// Von Mises stress at each node (Pa).
    pub von_mises: Vec<f64>,
    /// Whether any node has yielded.
    pub any_yielded: bool,
    /// Index of first yielded node, if any.
    pub first_yield_node: Option<usize>,
    /// Time at this step (s).
    pub time: f64,
}

/// Run a transient 1D thermal-mechanical simulation.
///
/// Steps the thermal grid through time using Crank-Nicolson, computing
/// the constrained thermal stress at each node after each step.
///
/// Returns the stress history at every `record_interval` steps.
pub fn transient_thermal_stress_1d(
    material: &Material,
    grid: &mut ushma::numerical::ThermalGrid1D,
    reference_temp: f64,
    dt: f64,
    n_steps: usize,
    record_interval: usize,
) -> Vec<ThermalStepResult> {
    let interval = record_interval.max(1);
    let mut history = Vec::with_capacity(n_steps / interval + 1);
    let mut time = 0.0;

    for step in 0..n_steps {
        let _ = grid.step_crank_nicolson(dt);
        time += dt;

        if (step + 1) % interval == 0 || step == n_steps - 1 {
            let mut von_mises = Vec::with_capacity(grid.nodes.len());
            let mut any_yielded = false;
            let mut first_yield = None;

            for (i, &t) in grid.nodes.iter().enumerate() {
                let stress = constrained_thermal_stress(material, t - reference_temp);
                let vm = stress.von_mises();
                if vm >= material.yield_strength && first_yield.is_none() {
                    any_yielded = true;
                    first_yield = Some(i);
                }
                von_mises.push(vm);
            }

            history.push(ThermalStepResult {
                temperatures: grid.nodes.clone(),
                von_mises,
                any_yielded,
                first_yield_node: first_yield,
                time,
            });
        }
    }

    history
}

/// Extract peak thermal stress history at a specific node from transient results.
///
/// Returns a vector of von Mises stress values over time for rainflow
/// or fatigue analysis.
#[must_use]
pub fn node_stress_history(results: &[ThermalStepResult], node_index: usize) -> Vec<f64> {
    results
        .iter()
        .filter_map(|r| r.von_mises.get(node_index).copied())
        .collect()
}

/// Compute thermal fatigue damage at a node from transient results.
///
/// Extracts the stress history at the given node, runs rainflow counting,
/// and accumulates Miner's rule damage using Basquin's law.
///
/// `fatigue_strength_coeff` and `fatigue_exponent` are the cycle-based
/// Basquin parameters.
///
/// Returns cumulative damage (>= 1.0 means fatigue failure predicted).
#[must_use]
pub fn thermal_fatigue_damage(
    results: &[ThermalStepResult],
    node_index: usize,
    fatigue_strength_coeff: f64,
    fatigue_exponent: f64,
) -> f64 {
    let history = node_stress_history(results, node_index);
    let turning = crate::fatigue::extract_turning_points(&history);
    let cycles = crate::fatigue::rainflow_count(&turning);

    cycles
        .iter()
        .map(|&(range, _mean, count)| {
            let amplitude = range / 2.0;
            let n_f =
                crate::fatigue::basquin_cycles(amplitude, fatigue_strength_coeff, fatigue_exponent);
            if n_f > 0.0 { count / n_f } else { 0.0 }
        })
        .sum()
}

/// Find the time at which yield first occurs in a transient simulation.
///
/// Returns `None` if no yielding occurred during the simulation.
#[must_use]
pub fn time_to_yield(results: &[ThermalStepResult]) -> Option<f64> {
    results.iter().find(|r| r.any_yielded).map(|r| r.time)
}

/// Maximum von Mises stress across all nodes and all time steps.
#[must_use]
pub fn peak_thermal_stress(results: &[ThermalStepResult]) -> f64 {
    results
        .iter()
        .flat_map(|r| r.von_mises.iter())
        .copied()
        .fold(0.0_f64, f64::max)
}

/// Compute thermal stress at each node of an ushma thermal network.
///
/// Solves the network for steady-state temperatures, then computes
/// constrained thermal stress at each node.
pub fn stress_from_thermal_network(
    material: &Material,
    network: &ushma::numerical::ThermalNetwork,
    reference_temp: f64,
) -> crate::Result<Vec<StressTensor>> {
    let temps = network.solve().map_err(|e| {
        crate::DravyaError::InvalidMaterial(format!("thermal network solve failed: {e}"))
    })?;
    Ok(temps
        .iter()
        .map(|&t| constrained_thermal_stress(material, t - reference_temp))
        .collect())
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

    // --- Transient analysis tests ---

    #[test]
    fn transient_thermal_stress_runs() {
        let steel = Material::steel();
        let mut grid = ushma::numerical::ThermalGrid1D::new(
            5,
            0.01,
            1.12e-5,
            400.0,
            ushma::numerical::BoundaryCondition::Fixed(400.0),
            ushma::numerical::BoundaryCondition::Fixed(400.0),
        )
        .expect("grid creation");

        let history = transient_thermal_stress_1d(&steel, &mut grid, 293.0, 0.001, 20, 5);
        assert_eq!(history.len(), 4, "20 steps / 5 interval = 4 records");
        // Verify structure: each record has correct node count and advancing time
        for (i, record) in history.iter().enumerate() {
            assert_eq!(record.temperatures.len(), 5);
            assert_eq!(record.von_mises.len(), 5);
            if i > 0 {
                assert!(record.time > history[i - 1].time, "time should advance");
            }
        }
    }

    #[test]
    fn transient_time_progresses() {
        let steel = Material::steel();
        let mut grid = ushma::numerical::ThermalGrid1D::new(
            5,
            0.1,
            1e-5,
            293.0,
            ushma::numerical::BoundaryCondition::Fixed(293.0),
            ushma::numerical::BoundaryCondition::Fixed(293.0),
        )
        .expect("grid creation");

        let history = transient_thermal_stress_1d(&steel, &mut grid, 293.0, 0.1, 50, 25);
        assert_eq!(history.len(), 2);
        assert!((history[0].time - 2.5).abs() < 0.01);
        assert!((history[1].time - 5.0).abs() < 0.01);
    }

    #[test]
    fn node_stress_history_extraction() {
        let steel = Material::steel();
        let mut grid = ushma::numerical::ThermalGrid1D::new(
            5,
            0.1,
            1e-5,
            350.0,
            ushma::numerical::BoundaryCondition::Fixed(350.0),
            ushma::numerical::BoundaryCondition::Fixed(350.0),
        )
        .expect("grid creation");

        let history = transient_thermal_stress_1d(&steel, &mut grid, 293.0, 0.01, 20, 5);
        let node_hist = node_stress_history(&history, 2);
        assert_eq!(node_hist.len(), history.len());
    }

    #[test]
    fn thermal_fatigue_damage_is_non_negative() {
        let steel = Material::steel();
        let mut grid = ushma::numerical::ThermalGrid1D::new(
            5,
            0.1,
            1e-5,
            350.0,
            ushma::numerical::BoundaryCondition::Fixed(350.0),
            ushma::numerical::BoundaryCondition::Fixed(293.0),
        )
        .expect("grid creation");

        let history = transient_thermal_stress_1d(&steel, &mut grid, 293.0, 0.01, 50, 5);
        let damage = thermal_fatigue_damage(&history, 2, 1000e6, -0.1);
        assert!(damage >= 0.0, "damage should be non-negative");
    }

    #[test]
    fn peak_thermal_stress_from_snapshot() {
        // Test peak_thermal_stress utility using manually constructed results
        let results = vec![
            ThermalStepResult {
                temperatures: vec![300.0, 350.0, 400.0],
                von_mises: vec![10e6, 50e6, 100e6],
                any_yielded: false,
                first_yield_node: None,
                time: 1.0,
            },
            ThermalStepResult {
                temperatures: vec![310.0, 360.0, 410.0],
                von_mises: vec![20e6, 80e6, 150e6],
                any_yielded: false,
                first_yield_node: None,
                time: 2.0,
            },
        ];
        let peak = peak_thermal_stress(&results);
        assert!((peak - 150e6).abs() < 1.0, "peak should be 150 MPa");
    }

    #[test]
    fn time_to_yield_none_for_small_gradient() {
        let steel = Material::steel();
        let mut grid = ushma::numerical::ThermalGrid1D::new(
            5,
            0.1,
            1e-5,
            295.0, // tiny ΔT
            ushma::numerical::BoundaryCondition::Fixed(295.0),
            ushma::numerical::BoundaryCondition::Fixed(293.0),
        )
        .expect("grid creation");

        let history = transient_thermal_stress_1d(&steel, &mut grid, 293.0, 0.01, 20, 5);
        assert!(
            time_to_yield(&history).is_none(),
            "tiny gradient should not yield"
        );
    }

    #[test]
    fn thermal_network_stress() {
        let steel = Material::steel();
        let mut network = ushma::numerical::ThermalNetwork::new(3);
        network.add_resistance(0, 1, 0.1).unwrap();
        network.add_resistance(1, 2, 0.1).unwrap();
        network.set_fixed_temperature(0, 400.0).unwrap();
        network.set_fixed_temperature(2, 293.0).unwrap();

        let stresses = stress_from_thermal_network(&steel, &network, 293.0);
        assert!(stresses.is_ok());
        let stresses = stresses.unwrap();
        assert_eq!(stresses.len(), 3, "should have stress for each node");
    }
}
