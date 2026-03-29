//! Soorat integration — visualization data structures for material science analysis.
//!
//! Provides structured types that soorat can render: stress/strain field heatmaps,
//! deformed geometry, fracture surfaces, and fatigue damage maps.

use serde::{Deserialize, Serialize};

// ── Stress/strain field visualization ──────────────────────────────────────

/// A 2D scalar field of stress or strain values for heatmap rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StressFieldVisualization {
    /// Scalar values at grid points (von Mises stress, principal stress, etc.).
    /// Flattened row-major: `values[y * nx + x]`.
    pub values: Vec<f64>,
    /// Grid dimensions (nx, ny).
    pub dimensions: [usize; 2],
    /// World-space origin `[x, y]` in metres.
    pub origin: [f64; 2],
    /// Grid spacing in metres.
    pub spacing: f64,
    /// Field name (e.g. "von_mises", "sigma_xx", "epsilon_eq").
    pub field_name: String,
    /// Minimum value in the field.
    pub min_value: f64,
    /// Maximum value in the field.
    pub max_value: f64,
    /// Unit label (e.g. "Pa", "MPa", "dimensionless").
    pub unit: String,
}

impl StressFieldVisualization {
    /// Create from a flat array of scalar values.
    #[must_use]
    pub fn from_values(
        values: Vec<f64>,
        dimensions: [usize; 2],
        origin: [f64; 2],
        spacing: f64,
        field_name: &str,
        unit: &str,
    ) -> Self {
        let min_value = values.iter().cloned().fold(f64::MAX, f64::min);
        let max_value = values.iter().cloned().fold(f64::MIN, f64::max);
        Self {
            values,
            dimensions,
            origin,
            spacing,
            field_name: field_name.to_string(),
            min_value,
            max_value,
            unit: unit.to_string(),
        }
    }
}

// ── Deformed geometry ──────────────────────────────────────────────────────

/// Displacement field for deformed shape rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DeformationVisualization {
    /// Original node positions `[x, y, z]`.
    pub positions: Vec<[f64; 3]>,
    /// Displacement at each node `[dx, dy, dz]` in metres.
    pub displacements: Vec<[f64; 3]>,
    /// Displacement magnitude at each node (for color mapping).
    pub magnitudes: Vec<f64>,
    /// Maximum displacement magnitude.
    pub max_displacement: f64,
    /// Scale factor applied to displacements for visualization (typically > 1).
    pub display_scale: f64,
}

impl DeformationVisualization {
    /// Create from positions and displacements.
    #[must_use]
    pub fn new(positions: Vec<[f64; 3]>, displacements: Vec<[f64; 3]>, display_scale: f64) -> Self {
        let magnitudes: Vec<f64> = displacements
            .iter()
            .map(|d| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt())
            .collect();
        let max_displacement = magnitudes.iter().cloned().fold(0.0_f64, f64::max);
        Self {
            positions,
            displacements,
            magnitudes,
            max_displacement,
            display_scale,
        }
    }

    /// Get the deformed position of node `i` (with display scaling).
    #[must_use]
    pub fn deformed_position(&self, i: usize) -> Option<[f64; 3]> {
        let pos = self.positions.get(i)?;
        let disp = self.displacements.get(i)?;
        Some([
            pos[0] + disp[0] * self.display_scale,
            pos[1] + disp[1] * self.display_scale,
            pos[2] + disp[2] * self.display_scale,
        ])
    }
}

// ── Fracture surface ───────────────────────────────────────────────────────

/// Crack path data for line/surface rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FractureVisualization {
    /// Crack tip positions over time `[x, y]`.
    pub crack_path: Vec<[f64; 2]>,
    /// Stress intensity factor at each point (for color mapping).
    pub stress_intensity: Vec<f64>,
    /// Critical SIF (K_Ic) for reference.
    pub k_ic: f64,
    /// Current crack length (m).
    pub crack_length: f64,
}

// ── Fatigue damage map ─────────────────────────────────────────────────────

/// Fatigue damage fraction (0–1) at grid points for color-coded surface rendering.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FatigueDamageMap {
    /// Damage fraction at each point (0 = undamaged, 1 = failed).
    /// Flattened row-major: `damage[y * nx + x]`.
    pub damage: Vec<f64>,
    /// Grid dimensions (nx, ny).
    pub dimensions: [usize; 2],
    /// Grid spacing in metres.
    pub spacing: f64,
    /// Number of cycles applied.
    pub cycles: u64,
    /// Maximum damage value in the field.
    pub max_damage: f64,
}

impl FatigueDamageMap {
    /// Create from a flat array of damage values.
    #[must_use]
    pub fn from_values(
        damage: Vec<f64>,
        dimensions: [usize; 2],
        spacing: f64,
        cycles: u64,
    ) -> Self {
        let max_damage = damage.iter().cloned().fold(0.0_f64, f64::max);
        Self {
            damage,
            dimensions,
            spacing,
            cycles,
            max_damage,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stress_field_basic() {
        let values = vec![100.0, 200.0, 150.0, 300.0];
        let field = StressFieldVisualization::from_values(
            values,
            [2, 2],
            [0.0, 0.0],
            1.0,
            "von_mises",
            "Pa",
        );
        assert_eq!(field.dimensions, [2, 2]);
        assert!((field.min_value - 100.0).abs() < 0.01);
        assert!((field.max_value - 300.0).abs() < 0.01);
    }

    #[test]
    fn deformation_basic() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let displacements = vec![[0.0, 0.001, 0.0], [0.0, 0.002, 0.0]];
        let viz = DeformationVisualization::new(positions, displacements, 100.0);
        assert_eq!(viz.magnitudes.len(), 2);
        assert!((viz.max_displacement - 0.002).abs() < 1e-6);

        let def = viz.deformed_position(0).unwrap();
        assert!((def[1] - 0.1).abs() < 1e-6); // 0.001 × 100
    }

    #[test]
    fn deformation_out_of_bounds() {
        let viz = DeformationVisualization::new(vec![], vec![], 1.0);
        assert!(viz.deformed_position(0).is_none());
    }

    #[test]
    fn fracture_serializes() {
        let frac = FractureVisualization {
            crack_path: vec![[0.0, 0.0], [0.01, 0.0], [0.02, 0.001]],
            stress_intensity: vec![10e6, 15e6, 20e6],
            k_ic: 30e6,
            crack_length: 0.02,
        };
        let json = serde_json::to_string(&frac);
        assert!(json.is_ok());
    }

    #[test]
    fn fatigue_map_basic() {
        let damage = vec![0.0, 0.1, 0.5, 0.9];
        let map = FatigueDamageMap::from_values(damage, [2, 2], 0.01, 100_000);
        assert_eq!(map.cycles, 100_000);
        assert!((map.max_damage - 0.9).abs() < 0.01);
    }

    #[test]
    fn fatigue_map_undamaged() {
        let map = FatigueDamageMap::from_values(vec![0.0; 4], [2, 2], 1.0, 0);
        assert_eq!(map.max_damage, 0.0);
    }

    #[test]
    fn stress_field_serializes() {
        let field = StressFieldVisualization::from_values(
            vec![1.0, 2.0],
            [2, 1],
            [0.0, 0.0],
            0.5,
            "sigma_xx",
            "MPa",
        );
        let json = serde_json::to_string(&field);
        assert!(json.is_ok());
    }
}
