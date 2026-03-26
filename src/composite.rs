//! Composite laminate mechanics: Classical Laminate Theory (CLT),
//! lamina stiffness, ply failure criteria.

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Lamina (single ply) properties and stiffness
// ---------------------------------------------------------------------------

/// Unidirectional lamina mechanical properties.
///
/// Properties are in the material coordinate system (1 = fiber, 2 = transverse).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Lamina {
    /// Longitudinal modulus E1 (Pa).
    pub e1: f64,
    /// Transverse modulus E2 (Pa).
    pub e2: f64,
    /// In-plane shear modulus G12 (Pa).
    pub g12: f64,
    /// Major Poisson's ratio v12.
    pub nu12: f64,
    /// Longitudinal tensile strength (Pa).
    pub xt: f64,
    /// Longitudinal compressive strength (Pa, positive value).
    pub xc: f64,
    /// Transverse tensile strength (Pa).
    pub yt: f64,
    /// Transverse compressive strength (Pa, positive value).
    pub yc: f64,
    /// In-plane shear strength (Pa).
    pub s12: f64,
}

impl Lamina {
    /// Minor Poisson's ratio: v21 = v12 * E2 / E1.
    #[must_use]
    #[inline]
    pub fn nu21(&self) -> f64 {
        if self.e1.abs() < hisab::EPSILON_F64 {
            return 0.0;
        }
        self.nu12 * self.e2 / self.e1
    }

    /// Reduced stiffness matrix Q (3x3) in material coordinates.
    ///
    /// Maps in-plane strain to stress for plane stress conditions:
    /// {σ1, σ2, τ12} = \[Q\] {ε1, ε2, γ12}
    #[must_use]
    pub fn stiffness_matrix(&self) -> [[f64; 3]; 3] {
        let nu21 = self.nu21();
        let denom = 1.0 - self.nu12 * nu21;
        if denom.abs() < hisab::EPSILON_F64 {
            return [[0.0; 3]; 3];
        }
        let q11 = self.e1 / denom;
        let q22 = self.e2 / denom;
        let q12 = self.nu12 * self.e2 / denom;
        let q66 = self.g12;

        [[q11, q12, 0.0], [q12, q22, 0.0], [0.0, 0.0, q66]]
    }

    /// Transformed stiffness matrix Q-bar (3x3) for a ply at angle `theta` (radians).
    ///
    /// Rotates the material stiffness into the laminate coordinate system.
    #[must_use]
    pub fn transformed_stiffness(&self, theta: f64) -> [[f64; 3]; 3] {
        let q = self.stiffness_matrix();
        let c = theta.cos();
        let s = theta.sin();
        let c2 = c * c;
        let s2 = s * s;
        let cs = c * s;
        let c4 = c2 * c2;
        let s4 = s2 * s2;

        let q11 = q[0][0];
        let q12 = q[0][1];
        let q22 = q[1][1];
        let q66 = q[2][2];

        let qb11 = q11 * c4 + 2.0 * (q12 + 2.0 * q66) * c2 * s2 + q22 * s4;
        let qb22 = q11 * s4 + 2.0 * (q12 + 2.0 * q66) * c2 * s2 + q22 * c4;
        let qb12 = (q11 + q22 - 4.0 * q66) * c2 * s2 + q12 * (c4 + s4);
        let qb66 = (q11 + q22 - 2.0 * q12 - 2.0 * q66) * c2 * s2 + q66 * (c4 + s4);
        let qb16 = (q11 - q12 - 2.0 * q66) * c2 * cs - (q22 - q12 - 2.0 * q66) * s2 * cs;
        let qb26 = (q11 - q12 - 2.0 * q66) * cs * s2 - (q22 - q12 - 2.0 * q66) * c2 * cs;

        [[qb11, qb12, qb16], [qb12, qb22, qb26], [qb16, qb26, qb66]]
    }

    /// Typical carbon/epoxy UD lamina (T300/914C class).
    #[must_use]
    pub fn carbon_epoxy() -> Self {
        Self {
            e1: 138e9,
            e2: 11e9,
            g12: 5.5e9,
            nu12: 0.28,
            xt: 1500e6,
            xc: 1200e6,
            yt: 50e6,
            yc: 250e6,
            s12: 70e6,
        }
    }

    /// Typical glass/epoxy UD lamina (E-glass/epoxy class).
    #[must_use]
    pub fn glass_epoxy() -> Self {
        Self {
            e1: 39e9,
            e2: 8.6e9,
            g12: 3.8e9,
            nu12: 0.28,
            xt: 1080e6,
            xc: 620e6,
            yt: 39e6,
            yc: 128e6,
            s12: 89e6,
        }
    }
}

// ---------------------------------------------------------------------------
// Ply definition and laminate layup
// ---------------------------------------------------------------------------

/// A single ply in a laminate.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Ply {
    /// Lamina material properties.
    pub lamina: Lamina,
    /// Ply angle (radians, measured from laminate x-axis).
    pub angle: f64,
    /// Ply thickness (m).
    pub thickness: f64,
}

/// ABD stiffness matrix result from Classical Laminate Theory.
///
/// The 6x6 matrix relates mid-plane forces/moments to strains/curvatures:
/// ```text
/// {N}   [A  B] {ε⁰}
/// {M} = [B  D] {κ}
/// ```
#[derive(Debug, Clone)]
pub struct AbdMatrix {
    /// Extensional stiffness (3x3), units: N/m.
    pub a: [[f64; 3]; 3],
    /// Coupling stiffness (3x3), units: N.
    pub b: [[f64; 3]; 3],
    /// Bending stiffness (3x3), units: N*m.
    pub d: [[f64; 3]; 3],
}

/// Compute the ABD stiffness matrix for a symmetric or general laminate.
///
/// `plies` defines the layup from bottom to top. Each ply has material,
/// angle, and thickness. The laminate mid-plane is at z = 0.
#[must_use]
pub fn abd_matrix(plies: &[Ply]) -> AbdMatrix {
    let mut a = [[0.0; 3]; 3];
    let mut b = [[0.0; 3]; 3];
    let mut d = [[0.0; 3]; 3];

    // Compute total thickness and z-coordinates
    let total_thickness: f64 = plies.iter().map(|p| p.thickness).sum();
    let mut z_bot = -total_thickness / 2.0;

    for ply in plies {
        let z_top = z_bot + ply.thickness;
        let qbar = ply.lamina.transformed_stiffness(ply.angle);

        let dz = z_top - z_bot;
        let dz2 = z_top * z_top - z_bot * z_bot;
        let dz3 = z_top * z_top * z_top - z_bot * z_bot * z_bot;

        for i in 0..3 {
            for j in 0..3 {
                a[i][j] += qbar[i][j] * dz;
                b[i][j] += qbar[i][j] * dz2 / 2.0;
                d[i][j] += qbar[i][j] * dz3 / 3.0;
            }
        }

        z_bot = z_top;
    }

    AbdMatrix { a, b, d }
}

// ---------------------------------------------------------------------------
// Ply failure criteria
// ---------------------------------------------------------------------------

/// Ply stress state in material coordinates (1-2 system).
#[derive(Debug, Clone, Copy)]
pub struct PlyStress {
    /// Stress in fiber direction (Pa).
    pub sigma1: f64,
    /// Stress transverse to fiber (Pa).
    pub sigma2: f64,
    /// In-plane shear stress (Pa).
    pub tau12: f64,
}

/// Maximum stress failure criterion.
///
/// Checks each stress component against its respective strength.
/// Returns failure index (>= 1.0 means failure).
#[must_use]
#[inline]
pub fn max_stress_failure_index(stress: &PlyStress, lamina: &Lamina) -> f64 {
    let f1 = if stress.sigma1 >= 0.0 {
        stress.sigma1 / lamina.xt
    } else {
        stress.sigma1.abs() / lamina.xc
    };
    let f2 = if stress.sigma2 >= 0.0 {
        stress.sigma2 / lamina.yt
    } else {
        stress.sigma2.abs() / lamina.yc
    };
    let f12 = stress.tau12.abs() / lamina.s12;

    f1.max(f2).max(f12)
}

/// Tsai-Hill failure criterion.
///
/// Returns failure index F (>= 1.0 means failure):
/// F = (σ1/X)² - σ1σ2/X² + (σ2/Y)² + (τ12/S)²
///
/// Uses tensile or compressive strength based on stress sign.
#[must_use]
#[inline]
pub fn tsai_hill_failure_index(stress: &PlyStress, lamina: &Lamina) -> f64 {
    let x = if stress.sigma1 >= 0.0 {
        lamina.xt
    } else {
        lamina.xc
    };
    let y = if stress.sigma2 >= 0.0 {
        lamina.yt
    } else {
        lamina.yc
    };

    if x.abs() < hisab::EPSILON_F64
        || y.abs() < hisab::EPSILON_F64
        || lamina.s12.abs() < hisab::EPSILON_F64
    {
        return f64::INFINITY;
    }

    let s1 = stress.sigma1;
    let s2 = stress.sigma2;
    let t = stress.tau12;

    (s1 / x).powi(2) - (s1 * s2) / (x * x) + (s2 / y).powi(2) + (t / lamina.s12).powi(2)
}

/// Tsai-Wu failure criterion.
///
/// Returns failure index F (>= 1.0 means failure):
/// F = F1*σ1 + F2*σ2 + F11*σ1² + F22*σ2² + F66*τ12² + 2*F12*σ1*σ2
///
/// Uses the default interaction term F12 = -0.5 * sqrt(F11 * F22).
#[must_use]
pub fn tsai_wu_failure_index(stress: &PlyStress, lamina: &Lamina) -> f64 {
    if lamina.xt.abs() < hisab::EPSILON_F64
        || lamina.xc.abs() < hisab::EPSILON_F64
        || lamina.yt.abs() < hisab::EPSILON_F64
        || lamina.yc.abs() < hisab::EPSILON_F64
        || lamina.s12.abs() < hisab::EPSILON_F64
    {
        return f64::INFINITY;
    }

    let f1 = 1.0 / lamina.xt - 1.0 / lamina.xc;
    let f2 = 1.0 / lamina.yt - 1.0 / lamina.yc;
    let f11 = 1.0 / (lamina.xt * lamina.xc);
    let f22 = 1.0 / (lamina.yt * lamina.yc);
    let f66 = 1.0 / (lamina.s12 * lamina.s12);
    let f12 = -0.5 * (f11 * f22).sqrt();

    let s1 = stress.sigma1;
    let s2 = stress.sigma2;
    let t = stress.tau12;

    f1 * s1 + f2 * s2 + f11 * s1 * s1 + f22 * s2 * s2 + f66 * t * t + 2.0 * f12 * s1 * s2
}

/// Transform laminate-coordinate stress to ply material coordinates.
///
/// Given stress in (x, y, xy) laminate frame and ply angle theta,
/// returns stress in (1, 2, 12) material frame.
#[must_use]
pub fn transform_stress_to_material(
    sigma_x: f64,
    sigma_y: f64,
    tau_xy: f64,
    theta: f64,
) -> PlyStress {
    let c = theta.cos();
    let s = theta.sin();
    let c2 = c * c;
    let s2 = s * s;
    let cs = c * s;

    PlyStress {
        sigma1: c2 * sigma_x + s2 * sigma_y + 2.0 * cs * tau_xy,
        sigma2: s2 * sigma_x + c2 * sigma_y - 2.0 * cs * tau_xy,
        tau12: -cs * sigma_x + cs * sigma_y + (c2 - s2) * tau_xy,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn carbon() -> Lamina {
        Lamina::carbon_epoxy()
    }

    // --- Lamina stiffness tests ---

    #[test]
    fn q_matrix_symmetry() {
        let q = carbon().stiffness_matrix();
        assert!((q[0][1] - q[1][0]).abs() < 1.0, "Q should be symmetric");
    }

    #[test]
    fn q_matrix_diagonal_positive() {
        let q = carbon().stiffness_matrix();
        for (i, row) in q.iter().enumerate() {
            assert!(row[i] > 0.0, "Q[{i}][{i}] should be positive");
        }
    }

    #[test]
    fn q11_dominated_by_e1() {
        let l = carbon();
        let q = l.stiffness_matrix();
        // Q11 ≈ E1 for small v12*v21
        assert!(
            (q[0][0] - l.e1).abs() / l.e1 < 0.05,
            "Q11 should be close to E1"
        );
    }

    #[test]
    fn nu21_reciprocity() {
        let l = carbon();
        // v21/E2 = v12/E1 (reciprocity)
        let lhs = l.nu21() / l.e2;
        let rhs = l.nu12 / l.e1;
        assert!(
            (lhs - rhs).abs() < hisab::EPSILON_F64,
            "reciprocity: v21/E2 = v12/E1"
        );
    }

    #[test]
    fn qbar_at_zero_equals_q() {
        let l = carbon();
        let q = l.stiffness_matrix();
        let qbar = l.transformed_stiffness(0.0);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (q[i][j] - qbar[i][j]).abs() < 1.0,
                    "Q-bar at 0° should equal Q at [{i}][{j}]"
                );
            }
        }
    }

    #[test]
    fn qbar_at_90_swaps_axes() {
        let l = carbon();
        let q = l.stiffness_matrix();
        let qbar = l.transformed_stiffness(PI / 2.0);
        // Q-bar_11 at 90° ≈ Q22 (fiber now along y)
        assert!(
            (qbar[0][0] - q[1][1]).abs() < 1e3,
            "Q-bar_11 at 90° should ≈ Q22"
        );
        assert!(
            (qbar[1][1] - q[0][0]).abs() < 1e3,
            "Q-bar_22 at 90° should ≈ Q11"
        );
    }

    #[test]
    fn qbar_symmetric() {
        let qbar = carbon().transformed_stiffness(PI / 4.0);
        for (i, row) in qbar.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    (val - qbar[j][i]).abs() < 1.0,
                    "Q-bar should be symmetric at [{i}][{j}]"
                );
            }
        }
    }

    #[test]
    fn qbar_45_has_coupling() {
        let qbar = carbon().transformed_stiffness(PI / 4.0);
        // At 45°, Q-bar_16 and Q-bar_26 should be non-zero
        assert!(
            qbar[0][2].abs() > 1e6,
            "Q-bar_16 at 45° should be significant"
        );
    }

    // --- ABD matrix tests ---

    #[test]
    fn abd_symmetric_laminate_zero_b() {
        // Symmetric [0/90]s layup: B should be approximately zero
        let l = carbon();
        let plies = vec![
            Ply {
                lamina: l.clone(),
                angle: 0.0,
                thickness: 0.000125,
            },
            Ply {
                lamina: l.clone(),
                angle: PI / 2.0,
                thickness: 0.000125,
            },
            Ply {
                lamina: l.clone(),
                angle: PI / 2.0,
                thickness: 0.000125,
            },
            Ply {
                lamina: l.clone(),
                angle: 0.0,
                thickness: 0.000125,
            },
        ];
        let abd = abd_matrix(&plies);
        for (i, row) in abd.b.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    val.abs() < 1e-3,
                    "B[{i}][{j}] should be ~0 for symmetric layup, got {val}"
                );
            }
        }
    }

    #[test]
    fn abd_a_matrix_positive_diagonal() {
        let l = carbon();
        let plies = vec![
            Ply {
                lamina: l.clone(),
                angle: 0.0,
                thickness: 0.000125,
            },
            Ply {
                lamina: l.clone(),
                angle: PI / 4.0,
                thickness: 0.000125,
            },
            Ply {
                lamina: l.clone(),
                angle: -PI / 4.0,
                thickness: 0.000125,
            },
            Ply {
                lamina: l.clone(),
                angle: PI / 2.0,
                thickness: 0.000125,
            },
        ];
        let abd = abd_matrix(&plies);
        for (i, row) in abd.a.iter().enumerate() {
            assert!(row[i] > 0.0, "A[{i}][{i}] should be positive");
        }
    }

    #[test]
    fn abd_d_matrix_positive_diagonal() {
        let l = carbon();
        let plies = vec![
            Ply {
                lamina: l.clone(),
                angle: 0.0,
                thickness: 0.000125,
            },
            Ply {
                lamina: l.clone(),
                angle: PI / 2.0,
                thickness: 0.000125,
            },
        ];
        let abd = abd_matrix(&plies);
        for (i, row) in abd.d.iter().enumerate() {
            assert!(row[i] > 0.0, "D[{i}][{i}] should be positive");
        }
    }

    #[test]
    fn abd_single_ply() {
        let l = carbon();
        let plies = vec![Ply {
            lamina: l.clone(),
            angle: 0.0,
            thickness: 0.001,
        }];
        let abd = abd_matrix(&plies);
        let q = l.stiffness_matrix();
        // A = Q * t for single ply
        assert!(
            (abd.a[0][0] - q[0][0] * 0.001).abs() < 1e3,
            "A11 should equal Q11*t for single ply"
        );
    }

    #[test]
    fn abd_thicker_laminate_stiffer() {
        let l = carbon();
        let thin = vec![Ply {
            lamina: l.clone(),
            angle: 0.0,
            thickness: 0.001,
        }];
        let thick = vec![Ply {
            lamina: l.clone(),
            angle: 0.0,
            thickness: 0.002,
        }];
        let abd_thin = abd_matrix(&thin);
        let abd_thick = abd_matrix(&thick);
        assert!(abd_thick.a[0][0] > abd_thin.a[0][0]);
    }

    // --- Failure criteria tests ---

    #[test]
    fn max_stress_no_failure() {
        let l = carbon();
        let s = PlyStress {
            sigma1: 100e6,
            sigma2: 10e6,
            tau12: 5e6,
        };
        let fi = max_stress_failure_index(&s, &l);
        assert!(fi < 1.0, "should not fail at low stress, got {fi}");
    }

    #[test]
    fn max_stress_fiber_failure() {
        let l = carbon();
        let s = PlyStress {
            sigma1: 2000e6,
            sigma2: 0.0,
            tau12: 0.0,
        };
        let fi = max_stress_failure_index(&s, &l);
        assert!(fi >= 1.0, "should fail in fiber tension");
    }

    #[test]
    fn max_stress_matrix_failure() {
        let l = carbon();
        let s = PlyStress {
            sigma1: 0.0,
            sigma2: 100e6,
            tau12: 0.0,
        };
        let fi = max_stress_failure_index(&s, &l);
        assert!(fi >= 1.0, "should fail in matrix tension (Yt=50 MPa)");
    }

    #[test]
    fn tsai_hill_no_failure() {
        let l = carbon();
        let s = PlyStress {
            sigma1: 100e6,
            sigma2: 10e6,
            tau12: 5e6,
        };
        let fi = tsai_hill_failure_index(&s, &l);
        assert!(fi < 1.0, "should not fail, got {fi}");
    }

    #[test]
    fn tsai_hill_fiber_failure() {
        let l = carbon();
        let s = PlyStress {
            sigma1: 1600e6,
            sigma2: 0.0,
            tau12: 0.0,
        };
        let fi = tsai_hill_failure_index(&s, &l);
        assert!(fi >= 1.0, "should fail in fiber direction");
    }

    #[test]
    fn tsai_wu_no_failure() {
        let l = carbon();
        let s = PlyStress {
            sigma1: 100e6,
            sigma2: 10e6,
            tau12: 5e6,
        };
        let fi = tsai_wu_failure_index(&s, &l);
        assert!(fi < 1.0, "should not fail, got {fi}");
    }

    #[test]
    fn tsai_wu_failure() {
        let l = carbon();
        let s = PlyStress {
            sigma1: 1600e6,
            sigma2: 0.0,
            tau12: 0.0,
        };
        let fi = tsai_wu_failure_index(&s, &l);
        assert!(fi >= 1.0, "should fail");
    }

    #[test]
    fn tsai_wu_accounts_for_sign() {
        let l = carbon();
        // Compressive strength is higher than tensile for transverse
        let tension = PlyStress {
            sigma1: 0.0,
            sigma2: 40e6,
            tau12: 0.0,
        };
        let compression = PlyStress {
            sigma1: 0.0,
            sigma2: -40e6,
            tau12: 0.0,
        };
        let fi_t = tsai_wu_failure_index(&tension, &l);
        let fi_c = tsai_wu_failure_index(&compression, &l);
        assert!(
            fi_t > fi_c,
            "tension should be closer to failure than compression (Yt < Yc)"
        );
    }

    // --- Stress transformation tests ---

    #[test]
    fn transform_zero_angle_identity() {
        let ps = transform_stress_to_material(100e6, 50e6, 20e6, 0.0);
        assert!((ps.sigma1 - 100e6).abs() < 1.0);
        assert!((ps.sigma2 - 50e6).abs() < 1.0);
        assert!((ps.tau12 - 20e6).abs() < 1.0);
    }

    #[test]
    fn transform_90_swaps() {
        let ps = transform_stress_to_material(100e6, 50e6, 0.0, PI / 2.0);
        assert!((ps.sigma1 - 50e6).abs() < 1e3, "σ1 at 90° should be σy");
        assert!((ps.sigma2 - 100e6).abs() < 1e3, "σ2 at 90° should be σx");
    }

    #[test]
    fn transform_preserves_invariant() {
        // σ1 + σ2 = σx + σy (first invariant of 2D stress)
        let ps = transform_stress_to_material(100e6, 50e6, 30e6, PI / 6.0);
        let sum_material = ps.sigma1 + ps.sigma2;
        let sum_laminate = 100e6 + 50e6;
        assert!(
            (sum_material - sum_laminate).abs() < 1e3,
            "stress invariant should be preserved"
        );
    }

    #[test]
    fn glass_epoxy_preset() {
        let g = Lamina::glass_epoxy();
        assert!(g.e1 > g.e2);
        assert!(g.xt > 0.0);
        assert!(g.s12 > 0.0);
    }
}
