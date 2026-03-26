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
#[derive(Debug, Clone, PartialEq)]
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

/// Invert the ABD matrix to get the abd compliance (lowercase).
///
/// Assembles the full 6x6 ABD, inverts via hisab, and splits back into
/// 3x3 sub-matrices. Returns the inverse such that:
/// ```text
/// {ε⁰}   [a  b] {N}
/// {κ}   = [b  d] {M}
/// ```
///
/// Returns `None` if the matrix is singular.
#[must_use]
pub fn abd_inverse(abd: &AbdMatrix) -> Option<AbdMatrix> {
    // Assemble 6x6
    let mut full: Vec<Vec<f64>> = vec![vec![0.0; 6]; 6];
    for i in 0..3 {
        for j in 0..3 {
            full[i][j] = abd.a[i][j];
            full[i][j + 3] = abd.b[i][j];
            full[i + 3][j] = abd.b[i][j];
            full[i + 3][j + 3] = abd.d[i][j];
        }
    }

    let inv = hisab::num::matrix_inverse(&full).ok()?;

    let mut a = [[0.0; 3]; 3];
    let mut b = [[0.0; 3]; 3];
    let mut d = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            a[i][j] = inv[i][j];
            b[i][j] = inv[i][j + 3];
            d[i][j] = inv[i + 3][j + 3];
        }
    }
    Some(AbdMatrix { a, b, d })
}

// ---------------------------------------------------------------------------
// Ply failure criteria
// ---------------------------------------------------------------------------

/// Ply stress state in material coordinates (1-2 system).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
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
    if lamina.xt.abs() < hisab::EPSILON_F64
        || lamina.xc.abs() < hisab::EPSILON_F64
        || lamina.yt.abs() < hisab::EPSILON_F64
        || lamina.yc.abs() < hisab::EPSILON_F64
        || lamina.s12.abs() < hisab::EPSILON_F64
    {
        return f64::INFINITY;
    }
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
/// For custom interaction, use [`tsai_wu_failure_index_custom`].
#[must_use]
#[inline]
pub fn tsai_wu_failure_index(stress: &PlyStress, lamina: &Lamina) -> f64 {
    tsai_wu_failure_index_custom(stress, lamina, -0.5)
}

/// Tsai-Wu failure criterion with custom interaction parameter f*.
///
/// F12 = f* * sqrt(F11 * F22). Standard values:
/// - f* = -0.5 (default, most common)
/// - f* = 0.0 (conservative, no interaction)
/// - Must satisfy |f*| <= 1.0 for a closed failure envelope.
#[must_use]
pub fn tsai_wu_failure_index_custom(stress: &PlyStress, lamina: &Lamina, f_star: f64) -> f64 {
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
    let f12 = f_star * (f11 * f22).sqrt();

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

// ---------------------------------------------------------------------------
// Additional failure criteria
// ---------------------------------------------------------------------------

/// Allowable strain limits for maximum strain criterion.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct StrainAllowables {
    /// Longitudinal tensile strain limit.
    pub eps1t: f64,
    /// Longitudinal compressive strain limit (positive value).
    pub eps1c: f64,
    /// Transverse tensile strain limit.
    pub eps2t: f64,
    /// Transverse compressive strain limit (positive value).
    pub eps2c: f64,
    /// Shear strain limit.
    pub gamma12_max: f64,
}

/// Maximum strain failure criterion.
///
/// Checks each strain component against its respective allowable.
/// Returns failure index (>= 1.0 means failure).
///
/// Strain components (ε1, ε2, γ12) are in material coordinates.
#[must_use]
#[inline]
pub fn max_strain_failure_index(
    eps1: f64,
    eps2: f64,
    gamma12: f64,
    allow: &StrainAllowables,
) -> f64 {
    let f1 = if eps1 >= 0.0 {
        if allow.eps1t.abs() < hisab::EPSILON_F64 {
            f64::INFINITY
        } else {
            eps1 / allow.eps1t
        }
    } else if allow.eps1c.abs() < hisab::EPSILON_F64 {
        f64::INFINITY
    } else {
        eps1.abs() / allow.eps1c
    };
    let f2 = if eps2 >= 0.0 {
        if allow.eps2t.abs() < hisab::EPSILON_F64 {
            f64::INFINITY
        } else {
            eps2 / allow.eps2t
        }
    } else if allow.eps2c.abs() < hisab::EPSILON_F64 {
        f64::INFINITY
    } else {
        eps2.abs() / allow.eps2c
    };
    let f12 = if allow.gamma12_max.abs() < hisab::EPSILON_F64 {
        f64::INFINITY
    } else {
        gamma12.abs() / allow.gamma12_max
    };
    f1.max(f2).max(f12)
}

/// Hashin 2D failure criterion.
///
/// Distinguishes four failure modes: fiber tension, fiber compression,
/// matrix tension, matrix compression. Returns per-mode failure indices
/// as `HashinResult`.
#[must_use]
pub fn hashin_failure(stress: &PlyStress, lamina: &Lamina) -> HashinResult {
    let s1 = stress.sigma1;
    let s2 = stress.sigma2;
    let t12 = stress.tau12;

    // Fiber tension (σ1 >= 0): (σ1/Xt)² + (τ12/S12)²
    let fiber_tension = if s1 >= 0.0 {
        (s1 / lamina.xt).powi(2) + (t12 / lamina.s12).powi(2)
    } else {
        0.0
    };

    // Fiber compression (σ1 < 0): (σ1/Xc)²
    let fiber_compression = if s1 < 0.0 {
        (s1 / lamina.xc).powi(2)
    } else {
        0.0
    };

    // Matrix tension (σ2 >= 0): (σ2/Yt)² + (τ12/S12)²
    let matrix_tension = if s2 >= 0.0 {
        (s2 / lamina.yt).powi(2) + (t12 / lamina.s12).powi(2)
    } else {
        0.0
    };

    // Matrix compression (σ2 < 0): (σ2/(2*S23))² + [(Yc/(2*S23))²-1]*(σ2/Yc) + (τ12/S12)²
    // Simplified (assuming S23 ≈ Yc/2): (σ2/Yc)² + (τ12/S12)²
    let matrix_compression = if s2 < 0.0 {
        (s2 / lamina.yc).powi(2) + (t12 / lamina.s12).powi(2)
    } else {
        0.0
    };

    HashinResult {
        fiber_tension,
        fiber_compression,
        matrix_tension,
        matrix_compression,
    }
}

/// Result of Hashin failure analysis with per-mode failure indices.
///
/// Each index >= 1.0 indicates failure in that mode.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HashinResult {
    /// Fiber tension failure index.
    pub fiber_tension: f64,
    /// Fiber compression failure index.
    pub fiber_compression: f64,
    /// Matrix tension failure index.
    pub matrix_tension: f64,
    /// Matrix compression failure index.
    pub matrix_compression: f64,
}

impl HashinResult {
    /// Maximum failure index across all modes.
    #[must_use]
    #[inline]
    pub fn max_index(&self) -> f64 {
        self.fiber_tension
            .max(self.fiber_compression)
            .max(self.matrix_tension)
            .max(self.matrix_compression)
    }

    /// Whether any mode has failed.
    #[must_use]
    #[inline]
    pub fn is_failed(&self) -> bool {
        self.max_index() >= 1.0
    }
}

// ---------------------------------------------------------------------------
// Progressive failure analysis
// ---------------------------------------------------------------------------

/// Ply degradation factors for progressive failure.
///
/// When a ply fails, its stiffness is degraded (multiplied by these factors).
#[derive(Debug, Clone, Copy)]
pub struct DegradationFactors {
    /// Factor for E1 after fiber failure (typically 0.0-0.1).
    pub fiber: f64,
    /// Factor for E2, G12 after matrix failure (typically 0.0-0.2).
    pub matrix: f64,
}

impl Default for DegradationFactors {
    fn default() -> Self {
        Self {
            fiber: 0.07,
            matrix: 0.14,
        }
    }
}

/// Perform progressive failure analysis on a laminate under given stress.
///
/// Iteratively checks each ply for Hashin failure, degrades failed plies,
/// recomputes the ABD matrix, and repeats until no new failures occur
/// or the laminate fails catastrophically.
///
/// `laminate_stress` is (Nx, Ny, Nxy) in N/m.
///
/// Returns the degraded plies and the number of iterations taken.
/// If all plies fail, returns the fully degraded state.
#[must_use]
pub fn progressive_failure(
    plies: &[Ply],
    laminate_stress: [f64; 3],
    degradation: &DegradationFactors,
    max_iterations: usize,
) -> (Vec<Ply>, usize) {
    let mut current_plies: Vec<Ply> = plies.to_vec();
    let mut failed: Vec<bool> = vec![false; plies.len()];

    for iteration in 1..=max_iterations {
        let abd = abd_matrix(&current_plies);

        // Invert ABD to get strains from forces (simplified: use A-inverse only for membrane)
        let a_inv = invert_3x3(&abd.a);
        let a_inv = match a_inv {
            Some(inv) => inv,
            None => return (current_plies, iteration),
        };

        // Mid-plane strains: ε = A⁻¹ N
        let eps0 = mat3_vec3_mul(&a_inv, &laminate_stress);

        let mut new_failure = false;
        for (i, ply) in current_plies.iter_mut().enumerate() {
            if failed[i] {
                continue;
            }
            // Ply stress in material coordinates
            let qbar = ply.lamina.transformed_stiffness(ply.angle);
            let sig_x = qbar[0][0] * eps0[0] + qbar[0][1] * eps0[1] + qbar[0][2] * eps0[2];
            let sig_y = qbar[1][0] * eps0[0] + qbar[1][1] * eps0[1] + qbar[1][2] * eps0[2];
            let tau_xy = qbar[2][0] * eps0[0] + qbar[2][1] * eps0[1] + qbar[2][2] * eps0[2];

            let ps = transform_stress_to_material(sig_x, sig_y, tau_xy, ply.angle);
            let hashin = hashin_failure(&ps, &ply.lamina);

            if hashin.is_failed() {
                failed[i] = true;
                new_failure = true;
                // Degrade the failed ply
                if hashin.fiber_tension >= 1.0 || hashin.fiber_compression >= 1.0 {
                    ply.lamina.e1 *= degradation.fiber;
                    ply.lamina.e2 *= degradation.fiber;
                    ply.lamina.g12 *= degradation.fiber;
                }
                if hashin.matrix_tension >= 1.0 || hashin.matrix_compression >= 1.0 {
                    ply.lamina.e2 *= degradation.matrix;
                    ply.lamina.g12 *= degradation.matrix;
                }
            }
        }

        if !new_failure {
            return (current_plies, iteration);
        }
    }

    (current_plies, max_iterations)
}

/// Invert a 3x3 matrix (for ABD sub-matrix inversion).
fn invert_3x3(m: &[[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

    if det.abs() < hisab::EPSILON_F64 {
        return None;
    }

    let inv_det = 1.0 / det;
    Some([
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det,
        ],
    ])
}

/// Multiply 3x3 matrix by 3-vector.
fn mat3_vec3_mul(m: &[[f64; 3]; 3], v: &[f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
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

    // --- ABD inverse test ---

    #[test]
    fn abd_inverse_roundtrip() {
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
        let inv = abd_inverse(&abd).expect("should be invertible");
        // A * a should approximate identity (for symmetric laminate, B≈0)
        for i in 0..3 {
            for j in 0..3 {
                let mut product = 0.0;
                for k in 0..3 {
                    product += abd.a[i][k] * inv.a[k][j];
                }
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (product - expected).abs() < 0.01,
                    "A*a_inv should be identity at [{i}][{j}], got {product}"
                );
            }
        }
    }

    // --- Max strain failure test ---

    #[test]
    fn max_strain_no_failure() {
        let allow = StrainAllowables {
            eps1t: 0.01,
            eps1c: 0.008,
            eps2t: 0.005,
            eps2c: 0.02,
            gamma12_max: 0.02,
        };
        let fi = max_strain_failure_index(0.001, 0.0005, 0.001, &allow);
        assert!(fi < 1.0, "should not fail at low strain, got {fi}");
    }

    #[test]
    fn max_strain_fiber_failure() {
        let allow = StrainAllowables {
            eps1t: 0.01,
            eps1c: 0.008,
            eps2t: 0.005,
            eps2c: 0.02,
            gamma12_max: 0.02,
        };
        let fi = max_strain_failure_index(0.02, 0.0, 0.0, &allow);
        assert!(fi >= 1.0, "should fail in fiber tension");
    }

    // --- Hashin failure tests ---

    #[test]
    fn hashin_no_failure() {
        let l = carbon();
        let s = PlyStress {
            sigma1: 100e6,
            sigma2: 10e6,
            tau12: 5e6,
        };
        let h = hashin_failure(&s, &l);
        assert!(!h.is_failed(), "should not fail at low stress");
    }

    #[test]
    fn hashin_fiber_tension() {
        let l = carbon();
        let s = PlyStress {
            sigma1: 2000e6,
            sigma2: 0.0,
            tau12: 0.0,
        };
        let h = hashin_failure(&s, &l);
        assert!(h.fiber_tension >= 1.0, "should fail in fiber tension");
        assert!(
            h.fiber_compression < hisab::EPSILON_F64,
            "no compression failure"
        );
    }

    #[test]
    fn hashin_matrix_tension() {
        let l = carbon();
        let s = PlyStress {
            sigma1: 0.0,
            sigma2: 100e6,
            tau12: 0.0,
        };
        let h = hashin_failure(&s, &l);
        assert!(
            h.matrix_tension >= 1.0,
            "should fail in matrix tension (Yt=50 MPa)"
        );
    }

    #[test]
    fn hashin_fiber_compression() {
        let l = carbon();
        let s = PlyStress {
            sigma1: -1500e6,
            sigma2: 0.0,
            tau12: 0.0,
        };
        let h = hashin_failure(&s, &l);
        assert!(
            h.fiber_compression >= 1.0,
            "should fail in fiber compression (Xc=1200 MPa)"
        );
    }

    #[test]
    fn hashin_distinguishes_modes() {
        let l = carbon();
        // Pure transverse tension — should fail matrix, not fiber
        let s = PlyStress {
            sigma1: 0.0,
            sigma2: 100e6,
            tau12: 0.0,
        };
        let h = hashin_failure(&s, &l);
        assert!(h.matrix_tension > h.fiber_tension);
        assert!(h.matrix_tension > h.fiber_compression);
    }

    // --- Progressive failure test ---

    #[test]
    fn progressive_failure_low_load_no_damage() {
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
        let deg = DegradationFactors::default();
        let (result, iters) = progressive_failure(&plies, [1000.0, 500.0, 0.0], &deg, 10);
        // Low load should produce no failures
        assert_eq!(iters, 1, "should converge in 1 iteration (no failures)");
        // E1 should be unchanged
        assert!(
            (result[0].lamina.e1 - l.e1).abs() < 1.0,
            "no degradation expected"
        );
    }

    #[test]
    fn progressive_failure_high_load_degrades() {
        let l = carbon();
        let plies = vec![
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
        let deg = DegradationFactors::default();
        // High Nx load — 90° plies should fail in matrix (transverse to their fiber)
        let (result, _) = progressive_failure(&plies, [500_000.0, 0.0, 0.0], &deg, 10);
        // 90° plies (index 0, 3) should have degraded E2
        assert!(
            result[0].lamina.e2 < l.e2,
            "90° ply should have degraded E2"
        );
    }

    // --- Tsai-Wu custom f* test ---

    #[test]
    fn tsai_wu_custom_f_star_zero_more_conservative() {
        let l = carbon();
        let s = PlyStress {
            sigma1: 500e6,
            sigma2: 20e6,
            tau12: 30e6,
        };
        let fi_default = tsai_wu_failure_index(&s, &l);
        let fi_zero = tsai_wu_failure_index_custom(&s, &l, 0.0);
        // f*=0 is more conservative than f*=-0.5 under biaxial tension
        assert!(
            (fi_default - fi_zero).abs() > hisab::EPSILON_F64,
            "different f* should give different results"
        );
    }
}
