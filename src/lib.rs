//! # Dravya
//!
//! **Dravya** (द्रव्य — Sanskrit for "substance, matter") — material science engine
//! for the AGNOS ecosystem.
//!
//! Provides material properties, stress/strain tensors, elastic analysis, yield criteria,
//! beam mechanics, and fatigue life prediction. Built on [`hisab`] for math.

pub mod beam;
pub mod constitutive;
pub mod elastic;
pub mod error;
pub mod fatigue;
pub mod fracture;
pub mod material;
pub mod strain;
pub mod stress;
pub mod yield_criteria;

#[cfg(feature = "logging")]
pub mod logging;

pub use error::{DravyaError, Result};
pub use material::Material;
pub use strain::StrainTensor;
pub use stress::StressTensor;

// Elastic
pub use elastic::{
    bulk_modulus, hookes_law, lame_lambda, plane_strain_modulus, plane_stress_modulus,
    shear_modulus, strain_from_stress, youngs_from_bulk_shear, youngs_from_shear,
};

// Beam
pub use beam::{
    cantilever_deflection, cantilever_deflection_udl, euler_buckling_load,
    simply_supported_deflection, simply_supported_deflection_udl,
};

// Yield criteria
pub use yield_criteria::{safety_factor, tresca_check, von_mises_check};

// Fatigue
pub use fatigue::{basquin_cycles, endurance_limit_estimate, is_fatigue_failure, miners_rule};

// Strain
pub use strain::{engineering_strain, true_strain};

// Constitutive
pub use constitutive::{compliance_matrix, stiffness_matrix};

// Fracture
pub use fracture::{fracture_check, ki_center_crack_infinite, ki_edge_crack, paris_law_rate};
