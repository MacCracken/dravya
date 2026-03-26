//! # Dravya
//!
//! **Dravya** (द्रव्य — Sanskrit for "substance, matter") — material science engine
//! for the AGNOS ecosystem.
//!
//! Provides material properties, stress/strain tensors, elastic analysis, yield criteria,
//! beam mechanics, and fatigue life prediction. Built on [`hisab`] for math.

pub mod error;
pub mod material;
pub mod stress;
pub mod strain;
pub mod elastic;
pub mod yield_criteria;
pub mod beam;
pub mod fatigue;

#[cfg(feature = "logging")]
pub mod logging;

pub use error::{DravyaError, Result};
pub use material::Material;
pub use stress::StressTensor;
pub use strain::StrainTensor;
pub use elastic::{hookes_law, bulk_modulus, shear_modulus};
pub use beam::{cantilever_deflection, simply_supported_deflection};
