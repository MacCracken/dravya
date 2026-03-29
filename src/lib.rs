//! # Dravya
//!
//! **Dravya** (द्रव्य — Sanskrit for "substance, matter") — material science engine
//! for the AGNOS ecosystem.
//!
//! Provides material properties, stress/strain tensors, elastic analysis, yield criteria,
//! beam mechanics, fatigue life prediction, constitutive models, and fracture mechanics.
//! Built on [`hisab`] for math.

pub mod beam;
pub mod bridge;
pub mod composite;
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

#[cfg(feature = "thermal")]
pub mod thermal;

// Core types
pub use error::{DravyaError, Result};
pub use material::{Material, TempDependentMaterial};
pub use strain::StrainTensor;
pub use stress::StressTensor;

// Elastic
pub use elastic::{
    bulk_modulus, hookes_law, lame_lambda, p_wave_modulus, plane_strain_modulus,
    plane_stress_modulus, poisson_from_bulk_shear, poisson_from_youngs_shear, shear_modulus,
    strain_from_stress, try_bulk_modulus, try_shear_modulus, try_strain_from_stress,
    youngs_from_bulk_shear, youngs_from_shear,
};

// Beam
pub use beam::{
    angle_of_twist, bending_stress, cantilever_deflection, cantilever_deflection_udl,
    euler_buckling_load, fixed_fixed_deflection, moment_of_inertia_circle,
    moment_of_inertia_hollow_circle, moment_of_inertia_hollow_rect, moment_of_inertia_rect,
    polar_moment_circle, polar_moment_hollow_circle, section_modulus_circle, section_modulus_rect,
    shear_stress_beam, simply_supported_deflection, simply_supported_deflection_udl,
    torsional_stress,
};

// Yield criteria
pub use yield_criteria::{
    drucker_prager_check, drucker_prager_from_mohr_coulomb, safety_factor, safety_factor_tresca,
    tresca_check, von_mises_check,
};

// Fatigue
pub use fatigue::{
    basquin_cycles, basquin_cycles_reversals, coffin_manson_strain, coffin_manson_transition_life,
    endurance_limit_estimate, extract_turning_points, gerber_correction, goodman_correction,
    is_fatigue_failure, marin_corrected_endurance, marin_reliability_factor, marin_size_factor,
    marin_surface_factor, miners_rule, neuber_product, neuber_ramberg_osgood, rainflow_count,
    rainflow_count_periodic, sn_interpolate, soderberg_correction, stress_amplitude_mean,
    stress_ratio,
};

// Strain
pub use strain::{engineering_strain, true_strain, try_engineering_strain, try_true_strain};

// Constitutive
pub use constitutive::{
    CombinedHardening, IsotropicHardening, JohnsonCook, KinematicHardening, NeoHookean,
    Orthotropic3D, bilinear_hardening, compliance_matrix, elastic_perfectly_plastic,
    elastic_perfectly_plastic_material, ramberg_osgood_strain, ramberg_osgood_stress,
    stiffness_matrix, strain_from_stress_3d, stress_from_strain_3d,
};

// Composite
pub use composite::{
    AbdMatrix, DegradationFactors, HashinResult, Lamina, Ply, PlyStress, StrainAllowables,
    abd_inverse, abd_matrix, hashin_failure, max_strain_failure_index, max_stress_failure_index,
    progressive_failure, transform_stress_to_material, tsai_hill_failure_index,
    tsai_wu_failure_index, tsai_wu_failure_index_custom,
};

// Bridge
pub use bridge::{
    bone_density_to_yield_strength, bone_density_to_youngs_modulus, bone_safety_factor,
    collision_to_contact_stress, e_field_to_piezo_stress, grain_size_to_toughness_scale,
    impact_to_fracture_volume, magnetic_to_magnetostrictive_strain, mineral_fractions_to_density,
    muscle_force_to_tendon_stress, temperature_to_thermal_strain, tendon_strain_to_force,
    thermal_gradient_to_stress, velocity_to_strain_rate,
};
/// Integration APIs for downstream consumers (soorat rendering).
pub mod integration;

// Fracture
pub use fracture::{
    critical_crack_length, fracture_check, fracture_stress, j_integral_from_sifs,
    j_integral_mode_i, k_from_j_integral, ki_center_crack_finite, ki_center_crack_infinite,
    ki_crack_at_hole, ki_edge_crack, ki_penny_crack, kic_from_energy_release, kii_center_crack,
    kii_edge_crack, kiii_through_crack, paris_law_life, paris_law_rate,
};
