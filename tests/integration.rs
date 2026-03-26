use dravya::material::Material;
use dravya::stress::StressTensor;
use dravya::yield_criteria;
use dravya::*;

#[test]
fn steel_below_yield() {
    let steel = Material::steel();
    let s = StressTensor::uniaxial(200e6);
    assert!(!yield_criteria::von_mises_check(&s, steel.yield_strength));
}

#[test]
fn steel_above_yield() {
    let steel = Material::steel();
    let s = StressTensor::uniaxial(300e6);
    assert!(yield_criteria::von_mises_check(&s, steel.yield_strength));
}

#[test]
fn hookes_law_roundtrip() {
    let e = Material::steel().youngs_modulus;
    let stress = elastic::hookes_law(e, 0.001);
    let strain_back = elastic::strain_from_stress(e, stress);
    assert!((strain_back - 0.001).abs() < 1e-12);
}

#[test]
fn cantilever_vs_simply_supported() {
    let i = beam::moment_of_inertia_rect(0.1, 0.01);
    let e = Material::steel().youngs_modulus;
    let cant = beam::cantilever_deflection(1000.0, 1.0, e, i);
    let ss = beam::simply_supported_deflection(1000.0, 1.0, e, i);
    assert!(ss < cant, "simply supported should deflect less");
}

// --- Constitutive workflow ---

#[test]
fn stress_strain_3d_roundtrip() {
    let steel = Material::steel();
    let original = StressTensor::new(100e6, 50e6, 30e6, 20e6, 10e6, 5e6);
    let strain = constitutive::strain_from_stress_3d(&steel, &original);
    let recovered = constitutive::stress_from_strain_3d(&steel, &strain);
    for i in 0..6 {
        assert!(
            (original.components[i] - recovered.components[i]).abs() < 1.0,
            "3D roundtrip failed at component {i}"
        );
    }
}

#[test]
fn hardening_state_evolution() {
    let mut state = IsotropicHardening::new(250e6, 10e9);
    // Load in stages
    for strain in [0.002, 0.004, 0.006, 0.008] {
        let (stress, new_state) = state.apply_uniaxial(200e9, strain);
        assert!(stress > 0.0);
        state = new_state;
    }
    assert!(
        state.current_yield() > 250e6,
        "yield should expand after plastic loading"
    );
}

// --- Composite workflow ---

#[test]
fn composite_layup_and_failure_check() {
    let lamina = Lamina::carbon_epoxy();
    let plies = vec![
        Ply {
            lamina: lamina.clone(),
            angle: 0.0,
            thickness: 0.000125,
        },
        Ply {
            lamina: lamina.clone(),
            angle: std::f64::consts::FRAC_PI_2,
            thickness: 0.000125,
        },
        Ply {
            lamina: lamina.clone(),
            angle: std::f64::consts::FRAC_PI_2,
            thickness: 0.000125,
        },
        Ply {
            lamina: lamina.clone(),
            angle: 0.0,
            thickness: 0.000125,
        },
    ];
    let abd = abd_matrix(&plies);
    // Symmetric layup: B should be ~0
    for row in &abd.b {
        for &val in row {
            assert!(val.abs() < 1e-3);
        }
    }
    // Check a ply under low load doesn't fail
    let ps = transform_stress_to_material(100e6, 10e6, 5e6, 0.0);
    let fi = tsai_wu_failure_index(&ps, &lamina);
    assert!(fi < 1.0, "low-stress ply should not fail");
}

// --- Fracture workflow ---

#[test]
fn fracture_assessment_workflow() {
    // Steel plate with 10mm center crack under 100 MPa
    let ki = ki_center_crack_infinite(100e6, 0.01);
    let kic = 50e6; // 50 MPa√m toughness

    // Check fracture
    assert!(
        !fracture_check(ki, kic),
        "should not fracture at low stress"
    );

    // Find critical crack length
    let a_cr = critical_crack_length(100e6, kic);
    assert!(a_cr > 0.01, "critical crack should be larger than current");

    // J-integral
    let j = fracture::j_integral_mode_i(ki, 200e9, 0.30);
    assert!(j > 0.0);

    // Paris law life prediction
    let life = paris_law_life(1e-11, 3.0, 100e6, 0.005, 0.02, 1000);
    assert!(life > 0.0 && life.is_finite());
}

// --- Fatigue workflow ---

#[test]
fn fatigue_assessment_workflow() {
    let steel = Material::steel();

    // Endurance limit with Marin factors
    let se_prime = endurance_limit_estimate(steel.ultimate_tensile_strength);
    let ka = fatigue::marin_surface_factor(4.51, -0.265, steel.ultimate_tensile_strength);
    let kb = fatigue::marin_size_factor(0.025);
    let se = fatigue::marin_corrected_endurance(se_prime, &[ka, kb]);
    assert!(se < se_prime, "Marin factors should reduce endurance limit");

    // Goodman correction
    let sigma_ar = fatigue::goodman_correction(100e6, 50e6, steel.ultimate_tensile_strength);
    assert!(
        sigma_ar > 100e6,
        "mean stress should increase equivalent amplitude"
    );

    // Rainflow counting
    let peaks = vec![0.0, 200e6, 50e6, 180e6, 20e6, 0.0];
    let cycles = fatigue::rainflow_count(&peaks);
    assert!(!cycles.is_empty());

    // Miner's rule damage
    let damage: f64 = cycles
        .iter()
        .map(|(range, _, count)| {
            let n_f = fatigue::basquin_cycles(range / 2.0, 1000e6, -0.1);
            if n_f > 0.0 { count / n_f } else { 0.0 }
        })
        .sum();
    assert!(damage >= 0.0);
}

// --- Temperature-dependent material ---

#[test]
fn temp_dependent_material_workflow() {
    let mut m_hot = Material::steel();
    m_hot.youngs_modulus = 150e9;
    m_hot.yield_strength = 180e6;
    m_hot.ultimate_tensile_strength = 300e6;

    let tdm = TempDependentMaterial::new("Steel", vec![(293.0, Material::steel()), (800.0, m_hot)])
        .unwrap();

    let m_500 = tdm.at_temperature(500.0);
    assert!(m_500.youngs_modulus < Material::steel().youngs_modulus);
    assert!(m_500.youngs_modulus > 150e9);
}
