use dravya::*;
use dravya::material::Material;
use dravya::stress::StressTensor;
use dravya::strain;
use dravya::yield_criteria;

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
