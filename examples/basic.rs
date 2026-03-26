use dravya::{material::Material, stress::StressTensor, elastic, beam, yield_criteria};

fn main() {
    let steel = Material::steel();
    println!("{}: E={:.0} GPa, σ_y={:.0} MPa", steel.name, steel.youngs_modulus / 1e9, steel.yield_strength / 1e6);

    // Uniaxial stress
    let s = StressTensor::uniaxial(150e6);
    println!("Applied: 150 MPa uniaxial");
    println!("Von Mises: {:.1} MPa", s.von_mises() / 1e6);
    println!("Safety factor: {:.2}", yield_criteria::safety_factor(&s, steel.yield_strength));

    // Elastic moduli
    let k = elastic::bulk_modulus(steel.youngs_modulus, steel.poisson_ratio);
    let g = elastic::shear_modulus(steel.youngs_modulus, steel.poisson_ratio);
    println!("Bulk modulus: {:.1} GPa", k / 1e9);
    println!("Shear modulus: {:.1} GPa", g / 1e9);

    // Beam deflection
    let i = beam::moment_of_inertia_rect(0.1, 0.01);
    let delta = beam::cantilever_deflection(1000.0, 1.0, steel.youngs_modulus, i);
    println!("Cantilever deflection (1kN, 1m, 100×10mm): {:.4} mm", delta * 1000.0);
}
