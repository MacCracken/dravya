use criterion::{Criterion, black_box, criterion_group, criterion_main};

fn bench_von_mises(c: &mut Criterion) {
    let s = dravya::stress::StressTensor::new(100e6, 50e6, 30e6, 20e6, 10e6, 5e6);
    c.bench_function("stress/von_mises", |b| {
        b.iter(|| black_box(s).von_mises());
    });
}

fn bench_principal_stresses(c: &mut Criterion) {
    let s = dravya::stress::StressTensor::new(100e6, 50e6, 30e6, 20e6, 10e6, 5e6);
    c.bench_function("stress/principal_stresses", |b| {
        b.iter(|| black_box(s).principal_stresses());
    });
}

fn bench_max_shear(c: &mut Criterion) {
    let s = dravya::stress::StressTensor::new(100e6, 50e6, 30e6, 20e6, 10e6, 5e6);
    c.bench_function("stress/max_shear", |b| {
        b.iter(|| black_box(s).max_shear());
    });
}

fn bench_j2(c: &mut Criterion) {
    let s = dravya::stress::StressTensor::new(100e6, 50e6, 30e6, 20e6, 10e6, 5e6);
    c.bench_function("stress/j2", |b| {
        b.iter(|| black_box(s).j2());
    });
}

fn bench_deviatoric(c: &mut Criterion) {
    let s = dravya::stress::StressTensor::new(100e6, 50e6, 30e6, 20e6, 10e6, 5e6);
    c.bench_function("stress/deviatoric", |b| {
        b.iter(|| black_box(s).deviatoric());
    });
}

fn bench_hookes_law(c: &mut Criterion) {
    c.bench_function("elastic/hookes_law", |b| {
        b.iter(|| dravya::elastic::hookes_law(black_box(200e9), black_box(0.001)));
    });
}

fn bench_bulk_modulus(c: &mut Criterion) {
    c.bench_function("elastic/bulk_modulus", |b| {
        b.iter(|| dravya::elastic::bulk_modulus(black_box(200e9), black_box(0.30)));
    });
}

fn bench_cantilever(c: &mut Criterion) {
    let i = dravya::beam::moment_of_inertia_rect(0.1, 0.01);
    c.bench_function("beam/cantilever_deflection", |b| {
        b.iter(|| {
            dravya::beam::cantilever_deflection(
                black_box(1000.0),
                black_box(1.0),
                black_box(200e9),
                black_box(i),
            )
        });
    });
}

fn bench_euler_buckling(c: &mut Criterion) {
    c.bench_function("beam/euler_buckling", |b| {
        let i = dravya::beam::moment_of_inertia_rect(0.05, 0.05);
        b.iter(|| {
            dravya::beam::euler_buckling_load(black_box(200e9), black_box(i), black_box(2.0))
        });
    });
}

fn bench_safety_factor(c: &mut Criterion) {
    let s = dravya::stress::StressTensor::new(100e6, 50e6, 30e6, 20e6, 10e6, 5e6);
    c.bench_function("yield/safety_factor", |b| {
        b.iter(|| dravya::yield_criteria::safety_factor(black_box(&s), black_box(250e6)));
    });
}

fn bench_basquin(c: &mut Criterion) {
    c.bench_function("fatigue/basquin_cycles", |b| {
        b.iter(|| {
            dravya::fatigue::basquin_cycles(black_box(200e6), black_box(1000e6), black_box(-0.1))
        });
    });
}

fn bench_miners_rule(c: &mut Criterion) {
    let loads: Vec<(f64, f64)> = (0..100).map(|i| (i as f64 * 10.0, 10000.0)).collect();
    c.bench_function("fatigue/miners_rule_100", |b| {
        b.iter(|| dravya::fatigue::miners_rule(black_box(&loads)));
    });
}

fn bench_goodman(c: &mut Criterion) {
    c.bench_function("fatigue/goodman_correction", |b| {
        b.iter(|| {
            dravya::fatigue::goodman_correction(
                black_box(100e6),
                black_box(100e6),
                black_box(400e6),
            )
        });
    });
}

fn bench_effective_strain(c: &mut Criterion) {
    let s = dravya::strain::StrainTensor::new(0.01, -0.003, -0.003, 0.005, 0.001, 0.002);
    c.bench_function("strain/effective_strain", |b| {
        b.iter(|| black_box(s).effective_strain());
    });
}

fn bench_stress_from_strain(c: &mut Criterion) {
    let steel = dravya::Material::steel();
    let eps = dravya::StrainTensor::new(0.001, -0.0003, -0.0003, 0.0005, 0.0, 0.0);
    c.bench_function("constitutive/stress_from_strain", |b| {
        b.iter(|| dravya::constitutive::stress_from_strain_3d(black_box(&steel), black_box(&eps)));
    });
}

fn bench_stiffness_matrix(c: &mut Criterion) {
    c.bench_function("constitutive/stiffness_matrix", |b| {
        b.iter(|| dravya::constitutive::stiffness_matrix(black_box(200e9), black_box(0.30)));
    });
}

fn bench_epp(c: &mut Criterion) {
    c.bench_function("constitutive/elastic_perfectly_plastic", |b| {
        b.iter(|| {
            dravya::constitutive::elastic_perfectly_plastic(
                black_box(200e9),
                black_box(250e6),
                black_box(0.005),
            )
        });
    });
}

fn bench_bilinear(c: &mut Criterion) {
    c.bench_function("constitutive/bilinear_hardening", |b| {
        b.iter(|| {
            dravya::constitutive::bilinear_hardening(
                black_box(200e9),
                black_box(250e6),
                black_box(20e9),
                black_box(0.005),
            )
        });
    });
}

fn bench_ramberg_osgood(c: &mut Criterion) {
    c.bench_function("constitutive/ramberg_osgood_strain", |b| {
        b.iter(|| {
            dravya::constitutive::ramberg_osgood_strain(
                black_box(200e9),
                black_box(1000e6),
                black_box(10.0),
                black_box(300e6),
            )
        });
    });
}

fn bench_ramberg_osgood_inverse(c: &mut Criterion) {
    c.bench_function("constitutive/ramberg_osgood_stress", |b| {
        b.iter(|| {
            dravya::constitutive::ramberg_osgood_stress(
                black_box(200e9),
                black_box(1000e6),
                black_box(10.0),
                black_box(0.003),
                1e-6,
                50,
            )
        });
    });
}

criterion_group!(
    benches,
    bench_von_mises,
    bench_principal_stresses,
    bench_max_shear,
    bench_j2,
    bench_deviatoric,
    bench_hookes_law,
    bench_bulk_modulus,
    bench_cantilever,
    bench_euler_buckling,
    bench_safety_factor,
    bench_basquin,
    bench_miners_rule,
    bench_goodman,
    bench_effective_strain,
    bench_stress_from_strain,
    bench_stiffness_matrix,
    bench_epp,
    bench_bilinear,
    bench_ramberg_osgood,
    bench_ramberg_osgood_inverse,
    bench_abd_matrix,
    bench_tsai_wu,
    bench_transform_ply_stress,
    bench_iso_hardening,
    bench_rainflow,
    bench_neuber_ro
);
criterion_main!(benches);

fn bench_abd_matrix(c: &mut Criterion) {
    let l = dravya::Lamina::carbon_epoxy();
    let plies = vec![
        dravya::Ply {
            lamina: l.clone(),
            angle: 0.0,
            thickness: 0.000125,
        },
        dravya::Ply {
            lamina: l.clone(),
            angle: std::f64::consts::FRAC_PI_4,
            thickness: 0.000125,
        },
        dravya::Ply {
            lamina: l.clone(),
            angle: -std::f64::consts::FRAC_PI_4,
            thickness: 0.000125,
        },
        dravya::Ply {
            lamina: l.clone(),
            angle: std::f64::consts::FRAC_PI_2,
            thickness: 0.000125,
        },
    ];
    c.bench_function("composite/abd_matrix_4ply", |b| {
        b.iter(|| dravya::abd_matrix(black_box(&plies)));
    });
}

fn bench_tsai_wu(c: &mut Criterion) {
    let l = dravya::Lamina::carbon_epoxy();
    let s = dravya::PlyStress {
        sigma1: 500e6,
        sigma2: 20e6,
        tau12: 30e6,
    };
    c.bench_function("composite/tsai_wu", |b| {
        b.iter(|| dravya::tsai_wu_failure_index(black_box(&s), black_box(&l)));
    });
}

fn bench_transform_ply_stress(c: &mut Criterion) {
    c.bench_function("composite/transform_stress", |b| {
        b.iter(|| {
            dravya::transform_stress_to_material(
                black_box(100e6),
                black_box(50e6),
                black_box(30e6),
                black_box(0.785),
            )
        });
    });
}

fn bench_iso_hardening(c: &mut Criterion) {
    let state = dravya::IsotropicHardening::new(250e6, 10e9);
    c.bench_function("constitutive/iso_hardening", |b| {
        b.iter(|| state.apply_uniaxial(black_box(200e9), black_box(0.005)));
    });
}

fn bench_rainflow(c: &mut Criterion) {
    let peaks: Vec<f64> = (0..100).map(|i| 100.0 * (i as f64 * 0.1).sin()).collect();
    c.bench_function("fatigue/rainflow_100", |b| {
        b.iter(|| dravya::fatigue::rainflow_count(black_box(&peaks)));
    });
}

fn bench_neuber_ro(c: &mut Criterion) {
    c.bench_function("fatigue/neuber_ramberg_osgood", |b| {
        b.iter(|| {
            dravya::fatigue::neuber_ramberg_osgood(
                black_box(2.5),
                black_box(150e6),
                black_box(200e9),
                black_box(1000e6),
                black_box(10.0),
                1e-6,
                50,
            )
        });
    });
}
