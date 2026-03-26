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
    bench_effective_strain
);
criterion_main!(benches);
