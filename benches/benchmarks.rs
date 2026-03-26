use criterion::{black_box, criterion_group, criterion_main, Criterion};

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

fn bench_hookes_law(c: &mut Criterion) {
    c.bench_function("elastic/hookes_law", |b| {
        b.iter(|| dravya::elastic::hookes_law(black_box(200e9), black_box(0.001)));
    });
}

fn bench_cantilever(c: &mut Criterion) {
    let i = dravya::beam::moment_of_inertia_rect(0.1, 0.01);
    c.bench_function("beam/cantilever_deflection", |b| {
        b.iter(|| dravya::beam::cantilever_deflection(black_box(1000.0), black_box(1.0), black_box(200e9), black_box(i)));
    });
}

fn bench_bulk_modulus(c: &mut Criterion) {
    c.bench_function("elastic/bulk_modulus", |b| {
        b.iter(|| dravya::elastic::bulk_modulus(black_box(200e9), black_box(0.30)));
    });
}

criterion_group!(benches, bench_von_mises, bench_principal_stresses, bench_hookes_law, bench_cantilever, bench_bulk_modulus);
criterion_main!(benches);
