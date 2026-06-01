use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

use criterion::{Criterion, criterion_group, criterion_main};

fn bench_pcoa(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-pcoa");
    let dm = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden/braycurtis_dm.tsv");
    c.bench_function("rsomics-pcoa braycurtis golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .arg(&dm)
                .args(["-t", "1"])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_pcoa);
criterion_main!(benches);
