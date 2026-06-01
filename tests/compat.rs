use std::collections::HashMap;
use std::path::PathBuf;
use std::process::Command;

const EPS: f64 = 1e-6;

fn ours_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-pcoa"))
}

fn golden(name: &str) -> String {
    format!("{}/tests/golden/{}", env!("CARGO_MANIFEST_DIR"), name)
}

fn oracle_script() -> String {
    format!("{}/tests/oracle_skbio.py", env!("CARGO_MANIFEST_DIR"))
}

/// scikit-bio is the named oracle; skip loudly if it (or python) is unavailable.
/// `RSOMICS_SKBIO_PYTHON` overrides the interpreter (e.g. an isolated venv).
fn skbio_python() -> Option<String> {
    let mut candidates = Vec::new();
    if let Ok(p) = std::env::var("RSOMICS_SKBIO_PYTHON") {
        candidates.push(p);
    }
    candidates.push("python3".into());
    candidates.push("python".into());
    for py in candidates {
        let probe = Command::new(&py)
            .args(["-c", "import skbio.stats.ordination"])
            .output();
        if let Ok(out) = probe
            && out.status.success()
        {
            return Some(py);
        }
    }
    eprintln!("SKIP: scikit-bio not importable — install `scikit-bio` to run the differential");
    None
}

/// Parsed PCoA: the eigenvalue and proportion rows plus per-sample coordinate
/// vectors, keyed by sample ID.
struct Pcoa {
    eigvals: Vec<f64>,
    proportion: Vec<f64>,
    coords: HashMap<String, Vec<f64>>,
}

fn parse_ours(text: &str) -> Pcoa {
    let mut eigvals = Vec::new();
    let mut proportion = Vec::new();
    let mut coords = HashMap::new();
    let mut in_coords = false;
    for line in text.lines() {
        if line == "# coordinates" {
            in_coords = true;
            continue;
        }
        if line.starts_with('#') || line.starts_with('\t') {
            continue;
        }
        let mut it = line.split('\t');
        let label = it.next().unwrap();
        let vals: Vec<f64> = it.map(|s| s.parse().unwrap()).collect();
        match label {
            "eigval" => eigvals = vals,
            "proportion_explained" if !in_coords => proportion = vals,
            _ if in_coords => {
                coords.insert(label.to_string(), vals);
            }
            _ => {}
        }
    }
    Pcoa {
        eigvals,
        proportion,
        coords,
    }
}

fn parse_oracle(text: &str) -> Pcoa {
    let mut eigvals = Vec::new();
    let mut proportion = Vec::new();
    let mut coords = HashMap::new();
    for line in text.lines() {
        let mut it = line.split('\t');
        let label = it.next().unwrap();
        let vals: Vec<f64> = it.map(|s| s.parse().unwrap()).collect();
        match label {
            "eigvals" => eigvals = vals,
            "proportion_explained" => proportion = vals,
            _ => {
                coords.insert(label.to_string(), vals);
            }
        }
    }
    Pcoa {
        eigvals,
        proportion,
        coords,
    }
}

fn ours_output(table: &str) -> String {
    let out = Command::new(ours_bin())
        .arg(golden(table))
        .output()
        .expect("run rsomics-pcoa");
    assert!(
        out.status.success(),
        "ours failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

fn oracle_output(py: &str, table: &str) -> String {
    let out = Command::new(py)
        .arg(oracle_script())
        .arg(golden(table))
        .output()
        .expect("run scikit-bio oracle");
    assert!(
        out.status.success(),
        "oracle failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

fn approx(a: f64, b: f64) -> bool {
    (a - b).abs() <= EPS + EPS * b.abs()
}

fn differential(table: &str) {
    let Some(py) = skbio_python() else { return };
    let ours = parse_ours(&ours_output(table));
    let theirs = parse_oracle(&oracle_output(&py, table));

    assert_eq!(
        ours.eigvals.len(),
        theirs.eigvals.len(),
        "axis count differs"
    );
    for (a, &o) in ours.eigvals.iter().enumerate() {
        assert!(
            approx(o, theirs.eigvals[a]),
            "{table} eigval PC{} {o} vs {}",
            a + 1,
            theirs.eigvals[a]
        );
    }
    for (a, &o) in ours.proportion.iter().enumerate() {
        assert!(
            approx(o, theirs.proportion[a]),
            "{table} proportion PC{} {o} vs {}",
            a + 1,
            theirs.proportion[a]
        );
    }

    // Eigenvector sign is arbitrary, so coordinates are compared up to a
    // per-axis sign flip: pick the sign from the largest-magnitude entry of
    // each axis (consistent across both sides) and align before diffing.
    let n_axes = ours.eigvals.len();
    let ids: Vec<&String> = ours.coords.keys().collect();
    for a in 0..n_axes {
        let sign_ours = axis_sign(&ours.coords, &ids, a);
        let sign_theirs = axis_sign(&theirs.coords, &ids, a);
        for id in &ids {
            let o = ours.coords[*id][a] * sign_ours;
            let t = theirs.coords[*id][a] * sign_theirs;
            assert!(
                approx(o, t),
                "{table} coord {id} PC{} {o} vs {t} (sign-aligned)",
                a + 1
            );
        }
    }
}

/// The sign of the largest-magnitude coordinate on axis `a`, giving a stable
/// per-axis orientation independent of the eigensolver's arbitrary sign.
fn axis_sign(coords: &HashMap<String, Vec<f64>>, ids: &[&String], a: usize) -> f64 {
    let mut best = 0.0_f64;
    let mut sign = 1.0_f64;
    for id in ids {
        let v = coords[*id][a];
        if v.abs() > best {
            best = v.abs();
            sign = if v < 0.0 { -1.0 } else { 1.0 };
        }
    }
    sign
}

#[test]
fn matches_skbio_braycurtis_dm() {
    differential("braycurtis_dm.tsv");
}

#[test]
fn matches_skbio_euclidean_dm() {
    differential("euclidean_dm.tsv");
}
