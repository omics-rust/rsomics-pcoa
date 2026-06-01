use std::io::{BufRead, Write};

use faer::Mat;
use faer::linalg::solvers::SelfAdjointEigen;
use rsomics_common::{Result, RsomicsError};

mod fmt;
use fmt::push_pyrepr;

/// A square symmetric distance matrix in scikit-bio's `DistanceMatrix` (LSMat)
/// TSV form: an empty top-left cell then the IDs as the header row, then one row
/// per ID (ID + tab-separated distances). This is exactly what
/// `rsomics-beta-diversity` emits.
pub struct DistanceMatrix {
    pub ids: Vec<String>,
    /// Row-major `n × n`.
    pub data: Vec<f64>,
}

impl DistanceMatrix {
    /// # Errors
    /// Errors on a missing header, a non-square / ragged body, a row count that
    /// disagrees with the header, or a non-numeric cell.
    pub fn parse<R: BufRead>(reader: R, delim: char) -> Result<DistanceMatrix> {
        let mut lines = reader.lines();
        let header = loop {
            match lines.next() {
                Some(line) => {
                    let line = line.map_err(RsomicsError::Io)?;
                    if line.trim().is_empty() || line.starts_with('#') {
                        continue;
                    }
                    break line;
                }
                None => return Err(RsomicsError::InvalidInput("empty distance matrix".into())),
            }
        };
        let ids: Vec<String> = header
            .split(delim)
            .skip(1)
            .map(|s| s.trim().to_string())
            .collect();
        let n = ids.len();
        if n == 0 {
            return Err(RsomicsError::InvalidInput(
                "header has no ID columns (need an empty top-left cell + ≥1 ID)".into(),
            ));
        }

        let mut data = vec![0.0_f64; n * n];
        let mut row = 0usize;
        for line in lines {
            let line = line.map_err(RsomicsError::Io)?;
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            if row >= n {
                return Err(RsomicsError::InvalidInput(format!(
                    "more data rows than the {n} IDs in the header"
                )));
            }
            let mut fields = line.split(delim);
            let label = fields.next().unwrap_or("").trim();
            if label != ids[row] {
                return Err(RsomicsError::InvalidInput(format!(
                    "row {} label '{label}' does not match header ID '{}'",
                    row + 1,
                    ids[row]
                )));
            }
            let mut col = 0usize;
            for field in fields {
                if col >= n {
                    return Err(RsomicsError::InvalidInput(format!(
                        "row {} ('{label}') has more values than the {n} IDs",
                        row + 1
                    )));
                }
                data[row * n + col] = field.trim().parse().map_err(|_| {
                    RsomicsError::InvalidInput(format!(
                        "row {} ('{label}'), column {}: '{}' is not numeric",
                        row + 1,
                        col + 1,
                        field.trim()
                    ))
                })?;
                col += 1;
            }
            if col != n {
                return Err(RsomicsError::InvalidInput(format!(
                    "row {} ('{label}') has {col} values, expected {n}",
                    row + 1
                )));
            }
            row += 1;
        }
        if row != n {
            return Err(RsomicsError::InvalidInput(format!(
                "{row} data rows but {n} IDs in the header"
            )));
        }
        Ok(DistanceMatrix { ids, data })
    }

    #[must_use]
    pub fn n(&self) -> usize {
        self.ids.len()
    }
}

/// The result of a PCoA: every axis (one per sample) with its eigenvalue,
/// per-sample coordinate, and proportion explained, ordered by descending
/// eigenvalue. Matches `skbio.stats.ordination.pcoa` with the default
/// `method='eigh'` and all dimensions retained.
pub struct Ordination {
    pub ids: Vec<String>,
    pub eigvals: Vec<f64>,
    pub proportion_explained: Vec<f64>,
    /// Row-major `n_samples × n_axes`; `coordinates[s * n_axes + a]`.
    pub coordinates: Vec<f64>,
}

/// scikit-bio sets eigenvalues within this tolerance of zero to exactly zero
/// (`numpy.isclose` defaults) before deciding which are positive.
fn close_to_zero(x: f64) -> bool {
    const RTOL: f64 = 1e-5;
    const ATOL: f64 = 1e-8;
    (x - 0.0).abs() <= ATOL + RTOL * 0.0_f64.abs() || x.abs() <= ATOL
}

impl Ordination {
    /// PCoA via Gower double-centering then symmetric eigendecomposition.
    #[must_use]
    pub fn compute(dm: &DistanceMatrix) -> Ordination {
        let n = dm.n();
        let centered = gower_center(&dm.data, n);

        let mat = Mat::from_fn(n, n, |i, j| centered[i * n + j]);
        let eig: SelfAdjointEigen<f64> = mat.self_adjoint_eigen(faer::Side::Lower).unwrap();
        let s = eig.S();
        let u = eig.U();

        // faer (like LAPACK) returns ascending; PCoA wants descending.
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| s[b].partial_cmp(&s[a]).unwrap());

        let mut eigvals: Vec<f64> = order
            .iter()
            .map(|&k| {
                let v = s[k];
                if close_to_zero(v) { 0.0 } else { v }
            })
            .collect();

        let num_positive = eigvals.iter().filter(|&&v| v >= 0.0).count();
        for v in &mut eigvals[num_positive..] {
            *v = 0.0;
        }

        let sum_eig: f64 = eigvals.iter().sum();
        let proportion_explained: Vec<f64> = eigvals.iter().map(|&v| v / sum_eig).collect();

        // coordinates = eigvecs · sqrt(eigvals); a zeroed eigenvalue gives a
        // zero column, so the corresponding (dropped) eigenvector never shows.
        let sqrt_eig: Vec<f64> = eigvals.iter().map(|&v| v.sqrt()).collect();
        let mut coordinates = vec![0.0_f64; n * n];
        for sample in 0..n {
            for (axis, &k) in order.iter().enumerate() {
                coordinates[sample * n + axis] = u[(sample, k)] * sqrt_eig[axis];
            }
        }

        Ordination {
            ids: dm.ids.clone(),
            eigvals,
            proportion_explained,
            coordinates,
        }
    }

    /// Write a flat ordination TSV: an `# eigenvalues` block, then a sample
    /// coordinate table (ID + per-axis coordinate), then a
    /// `# proportion_explained` block. Floats use Python `repr` so the
    /// sign-independent columns are byte-comparable against the oracle.
    ///
    /// # Errors
    /// Propagates write errors.
    pub fn write_tsv<W: Write>(&self, mut out: W) -> Result<()> {
        let n_axes = self.eigvals.len();
        let mut line = String::new();

        writeln!(out, "# eigenvalues").map_err(RsomicsError::Io)?;
        write_axis_header(&mut out, n_axes)?;
        line.clear();
        line.push_str("eigval");
        for &v in &self.eigvals {
            line.push('\t');
            push_pyrepr(&mut line, v);
        }
        writeln!(out, "{line}").map_err(RsomicsError::Io)?;

        line.clear();
        line.push_str("proportion_explained");
        for &v in &self.proportion_explained {
            line.push('\t');
            push_pyrepr(&mut line, v);
        }
        writeln!(out, "{line}").map_err(RsomicsError::Io)?;

        writeln!(out, "# coordinates").map_err(RsomicsError::Io)?;
        write_axis_header(&mut out, n_axes)?;
        for (s, id) in self.ids.iter().enumerate() {
            line.clear();
            line.push_str(id);
            for axis in 0..n_axes {
                line.push('\t');
                push_pyrepr(&mut line, self.coordinates[s * n_axes + axis]);
            }
            writeln!(out, "{line}").map_err(RsomicsError::Io)?;
        }
        Ok(())
    }
}

fn write_axis_header<W: Write>(out: &mut W, n_axes: usize) -> Result<()> {
    let mut header = String::new();
    for a in 1..=n_axes {
        header.push('\t');
        header.push_str("PC");
        header.push_str(&a.to_string());
    }
    writeln!(out, "{header}").map_err(RsomicsError::Io)
}

/// Gower double-centering of `-D²/2`: `F = E - rowmean - colmean + grandmean`
/// where `E = -D²/2`.
fn gower_center(d: &[f64], n: usize) -> Vec<f64> {
    let mut e = vec![0.0_f64; n * n];
    for (idx, &dij) in d.iter().enumerate() {
        e[idx] = dij * dij / -2.0;
    }
    let inv_n = 1.0 / n as f64;
    let row_means: Vec<f64> = (0..n)
        .map(|i| e[i * n..i * n + n].iter().sum::<f64>() * inv_n)
        .collect();
    let col_means: Vec<f64> = (0..n)
        .map(|j| (0..n).map(|i| e[i * n + j]).sum::<f64>() * inv_n)
        .collect();
    let grand: f64 = e.iter().sum::<f64>() * inv_n * inv_n;
    for i in 0..n {
        for j in 0..n {
            e[i * n + j] = e[i * n + j] - row_means[i] - col_means[j] + grand;
        }
    }
    e
}

/// # Errors
/// Propagates parse and write errors.
pub fn run<R: BufRead, W: Write>(reader: R, out: W, delim: char) -> Result<()> {
    let dm = DistanceMatrix::parse(reader, delim)?;
    let ord = Ordination::compute(&dm);
    ord.write_tsv(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dm() -> &'static str {
        "\tA\tB\tC\tD\n\
         A\t0.0\t0.5\t0.8\t0.9\n\
         B\t0.5\t0.0\t0.6\t0.7\n\
         C\t0.8\t0.6\t0.0\t0.4\n\
         D\t0.9\t0.7\t0.4\t0.0\n"
    }

    #[test]
    fn parses_square_matrix() {
        let m = DistanceMatrix::parse(dm().as_bytes(), '\t').unwrap();
        assert_eq!(m.ids, ["A", "B", "C", "D"]);
        assert_eq!(m.data[1], 0.5);
        assert_eq!(m.data[2 * 4 + 3], 0.4);
    }

    #[test]
    fn eigenvalues_sum_equals_proportion_one() {
        let m = DistanceMatrix::parse(dm().as_bytes(), '\t').unwrap();
        let o = Ordination::compute(&m);
        let p: f64 = o.proportion_explained.iter().sum();
        assert!((p - 1.0).abs() < 1e-12, "proportions sum to {p}");
        assert_eq!(o.eigvals.len(), 4);
        // first axis is the largest
        assert!(o.eigvals[0] >= o.eigvals[1]);
    }

    #[test]
    fn coordinates_recover_squared_distance() {
        // Euclidean recovery: for a PCoA of a Euclidean distance matrix the
        // squared coordinate distance reproduces the input squared distance.
        let m = DistanceMatrix::parse(dm().as_bytes(), '\t').unwrap();
        let o = Ordination::compute(&m);
        let n = m.n();
        for i in 0..n {
            for j in 0..n {
                let mut s = 0.0;
                for a in 0..n {
                    let di = o.coordinates[i * n + a] - o.coordinates[j * n + a];
                    s += di * di;
                }
                let want = m.data[i * n + j] * m.data[i * n + j];
                assert!((s - want).abs() < 1e-6, "d²[{i}][{j}]: {s} vs {want}");
            }
        }
    }

    #[test]
    fn ragged_matrix_errors() {
        let bad = "\tA\tB\nA\t0.0\nB\t0.0\t0.0\n";
        assert!(DistanceMatrix::parse(bad.as_bytes(), '\t').is_err());
    }

    #[test]
    fn mismatched_label_errors() {
        let bad = "\tA\tB\nA\t0.0\t0.5\nX\t0.5\t0.0\n";
        assert!(DistanceMatrix::parse(bad.as_bytes(), '\t').is_err());
    }
}
