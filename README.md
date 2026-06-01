# rsomics-pcoa

**Principal Coordinates Analysis (PCoA)** of a symmetric distance matrix —
ordination of the matrix `rsomics-beta-diversity` (or any source) emits.

Reads a square symmetric distance-matrix TSV — an empty top-left cell, the IDs
in the header row, then one row per ID (ID followed by tab-separated distances),
which is scikit-bio's `DistanceMatrix` form — and writes the eigenvalues,
per-sample coordinates, and proportion explained for every principal coordinate
axis.

```
rsomics-pcoa dm.tsv
rsomics-beta-diversity counts.tsv | rsomics-pcoa
```

## Method

Classical (metric) PCoA, matching `skbio.stats.ordination.pcoa` with its default
`method='eigh'`:

1. Gower double-centering of `-D²/2`: `F = E - rowmean - colmean + grandmean`,
   `E = -D²/2`.
2. Symmetric eigendecomposition of `F`.
3. Eigenvalues sorted descending; eigenvalues within `numpy.isclose`'s tolerance
   of zero are set to exactly zero.
4. Non-positive eigenvalues (and their eigenvector columns) are zeroed — a
   distance matrix that is not Euclidean yields small negative eigenvalues, and
   PCoA keeps only the positive-eigenvalue axes.
5. `coordinates = eigenvectors · √eigenvalues` (a zeroed eigenvalue gives a zero
   axis).
6. `proportion_explained = eigenvalue / Σ(positive eigenvalues)`.

All `n` axes are reported (`PC1 … PCn`), including the zeroed ones, exactly as
scikit-bio does.

### Eigenvector sign

The sign of an eigenvector — and therefore of a coordinate axis — is arbitrary;
flipping a whole axis is an equally valid PCoA solution. Eigenvalues and
proportion explained are sign-independent. The compat differential compares
eigenvalues and proportions directly and compares coordinates up to a per-axis
sign flip (orienting each axis by its largest-magnitude entry on both sides).

### Output

A flat TSV: an `# eigenvalues` block (eigenvalue + proportion_explained rows
over `PC1 … PCn`) then a `# coordinates` block (one row per sample, ID followed
by its per-axis coordinate). Floats use Python's shortest round-trip `repr`.

## Origin

This crate is an independent Rust reimplementation of the PCoA operation
provided by `scikit-bio` (`skbio.stats.ordination.pcoa`, default `method='eigh'`,
which delegates the eigendecomposition to `scipy.linalg.eigh`), based on:

- Gower, J. C. (1966), *Some distance properties of latent root and vector
  methods used in multivariate analysis*, Biometrika 53(3-4):325-338,
  <https://doi.org/10.1093/biomet/53.3-4.325> — the classical-scaling method.
- The black-box behaviour of `skbio.stats.ordination.pcoa`: the `-D²/2`
  double-centering, descending eigenvalue order, the near-zero clipping and
  non-positive-eigenvalue zeroing, all-axes retention, and
  `proportion_explained` over the sum of the (post-zeroing) eigenvalues.

scikit-bio is BSD-3-Clause and was read and cited. The symmetric
eigendecomposition uses [`faer`](https://crates.io/crates/faer) (pure Rust, SIMD
+ rayon — external-dependency quadrant ①). Test fixtures are deterministically
generated distance matrices.

License: MIT OR Apache-2.0.
Upstream credit: scikit-bio <https://scikit-bio.org> (BSD-3-Clause).

## Compatibility & performance

`tests/compat.rs` runs this binary and the scikit-bio oracle
(`tests/oracle_skbio.py`) on the golden distance matrices and asserts the
eigenvalues + proportion explained match directly and the coordinates match up
to a per-axis sign flip (epsilon `1e-6`). The differential is skipped loudly
when scikit-bio is not importable.

The PCoA hot path is the `O(n³)` symmetric eigendecomposition; `faer` carries it
single-threaded and scales across cores.
