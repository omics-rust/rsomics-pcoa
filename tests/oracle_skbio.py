#!/usr/bin/env python3
"""scikit-bio PCoA oracle for rsomics-pcoa compat tests.

Reads a square symmetric distance-matrix TSV (skbio DistanceMatrix form: an
empty top-left cell, IDs in the header, ID + tab-separated distances per row) on
argv[1], runs skbio.stats.ordination.pcoa with the default method='eigh', and
dumps eigenvalues, proportion_explained, and per-sample coordinates as
tab-separated rows the Rust compat harness parses for a value-level differential.
"""

import sys

import numpy as np
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa


def main():
    path = sys.argv[1]
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip() and not ln.startswith("#")]
    ids = lines[0].split("\t")[1:]
    rows = []
    for ln in lines[1:]:
        rows.append([float(v) for v in ln.split("\t")[1:]])
    dm = DistanceMatrix(np.array(rows, dtype=float), ids)
    res = pcoa(dm)

    def row(label, values):
        return label + "\t" + "\t".join(repr(float(v)) for v in values)

    out = [row("eigvals", res.eigvals.values),
           row("proportion_explained", res.proportion_explained.values)]
    for sid in res.samples.index:
        out.append(row(sid, res.samples.loc[sid].values))
    sys.stdout.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
