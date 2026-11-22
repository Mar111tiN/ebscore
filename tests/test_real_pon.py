"""Opt-in checks on a real PONmatrix (study-derived data, never part of the repo).

Run: EBSCORE_REAL_PONMATRIX=<PONmatrix folder> pytest -q -rA tests/test_real_pon.py
Optional: EBSCORE_REAL_CHROM (default chr21), EBSCORE_REAL_ROWS (rows sampled, default 400; every 50th of the first
50 x ROWS matrix rows, read as a stream, so a mounted PON is never read in full).
"""

import gzip
import os
from pathlib import Path

import numpy as np
import pytest

from scipy.optimize import fmin_l_bfgs_b

from ebscore.core import FIT_BOUNDS, FIT_START, STATIONARY_TOL, _fit_objective, eb_score, fit_bb, stationarity
from ebscore.counts import split_cell

REAL = os.environ.get("EBSCORE_REAL_PONMATRIX")
CHROM = os.environ.get("EBSCORE_REAL_CHROM", "chr21")
ROWS = int(os.environ.get("EBSCORE_REAL_ROWS", "400"))
pytestmark = pytest.mark.skipif(not REAL, reason="set EBSCORE_REAL_PONMATRIX to a PONmatrix folder to run")


@pytest.fixture(scope="module")
def vectors():
    """Distinct per-strand (depth, alt) vectors of the non-reference letters of the sampled rows."""
    out = {}
    with gzip.open(Path(REAL) / "matrix" / f"{CHROM}.pon.gz", "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for n, line in enumerate(fh):
            if n >= 50 * ROWS:
                break
            if n % 50:
                continue
            f = line.rstrip("\n").split("\t")
            depth = split_cell(f[idx["Depth"]])
            if not any(depth[0]) and not any(depth[1]):
                continue
            for letter in "ACGT":
                if letter == f[idx["Ref"]].upper():
                    continue
                alt = split_cell(f[idx[letter]])
                for s in (0, 1):
                    cm = np.transpose([depth[s], alt[s]]).astype(int)
                    out[cm.tobytes()] = cm
    return list(out.values())


def _legacy_style_fit(cm):
    """Reference only (the package has no legacy fit since plan D19): numerical-gradient L-BFGS-B with default
    tolerances, rounded to 5 decimals, as the main code and Genomon EBFilter did."""
    x = fmin_l_bfgs_b(_fit_objective, FIT_START, args=(cm, 0.5), approx_grad=True, bounds=FIT_BOUNDS)[0]
    return round(float(x[0]), 5), round(float(x[1]), 5)


def test_exact_fit_is_stationary_and_never_worse_than_legacy(vectors):
    for cm in vectors:
        ex = fit_bb(cm)
        assert stationarity(ex, cm) <= STATIONARY_TOL, cm.tolist()
        assert _fit_objective(ex, cm, 0.5) <= _fit_objective(_legacy_style_fit(cm), cm, 0.5) + 1e-6


def test_report_eb_delta_exact_minus_legacy(vectors, capsys):
    # documents why EB thresholds set on legacy values do not carry over (plan D19); nothing to assert beyond finiteness
    for alt, depth in [((2, 2), (60, 60)), ((3, 3), (150, 150))]:
        d = []
        for cm in vectors:
            ex, lg = fit_bb(cm), _legacy_style_fit(cm)
            d.append((eb_score(alt, depth, (ex, ex)), eb_score(alt, depth, (lg, lg))))
        d = np.array(d)
        delta = d[:, 0] - d[:, 1]
        cross = int(((d[:, 0] >= 4) != (d[:, 1] >= 4)).sum())
        with capsys.disabled():
            print(f"\n{CHROM} {len(vectors)} strand vectors, tumor alt {alt} depth {depth}: EB exact - legacy "
                  f"median {np.median(delta):+.3f} p90 {np.quantile(delta, 0.9):+.3f} min {delta.min():+.3f} "
                  f"max {delta.max():+.3f}; crossing EB 4: {cross}/{len(vectors)}")
        assert np.isfinite(delta).all()
