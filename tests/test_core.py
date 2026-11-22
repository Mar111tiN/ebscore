"""Beta-binomial math against the legacy implementation."""

import math

import numpy as np
import pandas as pd
import pytest
import scipy

from conftest import LEGACY
from ebscore.core import (
    _fit_objective,
    bb_loglikelihood,
    bb_upper_tail,
    eb_score,
    fisher_combination,
    fit_bb,
    p_to_eb,
)

FITS = pd.read_csv(LEGACY / "fits.tsv", sep="\t", dtype=str)
OLD_SCIPY = tuple(int(x) for x in scipy.__version__.split(".")[:2]) < (1, 15)


def matrix(depths, alts):
    return np.transpose([depths.split("|"), alts.split("|")]).astype(int)


@pytest.mark.skipif(not OLD_SCIPY, reason="legacy numerical-gradient fits are only reproducible with old scipy")
def test_legacy_fit_reproduces_old_values():
    for d, a, ab in FITS.itertuples(index=False):
        assert "|".join(map(str, fit_bb(matrix(d, a), method="legacy"))) == ab


def test_exact_fit_is_never_worse_than_legacy():
    for d, a, ab in FITS.itertuples(index=False):
        cm = matrix(d, a)
        old = tuple(float(x) for x in ab.split("|"))
        assert _fit_objective(fit_bb(cm), cm, 0.5) <= _fit_objective(old, cm, 0.5) + 1e-6


def test_exact_fit_respects_bounds_and_penalty():
    cm = matrix("30|40|50|60|70|80", "0|0|0|0|0|0")
    a, b = fit_bb(cm)
    assert a == pytest.approx(0.1)  # all-zero alts push a to its lower bound
    a2, b2 = fit_bb(cm, penalty=2.0)
    assert b2 < b  # a larger penalty on ln(a+b) shrinks b


def test_upper_tail_matches_legacy_loop():
    def legacy(alt, depth, a, b):
        obs = np.array([[depth, s] for s in range(alt, depth + 1)])
        return np.exp([bb_loglikelihood((a, b), o) for o in obs]).sum()

    for alt, depth, a, b in [(0, 0, 1, 2), (3, 40, 0.1, 13.1), (20, 400, 0.5, 28.3), (7, 7, 13.6, 94)]:
        assert bb_upper_tail(alt, depth, a, b) == pytest.approx(legacy(alt, depth, a, b), rel=1e-12, abs=1e-300)


def test_upper_tail_at_zero_alt_is_one():
    assert bb_upper_tail(0, 50, 0.3, 20) == pytest.approx(1.0)


def test_fisher_sf_equals_legacy_in_the_normal_range():
    for ps in [(0.5, 0.5), (0.01, 0.2), (1e-5, 1e-4)]:
        assert fisher_combination(ps, "sf") == pytest.approx(fisher_combination(ps, "legacy"), rel=1e-9)


def test_fisher_legacy_precision_ceiling():
    # 1 - cdf cannot go below ~2.2e-16, so legacy EB jumps from ~16 straight to 60
    ps = (1e-10, 1e-10)
    assert p_to_eb(fisher_combination(ps, "legacy")) == 60.0
    eb_sf = p_to_eb(fisher_combination(ps, "sf"))
    assert 16 < eb_sf < 60
    assert eb_sf == pytest.approx(-math.log10(fisher_combination(ps, "sf")), abs=1e-3)


def test_fisher_zero_p():
    assert fisher_combination((0.0, 0.3)) == 0.0


def test_p_to_eb_caps():
    assert p_to_eb(1e-70) == 60.0
    assert p_to_eb(1.0) == 0.0
    assert p_to_eb(0.001) == 3.0


def test_eb_score_no_alt_reads_is_zero():
    assert eb_score((0, 0), (30, 25), ((0.1, 10.0), (0.1, 12.0))) == 0.0
