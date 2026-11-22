"""Beta-binomial math against the legacy implementation and golden exact fits."""

import math

import numpy as np
import pandas as pd
import pytest
import scipy
from scipy.optimize import fmin_l_bfgs_b, minimize

import ebscore.core as core
from conftest import DATA, LEGACY
from ebscore.core import (
    EXACT_OPTIONS,
    FIT_BOUNDS,
    FIT_START,
    STATIONARY_TOL,
    _fit_exact,
    _fit_objective,
    _fit_objective_grad,
    _scaled_pg,
    bb_loglikelihood,
    bb_upper_tail,
    eb_score,
    fisher_combination,
    fit_bb,
    p_to_eb,
    round_sig,
    stationarity,
)

FITS = pd.read_csv(LEGACY / "fits.tsv", sep="\t", dtype=str)
HARD = pd.read_csv(DATA / "fits_hard.tsv", sep="\t", dtype=str)
GOLDEN = pd.read_csv(DATA / "fits_exact.tsv", sep="\t", dtype=str)
OLD_SCIPY = tuple(int(x) for x in scipy.__version__.split(".")[:2]) < (1, 15)


def matrix(depths, alts):
    return np.transpose([depths.split("|"), alts.split("|")]).astype(int)


def sparse_low_depth(seed, n):
    """Low-depth PON vectors with sparse alts (26-100 normals, Poisson depth 1-40, 0-3 alt samples)."""
    rng = np.random.default_rng(seed)
    for _ in range(n):
        ns = int(rng.integers(26, 101))
        depth = rng.poisson(rng.uniform(1, 40), size=ns)
        alt = np.zeros(ns, dtype=int)
        for i in rng.choice(ns, size=int(rng.integers(0, 4)), replace=False):
            if depth[i] > 0:
                alt[i] = int(rng.integers(1, min(3, depth[i]) + 1))
        yield np.transpose([depth, alt]).astype(int)


@pytest.mark.skipif(not OLD_SCIPY, reason="legacy numerical-gradient fits are only reproducible with old scipy")
def test_legacy_fit_reproduces_old_values():
    for d, a, ab in FITS.itertuples(index=False):
        assert "|".join(map(str, fit_bb(matrix(d, a), method="legacy"))) == ab


def test_legacy_fit_settings_reach_the_optimiser(monkeypatch):
    # scipy-independent check of --fit legacy: numerical gradient, legacy start/bounds, the given penalty
    calls = []

    def spy(*args, **kw):
        calls.append((args, kw))
        return fmin_l_bfgs_b(*args, **kw)

    monkeypatch.setattr(core, "fmin_l_bfgs_b", spy)
    cm = matrix("30|40|50|60|70|80", "0|0|1|0|0|0")
    fit_bb(cm, penalty=2.0, method="legacy")
    (args, kw), = calls
    assert args[0] is core._fit_objective and list(args[1]) == [20, 20]
    assert kw["approx_grad"] is True and kw["bounds"] == [(0.1, 10000000), (1, 10000000)]
    assert np.array_equal(kw["args"][0], cm)
    assert kw["args"][1] == 2.0
    assert set(kw) == {"args", "approx_grad", "bounds"}  # default tolerances, as in the old code


def test_legacy_fit_objective_close_to_stored_legacy_fits():
    # the legacy optimiser stops early at scipy-dependent points, but near the stored legacy objective
    for d, a, ab in FITS.itertuples(index=False):
        cm = matrix(d, a)
        old = tuple(float(x) for x in ab.split("|"))
        new = fit_bb(cm, method="legacy")
        assert new == tuple(round(x, 5) for x in new)  # legacy rounding: 5 decimals
        assert abs(_fit_objective(new, cm, 0.5) - _fit_objective(old, cm, 0.5)) < 0.15


def test_exact_fit_is_never_worse_than_legacy():
    for d, a, ab in FITS.itertuples(index=False):
        cm = matrix(d, a)
        old = tuple(float(x) for x in ab.split("|"))
        assert _fit_objective(fit_bb(cm), cm, 0.5) <= _fit_objective(old, cm, 0.5) + 1e-6


def test_exact_fit_matches_golden_fits():
    # fits_exact.tsv (tests/data/make_fits_exact.py, scipy 1.18.1): legacy fixture vectors, false-stop vectors,
    # large PONs (b in the thousands) and a germline-like vector (b at its lower bound). Across scipy 1.9.3/1.18 the
    # fits agree to 1.1e-5 relative; 1e-4 catches changed tolerances, bounds or penalty handling.
    for d, a, ab in GOLDEN.itertuples(index=False):
        cm = matrix(d, a)
        exp = tuple(float(x) for x in ab.split("|"))
        got = fit_bb(cm)
        assert got == pytest.approx(exp, rel=1e-4), (d, a)
    assert max(float(ab.split("|")[1]) for ab in GOLDEN["AB"]) > 5000  # the large-PON regime is covered
    assert "1577.08|1.0" in set(GOLDEN["AB"])  # b at its lower bound


def test_one_lbfgsb_run_can_stop_early_on_the_hard_vectors():
    # fits_hard.tsv holds vectors where a single exact L-BFGS-B run reported convergence at a non-stationary point
    # (3 with scipy 1.18.1, 6 others with scipy 1.9.3). This pins that the fixture still provokes the false stop.
    stops = 0
    for d, a in HARD.itertuples(index=False):
        cm = matrix(d, a)
        x, _, info = fmin_l_bfgs_b(_fit_objective_grad, FIT_START, args=(cm, 0.5), bounds=FIT_BOUNDS, **EXACT_OPTIONS)
        stops += _scaled_pg(x, info["grad"]) > 1.0
    if not stops:
        pytest.skip(f"no false stop with scipy {scipy.__version__}")
    assert stops >= 3


@pytest.mark.parametrize("source", ["hard", "seeded"])
def test_exact_fit_is_stationary_and_never_worse_than_legacy(source):
    vectors = [matrix(d, a) for d, a in HARD.itertuples(index=False)] if source == "hard" else sparse_low_depth(11, 600)
    for cm in vectors:
        ex = fit_bb(cm)
        assert stationarity(ex, cm) <= STATIONARY_TOL, cm.tolist()
        assert _fit_objective(ex, cm, 0.5) <= _fit_objective(fit_bb(cm, method="legacy"), cm, 0.5) + 1e-6


def test_exact_fit_falls_back_to_the_better_legacy_point(monkeypatch):
    # if restarts never reach a stationary point, the better of the exact and the legacy point is returned
    monkeypatch.setattr(core, "STATIONARY_TOL", -1.0)
    cm = matrix(*HARD.iloc[0])
    x = _fit_exact(cm, 0.5)
    legacy = core._fit_legacy(cm, 0.5)
    assert _fit_objective(x, cm, 0.5) <= _fit_objective(legacy, cm, 0.5)


def test_exact_fit_rounds_to_6_significant_digits():
    for d, a, _ in GOLDEN.head(20).itertuples(index=False):
        cm = matrix(d, a)
        assert fit_bb(cm) == tuple(round_sig(float(x), 6) for x in _fit_exact(cm, 0.5))
    assert round_sig(151896.3922) == 151896.0
    assert round_sig(7.278163) == 7.27816
    assert round_sig(0.000123456789) == 0.000123457


def test_exact_fit_respects_bounds_and_penalty():
    cm = matrix("30|40|50|60|70|80", "0|0|0|0|0|0")
    a, b = fit_bb(cm)
    assert a == pytest.approx(0.1)  # all-zero alts push a to its lower bound
    a2, b2 = fit_bb(cm, penalty=2.0)
    assert b2 < b  # a larger penalty on ln(a+b) shrinks b


PEN_DEPTHS = "70|34|38|41|39|70|73|59|31|34|46|51|61|53|43|37|64|66|31|35|52|49|74|55|51|51"  # 26 normals


@pytest.mark.parametrize(
    "alts, expected",
    [
        ("0|" * 25 + "0", {0.5: (0.1, 209.562), 2.0: (0.1, 16.8102)}),
        ("2|0|1|" + "0|" * 22 + "0", {0.5: (0.1, 31.8602), 2.0: (0.1, 8.08959)}),
        ("4|3|1|4|4|5|4|2|2|4|1|3|5|2|4|3|1|5|2|3|4|1|2|5|3|4", None),
    ],
)
def test_exact_fit_penalty_matches_an_independent_optimiser(alts, expected):
    # Nelder-Mead on the objective value alone (no gradient), started off the optimum, must not find a better point;
    # this pins the penalty in both the objective and its analytic gradient. Values: scipy 1.18.1.
    cm = matrix(PEN_DEPTHS, alts)
    for penalty in (0.5, 2.0):
        got = fit_bb(cm, penalty=penalty)
        for scale in (1.3, 0.7):
            nm = minimize(_fit_objective, [got[0] * scale + 0.05, got[1] * scale + 1], args=(cm, penalty),
                          method="Nelder-Mead", bounds=FIT_BOUNDS,
                          options=dict(xatol=1e-9, fatol=1e-13, maxiter=20000, maxfev=40000))
            assert _fit_objective(got, cm, penalty) <= nm.fun + 1e-7
            assert got == pytest.approx(tuple(nm.x), rel=1e-4)
        if expected:
            assert got == pytest.approx(expected[penalty], rel=1e-4)
    assert fit_bb(cm, penalty=2.0)[1] < fit_bb(cm, penalty=0.5)[1]


def test_upper_tail_matches_legacy_loop():
    def legacy(alt, depth, a, b):
        obs = np.array([[depth, s] for s in range(alt, depth + 1)])
        return np.exp([bb_loglikelihood((a, b), o) for o in obs]).sum()

    for alt, depth, a, b in [(0, 0, 1, 2), (3, 40, 0.1, 13.1), (20, 400, 0.5, 28.3), (7, 7, 13.6, 94)]:
        assert bb_upper_tail(alt, depth, a, b) == pytest.approx(legacy(alt, depth, a, b), rel=1e-12, abs=1e-300)


def test_upper_tail_at_zero_alt_is_one():
    assert bb_upper_tail(0, 50, 0.3, 20) == pytest.approx(1.0)


def test_fisher_sf_equals_legacy_in_the_normal_range():
    # identical below EB ~7 (relative difference < 1e-9) ...
    for ps in [(0.5, 0.5), (0.01, 0.2), (1e-5, 1e-4), (3e-4, 3e-4)]:
        assert fisher_combination(ps, "sf") == pytest.approx(fisher_combination(ps, "legacy"), rel=1e-9)
        assert p_to_eb(fisher_combination(ps, "sf")) < 8
    # ... within 0.005 EB below 15, diverging near the legacy ceiling (~15.95)
    for p in np.logspace(-3, -7.4, 60):
        eb_sf, eb_leg = p_to_eb(fisher_combination((p, p), "sf")), p_to_eb(fisher_combination((p, p), "legacy"))
        if eb_sf < 15:
            assert abs(eb_sf - eb_leg) <= 0.005, p


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
    # boundaries of the legacy caps: 60 below 1e-60, 0 above 1 - 1e-10
    assert p_to_eb(1e-55) == 55.0
    assert p_to_eb(1.1e-60) == pytest.approx(59.959)
    assert p_to_eb(9.9e-61) == 60.0
    assert math.copysign(1.0, p_to_eb(1 - 1e-11)) == 1.0  # the cap gives +0.0, the formula would give -0.0
    assert p_to_eb(1 - 1e-9) == -round(math.log10(1 - 1e-9), 3)


def test_eb_score_no_alt_reads_is_zero():
    assert eb_score((0, 0), (30, 25), ((0.1, 10.0), (0.1, 12.0))) == 0.0
