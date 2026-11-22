"""Beta-binomial core of the EB score (the Genomon EBFilter / EBCall math, github.com/friend1ws/EBCall).

Per strand, a beta-binomial (a, b) is fitted to the PON (depth, alt) pairs. The fit maximises the log-likelihood
minus a penalty ``penalty * ln(a + b)`` with L-BFGS-B (start [20, 20], a in [0.1, 1e7], b in [1, 1e7]); the default
``exact`` method uses the analytic gradient and restarts until the point is stationary (see ``fit_bb``).
The tumor p-value per strand is the upper tail P(X >= alt | depth, a, b). Both strands are combined with Fisher's
method (chi2, 4 degrees of freedom), and EB = -log10(p), capped at 60.
"""

from __future__ import annotations

import logging
import math

import numpy as np
from scipy.optimize import fmin_l_bfgs_b
from scipy.special import digamma, gammaln
from scipy.stats import chi2

log = logging.getLogger("ebscore")

# gammaln arguments [n+1, k+1, n-k+1, n+a+b, k+a, n-k+b, a+b, a, b] and their signs: the log of
# C(n, k) * B(k + a, n - k + b) / B(a, b) for depth n and alt count k.
KS_MATRIX = np.array([[1, 0, 1, 1, 0, 1, 0, 0, 0], [0, 1, -1, 0, 1, -1, 0, 0, 0]])
GAMMA_REDUCE = np.array([1, -1, -1, -1, 1, 1, 1, -1, -1])

FIT_START = [20, 20]
FIT_BOUNDS = [(0.1, 10000000), (1, 10000000)]
EB_MAX = 60.0


def bb_loglikelihood(params, count_matrix: np.ndarray) -> float:
    """Summed beta-binomial log-likelihood of an N x 2 array of [depth, alt] pairs (or one pair, shape (2,)).

    Kept exactly as in the legacy code so that fitted parameters stay comparable.
    """
    a, b = params
    ab_matrix = np.array([1, 1, 1, a + b, a, b, a + b, a, b])
    input_matrix = np.matmul(count_matrix, KS_MATRIX) + ab_matrix
    gamma_matrix = (
        gammaln(input_matrix) if count_matrix.shape == (2,) else np.sum(gammaln(input_matrix), axis=0)
    )
    return float(np.sum(gamma_matrix * GAMMA_REDUCE))


def _fit_objective(params, count_matrix, penalty):
    return penalty * math.log(sum(params)) - bb_loglikelihood(params, count_matrix)


def _fit_objective_grad(params, count_matrix, penalty):
    """Objective and its analytic gradient (digamma terms), used by the exact fit."""
    a, b = params
    n = count_matrix[:, 0].astype(float)
    k = count_matrix[:, 1].astype(float)
    value = penalty * math.log(a + b) - bb_loglikelihood(params, count_matrix)
    common = np.sum(digamma(a + b) - digamma(n + a + b))
    d_a = np.sum(digamma(k + a) - digamma(a)) + common
    d_b = np.sum(digamma(n - k + b) - digamma(b)) + common
    pen = penalty / (a + b)
    return value, np.array([pen - d_a, pen - d_b])


FIT_METHODS = ("exact", "legacy")
# exact fit: settings of one L-BFGS-B run, the stationarity test and the restart budget
EXACT_OPTIONS = dict(factr=10.0, pgtol=1e-10, maxiter=100000, maxfun=100000)
STATIONARY_TOL = 1e-3  # max |projected gradient_i| * max(|x_i|, 1); converged fits are <= 1e-4, false stops >= 1
MAX_RESTARTS = 5
AB_DIGITS = 6  # exact fits are rounded to 6 significant digits (legacy: 5 decimals)


def projected_gradient(params, grad) -> np.ndarray:
    """Gradient with the components zeroed where a parameter sits at its bound and the gradient points outward."""
    pg = np.array(grad, dtype=float)
    for i, (x, (lo, hi)) in enumerate(zip(params, FIT_BOUNDS)):
        if (x <= lo and pg[i] > 0) or (x >= hi and pg[i] < 0):
            pg[i] = 0.0
    return pg


def _scaled_pg(x, grad) -> float:
    return float(np.max(np.abs(projected_gradient(x, grad)) * np.maximum(np.abs(np.asarray(x, dtype=float)), 1.0)))


def stationarity(params, count_matrix: np.ndarray, penalty: float = 0.5) -> float:
    """Scaled projected gradient of the fit objective: max_i |pg_i| * max(|x_i|, 1) (about 0 at an optimum).

    For x > 1 this is the gradient with respect to ln(x), so the test does not depend on the size of b (which grows
    with the total PON depth for all-zero alt vectors).
    """
    x = np.asarray(params, dtype=float)
    return _scaled_pg(x, _fit_objective_grad(x, np.asarray(count_matrix), penalty)[1])


def round_sig(x: float, digits: int = AB_DIGITS) -> float:
    """Round to ``digits`` significant digits."""
    return float(f"{x:.{digits}g}")


def _fit_legacy(cm, penalty) -> np.ndarray:
    return fmin_l_bfgs_b(_fit_objective, FIT_START, args=(cm, penalty), approx_grad=True, bounds=FIT_BOUNDS)[0]


def _fit_exact(cm, penalty) -> np.ndarray:
    """L-BFGS-B with the analytic gradient, restarted from its own result until the point is stationary.

    One L-BFGS-B run can report convergence ("REL_REDUCTION_OF_F") at a point whose gradient is far from 0, mostly
    on low-depth vectors with sparse alts, with any scipy version. A restart (fresh curvature memory) moves on to the
    optimum. If the point is still not stationary after MAX_RESTARTS, the better of it and the legacy fit is used.
    """
    x = np.array(FIT_START, dtype=float)
    for _ in range(1 + MAX_RESTARTS):
        x, _, info = fmin_l_bfgs_b(_fit_objective_grad, x, args=(cm, penalty), bounds=FIT_BOUNDS, **EXACT_OPTIONS)
        if _scaled_pg(x, info["grad"]) <= STATIONARY_TOL:
            return x
    log.debug("exact fit not stationary after %d restarts; comparing with the legacy fit", MAX_RESTARTS)
    old = _fit_legacy(cm, penalty)
    return old if _fit_objective(old, cm, penalty) < _fit_objective(x, cm, penalty) else x


def fit_bb(count_matrix: np.ndarray, penalty: float = 0.5, method: str = "exact") -> tuple[float, float]:
    """Penalised maximum-likelihood (a, b) for an N x 2 array of PON [depth, alt] pairs.

    method="exact" (default): analytic gradient, tight tolerances and a stationarity check with restarts, so the
    result is the penalised optimum (projected gradient about 0), rounded to 6 significant digits. It is not
    bit-identical across scipy versions: for all-zero alt vectors the objective is very flat in b (b* is roughly
    proportional to the total PON depth), so the last digits of b depend on the optimiser's stopping point. The EB
    impact of that is in the 3rd decimal at most.

    method="legacy": numerical gradient and default tolerances like the old code (= Genomon EBFilter), rounded to 5
    decimals. For flat likelihoods (all-zero alt vectors, i.e. most PON cells) that stops early, far from the optimum
    on realistic PONs, at scipy-version-dependent points: with scipy 1.9.3 it reproduces the legacy fits exactly,
    with scipy 1.18 only 19 of 102 test fits (relative differences up to 0.6).

    Legacy note: the old code ignored the configured penalty and always used 0.5 (matrix2AB.py:25).
    """
    cm = np.asarray(count_matrix)
    if method == "legacy":
        params = _fit_legacy(cm, penalty)
        return round(float(params[0]), 5), round(float(params[1]), 5)
    if method == "exact":
        params = _fit_exact(cm, penalty)
        return round_sig(float(params[0])), round_sig(float(params[1]))
    raise ValueError(f"unknown fit method {method!r} (use one of {FIT_METHODS})")


def bb_upper_tail(alt: int, depth: int, a: float, b: float) -> float:
    """P(X >= alt) for X ~ BetaBinomial(depth, a, b), summed over s = alt..depth (vectorised legacy bb_pvalue)."""
    n = float(depth)
    s = np.arange(alt, depth + 1, dtype=float)
    loglik = (
        gammaln(n + 1)
        - gammaln(s + 1)
        - gammaln(n - s + 1)
        - gammaln(n + a + b)
        + gammaln(s + a)
        + gammaln(n - s + b)
        + gammaln(a + b)
        - gammaln(a)
        - gammaln(b)
    )
    return float(np.exp(loglik).sum())


def fisher_combination(p_values, method: str = "sf") -> float:
    """Fisher's combined p-value.

    method="sf" uses the chi2 survival function. method="legacy" uses ``1 - chi2.cdf`` like the old code, which
    cannot go below about 2.2e-16, so legacy EB values jump from about 15.95 straight to 60.
    """
    if 0 in p_values:
        return 0.0
    stat = sum(-2 * math.log(p) for p in p_values)
    if method == "sf":
        return float(chi2.sf(stat, 2 * len(p_values)))
    if method == "legacy":
        return float(1 - chi2.cdf(stat, 4))
    raise ValueError(f"unknown Fisher method {method!r} (use 'sf' or 'legacy')")


def p_to_eb(p: float) -> float:
    """EB = -log10(p), rounded to 3 decimals; 60 for p < 1e-60 and 0 for p > 1 - 1e-10 (legacy caps)."""
    if p < 1e-60:
        return EB_MAX
    if p > 1.0 - 1e-10:
        return 0.0
    return -round(math.log10(p), 3)


def eb_score(tumor_alt, tumor_depth, ab_params, fisher: str = "sf") -> float:
    """EB for one variant.

    tumor_alt, tumor_depth: (plus, minus) counts; ab_params: ((a+, b+), (a-, b-)).
    """
    p_values = [
        bb_upper_tail(int(k), int(n), float(a), float(b))
        for k, n, (a, b) in zip(tumor_alt, tumor_depth, ab_params)
    ]
    return p_to_eb(fisher_combination(p_values, fisher))
