"""Beta-binomial core of the EB score (the Genomon EBFilter / EBCall math, github.com/friend1ws/EBCall).

Per strand, a beta-binomial (a, b) is fitted to the PON (depth, alt) pairs. The fit maximises the log-likelihood
minus a penalty ``penalty * ln(a + b)`` with L-BFGS-B (start [20, 20], a in [0.1, 1e7], b in [1, 1e7]).
The tumor p-value per strand is the upper tail P(X >= alt | depth, a, b). Both strands are combined with Fisher's
method (chi2, 4 degrees of freedom), and EB = -log10(p), capped at 60.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import fmin_l_bfgs_b
from scipy.special import digamma, gammaln
from scipy.stats import chi2

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
    """Objective and its analytic gradient (digamma terms), for a converged, version-stable fit."""
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


def fit_bb(count_matrix: np.ndarray, penalty: float = 0.5, method: str = "exact") -> tuple[float, float]:
    """Penalised maximum-likelihood (a, b) for an N x 2 array of PON [depth, alt] pairs, rounded to 5 decimals.

    method="exact" (default): analytic gradient and tight tolerances, so the optimum does not depend on the scipy
    version. method="legacy": numerical gradient and default tolerances like the old code. For flat likelihoods
    (all-zero alt vectors) that stops early at implementation-dependent points: with scipy 1.9.3 it reproduces the
    legacy fits exactly, with scipy 1.18 only 19 of 102 test fits (relative differences up to 0.6).

    Legacy note: the old code ignored the configured penalty and always used 0.5 (matrix2AB.py:25).
    """
    cm = np.asarray(count_matrix)
    if method == "legacy":
        params = fmin_l_bfgs_b(_fit_objective, FIT_START, args=(cm, penalty), approx_grad=True, bounds=FIT_BOUNDS)[0]
    elif method == "exact":
        params = fmin_l_bfgs_b(
            _fit_objective_grad, FIT_START, args=(cm, penalty), bounds=FIT_BOUNDS, factr=10.0, pgtol=1e-10,
            maxiter=100000, maxfun=100000,
        )[0]
    else:
        raise ValueError(f"unknown fit method {method!r} (use one of {FIT_METHODS})")
    return round(float(params[0]), 5), round(float(params[1]), 5)


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
