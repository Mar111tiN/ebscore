"""String formats shared by the PONmatrix, the AB cache and the EB output.

- PONmatrix / pile2count cell: ``c1|c2|...|cN=c1|c2|...|cN`` (+ strand per sample ``=`` - strand per sample)
- AB cache cell: ``a+|b+=a-|b-``
- Tumor (EB debug column): ``alt+=alt-<depth+=depth-``
- PON+ / PON- (EB output): ``alt1|...|altN<depth1|...|depthN`` for one strand
"""

from __future__ import annotations

STRANDSEP = "="
ADSEP = "<"
SAMPLESEP = "|"
LETTERS = "AGCTID"
MATRIX_COLUMNS = ["Chr", "Start", "Ref", *LETTERS, "Depth"]
AB_COLUMNS = ["Chr", "Start", "Ref", *LETTERS]


def split_cell(cell: str) -> tuple[list[int], list[int]]:
    """``"1|0=2|3"`` -> ([1, 0], [2, 3])"""
    plus, minus = cell.split(STRANDSEP)
    return [int(x) for x in plus.split(SAMPLESEP)], [int(x) for x in minus.split(SAMPLESEP)]


def n_samples(cell: str) -> int:
    return cell.split(STRANDSEP, 1)[0].count(SAMPLESEP) + 1


def drop_samples(values: list[int], exclude: set[int]) -> list[int]:
    """Drop 0-based sample indices."""
    return [v for i, v in enumerate(values) if i not in exclude] if exclude else values


def pon_string(alts: list[int], depths: list[int]) -> str:
    return SAMPLESEP.join(map(str, alts)) + ADSEP + SAMPLESEP.join(map(str, depths))


def tumor_string(alt: tuple[int, int], depth: tuple[int, int]) -> str:
    return f"{alt[0]}{STRANDSEP}{alt[1]}{ADSEP}{depth[0]}{STRANDSEP}{depth[1]}"


def ab_string(ab_plus: tuple[float, float], ab_minus: tuple[float, float]) -> str:
    return f"{ab_plus[0]}{SAMPLESEP}{ab_plus[1]}{STRANDSEP}{ab_minus[0]}{SAMPLESEP}{ab_minus[1]}"


def parse_ab(cell: str) -> tuple[tuple[float, float], tuple[float, float]]:
    """``"a+|b+=a-|b-"`` -> ((a+, b+), (a-, b-))"""
    plus, minus = cell.split(STRANDSEP)
    ap, bp = plus.split(SAMPLESEP)
    am, bm = minus.split(SAMPLESEP)
    return (float(ap), float(bp)), (float(am), float(bm))


def fit_key(depths: list[int], alts: list[int]) -> tuple[tuple[int, int], ...]:
    """Cache key of one fit: the (depth, alt) pairs in PON sample order.

    Identical vectors are fitted once; this replaces the legacy shared zero cache (``zerocache.py``). The pairs are
    deliberately not sorted: the likelihood is order-invariant in exact arithmetic, but the fit is not. For all-zero
    alt vectors the objective is so flat that L-BFGS-B with a numerical gradient stops at order-dependent points
    (seen: a = 0.1 vs 0.134, about 5 % of AB cells). Keeping the order reproduces the legacy fits exactly.
    """
    return tuple(zip(depths, alts))
