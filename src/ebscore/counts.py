"""String formats shared by the PONmatrix and the EB output.

- PONmatrix / pile2count cell: ``c1|c2|...|cN=c1|c2|...|cN`` (+ strand per sample ``=`` - strand per sample)
- AB (EB debug column, ``--keep-matrix``): ``a+|b+=a-|b-``
- Tumor (EB debug column): ``alt+=alt-<depth+=depth-``
- PON+ / PON- (EB output): ``alt1|...|altN<depth1|...|depthN`` for one strand
"""

from __future__ import annotations

from itertools import chain

STRANDSEP = "="
ADSEP = "<"
SAMPLESEP = "|"
LETTERS = "AGCTID"
MATRIX_COLUMNS = ["Chr", "Start", "Ref", *LETTERS, "Depth"]


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


def fit_key(depths: list[int], alts: list[int]) -> tuple[int, ...]:
    """De-duplication key of one fit: the (depth, alt) pairs in PON sample order, flattened (d1, a1, d2, a2, ...).

    Identical vectors are fitted once per run; this replaces the legacy shared zero cache (``zerocache.py``). A flat
    tuple of ints takes about a fifth of the memory of a tuple of pairs (``fitting`` reshapes it to N x 2). The pairs
    are not sorted, so a fit's input is exactly the PON vector in sample order (the result is then deterministic for a
    given PON; sorting would only add order-dependent last digits on flat objectives).
    """
    return tuple(chain.from_iterable(zip(depths, alts)))
