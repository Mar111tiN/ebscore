"""Optional AB cache: beta-binomial parameters for every PONmatrix position, letter and strand.

Output ``AB/<chrom>.AB.gz``: ``Chr Start Ref A G C T I D``, each cell ``a+|b+=a-|b-``, sorted by Start (the legacy
``ebcache.PONmatrix2AB_multi`` format). The cache is only valid for the full PON; scoring with excluded PON samples
fits on the fly instead. Cost: 12 fits per covered position (legacy estimate ~600 M fits for a WES PON); exact
de-duplication of (depth, alt) vectors (``counts.fit_key``) replaces the legacy shared zero cache. Plan D-Q7 decides
whether production uses this cache at all.

Splitting (plan D16): the matrix is processed in row chunks (``chunk_rows``), each chunk fitted in parallel.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

import pandas as pd

from .counts import AB_COLUMNS, LETTERS, ab_string, fit_key, split_cell
from .fitting import fit_many
from .tools import ensure_dir

log = logging.getLogger("ebscore")


def _chunk_ab(chunk: pd.DataFrame, penalty: float, threads: int, fit_method: str) -> pd.DataFrame:
    depth = [split_cell(c) for c in chunk["Depth"]]
    keys = {}
    for letter in LETTERS:
        cells = [split_cell(c) for c in chunk[letter]]
        keys[letter] = [
            (fit_key(d[0], a[0]), fit_key(d[1], a[1])) for d, a in zip(depth, cells)
        ]
    fitted = fit_many([k for letter in LETTERS for pair in keys[letter] for k in pair], penalty, threads, fit_method)
    out = chunk[["Chr", "Start", "Ref"]].copy()
    for letter in LETTERS:
        out[letter] = [ab_string(fitted[kp], fitted[km]) for kp, km in keys[letter]]
    return out


def ab_cache(
    ponmatrix_dir, chrom: str, out=None, penalty: float = 0.5, threads: int = 1, chunk_rows: int = 50000,
    fit_method: str = "exact",
):
    """Build the AB cache of one chromosome from ``<ponmatrix_dir>/matrix/<chrom>.pon.gz``."""
    ponmatrix_dir = Path(ponmatrix_dir)
    matrix = ponmatrix_dir / "matrix" / f"{chrom}.pon.gz"
    if not matrix.is_file():
        raise FileNotFoundError(matrix)
    out = Path(out) if out else ponmatrix_dir / "AB" / f"{chrom}.AB.gz"
    ensure_dir(out.parent)
    tmp = out.with_suffix(".gz.part")
    first = True
    n_rows = 0
    if os.path.getsize(matrix) > 20:  # an empty gzip is about 20 bytes
        reader = pd.read_csv(matrix, sep="\t", dtype=str, keep_default_na=False, chunksize=chunk_rows)
        for chunk in reader:
            ab = _chunk_ab(chunk, penalty, threads, fit_method)
            ab.to_csv(tmp, sep="\t", index=False, header=first, mode="w" if first else "a", compression="gzip")
            first = False
            n_rows += len(ab)
            log.info("%s: %d positions done", chrom, n_rows)
    if first:  # no rows: write a header-only file
        pd.DataFrame(columns=AB_COLUMNS).to_csv(tmp, sep="\t", index=False, compression="gzip")
    os.replace(tmp, out)
    return out
