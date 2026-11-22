"""Optional AB cache: beta-binomial parameters for every PONmatrix position, letter and strand.

Output ``<ab_dir>/<chrom>.AB.gz``: ``Chr Start Ref A G C T I D``, each cell ``a+|b+=a-|b-``, sorted by Start (the
legacy ``ebcache.PONmatrix2AB_multi`` format), plus ``<ab_dir>/<chrom>.AB.json``, the manifest: penalty, fit method,
number of PON samples and the size of the matrix file it was fitted on. ``score`` uses a cache only if the manifest
matches the requested ``--penalty``/``--fit`` and the current matrix; a cache without manifest (legacy) counts as
penalty 0.5 + ``--fit legacy``. The cache is only valid for the full PON; scoring with excluded PON samples fits on
the fly instead. Cost: 12 fits per covered position (legacy estimate ~600 M fits for a WES PON); exact
de-duplication of (depth, alt) vectors (``counts.fit_key``) replaces the legacy shared zero cache. Plan D-Q7 decides
whether production uses this cache at all.

Splitting (plan D16): the matrix is processed in row chunks (``chunk_rows``), each chunk fitted in parallel. Keys
are flat int tuples, and the all-zero letters of a row share one key per strand: 5,000 rows x 100 normals take ~19 MB
of keys (the earlier tuple-of-pairs layout with 12 separate keys per row: ~440 MB), so the default 50,000-row chunk
needs ~0.2 GB for keys with zero alts, plus one key per strand for every non-zero letter.
"""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path

import pandas as pd

from . import __version__
from .counts import AB_COLUMNS, LETTERS, ab_string, fit_key, n_samples, split_cell
from .fitting import fit_many
from .tools import ensure_dir

log = logging.getLogger("ebscore")

MANIFEST_FORMAT = "ebscore-ab-cache/1"
LEGACY_MANIFEST = {"penalty": 0.5, "fit": "legacy"}  # what a cache without manifest was built with


def _chunk_ab(chunk: pd.DataFrame, penalty: float, threads: int, fit_method: str) -> pd.DataFrame:
    keys = {letter: [] for letter in LETTERS}
    for row in zip(chunk["Depth"], *(chunk[letter] for letter in LETTERS)):
        depth = split_cell(row[0])
        zero = [None, None]  # one shared key per strand for the all-zero letters of this row
        for letter, cell in zip(LETTERS, row[1:]):
            alt = split_cell(cell)
            pair = []
            for s in (0, 1):
                if any(alt[s]):
                    pair.append(fit_key(depth[s], alt[s]))
                else:
                    if zero[s] is None:
                        zero[s] = fit_key(depth[s], alt[s])
                    pair.append(zero[s])
            keys[letter].append(tuple(pair))
    fitted = fit_many([k for letter in LETTERS for pair in keys[letter] for k in pair], penalty, threads, fit_method)
    out = chunk[["Chr", "Start", "Ref"]].copy()
    for letter in LETTERS:
        out[letter] = [ab_string(fitted[kp], fitted[km]) for kp, km in keys[letter]]
    return out


def manifest_path(ab_dir, chrom: str) -> Path:
    return Path(ab_dir) / f"{chrom}.AB.json"


def read_manifest(ab_dir, chrom: str):
    """The manifest of ``<ab_dir>/<chrom>.AB.gz``, or None for a legacy cache without one."""
    p = manifest_path(ab_dir, chrom)
    return json.loads(p.read_text()) if p.is_file() else None


def cache_mismatch(manifest, penalty: float, fit_method: str, n_pon: int, matrix_file) -> list[str]:
    """Why a cache does not fit the requested scoring (empty list: it fits)."""
    m = manifest if manifest is not None else LEGACY_MANIFEST
    why = []
    if float(m["penalty"]) != float(penalty):
        why.append(f"penalty {m['penalty']} (requested {penalty})")
    if m["fit"] != fit_method:
        why.append(f"--fit {m['fit']} (requested {fit_method})")
    if manifest is not None:
        if n_pon and int(m["n_samples"]) != int(n_pon):
            why.append(f"{m['n_samples']} PON samples (matrix has {n_pon})")
        if matrix_file is not None and int(m["matrix_size"]) != os.path.getsize(matrix_file):
            why.append(f"matrix file of {m['matrix_size']} bytes (now {os.path.getsize(matrix_file)})")
    return why


def ab_cache(
    ponmatrix_dir, chrom: str, ab_dir, penalty: float = 0.5, threads: int = 1, chunk_rows: int = 50000,
    fit_method: str = "exact",
) -> Path:
    """Build ``<ab_dir>/<chrom>.AB.gz`` + manifest from ``<ponmatrix_dir>/matrix/<chrom>.pon.gz``."""
    matrix = Path(ponmatrix_dir) / "matrix" / f"{chrom}.pon.gz"
    if not matrix.is_file():
        raise FileNotFoundError(matrix)
    ab_dir = ensure_dir(ab_dir)
    out = ab_dir / f"{chrom}.AB.gz"
    tmp = ab_dir / f".{chrom}.AB.gz.part"
    first, n_rows, n_pon = True, 0, 0
    try:
        reader = pd.read_csv(matrix, sep="\t", dtype=str, keep_default_na=False, chunksize=chunk_rows)
        for chunk in reader:
            if chunk.empty:
                continue
            n_pon = n_pon or n_samples(chunk["Depth"].iloc[0])
            ab = _chunk_ab(chunk, penalty, threads, fit_method)
            ab.to_csv(tmp, sep="\t", index=False, header=first, mode="w" if first else "a", compression="gzip")
            first = False
            n_rows += len(ab)
            log.info("%s: %d positions done", chrom, n_rows)
    except pd.errors.EmptyDataError:
        pass  # a 0-byte matrix from an old build
    if first:  # no rows: write a header-only file
        pd.DataFrame(columns=AB_COLUMNS).to_csv(tmp, sep="\t", index=False, compression="gzip")
    os.replace(tmp, out)
    manifest = {
        "format": MANIFEST_FORMAT, "ebscore": __version__, "chrom": chrom, "penalty": penalty, "fit": fit_method,
        "n_samples": n_pon, "rows": n_rows, "matrix_file": matrix.name, "matrix_size": os.path.getsize(matrix),
    }
    manifest_path(ab_dir, chrom).write_text(json.dumps(manifest, indent=1) + "\n")
    return out
