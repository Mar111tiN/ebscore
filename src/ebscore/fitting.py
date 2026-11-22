"""Parallel, de-duplicated beta-binomial fits."""

from __future__ import annotations

import logging
from functools import partial
from multiprocessing import get_context

import numpy as np

from .core import fit_bb

log = logging.getLogger("ebscore")


def _fit_chunk(keys, penalty, method):
    return [fit_bb(np.array(key, dtype=int).reshape(-1, 2), penalty, method) for key in keys]


def fit_many(keys, penalty: float = 0.5, threads: int = 1, method: str = "exact", chunk: int = 200) -> dict:
    """Fit every distinct key ((depth, alt) pairs, see ``counts.fit_key``) once; return {key: (a, b)}."""
    unique = list(dict.fromkeys(keys))
    if not unique:
        return {}
    log.info("fitting %d distinct PON count vectors (%d requested) with %d threads", len(unique), len(keys), threads)
    chunks = [unique[i : i + chunk] for i in range(0, len(unique), chunk)]
    work = partial(_fit_chunk, penalty=penalty, method=method)
    if threads <= 1 or len(chunks) == 1:
        results = [ab for c in chunks for ab in work(c)]
    else:
        with get_context("spawn").Pool(threads) as pool:
            results = [ab for part in pool.map(work, chunks) for ab in part]
    return dict(zip(unique, results))
