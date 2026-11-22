"""Golden exact fits: tests/data/fits_exact.tsv (depths, alts, AB = fit_bb(method="exact") with the package).

Vectors: the 102 legacy fixture vectors (tests/data/legacy/fits.tsv), the false-stop vectors (tests/data/fits_hard.tsv)
and synthetic large-PON vectors (26-100 normals at depth 150-300, zero or sparse alts, where b is in the thousands),
plus a germline-like vector whose b sits at its lower bound.

Run: python tests/data/make_fits_exact.py   (env somvar; made with scipy 1.18.1)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ebscore.core import fit_bb

HERE = Path(__file__).parent
SEED = 20261122


def large_pon_vectors():
    rng = np.random.default_rng(SEED)
    out = []
    for ns, depth in [(26, 150), (30, 150), (60, 200), (100, 200), (100, 300)]:
        for n_alt in (0, 1, 3):
            d = rng.poisson(depth, size=ns)
            a = np.zeros(ns, dtype=int)
            for i in rng.choice(ns, size=n_alt, replace=False):
                a[i] = int(rng.integers(1, 4))
            out.append(("|".join(map(str, d)), "|".join(map(str, a))))
    # germline-like: nearly all reads alt, so a is large and b sits at its lower bound 1
    d = rng.poisson(40, size=20)
    out.append(("|".join(map(str, d)), "|".join(map(str, d))))
    return out


def main():
    vecs = pd.concat([
        pd.read_csv(HERE / "legacy" / "fits.tsv", sep="\t", dtype=str)[["depths", "alts"]],
        pd.read_csv(HERE / "fits_hard.tsv", sep="\t", dtype=str),
        pd.DataFrame(large_pon_vectors(), columns=["depths", "alts"]),
    ]).drop_duplicates()
    ab = []
    for d, a in vecs.itertuples(index=False):
        cm = np.transpose([d.split("|"), a.split("|")]).astype(int)
        ab.append("|".join(map(repr, fit_bb(cm))))
    vecs["AB"] = ab
    vecs.to_csv(HERE / "fits_exact.tsv", sep="\t", index=False)
    print(f"{len(vecs)} exact fits -> {HERE / 'fits_exact.tsv'}")


if __name__ == "__main__":
    main()
