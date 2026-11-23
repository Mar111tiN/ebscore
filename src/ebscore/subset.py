"""Run-level PONmatrix subset (plan T2): the PONmatrix rows at the positions that will be scored.

In the pipeline every EB task staged the whole PONmatrix (~20 GB at 100 normals). ``pon_subset`` reads the full set
once per run and keeps only the matrix rows at the union of the positions that ``score`` looks up for the given
variant files; each EB task then stages that subset (MBs). Scoring against the subset gives the same output as
scoring against the full set, bit for bit: the rows are the source lines unchanged (``filterpos``), and
``PON_list.txt`` is a byte copy.

A subset folder::

    matrix/<chrom>.pon.gz   header + the source rows at the listed positions (only chromosomes that have positions
                            and a source file; header-only if no listed position has a row)
    PON_list.txt            byte copy of the source's (if it has one; if not, a stale one in the folder is removed)
    positions.tsv           Chr Pos (tab, no header): every listed position, chromosomes in order of first appearance
                            in the inputs, positions numeric; Pos is the lookup position of ``score`` (``_pos``: the
                            Start of an SNV/insertion, Start - 1 for a deletion)
    subset.json             provenance: schema, ebscore version, source path + manifest, counts, md5 of positions.tsv

``score`` refuses to run against a subset (a folder with ``subset.json``) when a scorable position of its input is not
in ``positions.tsv`` (:func:`check_subset_positions`), so a wiring error cannot silently give EB NA. A listed position
without a source row is fine (EB NA, as with the full set). A subset can itself be the source of a smaller subset
(run level -> per pair); its positions must then be listed in the source subset.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import logging
import os
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__
from .counts import MATRIX_COLUMNS, n_samples
from .tools import cat_cmd, find_gz, first_line, first_lines, q, run, run_parallel, tool, unique
from .variants import lookup_positions, read_variants, restrict_chroms

log = logging.getLogger("ebscore")

SUBSET_SCHEMA = "ebscore-ponmatrix-subset/1"
SUBSET_JSON = "subset.json"
POSITIONS_TSV = "positions.tsv"


def md5_file(path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_positions(path) -> pd.DataFrame:
    """A positions file (``Chr Pos``, tab-separated, no header; the ``positions.tsv`` format) -> ``Chr`` str,
    ``Pos`` int64. An empty file gives no positions; anything else than two columns with a positive integer Pos fails.
    """
    empty = pd.DataFrame({"Chr": pd.Series(dtype=object), "Pos": pd.Series(dtype="int64")})
    if os.path.getsize(path) == 0:
        return empty
    try:
        df = pd.read_csv(path, sep="\t", header=None, dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:  # blank lines only
        return empty
    if df.shape[1] != 2:
        raise ValueError(f"{path}: expected 2 tab-separated columns Chr Pos (no header), found {df.shape[1]}")
    df.columns = ["Chr", "Pos"]
    bad = ~df["Pos"].str.fullmatch(r"[0-9]+") | df["Chr"].eq("")
    if bad.any():
        i = int(np.flatnonzero(bad.to_numpy())[0])
        raise ValueError(f"{path}: line {i + 1} ({df.iloc[i]['Chr']!r} {df.iloc[i]['Pos']!r}) is not Chr + a "
                         "positive integer Pos (the file has no header line)")
    pos = df["Pos"].astype("int64")
    if (pos < 1).any():
        raise ValueError(f"{path}: Pos must be >= 1")
    return pd.DataFrame({"Chr": df["Chr"].astype(str).astype(object), "Pos": pos})


def read_subset_info(ponmatrix) -> dict | None:
    """``subset.json`` of a PONmatrix subset (checked against its ``positions.tsv``), or None for a full set."""
    ponmatrix = Path(ponmatrix)
    p = ponmatrix / SUBSET_JSON
    if not p.is_file():
        return None
    try:
        info = json.loads(p.read_text())
    except json.JSONDecodeError as e:
        raise ValueError(f"{p}: not valid JSON ({e})") from None
    if not isinstance(info, dict) or info.get("schema") != SUBSET_SCHEMA:
        raise ValueError(f"{p}: unknown subset schema {info.get('schema') if isinstance(info, dict) else None!r} "
                         f"(this ebscore reads {SUBSET_SCHEMA})")
    pos_file = ponmatrix / POSITIONS_TSV
    if not pos_file.is_file():
        raise ValueError(f"{ponmatrix} has {SUBSET_JSON} but no {POSITIONS_TSV}: incomplete PONmatrix subset")
    if md5_file(pos_file) != info.get("positions_md5"):
        raise ValueError(f"{pos_file} does not match {SUBSET_JSON} (md5): incomplete or edited PONmatrix subset")
    return info


def check_subset_positions(ponmatrix, positions: pd.DataFrame, what: str = "scoring positions") -> dict | None:
    """Raise ValueError if ``ponmatrix`` is a subset that does not list every ``(Chr, _pos)`` of ``positions``.

    ``positions`` is :func:`variants.lookup_positions` output. Returns the subset info (None for a full set).
    """
    info = read_subset_info(ponmatrix)
    if info is None or positions.empty:
        return info
    listed = read_positions(Path(ponmatrix) / POSITIONS_TSV)
    want = pd.DataFrame({"Chr": positions["Chr"].astype(str).to_numpy(dtype=object),
                         "Pos": positions["_pos"].astype("int64").to_numpy()})
    m = want.merge(listed.drop_duplicates(), on=["Chr", "Pos"], how="left", indicator=True)
    missing = m.loc[m["_merge"] == "left_only", ["Chr", "Pos"]]
    if not missing.empty:
        first = " ".join(f"{c}:{p}" for c, p in missing.head(5).itertuples(index=False))
        more = " ..." if len(missing) > 5 else ""
        raise ValueError(
            f"{ponmatrix}: {len(missing)} {what} outside this PONmatrix subset (first: {first}{more}): rebuild it "
            "with these variants (ebscore pon-subset)"
        )
    return info


def _read_manifest(ponmatrix: Path):
    p = ponmatrix / "manifest.json"
    if not p.is_file():
        return None
    try:
        return json.loads(p.read_text())
    except json.JSONDecodeError as e:
        raise ValueError(f"{p}: not valid JSON ({e})") from None


def _add(union: dict, df: pd.DataFrame) -> int:
    """Add the (Chr, _pos) pairs of ``df`` to ``union`` (chrom -> list of int64 arrays; key order = first appearance)."""
    if df.empty:
        return 0
    for chrom, grp in df.groupby("Chr", sort=False):
        union.setdefault(str(chrom), []).append(grp["_pos"].astype("int64").to_numpy())
    return len(df)


def collect_positions(variants=(), positions=(), chroms=None) -> tuple[dict, list]:
    """Union of the lookup positions of the inputs: ``{chrom: sorted unique int64 array}`` in chromosome order of first
    appearance (variants files first, then positions files, each in the given order), and one record per input.

    Variant files go through the same reader and lookup as ``score`` (``read_variants``, ``restrict_chroms``,
    ``lookup_positions``); unscorable rows (no ``_pos``: MNV, complex, symbolic, or no valid position such as a
    deletion at Start 1) are left out and counted, so every position is one ``read_positions`` accepts.
    """
    union: dict = {}
    inputs = []
    for path in variants:
        vdf = restrict_chroms(read_variants(path), chroms)
        pos = lookup_positions(vdf)
        n_bad = int(vdf["_pos"].isna().sum())
        _add(union, pos)
        log.info("%s: %d variants, %d positions to look up, %d unscorable rows ignored", path, len(vdf), len(pos),
                 n_bad)
        inputs.append({"kind": "variants", "path": str(path), "rows": len(vdf), "unscorable": n_bad,
                       "positions": len(pos)})
    for path in positions:
        pdf = read_positions(path)
        if chroms:
            pdf = pdf[pdf["Chr"].isin(set(chroms))]
        pos = pdf.rename(columns={"Pos": "_pos"}).drop_duplicates()
        _add(union, pos)
        log.info("%s: %d positions (%d distinct)", path, len(pdf), len(pos))
        inputs.append({"kind": "positions", "path": str(path), "rows": len(pdf), "positions": len(pos)})
    return {c: np.unique(np.concatenate(a)) for c, a in union.items()}, inputs


def _check_matrix_header(src: Path) -> None:
    header = first_line(src).split("\t")
    if header[:2] != ["Chr", "Start"] or not set(MATRIX_COLUMNS) <= set(header):
        raise ValueError(f"{src}: not a PONmatrix file (header {' '.join(header[:10]) or 'missing'}; "
                         f"need {' '.join(MATRIX_COLUMNS)})")


def _matrix_stats(path: Path) -> tuple[int, set]:
    """Data rows of a (gzipped) PONmatrix file and the sample counts of their Depth cells."""
    rows, counts = 0, set()
    with gzip.open(path, "rt") as fh:
        d = fh.readline().rstrip("\n").split("\t").index("Depth")
        for line in fh:
            rows += 1
            counts.add(n_samples(line.rstrip("\n").split("\t")[d]))
    return rows, counts


def _first_row_samples(src: Path) -> set:
    lines = first_lines(src, 2)
    if len(lines) < 2:
        return set()
    header = lines[0].split("\t")
    return {n_samples(lines[1].split("\t")[header.index("Depth")])}


def _same_path(a: Path, b: Path) -> bool:
    return os.path.realpath(a) == os.path.realpath(b)


def _inside(child: Path, parent: Path) -> bool:
    c, p = os.path.realpath(child), os.path.realpath(parent)
    return c == p or c.startswith(p.rstrip(os.sep) + os.sep)


def _write_atomic_copy(src: Path, dst: Path) -> None:
    part = dst.with_name(f".{dst.name}.part")
    shutil.copyfile(src, part)
    os.replace(part, dst)


def pon_subset(
    ponmatrix,
    outdir,
    variants=(),
    positions=(),
    chroms=None,
    threads=1,
    allow_missing_chrom=False,
) -> dict:
    """Write the PONmatrix subset of ``ponmatrix`` at the lookup positions of ``variants`` + ``positions`` to
    ``outdir`` (see the module docstring); return the ``subset.json`` content.

    ``chroms`` restricts the positions as ``score --chrom`` does. Chromosomes are filtered in parallel (``threads``).
    Fails if ``outdir`` is the source or already has ``matrix/``, if the source is a subset that does not list every
    requested position, or if the source is a full set and no chromosome with positions has a file there (as ``score``
    against the full set would; ``allow_missing_chrom`` writes the subset anyway). A subset source is not checked for
    that (see the comment in the code), so one pair cannot fail the per-pair step of a whole run.
    ``PON_list.txt`` is copied if the source has one; a stale one in ``outdir`` is removed if it has none.
    The matrix is written in a hidden staging folder and moved to ``<outdir>/matrix`` last, after ``PON_list.txt``,
    ``positions.tsv`` and ``subset.json``; an interrupted run leaves no ``matrix/``, so nothing can score against it.
    """
    src, out = Path(ponmatrix), Path(outdir)
    if not (src / "matrix").is_dir():
        raise ValueError(f"{src} is not a PONmatrix folder (no matrix/ subfolder)")
    if _same_path(src, out) or _inside(out, src / "matrix"):
        raise ValueError(f"-o {out} is the source PONmatrix (or inside its matrix/): write the subset to a new folder")
    if (out / "matrix").exists() or (out / "matrix").is_symlink():
        raise ValueError(f"{out}/matrix already exists: pon-subset writes a new subset (remove it or use another -o)")
    if not variants and not positions:
        raise ValueError("pon-subset needs at least one variants or positions file")
    chroms = unique(chroms)

    union, inputs = collect_positions(variants, positions, chroms)
    n_pos = int(sum(len(a) for a in union.values()))
    want = pd.DataFrame({"Chr": np.repeat(np.array(list(union), dtype=object), [len(a) for a in union.values()]),
                         "_pos": np.concatenate(list(union.values())) if union else np.array([], dtype="int64")})
    source_subset = check_subset_positions(src, want, what="requested positions")
    if source_subset is not None:
        log.info("source %s is itself a PONmatrix subset (%d positions)", src, source_subset.get("n_positions", -1))

    sources = {c: find_gz(src / "matrix" / f"{c}.pon.gz") for c in union}
    no_file = [c for c, f in sources.items() if f is None]
    # Only a full set is checked: a subset has no file for a chromosome exactly when its full set had none, and the
    # run-level call from the full set has already checked that. A per-pair call (source = the run subset) whose
    # positions are all on such chromosomes writes its subset (no matrix files), and that pair's `score` decides
    # (it fails without --allow-missing-chrom, as against the full set), instead of failing the run's other pairs.
    if no_file and len(no_file) == len(union) and source_subset is None and not allow_missing_chrom:
        raise ValueError(
            f"{src}/matrix has no file for any chromosome with positions ({' '.join(no_file[:5])}"
            f"{' ...' if len(no_file) > 5 else ''}); wrong --ponmatrix? (--allow-missing-chrom writes it anyway)"
        )
    for c in no_file:
        log.warning("no PONmatrix file for %s: its %d positions get no PON data (as with the full set)", c,
                    len(union[c]))
    if no_file and len(no_file) == len(union) and source_subset is not None:
        log.warning("no chromosome with positions has a matrix file in the source subset: score against this subset "
                    "fails unless --allow-missing-chrom, as against the full set")
    todo = [c for c in union if sources[c] is not None]
    for c in todo:
        _check_matrix_header(sources[c])
    manifest = _read_manifest(src)
    log.info("PONmatrix subset: %d positions on %d chromosomes (%d with a matrix file), %d threads", n_pos,
             len(union), len(todo), threads)

    out.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".pon-subset.", dir=out))
    try:
        (stage / "matrix").mkdir()
        (stage / "pos").mkdir()
        pos_files = {}
        with open(stage / POSITIONS_TSV, "w") as all_fh:
            for i, (chrom, arr) in enumerate(union.items()):
                text = "".join(f"{chrom}\t{p}\n" for p in arr.tolist())
                all_fh.write(text)
                pos_files[chrom] = stage / "pos" / f"{i:05d}.tsv"
                pos_files[chrom].write_text(text)

        def one(chrom):
            dst = stage / "matrix" / f"{chrom}.pon.gz"
            run(f"{cat_cmd(sources[chrom])} | {tool('filterpos', pos_files[chrom])} | gzip -n -c > {q(dst)}")
            rows, counts = _matrix_stats(dst)
            counts |= _first_row_samples(sources[chrom])
            log.info("subset %s: %d positions, %d matrix rows", chrom, len(union[chrom]), rows)
            return chrom, rows, counts

        results = {c: (rows, counts) for c, rows, counts in run_parallel(one, todo, threads)}
        all_counts = set().union(*(counts for _, counts in results.values()))
        if len(all_counts) > 1:
            raise ValueError(f"{src}/matrix: rows with different sample counts {sorted(all_counts)}; corrupt matrix?")
        n_smp = next(iter(all_counts)) if all_counts else None

        info = {
            "schema": SUBSET_SCHEMA,
            "ebscore": __version__,
            "source": {"path": str(ponmatrix), "manifest": manifest, "subset": source_subset},
            "inputs": inputs,
            "chrom_filter": chroms or None,
            "n_positions": n_pos,
            "n_samples": n_smp,
            "chroms": {c: {"positions": int(len(union[c])), "rows": results[c][0] if c in results else 0,
                           "matrix": c in results} for c in union},
            "positions_md5": md5_file(stage / POSITIONS_TSV),
        }

        pon_list = src / "PON_list.txt"
        if pon_list.is_file():
            _write_atomic_copy(pon_list, out / "PON_list.txt")
        elif (out / "PON_list.txt").exists() or (out / "PON_list.txt").is_symlink():
            log.warning("%s has no PON_list.txt: removing the stale %s", src, out / "PON_list.txt")
            (out / "PON_list.txt").unlink()
        os.replace(stage / POSITIONS_TSV, out / POSITIONS_TSV)
        part = out / f".{SUBSET_JSON}.part"
        part.write_text(json.dumps(info, indent=2) + "\n")
        os.replace(part, out / SUBSET_JSON)
        os.replace(stage / "matrix", out / "matrix")
    finally:
        shutil.rmtree(stage, ignore_errors=True)
    rows = sum(r for r, _ in results.values())
    log.info("wrote %s: %d positions, %d matrix rows, %d samples", out, n_pos, rows, n_smp or 0)
    return info
