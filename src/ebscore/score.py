"""EB score of variants: tumor counts + PONmatrix -> PON+ PON- EB.

Replaces the legacy ``ebrun.run_ebscore`` / ``file2matrix.tumor2matrix`` / ``tumor2matrix.mawk`` chain:
- positions are looked up with a hash filter (``filterpos.mawk``) instead of a sorted merge, so any number of
  chromosomes works in one run, nothing is lost when the PON misses a position, and pure-indel sites are kept
- the beta-binomial parameters are always fitted on the fly at the called positions (de-duplicated, parallel); the
  legacy on-the-fly path crashed (tumor_matrix.py:80), and the AB cache was removed (plan D32)
- PON samples are excluded explicitly by name or index instead of by a file-name prefix convention
- a PONmatrix subset (``ebscore pon-subset``, plan T2) scores exactly like its source; positions it does not list are
  refused (``subset.check_subset_positions``) instead of silently getting EB NA
"""

from __future__ import annotations

import logging
import os
import tempfile
from pathlib import Path

import pandas as pd

from .core import eb_score
from .counts import ab_string, drop_samples, fit_key, n_samples, pon_string, split_cell, tumor_string
from .fitting import fit_many
from .pileup import MIN_BASEQ, MIN_MAPQ, mpileup_cmd
from .subset import check_subset_positions
from .tools import cat_cmd, cmd2df, find_gz, first_line, run_parallel, tool
from .variants import VCF_COLUMNS, lookup_positions, read_variants, restrict_chroms

log = logging.getLogger("ebscore")

MISSING_POLICIES = ("na", "zero")

_find = find_gz


def _rows_by_pos(df: pd.DataFrame) -> dict:
    if df.empty:
        return {}
    return {(c, int(s)): row for c, s, row in zip(df["Chr"], df["Start"], df.to_dict("records"))}


def read_pon_list(ponmatrix) -> list[str]:
    p = Path(ponmatrix) / "PON_list.txt"
    if not p.is_file():
        return []
    return [line.strip() for line in p.read_text().splitlines() if line.strip() and not line.startswith("#")]


def _sample_names(entry: str) -> set[str]:
    base = os.path.basename(entry)
    stem = base.split(".")[0]
    return {entry, base, stem}


def resolve_exclusions(pon_list: list[str], n_pon: int, names=(), indices=()) -> set[int]:
    """0-based PON column indices to drop. ``indices`` are 1-based; ``names`` match a PON_list entry, its file name,
    or the file name up to the first dot."""
    exclude = set()
    for i in indices:
        i = int(i)
        if not 1 <= i <= n_pon:
            raise ValueError(f"--exclude-index {i} is outside 1..{n_pon}")
        exclude.add(i - 1)
    if names:
        if len(pon_list) != n_pon:
            raise ValueError(
                f"PON_list.txt has {len(pon_list)} entries but the PONmatrix has {n_pon} samples; "
                "exclusion by name is unsafe (use --exclude-index)"
            )
        for name in names:
            hits = [i for i, entry in enumerate(pon_list) if name in _sample_names(entry)]
            if not hits:
                raise ValueError(f"--exclude-sample {name!r} matches no PON_list entry")
            exclude.update(hits)
    return exclude


def tumor_counts_from_pileup(pileup, positions_file, tumor_column=None) -> pd.DataFrame:
    """pile2count rows (one sample) at the given positions from a cleanpileup file (``-d`` or not).

    Default tumor column: ``Read2`` if the header has it (paired pileup, normal = ``Read1``), else the last ``Read``
    column, i.e. ``Read1`` of a tumor-only pileup (``ebscore pileup`` without ``--normal``).
    """
    line = first_line(pileup)
    if not line:
        raise ValueError(f"{pileup}: empty file, not even the cleanpileup header (make it with `ebscore pileup`)")
    header = line.split("\t")
    if header[0] != "Chr":
        raise ValueError(f"{pileup}: no cleanpileup header (Chr Start Ref ...)")
    reads = [c for c in header if c.startswith("Read")]
    col = tumor_column or ("Read2" if "Read2" in header else reads[-1] if reads else None)
    if col not in header:
        raise ValueError(f"{pileup}: column {col!r} not found (have {' '.join(reads)})")
    idx = header.index(col) + 1
    log.info("tumor counts from %s column %s", pileup, col)
    select = f"mawk -F'\\t' -v OFS='\\t' -v c={idx} '{{ print $1, $2, $3, $c }}'"
    return cmd2df(f"{cat_cmd(pileup)} | {tool('filterpos', positions_file)} | {select} | {tool('pile2count')}")


def tumor_counts_from_bam(bam, fasta, positions: pd.DataFrame, workdir, min_mapq, min_baseq) -> pd.DataFrame:
    bed = Path(workdir) / "positions.bed"
    with open(bed, "w") as fh:
        for chrom, pos in positions.itertuples(index=False):
            fh.write(f"{chrom}\t{int(pos) - 1}\t{int(pos)}\n")
    log.info("tumor counts from %s at %d positions", bam, len(positions))
    cmd = f"{mpileup_cmd([bam], fasta, bed, None, min_mapq, min_baseq)} | {tool('cleanpileup')} | {tool('pile2count')}"
    return cmd2df(cmd)


def _lookup_per_chrom(folder: Path, suffix: str, chroms, positions_file, threads) -> tuple[pd.DataFrame, set]:
    """Look the positions up in ``<folder>/<chrom><suffix>`` for each chromosome; return rows + chroms found."""

    def one(chrom):
        f = _find(folder / f"{chrom}{suffix}")
        if f is None:
            return chrom, None
        return chrom, cmd2df(f"{cat_cmd(f)} | {tool('filterpos', positions_file)}")

    results = run_parallel(one, chroms, threads)
    found = {c for c, df in results if df is not None}
    frames = [df for _, df in results if df is not None and not df.empty]
    return (pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()), found


def score_variants(
    variants,
    ponmatrix,
    pileup=None,
    tumor_column=None,
    tumor_bam=None,
    fasta=None,
    exclude_samples=(),
    exclude_index=(),
    penalty=0.5,
    fisher="sf",
    missing="na",
    threads=1,
    chroms=None,
    min_mapq=MIN_MAPQ,
    min_baseq=MIN_BASEQ,
    keep_matrix=False,
    tmpdir=None,
    allow_missing_chrom=False,
) -> pd.DataFrame:
    """Score a variant table/VCF (path) or a prepared DataFrame. Returns the input columns + PON+ PON- EB.

    The beta-binomial parameters are fitted on the fly for every scored position (identical PON vectors once).
    ``EB`` is float64 (NaN = no data). Fails if ``<ponmatrix>/matrix`` is missing, if no variant chromosome has a
    matrix file (unless ``allow_missing_chrom``), or if ``ponmatrix`` is a subset (``subset.json``, plan T2) whose
    ``positions.tsv`` lacks a position to score.
    """
    if missing not in MISSING_POLICIES:
        raise ValueError(f"missing must be one of {MISSING_POLICIES}")
    if (pileup is None) == (tumor_bam is None):
        raise ValueError("give exactly one of pileup or tumor_bam")
    if tumor_bam is not None and fasta is None:
        raise ValueError("tumor_bam needs fasta")
    ponmatrix = Path(ponmatrix)
    if not (ponmatrix / "matrix").is_dir():
        raise ValueError(f"{ponmatrix} is not a PONmatrix folder (no matrix/ subfolder)")
    vdf = read_variants(variants) if not isinstance(variants, pd.DataFrame) else variants.copy()
    vdf = restrict_chroms(vdf, chroms)
    positions = lookup_positions(vdf)
    var_chroms = list(dict.fromkeys(positions["Chr"]))
    log.info("%d variants (%d scorable) on %d chromosomes", len(vdf), int(vdf["_pos"].notna().sum()), len(var_chroms))
    subset = check_subset_positions(ponmatrix, positions)
    if subset is not None:
        log.info("PONmatrix subset of %s (%s positions); all %d positions to score are listed",
                 (subset.get("source") or {}).get("path"), subset.get("n_positions"), len(positions))
    no_pon = [c for c in var_chroms if _find(ponmatrix / "matrix" / f"{c}.pon.gz") is None]
    if no_pon and len(no_pon) == len(var_chroms) and not allow_missing_chrom:
        raise ValueError(
            f"{ponmatrix}/matrix has no file for any variant chromosome ({' '.join(no_pon[:5])}"
            f"{' ...' if len(no_pon) > 5 else ''}); wrong --ponmatrix? (--allow-missing-chrom gives EB NA instead)"
        )
    for c in no_pon:
        log.warning("no PONmatrix file for %s: its variants get no PON data", c)

    with tempfile.TemporaryDirectory(prefix="ebscore_", dir=tmpdir) as work:
        pos_file = Path(work) / "positions.tsv"
        positions.to_csv(pos_file, sep="\t", header=False, index=False)
        if positions.empty:
            tum = pon = pd.DataFrame()
        else:
            if pileup is not None:
                tum = tumor_counts_from_pileup(pileup, pos_file, tumor_column)
            else:
                tum = tumor_counts_from_bam(tumor_bam, fasta, positions, work, min_mapq, min_baseq)
            pon_chroms = [c for c in var_chroms if c not in set(no_pon)]
            pon, _ = _lookup_per_chrom(ponmatrix / "matrix", ".pon.gz", pon_chroms, pos_file, threads)

        n_pon = n_samples(pon.iloc[0]["Depth"]) if not pon.empty else 0
        pon_list = read_pon_list(ponmatrix)
        if pon_list and n_pon and len(pon_list) != n_pon:
            log.warning("PON_list.txt has %d entries but the PONmatrix has %d samples", len(pon_list), n_pon)
        exclude = resolve_exclusions(pon_list, n_pon, exclude_samples, exclude_index) if n_pon else set()
        if exclude:
            log.info("excluding PON samples (1-based) %s", sorted(i + 1 for i in exclude))

    tum_rows, pon_rows = _rows_by_pos(tum), _rows_by_pos(pon)

    records, to_fit = [], []
    for chrom, pos, letter in zip(vdf["Chr"], vdf["_pos"], vdf["_letter"]):
        rec = {"PON+": pd.NA, "PON-": pd.NA, "EB": float("nan"), "Tumor": pd.NA, "AB": pd.NA, "_ab": None,
               "_keys": None}
        records.append(rec)
        if pd.isna(pos):
            continue
        key = (chrom, int(pos))
        t, p = tum_rows.get(key), pon_rows.get(key)
        if t is not None:
            ta, td = split_cell(t[letter]), split_cell(t["Depth"])
            rec["_t"] = ((ta[0][0], ta[1][0]), (td[0][0], td[1][0]))
            rec["Tumor"] = tumor_string(*rec["_t"])
        if p is None:
            continue
        pa, pd_ = split_cell(p[letter]), split_cell(p["Depth"])
        plus = (drop_samples(pa[0], exclude), drop_samples(pd_[0], exclude))
        minus = (drop_samples(pa[1], exclude), drop_samples(pd_[1], exclude))
        rec["PON+"], rec["PON-"] = pon_string(*plus), pon_string(*minus)
        if sum(plus[1]) + sum(minus[1]) == 0:
            continue  # only excluded samples cover this position
        rec["_keys"] = (fit_key(plus[1], plus[0]), fit_key(minus[1], minus[0]))
        to_fit.extend(rec["_keys"])

    fitted = fit_many(to_fit, penalty, threads)
    n_missing = 0
    for rec in records:
        if rec["_keys"] is not None:
            rec["_ab"] = (fitted[rec["_keys"][0]], fitted[rec["_keys"][1]])
        t = rec.get("_t")
        if rec["_ab"] is None or t is None or sum(t[1]) == 0:
            n_missing += 1
            rec["EB"] = 0.0 if missing == "zero" else float("nan")
            continue
        rec["AB"] = ab_string(*rec["_ab"])
        rec["EB"] = eb_score(t[0], t[1], rec["_ab"], fisher)
    if n_missing:
        log.info("%d variants without tumor or PON data -> EB %s", n_missing, "0" if missing == "zero" else "NA")

    out = vdf.drop(columns=["_pos", "_letter"])
    res = pd.DataFrame(records, index=out.index, columns=["PON+", "PON-", "EB", "Tumor", "AB"])
    out["PON+"], out["PON-"] = res["PON+"].astype(object), res["PON-"].astype(object)
    out["EB"] = res["EB"].astype("float64")
    if keep_matrix:
        out["Tumor"], out["AB"] = res["Tumor"].astype(object), res["AB"].astype(object)
    return out


EB_INFO_HEADER = (
    '##INFO=<ID=EB,Number=A,Type=Float,Description="EB artifact score vs the PONmatrix (ebscore): -log10 p of the '
    'tumor alt reads under the PON beta-binomial, 0-60; low = artifact-like">'
)


def write_vcf_annotation(df: pd.DataFrame, path) -> int:
    """Write ``#CHROM POS REF ALT EB`` (VCF coordinates, scored rows only) for bgzip + tabix + bcftools annotate.

    Rows are sorted by chromosome (input order) and numeric POS, so ``tabix -s1 -b2 -e2`` accepts the file. Also
    writes ``<path>.hdr`` with the INFO header line (``Number=A``: one value per ALT). One line per ALT is not enough
    for multi-allelic VCF records: bcftools 1.21 keeps only one ALT's EB (``EB=1.5,.``). Merge the lines of one
    CHROM POS REF first (ALTs and EBs comma-joined in the same order; see README) before ``bcftools annotate``.
    """
    if not set(VCF_COLUMNS) <= set(df.columns):
        raise ValueError("VCF annotation needs a VCF input (vcf_* columns)")
    sub = df.loc[df["EB"].notna(), VCF_COLUMNS + ["EB"]].copy()
    order = {c: i for i, c in enumerate(dict.fromkeys(df["vcf_chrom"]))}
    sub["_o"] = sub["vcf_chrom"].map(order)
    sub["_p"] = pd.to_numeric(sub["vcf_pos"])
    sub = sub.sort_values(["_o", "_p"], kind="stable").drop(columns=["_o", "_p"])
    sub.columns = ["#CHROM", "POS", "REF", "ALT", "EB"]
    sub.to_csv(path, sep="\t", index=False)
    Path(f"{path}.hdr").write_text(EB_INFO_HEADER + "\n")
    return len(sub)
