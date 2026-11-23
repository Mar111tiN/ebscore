"""Variant input: ANNOVAR-style tables (``Chr Start End Ref Alt``) or VCF, plus the pileup lookup keys."""

from __future__ import annotations

import gzip
import logging

import pandas as pd

log = logging.getLogger("ebscore")

VARIANT_COLUMNS = ["Chr", "Start", "End", "Ref", "Alt"]
VCF_COLUMNS = ["vcf_chrom", "vcf_pos", "vcf_ref", "vcf_alt"]
SNV_BASES = set("ACGT")


def is_vcf(path) -> bool:
    return str(path).endswith((".vcf", ".vcf.gz", ".vcf.bgz"))


def vcf_record_to_annovar(chrom: str, pos: int, ref: str, alt: str):
    """One VCF allele (REF, one ALT) -> (Chr, Start, End, Ref, Alt) in ANNOVAR coordinates, or None if unsupported.

    SNV ``chr1 100 A G`` -> ``chr1 100 100 A G``; insertion ``chr1 100 A AT`` -> ``chr1 100 100 - T``;
    deletion ``chr1 100 AT A`` -> ``chr1 101 101 T -``. MNVs, complex and symbolic alleles return None.

    The shared trailing bases are trimmed first, then the shared leading bases (as ``bcftools norm`` does), so the
    alleles of a multi-allelic record, which keep the right-hand context of the longest REF (``ATT AT,ATTT``), map to
    the same events as their biallelic forms (``AT A``, ``A AT``). Without the reference an indel inside a repeat
    cannot be left-aligned: normalise the VCF upstream (``bcftools norm -f <fasta>``).
    """
    if alt in ("*", ".") or alt.startswith("<") or "[" in alt or "]" in alt:
        return None
    ref, alt = ref.upper(), alt.upper()
    # trim the shared trailing bases (right-hand context of multi-allelic records), keeping one base on each side
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    # trim the shared leading bases (VCF anchor)
    while len(ref) > 0 and len(alt) > 0 and ref[0] == alt[0] and (len(ref) > 1 or len(alt) > 1):
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    if not ref:  # insertion after base pos-1
        return chrom, pos - 1, pos - 1, "-", alt
    if not alt:  # deletion of pos .. pos+len-1
        return chrom, pos, pos + len(ref) - 1, ref, "-"
    if len(ref) == 1 and len(alt) == 1:
        return chrom, pos, pos, ref, alt
    return None


def read_vcf(path) -> pd.DataFrame:
    """Read the variants of a VCF (multi-allelic records are split per ALT) into ANNOVAR coordinates.

    Unsupported alleles are kept with empty ANNOVAR columns so that they appear in the output with EB = NA.
    The original VCF coordinates are kept in ``vcf_chrom vcf_pos vcf_ref vcf_alt``.
    """
    opener = gzip.open if str(path).endswith((".gz", ".bgz")) else open
    rows = []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, ref, alts = f[0], int(f[1]), f[3], f[4]
            for alt in alts.split(","):
                annovar = vcf_record_to_annovar(chrom, pos, ref, alt)
                rows.append((*(annovar or (chrom, pd.NA, pd.NA, pd.NA, pd.NA)), chrom, pos, ref, alt))
    df = pd.DataFrame(rows, columns=VARIANT_COLUMNS + VCF_COLUMNS)
    for col in ("Start", "End", "vcf_pos"):
        df[col] = pd.array(df[col], dtype="Int64")
    n_bad = int(df["Start"].isna().sum())
    if n_bad:
        log.warning("%d VCF alleles are MNV/complex/symbolic and get EB = NA", n_bad)
    return df


def read_table(path) -> pd.DataFrame:
    """Read a tab-separated variant table with at least ``Chr Start End Ref Alt`` (extra columns pass through).

    Every column is read as text, so values pass through unchanged (``vcf_pos`` of a ``vcf2table`` table too).
    """
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = [c for c in VARIANT_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing columns {missing} (need {' '.join(VARIANT_COLUMNS)})")
    return df


def read_variants(path) -> pd.DataFrame:
    df = read_vcf(path) if is_vcf(path) else read_table(path)
    return add_lookup(df)


def restrict_chroms(df: pd.DataFrame, chroms=None) -> pd.DataFrame:
    """The rows on ``chroms`` (index reset; all rows if ``chroms`` is empty): ``score``/``pon-subset --chrom``."""
    if not chroms:
        return df
    return df[df["Chr"].isin(set(chroms))].reset_index(drop=True)


def lookup_positions(df: pd.DataFrame) -> pd.DataFrame:
    """The distinct ``(Chr, _pos)`` pairs of an :func:`add_lookup` table, in first-appearance order: exactly the
    positions ``score`` looks up in the pileup and the PONmatrix, and those a ``pon-subset`` keeps. Unscorable rows
    (``_pos`` NA) are left out."""
    return df.loc[df["_pos"].notna(), ["Chr", "_pos"]].drop_duplicates()


def add_lookup(df: pd.DataFrame) -> pd.DataFrame:
    """Add the pileup lookup columns ``_pos`` (int, NA if unsupported) and ``_letter`` (A/G/C/T/I/D).

    SNV -> its Alt base at Start; insertion (Ref "-") -> I at Start; deletion (Alt "-") -> D at Start - 1, the anchor
    base where the pileup reports the deletion (legacy tumor2matrix.mawk logic).

    A row with an empty Chr, or whose lookup position would be < 1 (a deletion at Start 1, a Start <= 0), has no
    position either: no pileup or PONmatrix row can exist there, so EB is NA as before, and ``score`` and
    ``pon-subset`` (whose ``positions.tsv`` only takes Chr + a position >= 1) treat it alike, as unscorable.
    """
    df = df.reset_index(drop=True).copy()
    pos, letter = [], []
    n_type = n_where = 0
    for chrom, start, ref, alt in zip(df["Chr"], df["Start"], df["Ref"], df["Alt"]):
        if pd.isna(start) or start == "":
            n_type += 1  # VCF alleles that are MNV/complex/symbolic come with an empty Start
            pos.append(pd.NA)
            letter.append("")
            continue
        start = int(start)
        if ref == "-":
            p, let = start, "I"
        elif alt == "-":
            p, let = start - 1, "D"
        elif len(ref) == 1 and len(alt) == 1 and alt.upper() in SNV_BASES:
            p, let = start, alt.upper()
        else:
            n_type += 1
            pos.append(pd.NA)
            letter.append("")
            continue
        if p < 1 or pd.isna(chrom) or str(chrom) == "":
            n_where += 1
            pos.append(pd.NA)
            letter.append("")
            continue
        pos.append(p)
        letter.append(let)
    df["_pos"] = pd.array(pos, dtype="Int64")
    df["_letter"] = letter
    if n_type:
        log.warning("%d variants are not SNV/insertion/deletion and get EB = NA", n_type)
    if n_where:
        log.warning("%d variants have no lookup position (empty Chr, or position < 1, e.g. a deletion at Start 1) "
                    "and get EB = NA", n_where)
    return df
