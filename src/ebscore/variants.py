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
    """One biallelic VCF record -> (Chr, Start, End, Ref, Alt) in ANNOVAR coordinates, or None if unsupported.

    SNV ``chr1 100 A G`` -> ``chr1 100 100 A G``; insertion ``chr1 100 A AT`` -> ``chr1 100 100 - T``;
    deletion ``chr1 100 AT A`` -> ``chr1 101 101 T -``. MNVs, complex and symbolic alleles return None.
    """
    if alt in ("*", ".") or alt.startswith("<") or "[" in alt or "]" in alt:
        return None
    ref, alt = ref.upper(), alt.upper()
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
    n_bad = int(df["Start"].isna().sum())
    if n_bad:
        log.warning("%d VCF alleles are MNV/complex/symbolic and get EB = NA", n_bad)
    return df


def read_table(path) -> pd.DataFrame:
    """Read a tab-separated variant table with at least ``Chr Start End Ref Alt`` (extra columns pass through)."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = [c for c in VARIANT_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing columns {missing} (need {' '.join(VARIANT_COLUMNS)})")
    return df


def read_variants(path) -> pd.DataFrame:
    df = read_vcf(path) if is_vcf(path) else read_table(path)
    return add_lookup(df)


def add_lookup(df: pd.DataFrame) -> pd.DataFrame:
    """Add the pileup lookup columns ``_pos`` (int, NA if unsupported) and ``_letter`` (A/G/C/T/I/D).

    SNV -> its Alt base at Start; insertion (Ref "-") -> I at Start; deletion (Alt "-") -> D at Start - 1, the anchor
    base where the pileup reports the deletion (legacy tumor2matrix.mawk logic).
    """
    df = df.reset_index(drop=True).copy()
    pos, letter = [], []
    for start, ref, alt in zip(df["Start"], df["Ref"], df["Alt"]):
        if pd.isna(start) or start == "":
            pos.append(pd.NA)
            letter.append("")
            continue
        start = int(start)
        if ref == "-":
            pos.append(start)
            letter.append("I")
        elif alt == "-":
            pos.append(start - 1)
            letter.append("D")
        elif len(ref) == 1 and len(alt) == 1 and alt.upper() in SNV_BASES:
            pos.append(start)
            letter.append(alt.upper())
        else:
            pos.append(pd.NA)
            letter.append("")
    df["_pos"] = pd.array(pos, dtype="Int64")
    df["_letter"] = letter
    n_bad = int(df["_pos"].isna().sum())
    if n_bad:
        log.warning("%d variants are not SNV/insertion/deletion and get EB = NA", n_bad)
    return df
