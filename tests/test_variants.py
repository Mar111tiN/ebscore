"""VCF -> ANNOVAR coordinates and lookup keys."""

import pandas as pd
import pytest

from conftest import SYN
from ebscore.variants import add_lookup, read_table, read_vcf, vcf_record_to_annovar


@pytest.mark.parametrize(
    "rec, expected",
    [
        (("chr1", 100, "A", "G"), ("chr1", 100, 100, "A", "G")),
        (("chr1", 100, "A", "AT"), ("chr1", 100, 100, "-", "T")),
        (("chr1", 100, "A", "ATTG"), ("chr1", 100, 100, "-", "TTG")),
        (("chr1", 100, "AT", "A"), ("chr1", 101, 101, "T", "-")),
        (("chr1", 100, "ATGC", "A"), ("chr1", 101, 103, "TGC", "-")),
        (("chr1", 100, "AT", "GC"), None),  # MNV
        (("chr1", 100, "ATG", "AC"), None),  # complex
        (("chr1", 100, "A", "<DEL>"), None),
        (("chr1", 100, "A", "*"), None),
        # alleles of multi-allelic records keep the right-hand context of the longest REF
        (("chr1", 100, "ATT", "AT"), ("chr1", 101, 101, "T", "-")),
        (("chr1", 100, "ATT", "ATTT"), ("chr1", 100, 100, "-", "T")),
        (("chr1", 100, "AT", "GT"), ("chr1", 100, 100, "A", "G")),
        (("chr1", 100, "ATG", "AG"), ("chr1", 101, 101, "T", "-")),
        (("chr1", 100, "AT", "ACT"), ("chr1", 100, 100, "-", "C")),
        (("chr1", 103, "GCCC", "GCC"), ("chr1", 104, 104, "C", "-")),
        (("chr1", 103, "GCCC", "GC"), ("chr1", 104, 105, "CC", "-")),
        (("chr1", 100, "ACGT", "TCGT"), ("chr1", 100, 100, "A", "T")),
        (("chr1", 100, "ATGT", "ACGT"), ("chr1", 101, 101, "T", "C")),
        (("chr1", 100, "ATG", "ACT"), None),  # MNV with context
    ],
)
def test_vcf_record_to_annovar(rec, expected):
    assert vcf_record_to_annovar(*rec) == expected


def test_read_vcf_splits_multiallelic_records(tmp_path):
    f = tmp_path / "m.vcf"
    f.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                 "chr1\t100\t.\tATT\tAT,ATTT,GTT,<NON_REF>\t.\tPASS\t.\n"
                 "chr1\t200\t.\tC\tG,CA\t.\tPASS\t.\n")
    df = read_vcf(f)
    assert df["vcf_alt"].tolist() == ["AT", "ATTT", "GTT", "<NON_REF>", "G", "CA"]
    assert (df["vcf_pos"] == [100] * 4 + [200] * 2).all() and set(df["vcf_ref"]) == {"ATT", "C"}
    got = [tuple(r) if pd.notna(r[1]) else None for r in df[["Chr", "Start", "End", "Ref", "Alt"]].itertuples(index=False)]
    assert got == [("chr1", 101, 101, "T", "-"), ("chr1", 100, 100, "-", "T"), ("chr1", 100, 100, "A", "G"), None,
                   ("chr1", 200, 200, "C", "G"), ("chr1", 200, 200, "-", "A")]


def test_synthetic_vcf_matches_table():
    vcf = read_vcf(SYN / "variants.vcf")
    tab = read_table(SYN / "variants.tsv")
    ok = vcf["Start"].notna()
    assert int((~ok).sum()) == 2  # the two MNVs
    got = vcf.loc[ok, ["Chr", "Start", "End", "Ref", "Alt"]].astype(str).reset_index(drop=True)
    exp = tab.loc[tab["scenario"] != "mnv", ["Chr", "Start", "End", "Ref", "Alt"]].reset_index(drop=True)
    pd.testing.assert_frame_equal(got, exp)


def test_lookup_keys():
    df = pd.DataFrame(
        [("c", "10", "10", "A", "G"), ("c", "10", "10", "-", "TT"), ("c", "11", "12", "GA", "-"), ("c", "5", "6", "AT", "GC")],
        columns=["Chr", "Start", "End", "Ref", "Alt"],
    )
    out = add_lookup(df)
    assert out["_pos"].tolist()[:3] == [10, 10, 10]
    assert out["_letter"].tolist() == ["G", "I", "D", ""]
    assert pd.isna(out["_pos"].iloc[3])


def test_lookup_needs_chr_and_a_position_of_at_least_1(caplog):
    # no pileup or PONmatrix row exists at position < 1 or without a chromosome: unscorable (EB NA), never a
    # position that pon-subset would write to positions.tsv and read_positions reject
    df = pd.DataFrame(
        [("c", "1", "1", "A", "-"), ("c", "0", "0", "A", "G"), ("c", "-3", "-3", "-", "T"), ("", "10", "10", "A", "G"),
         ("c", "2", "2", "A", "-"), ("c", "1", "1", "-", "T"), ("c", "1", "1", "A", "G"), ("c", "5", "6", "AT", "GC")],
        columns=["Chr", "Start", "End", "Ref", "Alt"],
    )
    with caplog.at_level("WARNING", logger="ebscore"):
        out = add_lookup(df)
    assert out["_pos"].isna().tolist() == [True] * 4 + [False] * 3 + [True]
    assert out["_pos"].tolist()[4:7] == [1, 1, 1] and out["_letter"].tolist() == [""] * 4 + ["D", "I", "G", ""]
    assert "1 variants are not SNV/insertion/deletion" in caplog.text
    assert "4 variants have no lookup position" in caplog.text
