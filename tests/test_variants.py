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
    ],
)
def test_vcf_record_to_annovar(rec, expected):
    assert vcf_record_to_annovar(*rec) == expected


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
