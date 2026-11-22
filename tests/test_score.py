"""EB scoring against the legacy run_ebscore output and the documented deltas."""

import pandas as pd
import pytest

from conftest import LEGACY, SYN
from ebscore.score import resolve_exclusions, score_variants, write_vcf_annotation

KEY = ["Chr", "Start", "End", "Ref", "Alt"]
PAIR = LEGACY / "pair.cpile.tsv.gz"


@pytest.fixture(scope="module")
def legacy_eb():
    df = pd.concat([pd.read_csv(LEGACY / "eb" / f"{c}.EB", sep="\t", dtype=str, keep_default_na=False)
                    for c in ("chr7", "chr12")])
    return df.reset_index(drop=True)


@pytest.fixture(scope="module")
def scenarios():
    return pd.read_csv(SYN / "variants.tsv", sep="\t", dtype=str)


def scored(**kw):
    df = score_variants(SYN / "variants.tsv", kw.pop("ponmatrix"), pileup=PAIR, **kw)
    return df.astype({"Start": str, "End": str})


def test_parity_with_legacy_ab_cache(legacy_ponmatrix, legacy_eb, scenarios):
    new = scored(ponmatrix=legacy_ponmatrix, fisher="legacy", missing="zero", threads=2)
    m = legacy_eb.merge(new, on=KEY, suffixes=("_leg", ""), how="left")  # new passes `scenario` through
    assert m["EB"].notna().all(), "every legacy row is also scored"
    # legacy bug: the variant right after a position without PON coverage loses its PON line (tumor2matrix.mawk)
    after_gap = m["PON+_leg"].eq("") & m["scenario"].ne("pon_gap")
    assert after_gap.sum() == 2
    assert m.loc[after_gap, "PON+"].str.contains("<").all()
    same = m[~after_gap]
    assert (same["EB"].astype(float) == same["EB_leg"].astype(float)).all()
    covered = same["PON+_leg"].ne("")
    assert (same.loc[covered, "PON+"] == same.loc[covered, "PON+_leg"]).all()
    assert (same.loc[covered, "PON-"] == same.loc[covered, "PON-_leg"]).all()


def test_rows_legacy_dropped_are_scored(legacy_ponmatrix, legacy_eb, scenarios):
    new = scored(ponmatrix=legacy_ponmatrix)
    assert len(new) == len(scenarios)
    extra = new.merge(legacy_eb[KEY], on=KEY, how="left", indicator=True).query("_merge == 'left_only'")
    # pure-indel sites (no SNV mismatch in the tumor) were dropped by the legacy $7 filter
    assert {"ins", "del"} <= set(extra["scenario"])
    assert extra.query("scenario in ['ins', 'del']")["EB"].notna().all()
    # both alleles of a multi-allelic site are scored (legacy kept only the last one per position)
    assert new.query("scenario == 'multi'")["EB"].notna().all()


def test_on_the_fly_fits_match_ab_cache(legacy_ponmatrix, legacy_matrix_only):
    cache = scored(ponmatrix=legacy_ponmatrix)
    fly = scored(ponmatrix=legacy_matrix_only, threads=2)
    both = cache["EB"].notna()
    assert (fly["EB"].notna() == both).all()
    diff = (cache.loc[both, "EB"].astype(float) - fly.loc[both, "EB"].astype(float)).abs()
    assert diff.max() <= 0.005  # exact fits vs the legacy (early-stopped) AB cache


def test_missing_policies(legacy_ponmatrix):
    na = scored(ponmatrix=legacy_ponmatrix)
    zero = scored(ponmatrix=legacy_ponmatrix, missing="zero")
    for kind in ("pon_gap", "tumor_gap", "mnv"):
        assert na.query("scenario == @kind")["EB"].isna().all()
        assert (zero.query("scenario == @kind")["EB"].astype(float) == 0).all()


def test_exclusion_drops_pon_column_and_skips_cache(legacy_ponmatrix):
    by_index = scored(ponmatrix=legacy_ponmatrix, exclude_index=[2], keep_matrix=True)
    by_name = scored(ponmatrix=legacy_ponmatrix, exclude_samples=["N02"], keep_matrix=True)
    pd.testing.assert_frame_equal(by_index, by_name)
    pons = by_index["PON+"].dropna()
    assert (pons.str.split("<").str[0].str.count(r"\|") == 4).all()  # 5 of 6 samples left


def test_exclusion_by_name_needs_consistent_pon_list():
    with pytest.raises(ValueError, match="unsafe"):
        resolve_exclusions(["a.bam", "b.bam"], 3, names=["a"])
    with pytest.raises(ValueError, match="matches no"):
        resolve_exclusions(["a.bam", "b.bam"], 2, names=["zz"])
    assert resolve_exclusions(["x/N01.bam", "x/N02.bam"], 2, names=["N02"]) == {1}


def test_vcf_input_and_annotation(legacy_ponmatrix, tmp_path):
    tab = scored(ponmatrix=legacy_ponmatrix)
    vcf = score_variants(SYN / "variants.vcf", legacy_ponmatrix, pileup=PAIR)
    m = tab.merge(vcf.dropna(subset=["Start"]).astype({"Start": int}).astype({"Start": str, "End": str}),
                  on=KEY, suffixes=("", "_vcf"))
    assert len(m) == len(tab) - 2  # all but the MNVs
    assert (m["EB"].fillna(-1).astype(float) == m["EB_vcf"].fillna(-1).astype(float)).all()
    n = write_vcf_annotation(vcf, tmp_path / "eb.tsv")
    ann = pd.read_csv(tmp_path / "eb.tsv", sep="\t")
    assert list(ann.columns) == ["#CHROM", "POS", "REF", "ALT", "EB"] and len(ann) == n == vcf["EB"].notna().sum()


def test_chrom_restriction(legacy_ponmatrix):
    only7 = scored(ponmatrix=legacy_ponmatrix, chroms=["chr7"])
    full = scored(ponmatrix=legacy_ponmatrix)
    assert set(only7["Chr"]) == {"chr7"}
    pd.testing.assert_frame_equal(only7.reset_index(drop=True),
                                  full[full["Chr"] == "chr7"].reset_index(drop=True))
