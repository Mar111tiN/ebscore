"""EB scoring against the legacy run_ebscore output, the documented deltas and the scoring options."""

import shutil
import stat
import subprocess

import numpy as np
import pandas as pd
import pytest

import ebscore.fitting as fitting
from conftest import LEGACY, SYN
from ebscore.cli import main
from ebscore.core import eb_score, fit_bb
from ebscore.score import EB_INFO_HEADER, resolve_exclusions, score_variants
from ebscore.variants import add_lookup, read_vcf

KEY = ["Chr", "Start", "End", "Ref", "Alt"]
PAIR = LEGACY / "pair.cpile.tsv.gz"
LEGACY_AB = LEGACY / "ponmatrix" / "AB"  # the AB cache the main code built; test data only (plan D32)


def _parse_ab(cell):
    plus, minus = cell.split("=")
    return tuple(tuple(float(x) for x in s.split("|")) for s in (plus, minus))


def eb_from_legacy_ab(df, fisher="legacy", missing=0.0):
    """EB recomputed the legacy way: the port's tumor counts (``Tumor``) + the legacy AB cache values.

    This is what the legacy run_ebscore did with its AB cache; rows without an AB cache cell or tumor reads get
    ``missing`` (legacy: 0).
    """
    rows = {}
    for f in LEGACY_AB.glob("*.AB.gz"):
        ab = pd.read_csv(f, sep="\t", dtype=str)
        rows.update({(c, int(s)): r for c, s, r in zip(ab["Chr"], ab["Start"], ab.to_dict("records"))})
    look = add_lookup(df[KEY].astype(str))
    out = []
    for (_, v), tumor in zip(look.iterrows(), df["Tumor"]):
        cell = rows.get((v["Chr"], int(v["_pos"]))) if pd.notna(v["_pos"]) else None
        if cell is None or pd.isna(tumor):
            out.append(missing)
            continue
        (tp, tm), (dp, dm) = [x.split("=") for x in tumor.split("<")]
        if int(dp) + int(dm) == 0:
            out.append(missing)
            continue
        out.append(eb_score((int(tp), int(tm)), (int(dp), int(dm)), _parse_ab(cell[v["_letter"]]), fisher))
    return pd.Series(out, index=df.index, dtype="float64")


@pytest.fixture(scope="module")
def legacy_eb():
    df = pd.concat([pd.read_csv(LEGACY / "eb" / f"{c}.EB", sep="\t", dtype=str, keep_default_na=False)
                    for c in ("chr7", "chr12")])
    return df.reset_index(drop=True)


@pytest.fixture(scope="module")
def scenarios():
    return pd.read_csv(SYN / "variants.tsv", sep="\t", dtype=str)


def scored(**kw):
    df = score_variants(kw.pop("variants", SYN / "variants.tsv"), kw.pop("ponmatrix"), pileup=kw.pop("pileup", PAIR),
                        **kw)
    return df.astype({"Start": str, "End": str})


def cli(*args):
    assert main([*map(str, args), "--quiet"]) == 0


def read_out(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def eb_of(df):
    return pd.to_numeric(df["EB"].replace("NA", np.nan))


def drop_column(src, dst, col):
    """Copy a PONmatrix without PON sample ``col`` (0-based): a reference for exclusion."""
    (dst / "matrix").mkdir(parents=True)
    for f in (src / "matrix").iterdir():
        m = pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False)
        for c in [*"AGCTID", "Depth"]:
            m[c] = ["=".join("|".join(v for i, v in enumerate(s.split("|")) if i != col) for s in cell.split("="))
                    for cell in m[c]]
        m.to_csv(dst / "matrix" / f.name, sep="\t", index=False, compression="gzip")
    names = (src / "PON_list.txt").read_text().splitlines()
    (dst / "PON_list.txt").write_text("".join(f"{n}\n" for i, n in enumerate(names) if i != col))
    return dst


# ------------------------------------------------------------------ legacy parity and documented deltas
def test_parity_with_legacy_eb(legacy_ponmatrix, legacy_eb, scenarios):
    # the port's lookups (tumor counts, PON strings) + the legacy AB values reproduce the legacy EB exactly
    new = scored(ponmatrix=legacy_ponmatrix, missing="zero", threads=2, keep_matrix=True)
    new["EB_legacy_ab"] = eb_from_legacy_ab(new)
    m = legacy_eb.merge(new, on=KEY, suffixes=("_leg", ""), how="left")  # new passes `scenario` through
    assert m["EB"].notna().all(), "every legacy row is also scored"
    # legacy bug: the variant right after a position without PON coverage loses its PON line (tumor2matrix.mawk)
    after_gap = m["PON+_leg"].eq("") & m["scenario"].ne("pon_gap")
    assert after_gap.sum() == 2
    assert m.loc[after_gap, "PON+"].str.contains("<").all()
    same = m[~after_gap]
    assert len(same) == 33
    assert (same["EB_legacy_ab"] == same["EB_leg"].astype(float)).all()
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


def test_exact_fits_vs_legacy_ab_values(legacy_ponmatrix):
    fly = scored(ponmatrix=legacy_ponmatrix, threads=2, keep_matrix=True)
    leg = eb_from_legacy_ab(fly, fisher="sf", missing=np.nan)
    both = fly["EB"].notna() & leg.notna()
    assert both.sum() >= 40
    diff = (leg[both] - fly.loc[both, "EB"]).abs()
    assert 0 < diff.max() <= 0.005  # 6 normals at 20-80x: the legacy fit stops close to the optimum here


def test_missing_policies(legacy_ponmatrix):
    na = scored(ponmatrix=legacy_ponmatrix)
    zero = scored(ponmatrix=legacy_ponmatrix, missing="zero")
    for kind in ("pon_gap", "tumor_gap", "mnv"):
        assert na.query("scenario == @kind")["EB"].isna().all()
        assert (zero.query("scenario == @kind")["EB"] == 0).all()


def test_eb_column_is_float(legacy_ponmatrix):
    df = score_variants(SYN / "variants.vcf", legacy_ponmatrix, pileup=PAIR)
    assert df["EB"].dtype == np.float64 and df["EB"].isna().any()
    assert str(df["Start"].dtype) == "Int64" and str(df["End"].dtype) == "Int64"


# ------------------------------------------------------------------ exclusion of PON samples
def test_exclusion_equals_a_pon_without_that_sample(legacy_ponmatrix, tmp_path):
    reduced = drop_column(legacy_ponmatrix, tmp_path / "pon5", 1)  # N02
    ref = scored(ponmatrix=reduced, keep_matrix=True)
    full = scored(ponmatrix=legacy_ponmatrix, keep_matrix=True)
    by_index = scored(ponmatrix=legacy_ponmatrix, exclude_index=[2], keep_matrix=True)
    by_name = scored(ponmatrix=legacy_ponmatrix, exclude_samples=["N02"], keep_matrix=True)
    pd.testing.assert_frame_equal(by_index, by_name)
    pd.testing.assert_frame_equal(by_index, ref)  # same PON strings, same fits, same EB
    changed = (by_index["EB"].fillna(-1) != full["EB"].fillna(-1)).sum()
    assert changed >= 20, "excluding a normal changes most EB values (fit on 5 samples)"

    out = tmp_path / "cli.tsv"
    cli("score", "--variants", SYN / "variants.tsv", "--ponmatrix", legacy_ponmatrix, "--pileup", PAIR, "-o", out,
        "--exclude-sample", "N02")
    assert (eb_of(read_out(out)).fillna(-1).to_numpy() == ref["EB"].fillna(-1).to_numpy()).all()
    cli("score", "--variants", SYN / "variants.tsv", "--ponmatrix", legacy_ponmatrix, "--pileup", PAIR, "-o", out,
        "--exclude-index", 2)
    assert (eb_of(read_out(out)).fillna(-1).to_numpy() == ref["EB"].fillna(-1).to_numpy()).all()


def test_exclusion_by_name_needs_consistent_pon_list():
    with pytest.raises(ValueError, match="unsafe"):
        resolve_exclusions(["a.bam", "b.bam"], 3, names=["a"])
    with pytest.raises(ValueError, match="matches no"):
        resolve_exclusions(["a.bam", "b.bam"], 2, names=["zz"])
    assert resolve_exclusions(["x/N01.bam", "x/N02.bam"], 2, names=["N02"]) == {1}


def test_only_excluded_samples_cover_a_position(large_pon):
    lp = large_pon
    base = dict(variants=lp["variants"], ponmatrix=lp["ponmatrix"], pileup=lp["pileup"])
    full = scored(**base).set_index("Start")
    assert full.loc["1003", "EB"] > 0  # covered by sample 5 only
    excl = scored(**base, exclude_index=[5]).set_index("Start")
    assert np.isnan(excl.loc["1003", "EB"]) and excl.loc["1003", "PON+"].startswith("0|")
    assert excl.drop(index="1003")["EB"].notna().all()
    zero = scored(**base, exclude_index=[5], missing="zero").set_index("Start")
    assert zero.loc["1003", "EB"] == 0.0


# ------------------------------------------------------------------ large PON: Fisher sf vs legacy, penalty
def test_large_pon_reaches_eb_between_16_and_60(large_pon, tmp_path):
    lp = large_pon
    base = dict(variants=lp["variants"], ponmatrix=lp["ponmatrix"], pileup=lp["pileup"])
    sf = scored(**base).set_index("Start")
    leg = scored(**base, fisher="legacy").set_index("Start")
    assert 16 < sf.loc["1001", "EB"] < 60 and leg.loc["1001", "EB"] == 60.0
    assert sf.loc["1002", "EB"] < sf.loc["1001", "EB"]  # PON alt reads lower the EB
    assert sf.loc["1004", "EB"] == leg.loc["1004", "EB"] < 16

    out_sf, out_leg = tmp_path / "sf.tsv", tmp_path / "leg.tsv"
    args = ["score", "--variants", lp["variants"], "--ponmatrix", lp["ponmatrix"], "--pileup", lp["pileup"]]
    cli(*args, "-o", out_sf)
    cli(*args, "-o", out_leg, "--fisher", "legacy")
    assert (eb_of(read_out(out_sf)).to_numpy() == sf["EB"].to_numpy()).all()  # the CLI default is sf
    assert (eb_of(read_out(out_leg)).to_numpy() == leg["EB"].to_numpy()).all()


def test_penalty_reaches_the_fit(large_pon, tmp_path):
    lp = large_pon
    base = dict(variants=lp["variants"], ponmatrix=lp["ponmatrix"], pileup=lp["pileup"], keep_matrix=True)
    default = scored(**base)
    pen = scored(**base, penalty=2.0)
    assert (pen["EB"] != default["EB"]).sum() == 3  # all but 1003 (one PON sample: fits at the bounds)
    for _, row in pen.iterrows():  # EB recomputed from fit_bb(penalty=2.0) on the reported PON vectors
        cms = [np.transpose([[int(x) for x in row[s].split("<")[1].split("|")],
                             [int(x) for x in row[s].split("<")[0].split("|")]]) for s in ("PON+", "PON-")]
        ab = [fit_bb(cm, penalty=2.0) for cm in cms]
        (tp, tm), (dp, dm) = [x.split("=") for x in row["Tumor"].split("<")]
        assert row["EB"] == eb_score((int(tp), int(tm)), (int(dp), int(dm)), ab)
    out = tmp_path / "pen.tsv"
    cli("score", "--variants", lp["variants"], "--ponmatrix", lp["ponmatrix"], "--pileup", lp["pileup"], "-o", out,
        "--penalty", 2)
    assert (eb_of(read_out(out)).to_numpy() == pen["EB"].to_numpy()).all()


# ------------------------------------------------------------------ read-only PONmatrix
def test_read_only_ponmatrix(legacy_ponmatrix, tmp_path):
    pon = tmp_path / "ro"
    shutil.copytree(legacy_ponmatrix, pon)
    paths = [pon, *pon.rglob("*")]
    try:
        for p in paths:
            p.chmod(p.stat().st_mode & ~(stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH))
        a = scored(ponmatrix=pon, tmpdir=tmp_path)
        assert a["EB"].notna().sum() > 40
        pd.testing.assert_frame_equal(a, scored(ponmatrix=legacy_ponmatrix))
        cli("score", "--variants", SYN / "variants.tsv", "--ponmatrix", pon, "--pileup", PAIR,
            "-o", tmp_path / "ro.tsv", "--tmpdir", tmp_path)
    finally:
        for p in paths:
            p.chmod(p.stat().st_mode | stat.S_IWUSR)
    assert sorted(x.relative_to(pon) for x in pon.rglob("*")) == sorted(
        x.relative_to(legacy_ponmatrix) for x in legacy_ponmatrix.rglob("*"))


# ------------------------------------------------------------------ inputs, outputs, CLI defaults
def test_cli_defaults(legacy_ponmatrix, tmp_path):
    args = ["score", "--ponmatrix", legacy_ponmatrix, "--pileup", PAIR]
    tab, vcf = tmp_path / "tab.tsv", tmp_path / "vcf.tsv"
    cli(*args, "--variants", SYN / "variants.tsv", "-o", tab)
    cli(*args, "--variants", SYN / "variants.vcf", "-o", vcf)
    t = read_out(tab)
    for kind in ("pon_gap", "tumor_gap", "mnv"):
        assert (t.loc[t["scenario"] == kind, "EB"] == "NA").all()  # default --missing na
    api = scored(ponmatrix=legacy_ponmatrix)  # defaults: on-the-fly exact fit, sf, na (the AB/ folder is ignored)
    assert (eb_of(t).fillna(-1).to_numpy() == api["EB"].fillna(-1).to_numpy()).all()
    assert (read_out(vcf)["EB"] == "NA").sum() == (t["EB"] == "NA").sum() == 6  # pon_gap, tumor_gap, MNV x 2
    zero = tmp_path / "zero.tsv"
    cli(*args, "--variants", SYN / "variants.tsv", "-o", zero, "--missing", "zero")
    z = read_out(zero)
    assert (z.loc[z["scenario"].isin(["pon_gap", "tumor_gap", "mnv"]), "EB"] == "0.0").all()


@pytest.mark.parametrize("removed", [["--fit", "legacy"], ["--fit", "exact"], ["--no-ab-cache"], ["--ab-cache", "x"]])
def test_removed_score_options_are_rejected(legacy_ponmatrix, tmp_path, removed):
    # plan D19 (exact fit only) and D32 (no AB cache): the options are gone, argparse rejects them
    args = ["score", "--ponmatrix", str(legacy_ponmatrix), "--pileup", str(PAIR), "--variants",
            str(SYN / "variants.tsv"), "-o", str(tmp_path / "o.tsv"), *removed]
    with pytest.raises(SystemExit) as e:
        main(args)
    assert e.value.code == 2


def test_ab_cache_command_is_gone(tmp_path):
    with pytest.raises(SystemExit) as e:
        main(["ab-cache", "--ponmatrix", str(tmp_path), "--chrom", "chr1", "-o", str(tmp_path)])
    assert e.value.code == 2


def test_tumor_column(legacy_ponmatrix, tmp_path):
    out = tmp_path / "r1.tsv"
    cli("score", "--variants", SYN / "variants.tsv", "--ponmatrix", legacy_ponmatrix, "--pileup", PAIR, "-o", out,
        "--tumor-column", "Read1", "--keep-matrix")
    r1 = read_out(out)
    api = scored(ponmatrix=legacy_ponmatrix, tumor_column="Read1", keep_matrix=True)
    tum = scored(ponmatrix=legacy_ponmatrix, keep_matrix=True)
    assert (r1["Tumor"].replace("NA", np.nan).fillna("-").to_numpy() == api["Tumor"].fillna("-").to_numpy()).all()
    assert (api["Tumor"].fillna("-") != tum["Tumor"].fillna("-")).sum() > 40


def test_empty_inputs_give_header_only_outputs(legacy_ponmatrix, tmp_path):
    empty_vcf = tmp_path / "empty.vcf"
    empty_vcf.write_text("".join(line for line in open(SYN / "variants.vcf") if line.startswith("#")))
    header_only = tmp_path / "header.tsv"
    header_only.write_text("Chr\tStart\tEnd\tRef\tAlt\n")
    args = ["score", "--ponmatrix", legacy_ponmatrix, "--pileup", PAIR]
    for variants, extra in ((empty_vcf, []), (header_only, []), (SYN / "variants.tsv", ["--chrom", "chrX"])):
        out, ann = tmp_path / "o.tsv", tmp_path / "o.annot.tsv"
        more = ["--vcf-annotation", ann] if variants == empty_vcf else []
        cli(*args, "--variants", variants, "-o", out, *extra, *more)
        df = read_out(out)
        assert df.empty and list(df.columns[-3:]) == ["PON+", "PON-", "EB"]
        if more:
            assert open(ann).read() == "#CHROM\tPOS\tREF\tALT\tEB\n"


def test_wrong_ponmatrix_fails(legacy_ponmatrix, tmp_path):
    args = ["score", "--variants", str(SYN / "variants.tsv"), "--pileup", str(PAIR), "-o", str(tmp_path / "o.tsv")]
    assert main([*args, "--ponmatrix", str(tmp_path / "nope"), "--quiet"]) == 2
    assert main([*args, "--ponmatrix", str(legacy_ponmatrix / "matrix"), "--quiet"]) == 2  # the subfolder
    other = tmp_path / "other"
    (other / "matrix").mkdir(parents=True)
    assert main([*args, "--ponmatrix", str(other), "--quiet"]) == 2  # no variant chromosome has a file
    assert main([*args, "--ponmatrix", str(other), "--allow-missing-chrom", "--quiet"]) == 0
    assert (read_out(tmp_path / "o.tsv")["EB"] == "NA").all()
    shutil.copy(legacy_ponmatrix / "matrix" / "chr7.pon.gz", other / "matrix")  # partial miss: warning only
    assert main([*args, "--ponmatrix", str(other), "--quiet"]) == 0


def _vcf_records(path):
    """(CHROM, POS, REF, ALT) per ALT allele, read straight from the VCF text."""
    out = []
    for line in open(path):
        if not line.startswith("#"):
            f = line.rstrip("\n").split("\t")
            out += [(f[0], int(f[1]), f[3], alt) for alt in f[4].split(",")]
    return out


def test_vcf_input_and_annotation(legacy_ponmatrix, tmp_path):
    tab = scored(ponmatrix=legacy_ponmatrix)
    vcf = score_variants(SYN / "variants.vcf", legacy_ponmatrix, pileup=PAIR)
    assert (vcf["vcf_alt"].str.contains(",")).sum() == 0 and vcf["vcf_pos"].duplicated().sum() >= 2  # split
    m = tab.merge(vcf.dropna(subset=["Start"]).astype({"Start": str, "End": str}), on=KEY, suffixes=("", "_vcf"))
    assert len(m) == len(tab) - 2  # all but the MNVs
    assert (m["EB"].fillna(-1) == m["EB_vcf"].fillna(-1)).all()

    # shuffled input: the annotation is sorted by chromosome (input order) and numeric POS, VCF coordinates
    lines = open(SYN / "variants.vcf").read().splitlines()
    head, body = [x for x in lines if x.startswith("#")], [x for x in lines if not x.startswith("#")]
    rng = np.random.default_rng(1)
    shuffled = tmp_path / "shuffled.vcf"
    shuffled.write_text("\n".join(head + [body[i] for i in rng.permutation(len(body))]) + "\n")
    out, ann = tmp_path / "eb.tsv", tmp_path / "eb.annot.tsv"
    cli("score", "--variants", shuffled, "--ponmatrix", legacy_ponmatrix, "--pileup", PAIR, "-o", out,
        "--vcf-annotation", ann)
    a = pd.read_csv(ann, sep="\t", dtype={"POS": int})
    assert list(a.columns) == ["#CHROM", "POS", "REF", "ALT", "EB"]
    assert open(f"{ann}.hdr").read() == EB_INFO_HEADER + "\n" and "Number=A" in EB_INFO_HEADER
    first = list(dict.fromkeys(read_vcf(shuffled)["vcf_chrom"]))
    assert list(dict.fromkeys(a["#CHROM"])) == first
    for _, g in a.groupby("#CHROM", sort=False):
        assert g["POS"].is_monotonic_increasing
    scored_vcf = vcf.dropna(subset=["EB"])
    exp = sorted(zip(scored_vcf["vcf_chrom"], scored_vcf["vcf_pos"].astype(int), scored_vcf["vcf_ref"],
                     scored_vcf["vcf_alt"], scored_vcf["EB"]))
    got = sorted(zip(a["#CHROM"], a["POS"], a["REF"], a["ALT"], a["EB"]))
    assert got == exp
    assert set((c, p, r, x) for c, p, r, x, _ in got) <= set(_vcf_records(SYN / "variants.vcf"))
    indels = [g for g in got if len(g[2]) != len(g[3])]
    assert len(indels) >= 10 and all(g[2][0] == g[3][0] for g in indels)  # with their VCF anchor base

    if shutil.which("bcftools") and shutil.which("bgzip") and shutil.which("tabix"):
        sh = lambda c: subprocess.run(["bash", "-o", "pipefail", "-c", c], check=True, capture_output=True,  # noqa
                                      text=True).stdout
        sh(f"bgzip -c {ann} > {ann}.gz && tabix -s1 -b2 -e2 {ann}.gz")
        res = sh(f"bcftools annotate -a {ann}.gz -h {ann}.hdr -c CHROM,POS,REF,ALT,INFO/EB {SYN / 'variants.vcf'}"
                 " | grep -v '^#' | grep SC=multi")
        assert all(len(line.split("EB=")[1].split(",")) == 2 for line in res.splitlines())


def test_vcf2table_annotation_sorts_numerically(legacy_ponmatrix, tmp_path):
    vt = tmp_path / "vt.tsv"
    cli("vcf2table", "--vcf", SYN / "variants.vcf", "-o", vt)
    t = pd.read_csv(vt, sep="\t", dtype=str, keep_default_na=False)
    exp = read_vcf(SYN / "variants.vcf").astype(object)
    exp = exp.where(exp.notna(), "").astype(str)
    assert list(t.columns) == list(exp.columns)  # no index column, tab-separated
    pd.testing.assert_frame_equal(t, exp)
    ann = tmp_path / "vt.annot.tsv"
    cli("score", "--variants", vt, "--ponmatrix", legacy_ponmatrix, "--pileup", PAIR, "-o", tmp_path / "x.tsv",
        "--vcf-annotation", ann)
    a = pd.read_csv(ann, sep="\t")
    for _, g in a.groupby("#CHROM", sort=False):
        assert g["POS"].is_monotonic_increasing  # 5031 < 5073 < 5100, not as text


def test_multiallelic_records_with_right_hand_context(legacy_ponmatrix, scenarios, tmp_path):
    # Mutect2-style records: REF is the longest allele, other ALTs keep its right-hand context (ATT AT,ATTT)
    ref = pd.read_csv(SYN / "ref.tsv", sep="\t", header=None, names=["c", "p", "b"]).set_index(["c", "p"])["b"]
    lines, expect = [], []
    for _, v in scenarios.query("scenario in ['ins', 'del']").iterrows():
        c, s = v["Chr"], int(v["Start"])
        if v["Ref"] == "-":  # insertion after s: REF = anchor + next, ALT = anchor + ins + next, plus an SNV allele
            anchor, nxt = ref[(c, s)], ref[(c, s + 1)]
            snv = next(b for b in "ACGT" if b != anchor)
            lines.append(f"{c}\t{s}\t.\t{anchor}{nxt}\t{anchor}{v['Alt']}{nxt},{snv}{nxt}\t.\tPASS\t.")
        else:  # deletion of s: REF = anchor + del + next, ALT = anchor + next
            anchor, nxt = ref[(c, s - 1)], ref[(c, s + 1)]
            snv = next(b for b in "ACGT" if b != anchor)
            lines.append(f"{c}\t{s - 1}\t.\t{anchor}{v['Ref']}{nxt}\t{anchor}{nxt},{snv}{v['Ref']}{nxt}\t.\tPASS\t.")
        expect.append(tuple(v[KEY]))
    f = tmp_path / "multi.vcf"
    f.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + "\n".join(lines) + "\n")
    got = score_variants(f, legacy_ponmatrix, pileup=PAIR).astype({"Start": str, "End": str})
    tab = scored(ponmatrix=legacy_ponmatrix).set_index(KEY)
    indel = got[got["Ref"].eq("-") | got["Alt"].eq("-")]
    assert sorted(map(tuple, indel[KEY].to_numpy())) == sorted(expect)
    for _, row in indel.iterrows():
        assert row["EB"] == tab.loc[tuple(row[KEY]), "EB"] and not np.isnan(row["EB"])
    assert got["Start"].ne("<NA>").all()  # the SNV alleles with context are scorable too


def test_chrom_restriction(legacy_ponmatrix):
    only7 = scored(ponmatrix=legacy_ponmatrix, chroms=["chr7"])
    full = scored(ponmatrix=legacy_ponmatrix)
    assert set(only7["Chr"]) == {"chr7"}
    pd.testing.assert_frame_equal(only7.reset_index(drop=True),
                                  full[full["Chr"] == "chr7"].reset_index(drop=True))


def test_fit_many_dedups_and_keeps_order(monkeypatch):
    keys = [(10, 0, 20, 0), (10, 0, 20, 1), (10, 0, 20, 0)]
    calls = []
    monkeypatch.setattr(fitting, "fit_bb", lambda cm, pen: calls.append(cm.tolist()) or (float(cm[1, 1]), pen))
    res = fitting.fit_many(keys, 0.7)
    assert calls == [[[10, 0], [20, 0]], [[10, 0], [20, 1]]]
    assert res == {(10, 0, 20, 0): (0.0, 0.7), (10, 0, 20, 1): (1.0, 0.7)}
