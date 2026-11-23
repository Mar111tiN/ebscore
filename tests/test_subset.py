"""Run-level PONmatrix subset (plan T2): `ebscore pon-subset`, and `score` against a subset (bit-identical, guarded)."""

import gzip
import hashlib
import json
import logging
import shutil

import pandas as pd
import pytest

import ebscore.score as score_mod
import ebscore.subset as subset_mod
from conftest import LEGACY, SYN, read_text
from ebscore import __version__
from ebscore.cli import main
from ebscore.score import score_variants
from ebscore.subset import SUBSET_SCHEMA, pon_subset

PAIR = LEGACY / "pair.cpile.tsv.gz"
TSV, VCF = SYN / "variants.tsv", SYN / "variants.vcf"


def cli(*args):
    assert main([*map(str, args), "--quiet"]) == 0


def cli_fails(*args):
    assert main([*map(str, args), "--quiet"]) == 2


def score_cli(ponmatrix, variants, out, *extra, pileup=PAIR):
    """`ebscore score --keep-matrix` (Tumor and AB too); returns the output bytes."""
    cli("score", "--variants", variants, "--ponmatrix", ponmatrix, "--pileup", pileup, "-o", out, "--keep-matrix",
        *extra)
    return out.read_bytes()


def table(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def write_table(df, path):
    df.to_csv(path, sep="\t", index=False)
    return path


def listed(subset_dir):
    text = (subset_dir / "positions.tsv").read_text()
    return [tuple(line.split("\t")) for line in text.splitlines()]


def files_of(d):
    return sorted(str(p.relative_to(d)) for p in d.rglob("*"))


def all_bytes(d):
    return {str(p.relative_to(d)): p.read_bytes() for p in sorted(d.rglob("*")) if p.is_file()}


def full_rows_at(full, chrom, keys):
    """Header + the rows of the full matrix file at ``keys`` ((Chr, "Pos") strings), in source order."""
    lines = read_text(full / "matrix" / f"{chrom}.pon.gz").splitlines()
    return [lines[0]] + [x for x in lines[1:] if tuple(x.split("\t")[:2]) in keys]


@pytest.fixture(scope="module")
def split_variants(tmp_path_factory):
    """variants.tsv in two overlapping parts: a = chr7 + the first 3 chr12 rows, b = chr12 + the first 2 chr7 rows."""
    d = tmp_path_factory.mktemp("split_variants")
    v = table(TSV)
    c7, c12 = v[v["Chr"] == "chr7"], v[v["Chr"] == "chr12"]
    a = write_table(pd.concat([c7, c12.head(3)]), d / "a.tsv")
    b = write_table(pd.concat([c12, c7.head(2)]), d / "b.tsv")
    return a, b


# ------------------------------------------------------------------ bit-identical EB against the subset
@pytest.mark.parametrize("extra", [[], ["--exclude-sample", "N02"], ["--exclude-index", "2"], ["--missing", "zero"]],
                         ids=["default", "exclude-sample", "exclude-index", "missing-zero"])
def test_score_against_subset_is_identical(legacy_ponmatrix, tmp_path, extra):
    sub = tmp_path / "sub"
    pon_subset(legacy_ponmatrix, sub, variants=[TSV])
    full = score_cli(legacy_ponmatrix, TSV, tmp_path / "full.tsv", *extra)
    assert score_cli(sub, TSV, tmp_path / "sub.tsv", *extra) == full  # PON+ PON- EB Tumor AB, every byte
    df = table(tmp_path / "sub.tsv")
    assert (df["AB"] != "NA").sum() > 40 and (df["EB"] != "NA").sum() > 40
    n_full = sum(len(read_text(f).splitlines()) - 1 for f in (legacy_ponmatrix / "matrix").glob("*.pon.gz"))
    n_sub = sum(len(read_text(f).splitlines()) - 1 for f in (sub / "matrix").glob("*.pon.gz"))
    assert n_sub == 49 < n_full == 398


def test_score_api_frames_are_exactly_equal(legacy_ponmatrix, tmp_path):
    sub = tmp_path / "sub"
    pon_subset(legacy_ponmatrix, sub, variants=[VCF])
    for kw in ({}, {"exclude_samples": ["N02"]}, {"exclude_index": [3], "fisher": "legacy"}):
        a = score_variants(VCF, legacy_ponmatrix, pileup=PAIR, keep_matrix=True, **kw)
        b = score_variants(VCF, sub, pileup=PAIR, keep_matrix=True, threads=2, **kw)
        pd.testing.assert_frame_equal(a, b, check_exact=True)
    ann_full, ann_sub = tmp_path / "full.annot.tsv", tmp_path / "sub.annot.tsv"
    full = score_cli(legacy_ponmatrix, VCF, tmp_path / "full.tsv", "--vcf-annotation", ann_full)
    assert score_cli(sub, VCF, tmp_path / "sub.tsv", "--vcf-annotation", ann_sub) == full
    assert ann_sub.read_bytes() == ann_full.read_bytes()


def test_large_pon_subset_is_identical(large_pon, tmp_path):
    lp = large_pon
    sub = tmp_path / "sub"
    pon_subset(lp["ponmatrix"], sub, variants=[lp["variants"]])
    for kw in ({}, {"exclude_index": [5]}, {"fisher": "legacy", "penalty": 2.0}):
        a = score_variants(lp["variants"], lp["ponmatrix"], pileup=lp["pileup"], keep_matrix=True, **kw)
        b = score_variants(lp["variants"], sub, pileup=lp["pileup"], keep_matrix=True, **kw)
        pd.testing.assert_frame_equal(a, b, check_exact=True)
    assert json.loads((sub / "subset.json").read_text())["n_samples"] == 32


def test_positions_are_what_score_looks_up(legacy_ponmatrix, tmp_path, monkeypatch):
    # record the positions file score writes for its lookups; the subset lists exactly these
    seen = []
    orig = score_mod.tumor_counts_from_pileup

    def spy(pileup, positions_file, tumor_column=None):
        seen.append(set(tuple(x.split("\t")) for x in open(positions_file).read().splitlines()))
        return orig(pileup, positions_file, tumor_column)

    monkeypatch.setattr(score_mod, "tumor_counts_from_pileup", spy)
    for variants in (TSV, VCF):
        sub = tmp_path / variants.name
        pon_subset(legacy_ponmatrix, sub, variants=[variants])
        score_variants(variants, legacy_ponmatrix, pileup=PAIR)
        assert set(listed(sub)) == seen[-1] and len(listed(sub)) == len(seen[-1]) == 51


# ------------------------------------------------------------------ one run-level subset for several variant files
def test_one_subset_serves_two_variant_files(legacy_ponmatrix, split_variants, tmp_path):
    a, b = split_variants
    sub = tmp_path / "run"
    cli("pon-subset", "--ponmatrix", legacy_ponmatrix, "-o", sub, "--variants", a, "--variants", b, "-t", 2)
    for v in (a, b):
        full = score_cli(legacy_ponmatrix, v, tmp_path / f"{v.stem}.full.tsv")
        assert score_cli(sub, v, tmp_path / f"{v.stem}.sub.tsv", "-t", 2) == full
    single = tmp_path / "single"
    pon_subset(legacy_ponmatrix, single, variants=[TSV])  # a + b = variants.tsv
    assert sorted(listed(sub)) == sorted(listed(single))
    info = json.loads((sub / "subset.json").read_text())
    assert [(i["kind"], i["path"]) for i in info["inputs"]] == [("variants", str(a)), ("variants", str(b))]
    assert info["n_positions"] == 51 and sum(i["positions"] for i in info["inputs"]) > 51  # a and b overlap


# ------------------------------------------------------------------ the guard
def test_guard_refuses_positions_outside_the_subset(legacy_ponmatrix, split_variants, tmp_path):
    a, b = split_variants
    sub = tmp_path / "a_only"
    pon_subset(legacy_ponmatrix, sub, variants=[a])
    pon_subset(legacy_ponmatrix, tmp_path / "b_sub", variants=[b])
    in_a = set(listed(sub))
    b_pos = [p for p in listed(tmp_path / "b_sub") if p not in in_a]  # b's order: chr12 first
    assert len(b_pos) == 22
    with pytest.raises(ValueError, match="22 scoring positions outside this PONmatrix subset") as e:
        score_variants(b, sub, pileup=PAIR)
    assert f"{b_pos[0][0]}:{b_pos[0][1]}" in str(e.value) and "rebuild it with these variants" in str(e.value)
    cli_fails("score", "--variants", b, "--ponmatrix", sub, "--pileup", PAIR, "-o", tmp_path / "x.tsv")
    assert not (tmp_path / "x.tsv").exists()
    cli_fails("score", "--variants", TSV, "--ponmatrix", sub, "--tumor-bam", "nope.bam", "--fasta", "nope.fa",
              "-o", tmp_path / "x.tsv")  # before any tumor lookup


def test_listed_position_without_a_source_row_is_fine(legacy_ponmatrix, tmp_path):
    sub = tmp_path / "sub"
    pon_subset(legacy_ponmatrix, sub, variants=[TSV])
    gaps = table(TSV).query("scenario == 'pon_gap'")
    keys = {tuple(x) for x in gaps[["Chr", "Start"]].to_numpy()}
    assert keys <= set(listed(sub))  # listed ...
    rows = {tuple(x.split("\t")[:2]) for f in (sub / "matrix").glob("*.pon.gz") for x in read_text(f).splitlines()}
    assert not keys & rows  # ... but no matrix row, in the full set neither
    info = json.loads((sub / "subset.json").read_text())
    assert {c: (v["positions"], v["rows"]) for c, v in info["chroms"].items()} == {"chr7": (26, 25), "chr12": (25, 24)}
    score_cli(sub, TSV, tmp_path / "o.tsv")
    out = table(tmp_path / "o.tsv")
    assert (out.loc[out["scenario"] == "pon_gap", "EB"] == "NA").all()


def test_guard_checks_the_subset_files(legacy_ponmatrix, tmp_path):
    sub = tmp_path / "sub"
    pon_subset(legacy_ponmatrix, sub, variants=[TSV])
    pos = (sub / "positions.tsv").read_bytes()
    (sub / "positions.tsv").write_bytes(pos + b"chr7\t5\n")  # edited: md5 mismatch
    with pytest.raises(ValueError, match="does not match subset.json"):
        score_variants(TSV, sub, pileup=PAIR)
    (sub / "positions.tsv").unlink()
    with pytest.raises(ValueError, match="incomplete PONmatrix subset"):
        score_variants(TSV, sub, pileup=PAIR)
    (sub / "positions.tsv").write_bytes(pos)
    info = json.loads((sub / "subset.json").read_text())
    (sub / "subset.json").write_text(json.dumps({**info, "schema": "ebscore-ponmatrix-subset/9"}))
    with pytest.raises(ValueError, match="unknown subset schema"):
        score_variants(TSV, sub, pileup=PAIR)
    (sub / "subset.json").write_text("{")
    with pytest.raises(ValueError, match="not valid JSON"):
        score_variants(TSV, sub, pileup=PAIR)


# ------------------------------------------------------------------ inputs
def test_empty_and_header_only_inputs(legacy_ponmatrix, tmp_path):
    empty_vcf = tmp_path / "empty.vcf"
    empty_vcf.write_text("".join(x for x in open(VCF) if x.startswith("#")))
    header_only = tmp_path / "header.tsv"
    header_only.write_text("Chr\tStart\tEnd\tRef\tAlt\n")
    zero_bytes = tmp_path / "zero.pos"
    zero_bytes.write_text("")
    for i, (opt, path) in enumerate((("--variants", empty_vcf), ("--variants", header_only),
                                     ("--positions", zero_bytes))):
        sub = tmp_path / f"sub{i}"
        cli("pon-subset", "--ponmatrix", legacy_ponmatrix, "-o", sub, opt, path)
        assert files_of(sub) == ["PON_list.txt", "matrix", "positions.tsv", "subset.json"]
        assert (sub / "positions.tsv").read_bytes() == b""
        info = json.loads((sub / "subset.json").read_text())
        assert (info["n_positions"], info["n_samples"], info["chroms"]) == (0, None, {})
        if opt == "--variants":  # the same empty input scores to a header-only output
            score_cli(sub, path, tmp_path / "o.tsv")
            out = table(tmp_path / "o.tsv")
            assert out.empty and list(out.columns[-5:]) == ["PON+", "PON-", "EB", "Tumor", "AB"]
        cli_fails("score", "--variants", TSV, "--ponmatrix", sub, "--pileup", PAIR, "-o", tmp_path / "x.tsv")


def test_unscorable_rows_are_ignored_and_counted(legacy_ponmatrix, tmp_path, caplog):
    with caplog.at_level(logging.INFO, logger="ebscore"):
        info = pon_subset(legacy_ponmatrix, tmp_path / "sub", variants=[TSV, VCF])
    assert "55 variants, 51 positions to look up, 2 unscorable rows ignored" in caplog.text
    assert [(i["rows"], i["unscorable"], i["positions"]) for i in info["inputs"]] == [(55, 2, 51), (55, 2, 51)]
    assert info["n_positions"] == 51  # the MNVs have no position


def test_positions_file_input(legacy_ponmatrix, tmp_path):
    ref = tmp_path / "ref"
    pon_subset(legacy_ponmatrix, ref, variants=[TSV])
    again = tmp_path / "again"
    cli("pon-subset", "--ponmatrix", legacy_ponmatrix, "-o", again, "--positions", ref / "positions.tsv")
    for name in ("matrix/chr7.pon.gz", "matrix/chr12.pon.gz", "positions.tsv", "PON_list.txt"):
        assert (again / name).read_bytes() == (ref / name).read_bytes(), name
    assert score_cli(again, TSV, tmp_path / "a.tsv") == score_cli(legacy_ponmatrix, TSV, tmp_path / "f.tsv")
    bad = {"header": "Chr\tPos\nchr7\t1020\n", "bed": "chr7\t1019\t1020\n", "zero": "chr7\t0\n",
           "text": "chr7\tx\n", "nochr": "\t1020\n"}
    for name, text in bad.items():
        p = tmp_path / f"{name}.pos"
        p.write_text(text)
        with pytest.raises(ValueError):
            pon_subset(legacy_ponmatrix, tmp_path / f"bad_{name}", positions=[p])
        assert not (tmp_path / f"bad_{name}" / "matrix").exists()


# ------------------------------------------------------------------ outputs
def test_cli_end_to_end_and_subset_json(legacy_ponmatrix, split_variants, tmp_path):
    a, _ = split_variants
    src = tmp_path / "full"
    shutil.copytree(legacy_ponmatrix, src)
    manifest = {"pon_set": "synthetic", "n_normals_in_data": 6}
    (src / "manifest.json").write_text(json.dumps(manifest))
    extra = tmp_path / "extra.pos"  # chr12 first; 999 < 1020 numerically, not as text; 999 has no matrix row
    extra.write_text("chr12\t5130\nchr7\t1020\nchr7\t999\nchr12\t5011\n")
    sub = tmp_path / "sub"
    cli("pon-subset", "--ponmatrix", src, "-o", sub, "--variants", a, "--positions", extra, "-t", 2)
    assert files_of(sub) == ["PON_list.txt", "matrix", "matrix/chr12.pon.gz", "matrix/chr7.pon.gz",
                             "positions.tsv", "subset.json"]  # nothing else: no staging left-overs, no AB/
    assert (sub / "PON_list.txt").read_bytes() == (src / "PON_list.txt").read_bytes()

    pos = listed(sub)
    order = list(dict.fromkeys(c for c, _ in pos))
    assert order == ["chr7", "chr12"]  # first appearance: a (chr7 first), then the positions file
    for c in order:
        p = [int(x) for cc, x in pos if cc == c]
        assert p == sorted(set(p))
    assert ("chr7", "999") in pos and pos[0] == ("chr7", "999")

    info = json.loads((sub / "subset.json").read_text())
    assert info["schema"] == SUBSET_SCHEMA == "ebscore-ponmatrix-subset/1" and info["ebscore"] == __version__
    assert info["source"] == {"path": str(src), "manifest": manifest, "subset": None}
    assert info["chrom_filter"] is None
    assert info["n_positions"] == len(pos) and info["n_samples"] == 6
    assert info["positions_md5"] == hashlib.md5((sub / "positions.tsv").read_bytes()).hexdigest()
    keys = set(pos)
    for c in order:
        expect = full_rows_at(src, c, keys)
        got = read_text(sub / "matrix" / f"{c}.pon.gz").splitlines()
        assert got == expect  # header + the source lines unchanged, in source order
        assert info["chroms"][c] == {"positions": sum(1 for cc, _ in pos if cc == c), "rows": len(expect) - 1,
                                     "matrix": True}
    assert list(info["chroms"]) == order
    assert files_of(src) == sorted(files_of(legacy_ponmatrix) + ["manifest.json"])  # the source is untouched


def test_threads_give_identical_files(legacy_ponmatrix, split_variants, tmp_path):
    a, b = split_variants
    outs = []
    for t in (1, 2):
        out = tmp_path / f"t{t}"
        cli("pon-subset", "--ponmatrix", legacy_ponmatrix, "-o", out, "--variants", a, "--variants", b, "-t", t)
        outs.append(all_bytes(out))
    assert outs[0] == outs[1] and len(outs[0]) == 5


def test_outdir_and_source_checks(legacy_ponmatrix, tmp_path):
    before = files_of(legacy_ponmatrix)
    base = ["pon-subset", "--variants", TSV]
    cli_fails(*base, "--ponmatrix", legacy_ponmatrix, "-o", legacy_ponmatrix)
    link = tmp_path / "link"
    link.symlink_to(legacy_ponmatrix)
    cli_fails(*base, "--ponmatrix", legacy_ponmatrix, "-o", link)  # the source under another name
    cli_fails(*base, "--ponmatrix", legacy_ponmatrix, "-o", legacy_ponmatrix / "matrix")
    cli_fails(*base, "--ponmatrix", legacy_ponmatrix / "matrix", "-o", tmp_path / "x")  # not a PONmatrix folder
    cli_fails("pon-subset", "--ponmatrix", legacy_ponmatrix, "-o", tmp_path / "x")  # no input
    assert files_of(legacy_ponmatrix) == before
    sub = tmp_path / "sub"
    cli(*base, "--ponmatrix", legacy_ponmatrix, "-o", sub)
    cli_fails(*base, "--ponmatrix", legacy_ponmatrix, "-o", sub)  # matrix/ exists
    with pytest.raises(ValueError, match="already exists"):
        pon_subset(legacy_ponmatrix, sub, variants=[TSV])


def test_interrupted_run_leaves_no_matrix(legacy_ponmatrix, tmp_path, monkeypatch):
    sub = tmp_path / "sub"

    def boom(path):
        raise RuntimeError("killed")

    monkeypatch.setattr(subset_mod, "_matrix_stats", boom)
    with pytest.raises(RuntimeError):
        pon_subset(legacy_ponmatrix, sub, variants=[TSV], threads=2)
    assert files_of(sub) == []  # no matrix/, no staging folder: nothing can score against it
    monkeypatch.undo()
    pon_subset(legacy_ponmatrix, sub, variants=[TSV])
    assert score_cli(sub, TSV, tmp_path / "s.tsv") == score_cli(legacy_ponmatrix, TSV, tmp_path / "f.tsv")


# ------------------------------------------------------------------ run level -> per pair, --chrom, missing chromosomes
def test_subset_of_a_subset(legacy_ponmatrix, split_variants, tmp_path):
    a, b = split_variants
    run_dir, pair, direct = tmp_path / "run", tmp_path / "pair_a", tmp_path / "direct_a"
    run_info = pon_subset(legacy_ponmatrix, run_dir, variants=[a, b], threads=2)
    info = pon_subset(run_dir, pair, variants=[a])  # the cheap per-pair step: reads only the run subset
    pon_subset(legacy_ponmatrix, direct, variants=[a])
    for name in ("matrix/chr7.pon.gz", "matrix/chr12.pon.gz", "positions.tsv", "PON_list.txt"):
        assert (pair / name).read_bytes() == (direct / name).read_bytes(), name
    assert info["source"]["subset"] == run_info and info["source"]["path"] == str(run_dir)
    assert score_cli(pair, a, tmp_path / "p.tsv") == score_cli(legacy_ponmatrix, a, tmp_path / "f.tsv")
    with pytest.raises(ValueError, match="requested positions outside this PONmatrix subset"):
        pon_subset(direct, tmp_path / "pair_b", variants=[b])  # b has positions the a-only subset lacks
    assert not (tmp_path / "pair_b" / "matrix").exists()


def test_chrom_restriction(legacy_ponmatrix, tmp_path):
    sub = tmp_path / "sub"
    cli("pon-subset", "--ponmatrix", legacy_ponmatrix, "-o", sub, "--variants", TSV, "--chrom", "chr7")
    assert files_of(sub / "matrix") == ["chr7.pon.gz"]
    info = json.loads((sub / "subset.json").read_text())
    assert info["chrom_filter"] == ["chr7"] and list(info["chroms"]) == ["chr7"]
    full = score_cli(legacy_ponmatrix, TSV, tmp_path / "f.tsv", "--chrom", "chr7")
    assert score_cli(sub, TSV, tmp_path / "s.tsv", "--chrom", "chr7") == full
    with pytest.raises(ValueError, match="25 scoring positions outside"):  # the chr12 positions
        score_variants(TSV, sub, pileup=PAIR)


def test_chromosome_without_a_matrix_file(legacy_ponmatrix, tmp_path, caplog):
    v = table(TSV)
    on_x = v.head(4).assign(Chr="chrX")
    mixed = write_table(pd.concat([v, on_x]), tmp_path / "mixed.tsv")
    only_x = write_table(on_x, tmp_path / "x.tsv")
    sub = tmp_path / "sub"
    with caplog.at_level(logging.WARNING, logger="ebscore"):
        info = pon_subset(legacy_ponmatrix, sub, variants=[mixed])
    assert "no PONmatrix file for chrX" in caplog.text
    assert info["chroms"]["chrX"] == {"positions": 3, "rows": 0, "matrix": False}
    assert files_of(sub / "matrix") == ["chr12.pon.gz", "chr7.pon.gz"]
    assert score_cli(sub, mixed, tmp_path / "s.tsv") == score_cli(legacy_ponmatrix, mixed, tmp_path / "f.tsv")

    # only chromosomes without a matrix file: fails as score against the full set does
    cli_fails("pon-subset", "--ponmatrix", legacy_ponmatrix, "-o", tmp_path / "x", "--variants", only_x)
    cli_fails("score", "--variants", only_x, "--ponmatrix", legacy_ponmatrix, "--pileup", PAIR,
              "-o", tmp_path / "x.tsv")
    cli("pon-subset", "--ponmatrix", legacy_ponmatrix, "-o", tmp_path / "x", "--variants", only_x,
        "--allow-missing-chrom")
    assert files_of(tmp_path / "x" / "matrix") == []
    full = score_cli(legacy_ponmatrix, only_x, tmp_path / "xf.tsv", "--allow-missing-chrom")
    assert score_cli(tmp_path / "x", only_x, tmp_path / "xs.tsv", "--allow-missing-chrom") == full


def test_plain_matrix_source_and_corrupt_sources(legacy_ponmatrix, tmp_path):
    src = tmp_path / "plain"
    shutil.copytree(legacy_ponmatrix, src)
    (src / "matrix" / "chr7.pon").write_text(read_text(src / "matrix" / "chr7.pon.gz"))
    (src / "matrix" / "chr7.pon.gz").unlink()
    sub, ref = tmp_path / "sub", tmp_path / "ref"
    pon_subset(src, sub, variants=[TSV])
    pon_subset(legacy_ponmatrix, ref, variants=[TSV])
    assert (sub / "matrix" / "chr7.pon.gz").read_bytes() == (ref / "matrix" / "chr7.pon.gz").read_bytes()
    assert score_cli(sub, TSV, tmp_path / "s.tsv") == score_cli(src, TSV, tmp_path / "f.tsv")

    # a matrix file without the PONmatrix header, and one with a different sample count
    lines = read_text(src / "matrix" / "chr7.pon").splitlines()
    (src / "matrix" / "chr7.pon").write_text("\n".join(lines[1:]) + "\n")
    with pytest.raises(ValueError, match="not a PONmatrix file"):
        pon_subset(src, tmp_path / "bad1", variants=[TSV])
    five = [lines[0]] + ["\t".join(f[:3] + ["=".join("|".join(s.split("|")[:5]) for s in c.split("="))
                                             for c in f[3:]]) for f in (x.split("\t") for x in lines[1:])]
    (src / "matrix" / "chr7.pon").write_text("\n".join(five) + "\n")
    with pytest.raises(ValueError, match="different sample counts"):
        pon_subset(src, tmp_path / "bad2", variants=[TSV])
    for d in ("bad1", "bad2"):
        assert not (tmp_path / d / "matrix").exists()


def test_subset_files_are_gzip_with_header(legacy_ponmatrix, tmp_path):
    sub = tmp_path / "sub"
    pon_subset(legacy_ponmatrix, sub, positions=[_one_position(tmp_path)])
    with gzip.open(sub / "matrix" / "chr12.pon.gz", "rt") as fh:
        assert fh.read().splitlines() == full_rows_at(legacy_ponmatrix, "chr12", {("chr12", "5115")})
    assert json.loads((sub / "subset.json").read_text())["chroms"]["chr12"]["rows"] == 0  # pon_gap: header only


def _one_position(tmp_path):
    p = tmp_path / "one.pos"
    p.write_text("chr12\t5115\n")
    return p
