"""Tumor-only pileup (plan T1): `ebscore pileup` without --normal, and EB scoring from it."""

import gzip

import pandas as pd
import pytest

import ebscore.pileup as pileup
from conftest import sh
from ebscore.cli import main
from ebscore.pileup import PAIRED_HEADER, TUMOR_HEADER, pileup_header


def run(*args):
    assert main([*map(str, args), "--quiet"]) == 0


def lines(path):
    with gzip.open(path, "rt") as fh:
        return fh.read().splitlines()


def read_pileup(path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def tumor_args(b, out, *extra, tumor=None, bed=None):
    return ["pileup", "--tumor", tumor or b["tumor"], "--fasta", b["fasta"], "--bed", bed or b["bed"], "-o", out,
            *extra]


def pair_args(b, out, *extra, tumor=None, bed=None):
    return ["pileup", "--normal", b["normal"], *tumor_args(b, out, *extra, tumor=tumor, bed=bed)[1:]]


def compare(one_file, pair_file):
    """Join a tumor-only pileup with the tumor columns of a paired one; return the paired-only rows.

    Asserts: the tumor-only file has no position that the paired file lacks, its rows equal the paired tumor columns
    (Ref, Depth2, Read2) and come in the same order, and every paired-only row has tumor ``0 *``.
    """
    one = read_pileup(one_file)
    pair = read_pileup(pair_file)
    assert list(one.columns) == TUMOR_HEADER.split("\t")
    assert list(pair.columns) == PAIRED_HEADER.split("\t")
    tum = pair[["Chr", "Start", "Ref", "Depth2", "Read2"]].set_axis(list(one.columns), axis=1)
    m = tum.merge(one, on=["Chr", "Start"], how="outer", suffixes=("_pair", "_one"), indicator=True, sort=False)
    assert not (m["_merge"] == "right_only").any(), "the tumor-only file has positions the paired file lacks"
    both = m[m["_merge"] == "both"]
    for col in ("Ref", "Depth1", "Read1"):
        assert (both[f"{col}_pair"] == both[f"{col}_one"]).all(), col
    keys_one = list(zip(one["Chr"], one["Start"]))
    in_one = set(keys_one)
    assert keys_one == [k for k in zip(tum["Chr"], tum["Start"]) if k in in_one], "same row order"
    pair_only = m[m["_merge"] == "left_only"]
    assert (pair_only["Depth1_pair"] == "0").all() and (pair_only["Read1_pair"] == "*").all()
    return pair_only


@pytest.fixture(scope="module")
def tumor_part(bam_set, tmp_path_factory):
    """The tumor BAM cut to chrA:1-600 + chrB: chrA from ~660 on is covered by the normal only."""
    b = bam_set
    out = tmp_path_factory.mktemp("tumor_part") / "tumor_part.bam"
    sh(f"samtools view -b -o {out} {b['tumor']} chrA:1-600 chrB")
    sh(f"samtools index {out}")
    return out


@pytest.fixture(scope="module")
def ponmatrix(bam_set, tmp_path_factory):
    b = bam_set
    pon = tmp_path_factory.mktemp("tumor_only_pon") / "ponmatrix"
    run("pon-matrix", "--bam-list", b["pon_list"], "--fasta", b["fasta"], "--bed", b["bed"], "-o", pon, "-t", 2)
    return pon


# ------------------------------------------------------------------ unit level
def test_headers():
    assert pileup_header(1) == TUMOR_HEADER == "Chr\tStart\tRef\tDepth1\tRead1"
    assert pileup_header(2) == PAIRED_HEADER


def test_cli_without_normal_calls_the_tumor_pileup(monkeypatch):
    seen = {}

    def fake_tumor(tumor, fasta, bed, out, chroms, min_mapq, min_baseq, threads, tmpdir):
        seen.update(tumor=tumor, fasta=fasta, bed=bed, out=out, chroms=chroms, mapq=min_mapq, baseq=min_baseq,
                    threads=threads, tmpdir=tmpdir)

    def fake_paired(*a, **kw):
        raise AssertionError("paired_pileup called without --normal")

    monkeypatch.setattr(pileup, "tumor_pileup", fake_tumor)
    monkeypatch.setattr(pileup, "paired_pileup", fake_paired)
    assert main(["pileup", "--tumor", "t", "--fasta", "f", "--bed", "b", "-o", "o", "--min-mapq", "7",
                 "--min-baseq", "9", "--chrom", "c2", "--chrom", "c1", "--chrom", "c2", "-t", "3",
                 "--tmpdir", "s"]) == 0
    assert seen == {"tumor": "t", "fasta": "f", "bed": "b", "out": "o", "chroms": ["c2", "c1"], "mapq": 7, "baseq": 9,
                    "threads": 3, "tmpdir": "s"}
    assert main(["pileup", "--normal", "", "--tumor", "t", "--fasta", "f", "--bed", "b", "-o", "o", "--quiet"]) == 2


# ------------------------------------------------------------------ with samtools
def test_samtools_row_semantics(bam_set, tmp_path):
    """What samtools prints: a position is a row when a read of any input overlaps it; depth counts bases >= -Q.

    Normal reads at chrA 201-260 and 401-460 (1-based), one tumor read at 201-260 whose bases 221-230 have quality 2.
    Tumor-only: 201-260, with 221-230 as ``0 *``. Paired: also 401-460, there with tumor ``0 *``.
    """
    b = bam_set
    seq = b["contigs"]["chrA"]

    def bam(name, reads):
        sam = tmp_path / f"{name}.sam"
        rows = ["@HD\tVN:1.6\tSO:unsorted", *(f"@SQ\tSN:{c}\tLN:{len(s)}" for c, s in b["contigs"].items()),
                f"@RG\tID:{name}\tSM:{name}"]
        for i, (start, qual) in enumerate(reads):
            rows.append(f"r{i}\t0\tchrA\t{start + 1}\t60\t60M\t*\t0\t0\t{seq[start:start + 60]}\t{qual}\tRG:Z:{name}")
        sam.write_text("\n".join(rows) + "\n")
        sh(f"samtools sort -o {tmp_path / name}.bam {sam}")
        sh(f"samtools index {tmp_path / name}.bam")
        return tmp_path / f"{name}.bam"

    normal = bam("n", [(200, "I" * 60), (400, "I" * 60)])
    tumor = bam("t", [(200, "I" * 20 + "#" * 10 + "I" * 30)])
    bed = tmp_path / "w.bed"
    bed.write_text("chrA\t150\t500\n")
    one, pair = tmp_path / "one.gz", tmp_path / "pair.gz"
    run("pileup", "--tumor", tumor, "--fasta", b["fasta"], "--bed", bed, "-o", one)
    run("pileup", "--normal", normal, "--tumor", tumor, "--fasta", b["fasta"], "--bed", bed, "-o", pair)
    o = read_pileup(one)
    assert o["Start"].astype(int).tolist() == list(range(201, 261))
    low = o["Start"].astype(int).between(221, 230)
    assert (o.loc[low, "Depth1"] == "0").all() and (o.loc[low, "Read1"] == "*").all()
    assert (o.loc[~low, "Depth1"] == "1").all()
    pair_only = compare(one, pair)
    assert sorted(pair_only["Start"].astype(int)) == list(range(401, 461))


def test_tumor_only_equals_the_paired_tumor_columns(bam_set, tumor_part, tmp_path):
    b = bam_set
    for name, tumor in (("full", b["tumor"]), ("part", tumor_part)):
        one, pair = tmp_path / f"{name}.one.gz", tmp_path / f"{name}.pair.gz"
        run(*tumor_args(b, one, tumor=tumor))
        run(*pair_args(b, pair, tumor=tumor))
        pair_only = compare(one, pair)
        if name == "part":
            # the normal-only region: in the paired file with tumor 0 *, absent from the tumor-only file
            assert len(pair_only) > 500 and set(pair_only["Chr"]) == {"chrA"}
            assert pair_only["Start"].astype(int).min() > 600
        assert len(read_pileup(one)) > (2500 if name == "full" else 1800)  # the BED has 2 x 1300 positions


def test_cli_end_to_end(bam_set, tmp_path):
    b = bam_set
    t1, t2 = tmp_path / "t1.pileup.gz", tmp_path / "t2.pileup.gz"
    run(*tumor_args(b, t1))
    run(*tumor_args(b, t2, "-t", 2))
    for f in (t1, t2):
        assert f.read_bytes()[:2] == b"\x1f\x8b", "gzip"
    assert gzip.decompress(t1.read_bytes()) == gzip.decompress(t2.read_bytes()), "-t 2 gives the same pileup as -t 1"
    rows = lines(t1)
    assert rows[0] == TUMOR_HEADER
    assert list(dict.fromkeys(r.split("\t")[0] for r in rows[1:])) == ["chrA", "chrB"]
    assert all(len(r.split("\t")) == 5 for r in rows)

    only_b = tmp_path / "chrB.pileup.gz"
    run(*tumor_args(b, only_b, "--chrom", "chrB", "--chrom", "chrB"))
    assert lines(only_b) == rows[:1] + [r for r in rows[1:] if r.startswith("chrB\t")]

    bed_a = tmp_path / "a.bed"
    bed_a.write_text("chrA\t100\t1400\n")
    empty = tmp_path / "empty.pileup.gz"
    run(*tumor_args(b, empty, "--chrom", "chrB", bed=bed_a))
    assert lines(empty) == [TUMOR_HEADER]
    assert main([*map(str, tumor_args(b, tmp_path / "x.gz", bed=tmp_path / "none.bed")), "--quiet"]) == 2


def test_eb_same_via_tumor_only_pileup_paired_pileup_and_tumor_bam(bam_set, tumor_part, ponmatrix, tmp_path):
    b = bam_set
    for name, tumor in (("full", b["tumor"]), ("part", tumor_part)):
        one, pair = tmp_path / f"{name}.one.gz", tmp_path / f"{name}.pair.gz"
        run(*tumor_args(b, one, tumor=tumor))
        run(*pair_args(b, pair, tumor=tumor))
        outs = {k: tmp_path / f"{name}.{k}.tsv" for k in ("one", "one_read1", "pair", "bam")}
        common = ["score", "--variants", b["variants"], "--ponmatrix", ponmatrix, "--keep-matrix"]
        run(*common, "--pileup", one, "-o", outs["one"], "-t", 2)
        run(*common, "--pileup", one, "--tumor-column", "Read1", "-o", outs["one_read1"])
        run(*common, "--pileup", pair, "-o", outs["pair"])
        run(*common, "--tumor-bam", tumor, "--fasta", b["fasta"], "-o", outs["bam"])
        res = {k: pd.read_csv(v, sep="\t") for k, v in outs.items()}
        for k in ("one_read1", "bam"):
            pd.testing.assert_frame_equal(res["one"], res[k], obj=f"{name}: tumor-only pileup vs {k}")
        # where only the normal covers a variant, the paired route sees tumor depth 0 (Tumor 0=0<0=0) and the other
        # two see no tumor row (Tumor NA); PON+, PON-, AB and EB are the same
        pd.testing.assert_frame_equal(res["one"].drop(columns="Tumor"), res["pair"].drop(columns="Tumor"),
                                      obj=f"{name}: tumor-only vs paired pileup")
        has = res["one"]["Tumor"].notna()
        assert res["one"].loc[has, "Tumor"].equals(res["pair"].loc[has, "Tumor"])
        assert res["pair"].loc[~has, "Tumor"].tolist() == ([] if name == "full" else ["0=0<0=0"])
        eb = res["one"]
        n_real = len(b["real"])
        if name == "full":
            assert (eb["EB"].iloc[:n_real] > 1).all() and eb["EB"].notna().all()
        else:
            # chrA:901 (real SNV 2) lies in the normal-only region: no tumor data on every route -> EB NA
            assert eb["EB"].isna().tolist() == [False, True, False, False]
        bad = [*common, "--pileup", one, "--tumor-column", "Read2", "-o", tmp_path / "x.tsv", "--quiet"]
        assert main([str(a) for a in bad]) == 2  # a tumor-only pileup has no Read2
