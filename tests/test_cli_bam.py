"""End to end through samtools with tiny simulated BAMs: pileup, pon-matrix, ab-cache, score."""

import gzip

import pandas as pd

from ebscore.cli import main


def run(*args):
    assert main([*map(str, args), "--quiet"]) == 0


def test_end_to_end(bam_set, tmp_path):
    b = bam_set
    pon = tmp_path / "ponmatrix"
    run("pon-matrix", "--bam-list", b["pon_list"], "--fasta", b["fasta"], "--bed", b["bed"], "-o", pon, "-t", 2)
    assert sorted(p.name for p in (pon / "matrix").iterdir()) == ["chrA.pon.gz", "chrB.pon.gz"]
    assert (pon / "PON_list.txt").read_text() == b["pon_list"].read_text()  # entries as given, no build paths
    with gzip.open(pon / "matrix" / "chrA.pon.gz", "rt") as fh:
        assert fh.readline().rstrip("\n").split("\t") == ["Chr", "Start", "Ref", *"AGCTID", "Depth"]
        assert fh.readline().split("\t")[3].count("|") == 6  # 4 samples per strand: 3 separators each

    pair1, pair4 = tmp_path / "pair1.pileup.gz", tmp_path / "pair4.pileup.gz"
    run("pileup", "--normal", b["normal"], "--tumor", b["tumor"], "--fasta", b["fasta"], "--bed", b["bed"], "-o", pair1)
    run("pileup", "--normal", b["normal"], "--tumor", b["tumor"], "--fasta", b["fasta"], "--bed", b["bed"],
        "-o", pair4, "-t", 4)
    with gzip.open(pair1, "rt") as f1, gzip.open(pair4, "rt") as f4:
        t1 = f1.read()
        assert t1 == f4.read(), "per-chromosome threads give the same pileup"
    assert t1.split("\n", 1)[0] == "Chr\tStart\tRef\tDepth1\tRead1\tDepth2\tRead2"

    eb_pileup, eb_bam = tmp_path / "eb_pileup.tsv", tmp_path / "eb_bam.tsv"
    run("score", "--variants", b["variants"], "--ponmatrix", pon, "--pileup", pair1, "-o", eb_pileup, "-t", 2)
    run("score", "--variants", b["variants"], "--ponmatrix", pon, "--tumor-bam", b["tumor"], "--fasta", b["fasta"],
        "-o", eb_bam)
    a, c = pd.read_csv(eb_pileup, sep="\t"), pd.read_csv(eb_bam, sep="\t")
    pd.testing.assert_frame_equal(a, c)
    real = a.iloc[: len(b["real"])]
    artifact = a.iloc[len(b["real"]) :]
    # with 4 PON normals at ~20x, EB of a clean 30 % SNV is ~1.4-1.5; the PON-supported artifact must score lower
    assert (real["EB"] > 1).all(), real
    assert (artifact["EB"] < real["EB"].min()).all(), artifact

    run("ab-cache", "--ponmatrix", pon, "-o", pon, "--chrom", "chrA", "--chrom", "chrB", "-t", 2)
    eb_cache = tmp_path / "eb_cache.tsv"
    run("score", "--variants", b["variants"], "--ponmatrix", pon, "--pileup", pair1, "-o", eb_cache)
    pd.testing.assert_frame_equal(pd.read_csv(eb_cache, sep="\t"), a)


def test_cli_errors(tmp_path, capsys):
    assert main(["score", "--variants", str(tmp_path / "nope.tsv"), "--ponmatrix", str(tmp_path),
                 "--pileup", str(tmp_path / "p.gz"), "-o", str(tmp_path / "o.tsv")]) == 2
