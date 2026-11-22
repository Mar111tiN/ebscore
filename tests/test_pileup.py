"""samtools front ends: chromosome planning, order, headers, flags, BED formats, PON list handling."""

import gzip
import logging
import os

import pandas as pd
import pytest

import ebscore.pileup as pileup
from ebscore.cli import main
from ebscore.pileup import MATRIX_HEADER, PAIRED_HEADER, mpileup_cmd, plan_chroms, read_bam_list
from ebscore.tools import bed_chroms


def run(*args):
    assert main([*map(str, args), "--quiet"]) == 0


def lines(path):
    with gzip.open(path, "rt") as fh:
        return fh.read().splitlines()


def pair_args(b, bed, out, *extra):
    return ["pileup", "--normal", b["normal"], "--tumor", b["tumor"], "--fasta", b["fasta"], "--bed", bed, "-o", out,
            *extra]


# ------------------------------------------------------------------ unit level
def test_mpileup_flags():
    cmd = mpileup_cmd(["n.bam", "t.bam"], "ref.fa", "kit.bed", "chr7", 21, 26)
    assert cmd == "samtools mpileup -f ref.fa -q 21 -Q 26 -l kit.bed -r chr7 n.bam t.bam"
    assert mpileup_cmd(None, "ref.fa", "kit.bed", bam_list="l.txt") == \
        "samtools mpileup -f ref.fa -q 20 -Q 25 -l kit.bed -b l.txt"


def test_cli_passes_mapq_baseq_and_chrom(monkeypatch, tmp_path):
    seen = {}

    def fake(normal, tumor, fasta, bed, out, chroms, min_mapq, min_baseq, threads, tmpdir):
        seen.update(chroms=chroms, mapq=min_mapq, baseq=min_baseq, threads=threads)

    monkeypatch.setattr(pileup, "paired_pileup", fake)
    assert main(["pileup", "--normal", "n", "--tumor", "t", "--fasta", "f", "--bed", "b", "-o", "o",
                 "--min-mapq", "7", "--min-baseq", "9", "--chrom", "c2", "--chrom", "c1", "--chrom", "c2",
                 "-t", "3"]) == 0
    assert seen == {"chroms": ["c2", "c1"], "mapq": 7, "baseq": 9, "threads": 3}


def test_bed_chroms_reads_gzip_and_whitespace(tmp_path):
    (tmp_path / "a.bed").write_text("track name=x\n#c\nchr2 1 5\nchr1\t1\t5\nchr2\t7\t9\n")
    with gzip.open(tmp_path / "a.bed.gz", "wt") as fh:
        fh.write("browser x\nchr3\t1\t5\nchr1\t1\t5\n")
    assert bed_chroms(tmp_path / "a.bed") == ["chr2", "chr1"]
    assert bed_chroms(tmp_path / "a.bed.gz") == ["chr3", "chr1"]


def test_plan_chroms_orders_by_header_and_skips_missing(monkeypatch, tmp_path, caplog):
    headers = {"n.bam": ["chr1", "chr2", "chr10", "chrX"], "t.bam": ["chr1", "chr2", "chr10"]}
    monkeypatch.setattr(pileup, "bam_contigs", lambda bam: headers[bam])
    bed = tmp_path / "kit.bed"
    bed.write_text("chr1\t1\t2\nchr10\t1\t2\nchr2\t1\t2\nchrX\t1\t2\nchrZ\t1\t2\n")  # lexical vendor order
    with caplog.at_level(logging.WARNING, logger="ebscore"):
        assert plan_chroms(bed, ["n.bam", "t.bam"]) == ["chr1", "chr2", "chr10"]
    assert "chrX" in caplog.text and "chrZ" in caplog.text
    assert plan_chroms(bed, ["n.bam", "t.bam"], ["chr10", "chr1", "chr10"]) == ["chr1", "chr10"]
    with pytest.raises(ValueError, match="naming"):
        plan_chroms(bed, ["n.bam", "t.bam"], ["1", "2"])
    (tmp_path / "empty.bed").write_text("track name=empty\n")
    with pytest.raises(ValueError, match="none of the chromosomes"):
        plan_chroms(tmp_path / "empty.bed", ["n.bam"])


def test_read_bam_list_does_not_follow_the_list_symlink(tmp_path):
    src, work = tmp_path / "src", tmp_path / "work"
    src.mkdir()
    work.mkdir()
    (src / "list.txt").write_text("# normals\nA.bam\n/abs/B.bam\n")
    os.symlink(src / "list.txt", work / "list.txt")
    entries, resolved = read_bam_list(work / "list.txt")
    assert entries == ["A.bam", "/abs/B.bam"]
    assert resolved == [str(work / "A.bam"), "/abs/B.bam"]
    (src / "none.txt").write_text("# nothing\n")
    with pytest.raises(ValueError, match="no BAM"):
        read_bam_list(src / "none.txt")


# ------------------------------------------------------------------ with samtools
def test_pileup_same_rows_and_order_for_any_thread_count(bam_set, tmp_path):
    b = bam_set
    bed = tmp_path / "lexi.bed"  # chrB first, plus a contig that no BAM has
    bed.write_text("chrB\t100\t1400\nchrZ\t5\t50\nchrA\t100\t1400\n")
    one, four = tmp_path / "p1.gz", tmp_path / "p4.gz"
    run(*pair_args(b, bed, one))
    run(*pair_args(b, bed, four, "-t", 4))
    assert lines(one) == lines(four)
    chroms = list(dict.fromkeys(row.split("\t")[0] for row in lines(one)[1:]))
    assert chroms == ["chrA", "chrB"]  # BAM header order, chrZ skipped
    assert not list(tmp_path.glob(".ebscore_pileup_*")), "scratch removed"


def test_pileup_bed_gz_and_space_delimited(bam_set, tmp_path):
    b = bam_set
    with gzip.open(tmp_path / "t.bed.gz", "wt") as fh:
        fh.write("chrA\t100\t1400\nchrB\t100\t1400\n")
    (tmp_path / "s.bed").write_text("chrA 100 1400\nchrB 100 1400\n")
    ref = tmp_path / "ref.gz"
    run(*pair_args(b, b["bed"], ref))
    for bed in ("t.bed.gz", "s.bed"):
        out = tmp_path / f"{bed}.pileup.gz"
        run(*pair_args(b, tmp_path / bed, out, "-t", 2))
        assert lines(out) == lines(ref)


def test_pileup_respects_bed_mapq_and_chrom(bam_set, tmp_path):
    b = bam_set
    out = tmp_path / "p.gz"
    run(*pair_args(b, b["bed"], out, "-t", 2))
    df = pd.read_csv(out, sep="\t", dtype={"Read2": str})
    assert list(df.columns) == PAIRED_HEADER.split("\t")
    assert df["Start"].between(101, 1400).all()  # BED 100-1400 (0-based) only; reads cover 1-1500
    c, p = b["low_mapq"]
    alt = b["alt_of"](c, p)
    row = df[(df["Chr"] == c) & (df["Start"] == p + 1)].iloc[0]
    assert row["Read2"].upper().count(alt) == 0  # the 12 MAPQ-10 alt reads are filtered
    out0 = tmp_path / "p0.gz"
    run(*pair_args(b, b["bed"], out0, "--min-mapq", 0, "--chrom", "chrB"))
    df0 = pd.read_csv(out0, sep="\t", dtype={"Read2": str})
    assert set(df0["Chr"]) == {"chrB"}
    row0 = df0[df0["Start"] == p + 1].iloc[0]
    assert row0["Read2"].upper().count(alt) == 12
    dup = tmp_path / "dup.gz"
    run(*pair_args(b, b["bed"], dup, "--chrom", "chrB", "--chrom", "chrB"))
    assert lines(dup) == lines(tmp_path / "p.gz")[:1] + [r for r in lines(tmp_path / "p.gz")[1:] if r.startswith("chrB\t")]


def test_pileup_without_covered_positions_has_a_header(bam_set, tmp_path):
    b = bam_set
    bed = tmp_path / "a.bed"
    bed.write_text("chrA\t100\t1400\n")
    out = tmp_path / "empty.gz"
    run(*pair_args(b, bed, out, "--chrom", "chrB"))
    assert lines(out) == [PAIRED_HEADER]
    eb = tmp_path / "eb.tsv"
    pon = tmp_path / "pon"
    run("pon-matrix", "--bam-list", b["pon_list"], "--fasta", b["fasta"], "--bed", bed, "-o", pon, "--chrom", "chrB")
    assert lines(pon / "matrix" / "chrB.pon.gz") == [MATRIX_HEADER]
    run("score", "--variants", b["variants"], "--ponmatrix", pon, "--pileup", out, "--chrom", "chrB", "-o", eb)
    res = pd.read_csv(eb, sep="\t")
    assert len(res) == 2 and res["EB"].isna().all()


def test_pileup_fails_without_common_chromosome(bam_set, tmp_path):
    b = bam_set
    bed = tmp_path / "z.bed"
    bed.write_text("1\t100\t1400\n")
    for extra in ([], ["-t", "2"]):
        assert main([*map(str, pair_args(b, bed, tmp_path / "x.gz", *extra)), "--quiet"]) == 2
    assert main(["pon-matrix", "--bam-list", str(b["pon_list"]), "--fasta", str(b["fasta"]), "--bed", str(bed),
                 "-o", str(tmp_path / "pon"), "--quiet"]) == 2


def test_pon_matrix_column_order_follows_the_list(bam_set, tmp_path):
    b = bam_set
    pon = tmp_path / "pon"
    run("pon-matrix", "--bam-list", b["pon_list_mixed"], "--fasta", b["fasta"], "--bed", b["bed"], "-o", pon,
        "--chrom", "chrA", "--chrom", "chrA", "-t", 2)
    assert (pon / "PON_list.txt").read_text() == "PON3.bam\nPONX.bam\nPON1.bam\n"
    assert sorted(p.name for p in (pon / "matrix").iterdir()) == ["chrA.pon.gz"]
    m = pd.read_csv(pon / "matrix" / "chrA.pon.gz", sep="\t", dtype=str)
    c, p = b["private"]
    cell = m.loc[m["Start"] == str(p + 1), b["alt_of"](c, p)].iloc[0]
    plus, minus = ([int(x) for x in s.split("|")] for s in cell.split("="))
    alts = [x + y for x, y in zip(plus, minus)]
    assert alts[1] >= 5 and alts[0] == alts[2] == 0  # the private SNV sits in column 2 (PONX)

    # --exclude-sample PONX drops exactly that column
    var = tmp_path / "v.tsv"
    var.write_text(f"Chr\tStart\tEnd\tRef\tAlt\n{c}\t{p + 1}\t{p + 1}\t{b['contigs'][c][p]}\t{b['alt_of'](c, p)}\n")
    pair = tmp_path / "pair.gz"
    run(*pair_args(b, b["bed"], pair, "--chrom", "chrA"))
    full, excl = tmp_path / "full.tsv", tmp_path / "excl.tsv"
    run("score", "--variants", var, "--ponmatrix", pon, "--pileup", pair, "-o", full)
    run("score", "--variants", var, "--ponmatrix", pon, "--pileup", pair, "-o", excl, "--exclude-sample", "PONX")
    f, e = pd.read_csv(full, sep="\t").iloc[0], pd.read_csv(excl, sep="\t").iloc[0]
    fa, fd = f["PON+"].split("<")
    ea, ed = e["PON+"].split("<")
    keep = [0, 2]
    assert ea.split("|") == [fa.split("|")[i] for i in keep] and ed.split("|") == [fd.split("|")[i] for i in keep]


def test_pon_matrix_with_a_staged_symlinked_list(bam_set, tmp_path):
    # Nextflow-style staging: the list and the BAMs are symlinks in the task folder, the list's target lives elsewhere
    b = bam_set
    src, work = tmp_path / "src", tmp_path / "work"
    src.mkdir()
    work.mkdir()
    (src / "PON_list.txt").write_text("# PON normals\nPON1.bam\nPON2.bam\n")
    for name in ("PON1", "PON2"):
        os.symlink(b["dir"] / f"{name}.bam", work / f"{name}.bam")
        os.symlink(b["dir"] / f"{name}.bam.bai", work / f"{name}.bam.bai")
    os.symlink(src / "PON_list.txt", work / "PON_list.txt")
    cwd = os.getcwd()
    try:
        os.chdir(work)
        run("pon-matrix", "--bam-list", "PON_list.txt", "--fasta", b["fasta"], "--bed", b["bed"], "-o", ".",
            "--chrom", "chrB")
    finally:
        os.chdir(cwd)
    assert (src / "PON_list.txt").read_text() == "# PON normals\nPON1.bam\nPON2.bam\n", "input list untouched"
    assert not (work / "PON_list.txt").is_symlink()
    assert (work / "PON_list.txt").read_text() == "PON1.bam\nPON2.bam\n"
    assert (work / "matrix" / "chrB.pon.gz").is_file()
