import gzip
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

DATA = Path(__file__).parent / "data"
SYN = DATA / "synthetic"
LEGACY = DATA / "legacy"


def sh(cmd: str, **kw) -> subprocess.CompletedProcess:
    return subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True, capture_output=True, text=True, **kw)


def read_text(path) -> str:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        return fh.read()


@pytest.fixture(scope="session")
def legacy_ponmatrix(tmp_path_factory):
    """A writable copy of the legacy PONmatrix (matrix/ + AB/ + PON_list.txt)."""
    dst = tmp_path_factory.mktemp("ponmatrix") / "pon"
    shutil.copytree(LEGACY / "ponmatrix", dst)
    return dst


@pytest.fixture(scope="session")
def legacy_matrix_only(tmp_path_factory):
    """The legacy PONmatrix without its AB cache."""
    dst = tmp_path_factory.mktemp("ponmatrix_noab") / "pon"
    shutil.copytree(LEGACY / "ponmatrix", dst, ignore=shutil.ignore_patterns("AB"))
    return dst


# ---------------------------------------------------------------- a large synthetic PONmatrix (WES-like PON size)
LARGE_N = 32


@pytest.fixture(scope="session")
def large_pon(tmp_path_factory):
    """PONmatrix with 32 normals at ~200x (Poisson 100 per strand) on chrL, a paired pileup and a variant table.

    The 6-normal legacy fixture cannot reach EB > 16 (its fits keep the tail heavy); here a 20+20 of 200+200 tumor SNV
    over an all-zero PON scores EB ~26 with --fisher sf and 60 with --fisher legacy. Positions (Ref A):
    1001 strong SNV (C 20+20), PON all zero    1002 strong SNV (C 20+20), PON sparse alts (sample 2: 2+1 C)
    1003 SNV (C 3+3), covered by PON sample 5 only    1004 weak SNV (G 3+3), PON all zero
    """
    rng = np.random.default_rng(20261122)
    d = tmp_path_factory.mktemp("large_pon")
    (d / "pon" / "matrix").mkdir(parents=True)
    (d / "pon" / "PON_list.txt").write_text("".join(f"L{i:02d}.bam\n" for i in range(1, LARGE_N + 1)))
    rows = []
    for pos in (1001, 1002, 1003, 1004):
        dp, dm = rng.poisson(100, LARGE_N), rng.poisson(100, LARGE_N)
        if pos == 1003:
            keep = np.arange(LARGE_N) == 4
            dp, dm = dp * keep, dm * keep
        zero = "=".join(["|".join(["0"] * LARGE_N)] * 2)
        cells = {b: zero for b in "AGCTID"}
        cells["A"] = "|".join(map(str, dp)) + "=" + "|".join(map(str, dm))
        if pos == 1002:
            cp, cm = [0] * LARGE_N, [0] * LARGE_N
            cp[1], cm[1] = 2, 1
            cells["C"] = "|".join(map(str, cp)) + "=" + "|".join(map(str, cm))
        depth = "|".join(map(str, dp)) + "=" + "|".join(map(str, dm))
        rows.append("\t".join(["chrL", str(pos), "A", *(cells[b] for b in "AGCTID"), depth]))
    with gzip.open(d / "pon" / "matrix" / "chrL.pon.gz", "wt") as fh:
        fh.write("\t".join(["Chr", "Start", "Ref", *"AGCTID", "Depth"]) + "\n" + "\n".join(rows) + "\n")
    tumor = {1001: ("C", 20), 1002: ("C", 20), 1003: ("C", 3), 1004: ("G", 3)}
    with gzip.open(d / "pair.pileup.gz", "wt") as fh:
        fh.write("Chr\tStart\tRef\tDepth1\tRead1\tDepth2\tRead2\n")
        for pos, (alt, k) in tumor.items():
            read2 = alt * k + alt.lower() * k + "." * (200 - k) + "," * (200 - k)
            fh.write(f"chrL\t{pos}\tA\t200\t{'.' * 100 + ',' * 100}\t400\t{read2}\n")
    with open(d / "variants.tsv", "w") as fh:
        fh.write("Chr\tStart\tEnd\tRef\tAlt\n")
        for pos, (alt, _) in tumor.items():
            fh.write(f"chrL\t{pos}\t{pos}\tA\t{alt}\n")
    return {"ponmatrix": d / "pon", "pileup": d / "pair.pileup.gz", "variants": d / "variants.tsv"}


# ---------------------------------------------------------------- tiny BAMs for the samtools integration tests
READ_LEN = 60


def _write_sam(path, contigs, reads, sample):
    with open(path, "w") as fh:
        fh.write("@HD\tVN:1.6\tSO:unsorted\n")
        for name, seq in contigs.items():
            fh.write(f"@SQ\tSN:{name}\tLN:{len(seq)}\n")
        fh.write(f"@RG\tID:{sample}\tSM:{sample}\n")
        for i, (contig, start, seq, rev, *mapq) in enumerate(reads):
            flag = 16 if rev else 0
            qual = "I" * len(seq)
            mq = mapq[0] if mapq else 60
            fh.write(f"r{i}\t{flag}\t{contig}\t{start + 1}\t{mq}\t{len(seq)}M\t*\t0\t0\t{seq}\t{qual}\tRG:Z:{sample}\n")


def _simulate(rng, contigs, depth, snvs):
    """Single-end reads tiling the contigs; snvs: {(contig, 0-based pos): (alt, vaf)}."""
    reads = []
    for name, seq in contigs.items():
        n_reads = depth * len(seq) // READ_LEN
        for _ in range(n_reads):
            start = int(rng.integers(0, len(seq) - READ_LEN))
            read = list(seq[start : start + READ_LEN])
            for (c, pos), (alt, vaf) in snvs.items():
                if c == name and start <= pos < start + READ_LEN and rng.random() < vaf:
                    read[pos - start] = alt
            reads.append((name, start, "".join(read), bool(rng.random() < 0.5)))
    return reads


@pytest.fixture(scope="session")
def bam_set(tmp_path_factory):
    """Reference + 4 PON normals + normal + tumor BAMs, a BED and a variant table (needs samtools)."""
    if shutil.which("samtools") is None:
        pytest.skip("samtools not available")
    d = tmp_path_factory.mktemp("bams")
    rng = np.random.default_rng(7)
    contigs = {c: "".join(rng.choice(list("ACGT"), size=1500)) for c in ("chrA", "chrB")}
    with open(d / "ref.fa", "w") as fh:
        for name, seq in contigs.items():
            fh.write(f">{name}\n")
            fh.write("\n".join(seq[i : i + 60] for i in range(0, len(seq), 60)) + "\n")
    sh(f"samtools faidx {d / 'ref.fa'}")
    with open(d / "target.bed", "w") as fh:
        fh.write("chrA\t100\t1400\nchrB\t100\t1400\n")

    def alt_of(c, p):
        return next(b for b in "ACGT" if b != contigs[c][p])

    real = [("chrA", 300), ("chrA", 900), ("chrB", 500)]
    artifact = [("chrB", 1100)]
    pon_snvs = {k: (alt_of(*k), 0.05) for k in artifact}
    tum_snvs = {k: (alt_of(*k), 0.3) for k in real} | {k: (alt_of(*k), 0.05) for k in artifact}

    def make(name, depth, snvs, extra=()):
        sam = d / f"{name}.sam"
        _write_sam(sam, contigs, _simulate(rng, contigs, depth, snvs) + list(extra), name)
        sh(f"samtools sort -o {d / name}.bam {sam} && samtools index {d / name}.bam")
        sam.unlink()
        return d / f"{name}.bam"

    # 12 tumor reads with MAPQ 10 carrying an alt base: dropped by the default --min-mapq 20
    low_mapq = ("chrB", 1300)
    lm_reads = []
    for s in range(low_mapq[1] - 40, low_mapq[1] - 28):
        read = list(contigs["chrB"][s : s + READ_LEN])
        read[low_mapq[1] - s] = alt_of(*low_mapq)
        lm_reads.append(("chrB", s, "".join(read), s % 2 == 0, 10))
    pon = [make(f"PON{i}", 30, pon_snvs) for i in range(1, 5)]
    normal = make("normal", 30, {})
    tumor = make("tumor", 60, tum_snvs, lm_reads)
    (d / "pon_list.txt").write_text("".join(f"{p.name}\n" for p in pon))  # relative to the list folder
    # a PON normal with a private 50 % SNV, and an unsorted list that puts it in matrix column 2
    private = ("chrA", 700)
    ponx = make("PONX", 30, {private: (alt_of(*private), 0.5)})
    (d / "pon_list_mixed.txt").write_text("PON3.bam\nPONX.bam\nPON1.bam\n")
    with open(d / "variants.tsv", "w") as fh:
        fh.write("Chr\tStart\tEnd\tRef\tAlt\n")
        for c, p in real + artifact:
            fh.write(f"{c}\t{p + 1}\t{p + 1}\t{contigs[c][p]}\t{alt_of(c, p)}\n")
    return {"dir": d, "fasta": d / "ref.fa", "bed": d / "target.bed", "pon_list": d / "pon_list.txt",
            "normal": normal, "tumor": tumor, "variants": d / "variants.tsv", "real": real, "artifact": artifact,
            "low_mapq": low_mapq, "ponx": ponx, "private": private, "pon_list_mixed": d / "pon_list_mixed.txt",
            "alt_of": alt_of, "contigs": contigs}
