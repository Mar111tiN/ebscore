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


# ---------------------------------------------------------------- tiny BAMs for the samtools integration tests
READ_LEN = 60


def _write_sam(path, contigs, reads, sample):
    with open(path, "w") as fh:
        fh.write("@HD\tVN:1.6\tSO:unsorted\n")
        for name, seq in contigs.items():
            fh.write(f"@SQ\tSN:{name}\tLN:{len(seq)}\n")
        fh.write(f"@RG\tID:{sample}\tSM:{sample}\n")
        for i, (contig, start, seq, rev) in enumerate(reads):
            flag = 16 if rev else 0
            qual = "I" * len(seq)
            fh.write(f"r{i}\t{flag}\t{contig}\t{start + 1}\t60\t{len(seq)}M\t*\t0\t0\t{seq}\t{qual}\tRG:Z:{sample}\n")


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

    def make(name, depth, snvs):
        sam = d / f"{name}.sam"
        _write_sam(sam, contigs, _simulate(rng, contigs, depth, snvs), name)
        sh(f"samtools sort -o {d / name}.bam {sam} && samtools index {d / name}.bam")
        sam.unlink()
        return d / f"{name}.bam"

    pon = [make(f"PON{i}", 30, pon_snvs) for i in range(1, 5)]
    normal = make("normal", 30, {})
    tumor = make("tumor", 60, tum_snvs)
    (d / "pon_list.txt").write_text("".join(f"{p.name}\n" for p in pon))  # relative to the list folder
    with open(d / "variants.tsv", "w") as fh:
        fh.write("Chr\tStart\tEnd\tRef\tAlt\n")
        for c, p in real + artifact:
            fh.write(f"{c}\t{p + 1}\t{p + 1}\t{contigs[c][p]}\t{alt_of(c, p)}\n")
    return {"dir": d, "fasta": d / "ref.fa", "bed": d / "target.bed", "pon_list": d / "pon_list.txt",
            "normal": normal, "tumor": tumor, "variants": d / "variants.tsv", "real": real, "artifact": artifact}
