"""samtools mpileup front ends: the paired normal/tumor pileup and the PONmatrix.

Both use the legacy flag set: ``samtools mpileup -f <fasta> -q <min_mapq=20> -Q <min_baseq=25> -l <bed>`` (samtools
defaults otherwise: BAQ on, anomalous pairs skipped, overlap detection on, ``--ff UNMAP,SECONDARY,QCFAIL,DUP``,
max depth 8000 per file). The bundle FASTA with ``-r <chrom>`` replaces the legacy per-chromosome split FASTAs.

Splitting (plan D16): the paired pileup runs as one task over the whole BED unless ``chroms`` is given or
``threads > 1`` (then one samtools process per chromosome, concatenated in BED order). The PONmatrix is always
written per chromosome (``matrix/<chrom>.pon.gz``) and builds all BED chromosomes by default.
"""

from __future__ import annotations

import logging
import os
import shutil
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from .tools import bed_chroms, ensure_dir, q, run, tool

log = logging.getLogger("ebscore")

MIN_MAPQ = 20
MIN_BASEQ = 25


def mpileup_cmd(bams, fasta, bed, chrom=None, min_mapq=MIN_MAPQ, min_baseq=MIN_BASEQ, bam_list=None) -> str:
    cmd = f"samtools mpileup -f {q(fasta)} -q {int(min_mapq)} -Q {int(min_baseq)} -l {q(bed)}"
    if chrom:
        cmd += f" -r {q(chrom)}"
    if bam_list:
        cmd += f" -b {q(bam_list)}"
    for bam in bams or []:
        cmd += f" {q(bam)}"
    return cmd


def _run_chunks(make_cmd, chroms, out_gz: Path, threads: int, tmpdir) -> None:
    """Run one pipeline per chromosome (each prints its own header), then concatenate with one header and gzip."""
    tmp = Path(tempfile.mkdtemp(prefix="ebscore_pileup_", dir=tmpdir))
    try:
        parts = [tmp / f"{i:04d}.tsv" for i in range(len(chroms))]
        with ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
            list(pool.map(lambda ip: run(f"{make_cmd(ip[1])} > {q(ip[0])}"), zip(parts, chroms)))
        header_cmd = " ".join(q(p) for p in parts)
        # header from the first non-empty part, then all data lines
        run(
            f"mawk 'FNR == 1 {{ if (!h++) print; next }} {{ print }}' {header_cmd} | gzip -c > {q(out_gz)}"
        )
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def paired_pileup(
    normal,
    tumor,
    fasta,
    bed,
    out,
    chroms=None,
    min_mapq=MIN_MAPQ,
    min_baseq=MIN_BASEQ,
    threads=1,
    tmpdir=None,
) -> Path:
    """``samtools mpileup NORMAL TUMOR | cleanpileup -d | gzip`` -> ``Chr Start Ref Depth1 Read1 Depth2 Read2``.

    Normal = sample 1, tumor = sample 2 (legacy order, varscan.snk:24). EB and PON-CNV both read this file.
    """
    out = Path(out)
    ensure_dir(out.parent)

    def make_cmd(chrom):
        return f"{mpileup_cmd([normal, tumor], fasta, bed, chrom, min_mapq, min_baseq)} | {tool('cleanpileup', '-d')}"

    if threads > 1 or (chroms and len(chroms) > 1):
        chroms = chroms or bed_chroms(bed)
        log.info("paired pileup over %d chromosomes with %d threads", len(chroms), threads)
        _run_chunks(make_cmd, chroms, out, threads, tmpdir)
    else:
        chrom = chroms[0] if chroms else None
        log.info("paired pileup over %s", chrom or "the whole BED")
        run(f"{make_cmd(chrom)} | gzip -c > {q(out)}")
    return out


def read_bam_list(bam_list) -> tuple[list[str], list[str]]:
    """PON list entries as given, and as absolute paths (relative entries resolve against the list's folder)."""
    base = Path(bam_list).resolve().parent
    entries = []
    with open(bam_list) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                entries.append(line)
    resolved = [str(Path(e) if Path(e).is_absolute() else base / e) for e in entries]
    return entries, resolved


def pon_matrix(
    bam_list,
    fasta,
    bed,
    outdir,
    chroms=None,
    min_mapq=MIN_MAPQ,
    min_baseq=MIN_BASEQ,
    threads=1,
) -> list[Path]:
    """Build ``<outdir>/matrix/<chrom>.pon.gz`` (PONmatrix format) for every chromosome of the BED or ``chroms``.

    Also writes ``<outdir>/PON_list.txt`` (the list entries as given, in matrix column order) so that the list and
    the matrix always agree (the legacy HAEv7 PON has 28 matrix columns but 26 list entries). samtools gets a
    resolved copy in a temp folder, so no build-time paths end up in the PON set.
    """
    outdir = ensure_dir(outdir)
    matrix_dir = ensure_dir(outdir / "matrix")
    entries, resolved = read_bam_list(bam_list)
    (outdir / "PON_list.txt").write_text("".join(f"{e}\n" for e in entries))
    chroms = list(chroms) if chroms else bed_chroms(bed)
    log.info("PONmatrix: %d normals, %d chromosomes, %d threads", len(entries), len(chroms), threads)

    with tempfile.TemporaryDirectory(prefix="ebscore_pon_") as tmpdir:
        list_file = Path(tmpdir) / "bams.txt"
        list_file.write_text("".join(f"{b}\n" for b in resolved))

        def build(chrom):
            out = matrix_dir / f"{chrom}.pon.gz"
            tmp = out.with_suffix(".gz.part")
            run(
                f"{mpileup_cmd(None, fasta, bed, chrom, min_mapq, min_baseq, bam_list=list_file)}"
                f" | {tool('cleanpileup')} | {tool('pile2count')} | gzip -c > {q(tmp)}"
            )
            os.replace(tmp, out)
            log.info("PONmatrix %s done", out.name)
            return out

        with ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
            return list(pool.map(build, chroms))
