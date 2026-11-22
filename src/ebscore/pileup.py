"""samtools mpileup front ends: the paired normal/tumor pileup and the PONmatrix.

Both use the legacy flag set: ``samtools mpileup -f <fasta> -q <min_mapq=20> -Q <min_baseq=25> -l <bed>`` (samtools
defaults otherwise: BAQ on, anomalous pairs skipped, overlap detection on, ``--ff UNMAP,SECONDARY,QCFAIL,DUP``,
max depth 8000 per file). The bundle FASTA with ``-r <chrom>`` replaces the legacy per-chromosome split FASTAs.

Chromosomes: the BED's chromosomes (or the requested ``chroms``, de-duplicated) that every BAM header has, in the
order of the (first) BAM header. Chromosomes missing from a BAM header are skipped with a warning; if none is left,
the run fails (e.g. a chr1 vs 1 naming mismatch). The output order is therefore the same with any thread count.

Splitting (plan D16): the paired pileup runs as one samtools process over the whole BED unless ``chroms`` is given
or ``threads > 1`` (then one samtools process per chromosome, gzipped parts concatenated in BAM header order). The
PONmatrix is always written per chromosome (``matrix/<chrom>.pon.gz``) and builds all BED chromosomes by default.
Every output has its header line, also when no position is covered.
"""

from __future__ import annotations

import gzip
import logging
import os
import shutil
import tempfile
from pathlib import Path

from .counts import MATRIX_COLUMNS
from .tools import bam_contigs, bed_chroms, ensure_dir, first_lines, q, run, run_parallel, tool, unique
from .tools import write_text_atomic

log = logging.getLogger("ebscore")

MIN_MAPQ = 20
MIN_BASEQ = 25
PAIRED_HEADER = "Chr\tStart\tRef\tDepth1\tRead1\tDepth2\tRead2"
MATRIX_HEADER = "\t".join(MATRIX_COLUMNS)


def mpileup_cmd(bams, fasta, bed, chrom=None, min_mapq=MIN_MAPQ, min_baseq=MIN_BASEQ, bam_list=None) -> str:
    cmd = f"samtools mpileup -f {q(fasta)} -q {int(min_mapq)} -Q {int(min_baseq)} -l {q(bed)}"
    if chrom:
        cmd += f" -r {q(chrom)}"
    if bam_list:
        cmd += f" -b {q(bam_list)}"
    for bam in bams or []:
        cmd += f" {q(bam)}"
    return cmd


def plan_chroms(bed, bams, chroms=None, threads: int = 1) -> list[str]:
    """Chromosomes to pile up, in the order of the first BAM header (see the module docstring)."""
    wanted = unique(chroms) if chroms else bed_chroms(bed)
    headers = run_parallel(bam_contigs, bams, threads)
    common = set(headers[0]).intersection(*headers[1:])
    for bam, contigs in zip(bams, headers):
        lost = [c for c in wanted if c not in set(contigs)]
        if lost:
            log.warning("%s lacks %d %s chromosome(s), skipped: %s", bam, len(lost),
                        "requested" if chroms else "BED", " ".join(lost[:10]) + (" ..." if len(lost) > 10 else ""))
    wanted_set = set(wanted)
    keep = [c for c in headers[0] if c in wanted_set and c in common]
    if not keep:
        source = "requested chromosomes" if chroms else f"chromosomes of {bed}"
        raise ValueError(
            f"none of the {source} ({' '.join(wanted[:5]) or 'none'}) is in the BAM header(s) "
            f"(first BAM has {' '.join(headers[0][:5])} ...); check the chromosome naming"
        )
    return keep


def _check_rows(out: Path) -> None:
    if len(first_lines(out, 2)) < 2:
        log.warning("%s has no covered position (header only)", out)


def _write_parts(make_cmd, chroms, out: Path, header: str, threads: int, tmpdir) -> None:
    """One pipeline per chromosome (its header line dropped, gzipped), then header + parts concatenated.

    The result is a multi-member gzip (valid gzip; gzip -dc, zcat, Python and pandas read it). ``chroms == [None]``
    runs one pipeline over the whole BED. Scratch goes to ``tmpdir`` (default: next to the output).
    """
    tmp = Path(tempfile.mkdtemp(prefix=".ebscore_pileup_", dir=tmpdir or out.parent))
    try:
        head = tmp / "header.gz"
        head.write_bytes(gzip.compress(f"{header}\n".encode(), mtime=0))
        parts = [tmp / f"{i:04d}.tsv.gz" for i in range(len(chroms))]
        run_parallel(
            lambda i: run(f"{make_cmd(chroms[i])} | tail -n +2 | gzip -c > {q(parts[i])}"), range(len(chroms)), threads
        )
        part = out.with_name(f".{out.name}.part")
        run(f"cat {q(head)} {' '.join(q(p) for p in parts)} > {q(part)}")
        os.replace(part, out)
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

    Normal = sample 1, tumor = sample 2 (legacy order, varscan.snk:24). EB and PON-CNV both read this file. Rows are
    in BAM header order; the header line is always written.
    """
    out = Path(out)
    ensure_dir(out.parent)
    chroms = unique(chroms)
    planned = plan_chroms(bed, [normal, tumor], chroms, threads)

    def make_cmd(chrom):
        return f"{mpileup_cmd([normal, tumor], fasta, bed, chrom, min_mapq, min_baseq)} | {tool('cleanpileup', '-d')}"

    if threads > 1 or chroms:
        log.info("paired pileup over %d chromosomes with %d threads", len(planned), threads)
        _write_parts(make_cmd, planned, out, PAIRED_HEADER, threads, tmpdir)
    else:
        log.info("paired pileup over the whole BED (%d chromosomes)", len(planned))
        _write_parts(make_cmd, [None], out, PAIRED_HEADER, 1, tmpdir)
    _check_rows(out)
    return out


def read_bam_list(bam_list) -> tuple[list[str], list[str]]:
    """PON list entries as given, and as absolute paths.

    Relative entries resolve against the folder the list is seen in (``os.path.abspath``, symlinks not followed), so
    a list staged as a symlink next to staged BAMs (Nextflow) finds them.
    """
    base = Path(os.path.abspath(bam_list)).parent
    entries = []
    with open(bam_list) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                entries.append(line)
    if not entries:
        raise ValueError(f"{bam_list}: no BAM/CRAM entries")
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
    resolved copy in a temp folder, so no build-time paths end up in the PON set. Chromosomes missing from any PON
    BAM header are skipped with a warning.
    """
    outdir = ensure_dir(outdir)
    matrix_dir = ensure_dir(outdir / "matrix")
    entries, resolved = read_bam_list(bam_list)
    planned = plan_chroms(bed, resolved, unique(chroms), threads)
    write_text_atomic(outdir / "PON_list.txt", "".join(f"{e}\n" for e in entries))
    log.info("PONmatrix: %d normals, %d chromosomes, %d threads", len(entries), len(planned), threads)

    with tempfile.TemporaryDirectory(prefix="ebscore_pon_") as tmpdir:
        list_file = Path(tmpdir) / "bams.txt"
        list_file.write_text("".join(f"{b}\n" for b in resolved))

        def build(chrom):
            out = matrix_dir / f"{chrom}.pon.gz"
            tmp = matrix_dir / f".{chrom}.pon.gz.part"
            pipeline = (
                f"{mpileup_cmd(None, fasta, bed, chrom, min_mapq, min_baseq, bam_list=list_file)}"
                f" | {tool('cleanpileup')} | {tool('pile2count')} | tail -n +2"
            )
            run(f"{{ printf '%s\\n' {q(MATRIX_HEADER)}; {pipeline}; }} | gzip -c > {q(tmp)}")
            os.replace(tmp, out)
            _check_rows(out)
            log.info("PONmatrix %s done", out.name)
            return out

        return run_parallel(build, planned, threads)
