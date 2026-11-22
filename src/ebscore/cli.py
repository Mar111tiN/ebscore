"""ebscore command line: pileup, pon-matrix, ab-cache, score, vcf2table."""

from __future__ import annotations

import argparse
import logging
import sys

from . import __version__
from .pileup import MIN_BASEQ, MIN_MAPQ

log = logging.getLogger("ebscore")


def _common(p: argparse.ArgumentParser, threads=True) -> None:
    if threads:
        p.add_argument("-t", "--threads", type=int, default=1, help="parallel processes/threads (default 1)")
    p.add_argument("-v", "--verbose", action="store_true", help="debug logging (shows every shell command)")
    p.add_argument("--quiet", action="store_true", help="warnings and errors only")


def _fit_opt(p: argparse.ArgumentParser) -> None:
    p.add_argument(
        "--fit", choices=["exact", "legacy"], default="exact",
        help="beta-binomial fit: exact (analytic gradient, restarted until stationary: the penalised optimum; "
        "default) or legacy (numerical gradient, stops early like the old code and Genomon EBFilter; old values only "
        "with scipy < 1.15). On realistic PONs the two differ by several EB units (see CHANGELOG.md)",
    )


def _pileup_opts(p: argparse.ArgumentParser) -> None:
    p.add_argument("--min-mapq", type=int, default=MIN_MAPQ, help=f"samtools mpileup -q (default {MIN_MAPQ})")
    p.add_argument("--min-baseq", type=int, default=MIN_BASEQ, help=f"samtools mpileup -Q (default {MIN_BASEQ})")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ebscore",
        description="EB artifact score of somatic variants against a pileup Panel of Normals (PONmatrix).",
    )
    parser.add_argument("--version", action="version", version=f"ebscore {__version__}")
    sub = parser.add_subparsers(dest="command", required=True, metavar="COMMAND")

    p = sub.add_parser("pileup", help="paired normal/tumor pileup (cleanpileup -d), shared by EB and PON-CNV")
    p.add_argument("--normal", required=True, help="normal BAM/CRAM (sample 1)")
    p.add_argument("--tumor", required=True, help="tumor BAM/CRAM (sample 2)")
    p.add_argument("--fasta", required=True, help="reference FASTA (+ .fai)")
    p.add_argument("--bed", required=True, help="target BED (the kit's covered flavour, plan L9)")
    p.add_argument("-o", "--output", required=True, help="output .pileup.gz")
    p.add_argument("--chrom", action="append", help="only these chromosomes (repeatable; default: whole BED)")
    p.add_argument("--tmpdir", help="scratch folder (default: system temp)")
    _pileup_opts(p)
    _common(p)

    p = sub.add_parser("pon-matrix", help="build the PONmatrix (matrix/<chrom>.pon.gz) from normal BAMs")
    p.add_argument("--bam-list", required=True, help="one BAM/CRAM per line; relative paths resolve against the list")
    p.add_argument("--fasta", required=True)
    p.add_argument("--bed", required=True, help="target BED (same flavour as the samples, plan L9)")
    p.add_argument("-o", "--outdir", required=True, help="PONmatrix folder (gets matrix/ and PON_list.txt)")
    p.add_argument("--chrom", action="append", help="only these chromosomes (repeatable; default: all BED chromosomes)")
    _pileup_opts(p)
    _common(p)

    p = sub.add_parser("ab-cache", help="optional: fit the AB cache (AB/<chrom>.AB.gz) for a PONmatrix")
    p.add_argument("--ponmatrix", required=True, help="PONmatrix folder")
    p.add_argument("--chrom", action="append", required=True, help="chromosome(s) to fit (repeatable)")
    p.add_argument("-o", "--outdir", help="output folder for AB/ (default: the PONmatrix folder)")
    p.add_argument("--penalty", type=float, default=0.5, help="fit penalty on ln(a+b) (default 0.5)")
    p.add_argument("--chunk-rows", type=int, default=50000, help="PONmatrix rows per chunk (default 50000)")
    _fit_opt(p)
    _common(p)

    p = sub.add_parser("score", help="EB score of variants against the PONmatrix")
    p.add_argument("--variants", required=True, help="variant table (Chr Start End Ref Alt, ANNOVAR coords) or VCF")
    p.add_argument("--ponmatrix", required=True, help="PONmatrix folder (matrix/, optional AB/, PON_list.txt)")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--pileup", help="cleanpileup file, e.g. from `ebscore pileup` (tumor = Read2)")
    src.add_argument("--tumor-bam", help="tumor BAM/CRAM; pileup at the variant positions (needs --fasta)")
    p.add_argument("--tumor-column", help="pileup read column of the tumor (default Read2, else the last Read column)")
    p.add_argument("--fasta", help="reference FASTA for --tumor-bam")
    p.add_argument("-o", "--output", required=True, help="output TSV: input columns + PON+ PON- EB")
    p.add_argument("--vcf-annotation", help="also write #CHROM POS REF ALT EB (VCF input only) for bcftools annotate")
    p.add_argument("--exclude-sample", action="append", default=[], help="drop this PON sample (e.g. the matched normal; repeatable)")
    p.add_argument("--exclude-index", action="append", type=int, default=[], help="drop PON column N (1-based; repeatable)")
    p.add_argument("--no-ab-cache", action="store_true", help="fit on the fly even if AB/ exists")
    p.add_argument("--penalty", type=float, default=0.5, help="fit penalty on ln(a+b) (default 0.5)")
    p.add_argument("--fisher", choices=["sf", "legacy"], default="sf", help="sf (default) or legacy 1-cdf (EB jumps 16 -> 60)")
    _fit_opt(p)
    p.add_argument("--missing", choices=["na", "zero"], default="na", help="EB without tumor/PON data: NA (default) or 0 (legacy)")
    p.add_argument("--chrom", action="append", help="only variants on these chromosomes (repeatable)")
    p.add_argument("--keep-matrix", action="store_true", help="also output the Tumor and AB columns")
    p.add_argument("--tmpdir", help="scratch folder (default: system temp)")
    _pileup_opts(p)
    _common(p)

    p = sub.add_parser("vcf2table", help="VCF -> variant table in ANNOVAR coordinates (with vcf_* columns)")
    p.add_argument("--vcf", required=True)
    p.add_argument("-o", "--output", required=True)
    _common(p, threads=False)
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    from .tools import setup_logging

    setup_logging(args.verbose, args.quiet)
    try:
        if args.command == "pileup":
            from .pileup import paired_pileup

            paired_pileup(
                args.normal, args.tumor, args.fasta, args.bed, args.output, args.chrom,
                args.min_mapq, args.min_baseq, args.threads, args.tmpdir,
            )
        elif args.command == "pon-matrix":
            from .pileup import pon_matrix

            pon_matrix(args.bam_list, args.fasta, args.bed, args.outdir, args.chrom, args.min_mapq, args.min_baseq, args.threads)
        elif args.command == "ab-cache":
            from pathlib import Path

            from .abcache import ab_cache

            for chrom in args.chrom:
                out = Path(args.outdir) / "AB" / f"{chrom}.AB.gz" if args.outdir else None
                ab_cache(args.ponmatrix, chrom, out, args.penalty, args.threads, args.chunk_rows, args.fit)
        elif args.command == "score":
            from .score import score_variants, write_vcf_annotation

            df = score_variants(
                args.variants, args.ponmatrix, pileup=args.pileup, tumor_column=args.tumor_column,
                tumor_bam=args.tumor_bam, fasta=args.fasta, use_ab_cache=not args.no_ab_cache,
                exclude_samples=args.exclude_sample, exclude_index=args.exclude_index, penalty=args.penalty,
                fisher=args.fisher, missing=args.missing, threads=args.threads, chroms=args.chrom,
                min_mapq=args.min_mapq, min_baseq=args.min_baseq, keep_matrix=args.keep_matrix, tmpdir=args.tmpdir,
                fit_method=args.fit,
            )
            df.to_csv(args.output, sep="\t", index=False, na_rep="NA")
            log.info("wrote %d variants to %s", len(df), args.output)
            if args.vcf_annotation:
                n = write_vcf_annotation(df, args.vcf_annotation)
                log.info("wrote %d EB annotations to %s", n, args.vcf_annotation)
        elif args.command == "vcf2table":
            from .variants import read_vcf

            read_vcf(args.vcf).to_csv(args.output, sep="\t", index=False, na_rep="")
    except (ValueError, FileNotFoundError) as e:
        log.error("%s", e)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
