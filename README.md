# ebscore

EB score: how likely a somatic variant is a sequencing artifact, judged against a pileup Panel of Normals (the
**PONmatrix**). Per strand, a beta-binomial is fitted to the PON alt/depth counts at the variant position. The
tumor's upper-tail p-values of both strands are combined with Fisher's method, and `EB = -log10(p)` (0–60). High EB
means the PON does not explain the tumor's alt reads; low EB means artifact-like. The math is the Genomon
EBFilter / [EBCall](https://github.com/friend1ws/EBCall) approach (Shiraishi et al. 2013).

Branch `kumodev` is the port for `kumo-variant-calling-nf`: an installable package with a CLI, current Python/pandas,
the legacy defects fixed (see `CHANGELOG.md`) and a Dockerfile. `main` keeps the original Snakemake-era code.

## Install

```bash
conda env create -f environment.yml      # or any env with python>=3.10, numpy, pandas, scipy, samtools, mawk
conda activate ebscore
pip install --no-deps --no-build-isolation .
```

## Commands

All commands take explicit paths, log to stderr and write files. `-t/--threads` parallelises, `-v` shows every shell
command.

```bash
# 1. PONmatrix from normal BAMs/CRAMs (one-off per kit + normal set). Writes <dir>/matrix/<chrom>.pon.gz + PON_list.txt.
#    Builds every BED chromosome by default; --chrom chr7 builds one (scatter per chromosome).
ebscore pon-matrix --bam-list PON_list.txt --fasta genome.fa --bed kit.covered.bed -o PONmatrix/ -t 8

# 2. Optional AB cache (12 fits per PON position; see "AB cache" below). Writes <dir>/AB/<chrom>.AB.gz.
ebscore ab-cache --ponmatrix PONmatrix/ --chrom chr7 -t 8

# 3. Paired pileup per tumor/normal pair (normal first). Shared by EB and PON-CNV (myCNV).
ebscore pileup --normal N.bam --tumor T.bam --fasta genome.fa --bed kit.covered.bed -o pair.pileup.gz -t 4

# 4. EB per variant (table in ANNOVAR coordinates, or a VCF).
ebscore score --variants calls.vcf.gz --ponmatrix PONmatrix/ --pileup pair.pileup.gz -o calls.EB.tsv \
    --vcf-annotation calls.EB.annot.tsv --exclude-sample <matched normal, if it is in the PON> -t 4
#    or straight from the tumor BAM: --tumor-bam T.bam --fasta genome.fa

# VCF -> ANNOVAR-coordinate table (with the original VCF coordinates as vcf_* columns)
ebscore vcf2table --vcf calls.vcf.gz -o calls.tsv
```

`score` output: the input columns, then `PON+ PON- EB` (`PON±` = `alt1|…|altN<depth1|…|depthN` per strand;
`--keep-matrix` adds `Tumor` and `AB`). `--vcf-annotation` writes `#CHROM POS REF ALT EB` for
`bgzip` + `tabix -s1 -b2 -e2` + `bcftools annotate -c CHROM,POS,REF,ALT,INFO/EB`.

### Options that change values

| Option | Default | Legacy behaviour |
|---|---|---|
| `--fisher sf\|legacy` | `sf`: chi2 survival function | `legacy`: `1 - chi2.cdf`; no EB between ~16 and 60 |
| `--fit exact\|legacy` | `exact`: analytic gradient, tight tolerances, same optimum on any scipy | `legacy`: numerical gradient, early stop; old values only on scipy < 1.15 |
| `--missing na\|zero` | `na`: no tumor or PON data → EB `NA` | `zero`: EB 0 (looks like an artifact) |
| `--penalty` | 0.5 | was hard-coded 0.5 (config value ignored) |
| `--min-mapq/--min-baseq` | 20 / 25 | same |

### Pileup flags

`samtools mpileup -f <fasta> -q 20 -Q 25 -l <bed> [-r <chrom>]`, samtools defaults otherwise (BAQ on, anomalous pairs
skipped, overlap detection on, `--ff UNMAP,SECONDARY,QCFAIL,DUP`, max depth 8000). Build the PONmatrix and the
sample pileups with the same samtools version, flags and BED flavour.

### PONmatrix folder

```
PONmatrix/
  PON_list.txt          BAM paths in matrix column order
  matrix/<chrom>.pon.gz Chr Start Ref A G C T I D Depth; cells c1|…|cN=c1|…|cN (+ strand = - strand)
  AB/<chrom>.AB.gz      optional: Chr Start Ref A G C T I D; cells a+|b+=a-|b-
```

The same `matrix/` feeds the PON-CNV of myCNV. Scoring reads the folder only (temp files go to `--tmpdir`).

### AB cache

Fitting on the fly at the called positions is the default and needs no cache. The AB cache precomputes every
position × letter × strand (12 fits per position; a WES PON is ~50 M positions), which only pays off for very large
call sets. The cache is valid for the full PON only: with `--exclude-sample` the fits always run on the fly.

## Tests

```bash
pytest -q          # mawk tools, core math, legacy parity, VCF conversion, end-to-end with simulated BAMs
```

`tests/data/legacy/` holds outputs of the unmodified `main` code (py 3.9, pandas 1.5.2, scipy 1.9.3), made by
`tests/legacy/make_legacy_outputs.py`; `tests/data/make_synthetic.py` makes the synthetic pileups.

## Docker

```bash
docker build --target test .          # pytest inside the build env (not shipped)
docker build -t ebscore:<version> .   # runtime: debian slim + /opt/conda, non-root, no entrypoint
```
