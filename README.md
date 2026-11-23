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
conda env create -f environment.yml      # runtime deps + pytest; or any env with python>=3.10, numpy, pandas, scipy,
conda activate ebscore                   # samtools, mawk (and pytest for the tests)
pip install --no-deps --no-build-isolation .
```

## Commands

All commands take explicit paths, log to stderr and write files. `-t/--threads` parallelises, `-v/--verbose` shows
every shell command, `--quiet` keeps warnings and errors only. Repeated `--chrom` values count once.

```bash
# The BED is always the kit's padded BED (plan D20): the same file for the PONmatrix, the pileups and myCNV.

# 1. PONmatrix from normal BAMs/CRAMs (one-off per kit + normal set). Writes <dir>/matrix/<chrom>.pon.gz + PON_list.txt.
#    Builds every BED chromosome by default; --chrom chr7 builds one (scatter per chromosome).
ebscore pon-matrix --bam-list PON_list.txt --fasta genome.fa --bed kit.padded.bed -o PONmatrix/ -t 8

# 2. Paired pileup per tumor/normal pair (normal first). Shared by EB and PON-CNV (myCNV).
ebscore pileup --normal N.bam --tumor T.bam --fasta genome.fa --bed kit.padded.bed -o pair.pileup.gz -t 4

# 3. EB per variant (table in ANNOVAR coordinates, or a VCF; left-align indels first: bcftools norm -f genome.fa).
ebscore score --variants calls.vcf.gz --ponmatrix PONmatrix/ --pileup pair.pileup.gz -o calls.EB.tsv \
    --vcf-annotation calls.EB.annot.tsv --exclude-sample <matched normal, if it is in the PON> -t 4
#    or straight from the tumor BAM: --tumor-bam T.bam --fasta genome.fa
#    PON_list.txt and matrix disagree in size (legacy HAEv7 PON: 26 vs 28)? exclude by column: --exclude-index 7

# VCF -> ANNOVAR-coordinate table (with the original VCF coordinates as vcf_* columns)
ebscore vcf2table --vcf calls.vcf.gz -o calls.tsv
```

`score` output: the input columns, then `PON+ PON- EB` (`PON±` = `alt1|…|altN<depth1|…|depthN` per strand;
`--keep-matrix` adds `Tumor` and `AB`). No variants give a header-only output. `--vcf-annotation` writes
`#CHROM POS REF ALT EB` (VCF coordinates, one line per scored ALT, sorted by chromosome and POS) plus `<file>.hdr`,
the INFO header with `Number=A`, so multi-allelic records keep one EB per ALT:

```bash
bgzip -c calls.EB.annot.tsv > ann.tsv.gz && tabix -s1 -b2 -e2 ann.tsv.gz
bcftools annotate -a ann.tsv.gz -h calls.EB.annot.tsv.hdr -c CHROM,POS,REF,ALT,INFO/EB calls.vcf.gz -Oz -o calls.EB.vcf.gz
```

### Options

| Command | Option | Meaning |
|---|---|---|
| all | `-t/--threads N`, `-v/--verbose`, `--quiet` | parallelism (not `vcf2table`); debug logging with shell commands; warnings only |
| `pileup` | `--normal`, `--tumor`, `--fasta`, `--bed`, `-o/--output` | inputs (normal = sample 1) and the `.pileup.gz` |
| | `--chrom C` | only these chromosomes (repeatable; default: whole BED) |
| | `--min-mapq`, `--min-baseq` | samtools `-q`/`-Q` (20/25) |
| | `--tmpdir DIR` | scratch for the per-chromosome parts (default: the output's folder) |
| `pon-matrix` | `--bam-list`, `--fasta`, `--bed`, `-o/--outdir` | one BAM/CRAM per line (relative paths: against the list's folder, symlinks not followed) |
| | `--chrom`, `--min-mapq`, `--min-baseq` | as for `pileup` |
| `score` | `--variants`, `--ponmatrix`, `-o/--output` | inputs and the output TSV |
| | `--pileup` or `--tumor-bam` + `--fasta` | tumor counts from the paired pileup or straight from the BAM |
| | `--tumor-column ReadN` | pileup column of the tumor (default `Read2`, else the last `Read` column) |
| | `--exclude-sample NAME`, `--exclude-index N` | drop PON samples (matched normal) by list entry/file name, or by 1-based column |
| | `--penalty`, `--fisher`, `--missing` | see the next table |
| | `--chrom C` | only variants on these chromosomes |
| | `--allow-missing-chrom` | EB NA instead of exit 2 when no variant chromosome has a matrix file |
| | `--keep-matrix`, `--vcf-annotation FILE`, `--tmpdir DIR` | extra columns; bcftools annotation (+ `.hdr`); scratch (default: system temp) |
| | `--min-mapq`, `--min-baseq` | for `--tumor-bam` |
| `vcf2table` | `--vcf`, `-o/--output` | VCF -> ANNOVAR-coordinate table |

### Options that change values

| Option | Default | Legacy behaviour |
|---|---|---|
| `--fisher sf\|legacy` | `sf`: chi2 survival function | `legacy`: `1 - chi2.cdf`; no EB between ~16 and 60 |
| `--missing na\|zero` | `na`: no tumor or PON data → EB `NA` | `zero`: EB 0 (looks like an artifact) |
| `--penalty` | 0.5 | was hard-coded 0.5 (config value ignored) |
| `--min-mapq/--min-baseq` | 20 / 25 | same |

### Pileup flags

`samtools mpileup -f <fasta> -q 20 -Q 25 -l <bed> [-r <chrom>]`, samtools defaults otherwise (BAQ on, anomalous pairs
skipped, overlap detection on, `--ff UNMAP,SECONDARY,QCFAIL,DUP`, max depth 8000). Build the PONmatrix and the
sample pileups with the same samtools version, flags and BED flavour.

Chromosomes are those of the BED (plain or gzipped, tab- or blank-delimited), or of `--chrom`, that every BAM header
has; missing ones are skipped with a warning, and the run fails if none is left. Rows come in **BAM header order**
with any `-t` (myCNV's `filterBed.mawk` walks the pileup in this order). Every output has its header line, also
when nothing is covered. The threaded pileup is a multi-member gzip (read it with `gzip -dc`/`zcat`, Python or
pandas); its scratch parts go next to the output unless `--tmpdir` is set.

### PONmatrix folder

```
PONmatrix/
  PON_list.txt          BAM paths in matrix column order
  matrix/<chrom>.pon.gz Chr Start Ref A G C T I D Depth; cells c1|…|cN=c1|…|cN (+ strand = - strand)
```

The same `matrix/` feeds the PON-CNV of myCNV. Scoring reads the folder only (temp files go to `--tmpdir`), so it
can be mounted read-only. A `--ponmatrix` without `matrix/`, or without a file for any variant chromosome, fails
(exit 2).

### The fit

`score` fits the per-strand beta-binomial on the fly at the called positions (identical PON vectors once per run):
L-BFGS-B with the analytic gradient, restarted until the point is stationary, i.e. the penalised optimum; AB values
to 6 significant digits. There is no other fit method (plan D19) and no AB cache (plan D32): at ~14 ms per fit a WES
tumor needs minutes, a precomputed cache ~800 CPU-hours per PON, and the matched-normal exclusion needs on-the-fly
fits anyway. On realistic PONs EB differs from the legacy (Genomon EBFilter) fit by several units, so EB thresholds
set on legacy values do not carry over (see `CHANGELOG.md`). Across scipy versions EB agrees to ~0.001, not bit for
bit.

## Tests

```bash
pytest -q          # mawk tools, core math, legacy parity, VCF conversion, scoring options, simulated BAMs
```

`tests/data/legacy/` holds outputs of the unmodified `main` code (py 3.9, pandas 1.5.2, scipy 1.9.3), made by
`tests/legacy/make_legacy_outputs.py`; `tests/data/make_synthetic.py` makes the synthetic pileups and
`tests/data/make_fits_exact.py` the golden exact fits (`fits_exact.tsv`, incl. large-PON vectors). `fits_hard.tsv`
holds PON vectors where a single L-BFGS-B run stops at a non-stationary point. The suite also runs in
`somvar-legacy` (`PYTHONPATH=src`). The parity tests rebuild the legacy EB from the port's lookups and the stored
legacy AB values (`tests/data/legacy/ponmatrix/AB/`, reference data only). `tests/test_real_pon.py` is opt-in
(`EBSCORE_REAL_PONMATRIX=<PONmatrix folder>`): it checks the fit on real PON rows and prints the EB shift against a
legacy-style reference fit; real PON data never goes into the repo.

## Docker

```bash
docker build --target test .          # pytest inside the build env (not shipped)
docker build -t ebscore:<version> .   # runtime: debian slim + /opt/conda, non-root, no entrypoint
```

CI (`.github/workflows/ci.yml`, plan D24): on push/PR to `kumodev`, pytest in a micromamba env from `environment.yml`,
then the Docker test stage, the runtime image, a runtime smoke check and Trivy 0.69.3: a Dockerfile config scan,
the vulnerability gate (library packages, fixed vulnerabilities, HIGH/CRITICAL; accepted findings go into `.trivyignore`
with a reason) and an informational report that includes the OS packages.
