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
#    Chr Start Ref Depth1 Read1 Depth2 Read2 (normal = 1, tumor = 2)
ebscore pileup --normal N.bam --tumor T.bam --fasta genome.fa --bed kit.padded.bed -o pair.pileup.gz -t 4
#    Tumor-only sample (plan T1): leave out --normal -> Chr Start Ref Depth1 Read1 (tumor = 1)
ebscore pileup --tumor T.bam --fasta genome.fa --bed kit.padded.bed -o T.pileup.gz -t 4

# 3. EB per variant (table in ANNOVAR coordinates, or a VCF; left-align indels first: bcftools norm -f genome.fa).
ebscore score --variants calls.vcf.gz --ponmatrix PONmatrix/ --pileup pair.pileup.gz -o calls.EB.tsv \
    --vcf-annotation calls.EB.annot.tsv --exclude-sample <matched normal, if it is in the PON> -t 4
#    or straight from the tumor BAM: --tumor-bam T.bam --fasta genome.fa
#    PON_list.txt and matrix disagree in size (legacy HAEv7 PON: 26 vs 28)? exclude by column: --exclude-index 7
#    Tumor-only: --pileup T.pileup.gz (the tumor column defaults to Read1 there) or --tumor-bam T.bam --fasta genome.fa;
#    both give the same EB as the paired pileup of the same tumor.

# 3a. Run-level PONmatrix subset (plan T2): the matrix rows at the positions of all variant files of a run, once;
#     every score task then stages the subset (MBs) instead of the full matrix (~20 GB at 100 normals). Same EB.
ebscore pon-subset --ponmatrix PONmatrix/ -o run_subset/ --variants p1.union.tsv --variants p2.union.tsv -t 8
ebscore pon-subset --ponmatrix run_subset/ -o p1_subset/ --variants p1.union.tsv -t 8  # per pair, from the run subset
ebscore score --variants p1.union.tsv --ponmatrix p1_subset/ --pileup p1.pileup.gz -o p1.EB.tsv ...

# 4. PON-CNV from the same pileup: myCNV `mycnv sample --pileup pair.pileup.gz ...`.
#    Tumor-only is PLANNED (myCNV plan M1, not in mycnv yet; flag names may change):
#    `mycnv sample --tumor-only --pileup T.pileup.gz --snp-sites <kit common SNP sites> ...`

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
| `pileup` | `--normal`, `--tumor`, `--fasta`, `--bed`, `-o/--output` | inputs (normal = sample 1, tumor = sample 2) and the `.pileup.gz` |
| | no `--normal` | tumor-only pileup: `Chr Start Ref Depth1 Read1`, tumor = sample 1 (plan T1) |
| | `--chrom C` | only these chromosomes (repeatable; default: whole BED) |
| | `--min-mapq`, `--min-baseq` | samtools `-q`/`-Q` (20/25) |
| | `--tmpdir DIR` | scratch for the per-chromosome parts (default: the output's folder) |
| `pon-matrix` | `--bam-list`, `--fasta`, `--bed`, `-o/--outdir` | one BAM/CRAM per line (relative paths: against the list's folder, symlinks not followed) |
| | `--chrom`, `--min-mapq`, `--min-baseq` | as for `pileup` |
| `pon-subset` | `--ponmatrix`, `-o/--outdir` | source PONmatrix (full set or subset; read only) and a new folder (not the source, no `matrix/` yet) |
| | `--variants FILE`, `--positions FILE` | repeatable, combinable, at least one: variant tables/VCFs as `score` reads them; `Chr Pos` files (tab, no header, the `positions.tsv` format) |
| | `--chrom C` | only positions on these chromosomes (as `score --chrom`) |
| | `--allow-missing-chrom` | write the subset even when no chromosome with positions has a matrix file (checked for a full-set source only) |
| `score` | `--variants`, `--ponmatrix`, `-o/--output` | inputs and the output TSV |
| | `--pileup` or `--tumor-bam` + `--fasta` | tumor counts from an `ebscore pileup` file (paired or tumor-only) or straight from the BAM |
| | `--tumor-column ReadN` | pileup column of the tumor (default `Read2` if the header has it, i.e. paired; else the last `Read` column, i.e. `Read1` of a tumor-only pileup) |
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
pandas); its scratch parts go next to the output unless `--tmpdir` is set. The tumor-only pileup plans its
chromosomes on the tumor's header alone.

Which positions are rows: samtools prints a BED position when at least one read of any input BAM passes the read
filters (`-q`, `--ff`, anomalous pairs skipped) and overlaps it; the depth column counts only the bases that pass `-Q`
(after BAQ), so a row can read `0 *`. Reads that fail a read filter count for nothing, not even for the row. Hence the
tumor-only pileup has a row exactly where a filter-passing tumor read overlaps, and there it equals the paired
pileup's tumor columns (`Depth2 Read2`, `0 *` included). The paired-only rows are the positions that no
filter-passing tumor read overlaps but a filter-passing normal read does; the paired file has tumor `0 *` there. They
are the normal-only positions, and also those where every tumor read fails a read filter (e.g. MAPQ < 20, duplicates,
anomalous pairs). EB is the same either way: no tumor depth gives `NA` (`--missing`).

### PONmatrix folder

```
PONmatrix/
  PON_list.txt          BAM paths in matrix column order
  matrix/<chrom>.pon.gz Chr Start Ref A G C T I D Depth; cells c1|…|cN=c1|…|cN (+ strand = - strand)
```

The same `matrix/` feeds the PON-CNV of myCNV. Scoring reads the folder only (temp files go to `--tmpdir`), so it
can be mounted read-only. A `--ponmatrix` without `matrix/`, or without a file for any variant chromosome, fails
(exit 2).

### Run-level PONmatrix subset (`pon-subset`)

**Purpose.** `score` reads only the matrix rows at its variant positions, but a pipeline task has to stage the whole
`matrix/` folder first: ≈ 6 GB at 28 normals, ≈ 20 GB at 100, once per pair (≈ 2 TB per 100 tumors). `pon-subset`
reads the full matrix once and keeps the rows at the positions that will be scored, so each EB task stages a few MB.

**What it writes** (`-o DIR`):

```
DIR/
  matrix/<chrom>.pon.gz  header + the source rows at the listed positions, unchanged and in source order (filterpos);
                         one file per chromosome that has positions and a source file (header-only if none matches)
  PON_list.txt           byte copy of the source's, so --exclude-sample works as before; if the source has none, the
                         subset has none either (a stale PON_list.txt in DIR is removed, with a warning)
  positions.tsv          Chr Pos, tab, no header: the union of the positions, chromosomes in order of first appearance
                         in the inputs (all --variants, then all --positions, as given), positions numeric
  subset.json            schema ebscore-ponmatrix-subset/1, ebscore version, source (path as given, its parsed
                         manifest.json or null, its subset.json if the source is a subset), inputs (path, rows,
                         unscorable rows, positions), chrom_filter, n_positions, n_samples (from the Depth cells),
                         chroms {chrom: positions, rows, matrix}, positions_md5 (of positions.tsv)
```

The positions are exactly the `(Chr, Pos)` pairs `score` looks up for the same files (the same reader and lookup
code): the Start of an SNV or insertion, Start − 1 for a deletion (the anchor base of the pileup). Unscorable rows (MNV,
complex, symbolic, and rows without a valid position: empty Chr, or a position < 1 such as a deletion at Start 1) have
no position and get EB NA in `score` too; they are counted in the log and in `inputs`. A `--positions` file uses the
same `Pos`, so `--positions DIR/positions.tsv` rebuilds the same subset. Chromosomes are filtered in parallel (`-t`); the
output is byte-identical with any `-t`. The matrix is written in a hidden staging folder and moved to `DIR/matrix`
last, after the other files, so an interrupted run leaves no `matrix/` and nothing can score against it. Only
`matrix/` and `PON_list.txt` are carried over: the subset is for `score`; myCNV needs the full set (`cov/ snp/ CNV/`).

**Same EB.** Scoring against the subset gives the same output, byte for byte, as against the full set, also with
`--exclude-sample`/`--exclude-index`, `--keep-matrix` (`Tumor`, `AB`) and `--vcf-annotation`: the PON rows and
`PON_list.txt` are identical, so the fits are too. A listed position without a source row stays without one (EB NA, as
with the full set); a chromosome with positions but no source file gets no file (a warning, as `score` gives). If the
source is a full set and no chromosome with positions has a file there, `pon-subset` fails as `score` against the full
set would (a wrong `--ponmatrix`; `--allow-missing-chrom` writes it anyway). A subset source is not checked for this:
its chromosomes without a file are those the full set lacks, which the run-level call has already checked, so a
per-pair call writes the subset (with a warning) and that pair's `score` decides, as against the full set.

**The guard.** A folder with `subset.json` is a subset. `score` checks `positions.tsv` against its md5 in
`subset.json`, then refuses to run (exit 2) if any scorable position of its input (after `--chrom`) is not listed:
`N scoring positions outside this PONmatrix subset (first: chr12:5035 ...): rebuild it with these variants`. So a
wiring error (the wrong pair's subset, a subset built before a VCF changed) fails instead of silently giving EB NA. A
subset can be the source of a smaller one; the same check then refuses requested positions the source does not list.

**In the pipeline** (`EBSCORE_PONSUBSET`, design §3, N-Q8 option A + T2): one task per run, after the last pair's
caller VCFs exist and `EBSCORE_UNION` has written every union. It stages the full matrix once and runs
`pon-subset --ponmatrix ponmatrix -o run --variants <union> ... -t $task.cpus` over all unions of the run, then, per
union, `pon-subset --ponmatrix run -o pairs/<id> --variants <id>.union.tsv -t $task.cpus`. Each `EBSCORE_SCORE` task
stages only `pairs/<id>/` and passes it as `--ponmatrix`; `run/subset.json` is published as provenance. The run-level
call fails on a wrong `--ponmatrix` (no matrix file for any chromosome of the run). A per-pair call does not fail when
all of its pair's chromosomes lack a matrix file, so such a pair fails only its own `EBSCORE_SCORE` (as with the full
set), not the EB of the whole run. Pass the unions in a fixed order (e.g. `collect(sort: true)`), so the run's
`positions.tsv` and its md5 are reproducible; the matrix files do not depend on the order.

**Cost.** *Staged:* the full matrix once per run (≈ 20 GB at 100 normals, also the disk of that one task), then
10–30 MB per EB task (≈ 250 B per position at 100 normals × 50–100k union positions per WES pair). *Read*, all in the
subset task: one pass over the full matrix (CPU-bound on `gzip -dc` + mawk, parallel over chromosomes with `-t`), and
in the per-pair step the run subset once per pair. Each per-pair call streams the run subset's file of every chromosome
its pair has positions on and reads the run's whole `positions.tsv` (md5 + check). That is N_pairs × the run subset,
and the run subset grows with N_pairs, so this term grows with the square of the run size. At 10 pairs the run subset
is ≤ 250 MB and the step reads ≤ 2.5 GB, far below the full pass. At 100 pairs with mostly private positions the run
subset is ≈ 1.2–2.5 GB and the step reads ≈ 125–250 GB, 6–12 times the full pass (less where pairs share positions,
e.g. recurrent artifacts). The per-pair calls are independent and can run side by side. Staging the run subset in
every EB task instead (no per-pair step) moves the same N_pairs × run-subset volume into staging. A one-pass mode that
writes the per-pair subsets straight from the full matrix would make the step linear; it is not implemented. Scoring
waits for the last pair's VCFs (one barrier); EB itself is minutes per pair, so the barrier costs little next to
staging the full matrix per pair.

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
pytest -q          # mawk tools, core math, legacy parity, VCF conversion, scoring options, simulated BAMs, subsets,
                   # CLI contract
```

`tests/data/legacy/` holds outputs of the unmodified `main` code (py 3.9, pandas 1.5.2, scipy 1.9.3), made by
`tests/legacy/make_legacy_outputs.py`; `tests/data/make_synthetic.py` makes the synthetic pileups and
`tests/data/make_fits_exact.py` the golden exact fits (`fits_exact.tsv`, incl. large-PON vectors). `fits_hard.tsv`
holds PON vectors where a single L-BFGS-B run stops at a non-stationary point. The suite also runs in
`somvar-legacy` (`PYTHONPATH=src`). The parity tests rebuild the legacy EB from the port's lookups and the stored
legacy AB values (`tests/data/legacy/ponmatrix/AB/`, reference data only). `tests/test_real_pon.py` is opt-in
(`EBSCORE_REAL_PONMATRIX=<PONmatrix folder>`): it checks the fit on real PON rows and prints the EB shift against a
legacy-style reference fit; real PON data never goes into the repo.

**CLI contract (plan C1).** The kumo-variant-calling-nf modules call these subcommands with fixed flags (design §3
process table), so `tests/test_cli_contract.py` compares the parser with the snapshot `tests/data/cli_contract.json`.
Per command it holds the sorted option strings, and per option the flags, required, the action kind (`store`,
`store_true`, `append`, …), `nargs`, `type`, `choices` and the default when it is a constant; also the mutually
exclusive groups and positionals, and the `--version` output format with the version masked (`ebscore <version>`;
versions.yml parses it). Help texts are not part of it (argparse formats `--help` differently between Python
versions). The comparison is type-exact (a default `1` → `1.0` or `0` → `false` is a change). The test also parses
the module command lines of design §3 with prefix matching off, so a renamed flag fails even when the old one is a
prefix of the new one. After a deliberate CLI change:

```bash
python tests/make_cli_contract.py   # rewrites the JSON from this checkout's src/, prints the changes and the modules
```

Then commit the JSON with the change and review the kumo-variant-calling-nf modules that call the changed flags. The
failure message names the ebscore the test imported. If it warns about a stale install (another copy than this
checkout's `src/`, with other code), do not regenerate: reinstall (`pip install -e .`) or set `PYTHONPATH=src`.

## Docker

```bash
docker build --target test .          # pytest inside the build env (not shipped)
docker build -t ebscore:<version> .   # runtime: debian slim + /opt/conda, non-root, no entrypoint
```

CI (`.github/workflows/ci.yml`, plan D24): on push/PR to `kumodev`, pytest in a micromamba env from `environment.yml`,
then the Docker test stage, the runtime image, a runtime smoke check, the tier-A fixture smoke inside the image
(`tests/smoke/smoke_ebscore.sh` on `tests/data/plumbing/`, as the image user and as an arbitrary uid) and Trivy 0.69.3: a Dockerfile config scan,
the vulnerability gate (library packages, fixed vulnerabilities, HIGH/CRITICAL; accepted findings go into `.trivyignore`
with a reason) and an informational report that includes the OS packages.
