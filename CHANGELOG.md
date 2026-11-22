# Changelog

## 1.0.0.dev0 — branch `kumodev` (2026-11-22)

Port of `main` (HEAD `ef67010`) for kumo-variant-calling-nf. Values were checked against the unmodified `main` code
(py 3.9, pandas 1.5.2, scipy 1.9.3) on synthetic data (`tests/data/legacy/`).

### Structure
- Installable package `src/ebscore/` with the console script `ebscore` (`pileup`, `pon-matrix`, `score`, `vcf2table`).
  No config dict needed; every input is a path or option.
- mawk tools are package data, start with `#!/usr/bin/env bash` and always run through `bash`. Under `#!/bin/sh` (dash
  on Debian/Ubuntu) the bash-only option parsing silently failed, so every option fell back to its default.
- Current Python (≥ 3.10), pandas 2/3, numpy 2, scipy 1.18; unpinned `environment.yml`; Dockerfile with named
  `build`/`test`/`runtime` stages.
- Removed: `EBFilter.sh`, Mac envs, `.vscode`, notebook checkpoints, `tumor2matrix.mawk`, `csv2bed.mawk`, the zero
  cache (`zerocache.py`), the Snakemake-oriented `ebrun`/`file2matrix` entry points.
- The bundle FASTA with `-r <chrom>` replaces the per-chromosome split FASTAs.
- Splitting (plan D16): `pon-matrix` works per chromosome; `pileup` and `score` run on the whole sample by default,
  `--chrom` restricts them.
- One BED flavour (plan D20): the kit's **padded** BED for the PONmatrix, the paired pileup and myCNV alike (legacy:
  covered for the samples, padded for the PON matrix). Help texts and examples say so.
- GitHub Actions CI (`.github/workflows/ci.yml`, plan D24): tests in a micromamba env, then the Docker stages.

### Unchanged (verified)
- `cleanpileup` and `pile2count` output is byte-identical to the legacy tools on the legacy fixtures and on synthetic
  pileups with indels and start markers carrying pileup-like MAPQ characters (`^$`, `^+`, `^-`).
- Beta-binomial likelihood, penalty term, start values, bounds, upper-tail p-value, Fisher combination (with
  `--fisher legacy`), EB rounding and the 0/60 caps.
- The port's lookups (tumor counts, PON strings) plus the legacy AB values (the cache the `main` code built, kept as
  test data) give exactly the legacy `run_ebscore` EB for 33 of the 35 legacy rows. The other 2 are the sorted-merge
  bug below.

### Fixed (deliberate behaviour changes)
- **Lost PON data after a PON gap.** `tumor2matrix.mawk` consumed and dropped the PON line that follows a variant
  position missing from the PONmatrix, so the next variant got no PON data and EB 0 (in the test set a real SNV went
  from EB 3.036 to 0). Positions are now looked up with a hash filter (`filterpos.mawk`).
- **Variants silently dropped.** The legacy chain kept only tumor pileup lines with an SNV mismatch (`$7~/[ATCGatcg]/`),
  so pure-indel sites vanished. It collapsed multi-allelic positions to the last allele and lost rows on the inner
  merge. Now every input variant is in the output (20 of 55 test variants were missing before).
- **On-the-fly AB fitting crashed** (`tumor_matrix.py:80` unpacked 3 values from `zero2AB`, which returns 2). It is now
  the only path.
- **AB cache removed** (plan D32, owner 2026-11-22): no `ab-cache` command, no `--ab-cache`/`--no-ab-cache`, no
  `AB/` folder in the PONmatrix. One exact fit takes ~14 ms on one core, so a WES tumor (10,000 merged calls, ≤ 20,000
  fits) needs ~5 CPU-minutes; a precomputed cache needs 12 fits per PON position (~200 M fits, ~800 CPU-hours per padded
  HAEv7 PON), several GB, and invalidation rules. Excluding the matched normal needs on-the-fly fits anyway. An `AB/`
  folder left in an old PONmatrix is ignored.
- **Fit stopped far from the optimum.** The legacy L-BFGS-B with a numerical gradient (also what Genomon EBFilter
  does) stops early on flat likelihoods, i.e. all-zero alt vectors, which are most PON cells. Where it stops depends
  on the scipy version: with scipy 1.18 only 19 of 102 test fits match scipy 1.9.3 (relative differences up to 0.6).
  For an all-zero vector the penalised optimum is a = 0.1 with b growing with the total PON depth, while the legacy
  fit stops at b of a few hundred. The port uses the analytic (digamma) gradient, tight tolerances
  and a stationarity check. A single L-BFGS-B run can report convergence at a non-stationary point (3 of 4,000
  low-depth sparse-alt vectors with scipy 1.18.1 and 6 others with 1.9.3, objective up to 1.9 worse than legacy), so
  the fit restarts from its own result until the scaled projected gradient is ≤ 1e-3 (at most 5 restarts; if that
  never happens, the best point over three more start values, each with restarts). AB values are rounded to 6
  significant digits (legacy: 5 decimals). (EBCall's R code is a third variant: it uses fixed (0.1, 1) for all-zero
  vectors.)
  - **The legacy fit is gone** (plan D19, owner 2026-11-22): there is no `--fit` option and no `method=` argument; the
    exact fit is the only fit. The tests compare it with the stored legacy fits and with a test-only legacy-style
    reference optimiser.
  - **EB impact on realistic PONs is several units, not small.** On the 6-normal test set (depth 20–80) EB changes by
    ≤ 0.004, because the legacy fit happens to stop close to the optimum there. On the 28-normal HAEv7 PON
    (`somvar_HAEv7_NovaSeq_v1`, 1,200 sampled chr21 strand vectors), exact minus legacy EB for a tumor with 2+2 alt
    reads at 60+60× ranges from −3.98 to +1.52, and 142 of 1,200 sites cross EB 4. With 100 all-zero normals at
    ~300×, EB for 2 alt reads goes from 4.1 (legacy) to 8.2 (exact). EB thresholds set on legacy values
    (`filter_settings.xlsx`) therefore do not carry over; they are re-set on the new PONs (plan D19).
  - **Not bit-identical across scipy versions.** b of all-zero vectors sits on a very flat objective, so its last
    digits depend on the optimiser's stopping point: between scipy 1.9.3 and 1.18.1, 121 of 420 WES-like exact fits
    differ after rounding (relative ≤ 2e-4), which changes EB by ≤ 0.001. Compare AB values across environments with
    a tolerance, not byte for byte.
- **EB precision ceiling.** `1 - chi2.cdf` cannot go below ~2.2e-16, so no EB fell between ~16 and 60. The default
  `--fisher sf` uses the survival function. The two are identical below EB ~7 (relative difference < 1e-9), within
  0.005 EB below 15, and diverge only near the legacy ceiling of ~15.95, so thresholds in the usual range (EB 3–5)
  are unaffected (plan D-Q4, now D-Q15).
- **No data → EB 0.** A variant without PON or tumor coverage got EB 0, i.e. it looked like an artifact. The default
  is now `NA` (`--missing zero` restores 0). MNVs/complex alleles also get `NA` (the legacy code returned garbage for
  them).
- **Penalty ignored.** The configured `fit_pen` was never passed to the fit; `--penalty` now is.
- **Matched-normal exclusion** by file-name prefix (`<patient>_…`, zeros stripped) is replaced by explicit
  `--exclude-sample NAME` / `--exclude-index N`. Name matching refuses to run when `PON_list.txt` and the matrix
  disagree in size (the legacy HAEv7 PON has 28 matrix columns and 26 list entries).
- **PONmatrix needed write access** (`temp_dir` defaulted to `<pon>/temp`). Temp files now go to `--tmpdir`, and
  `pon-matrix` writes the `PON_list.txt` it used, so list and matrix always agree.
- **Zero cache** (shared mutable files, race-prone, `flatten_df` a silent no-op on pandas 1 and a crash on pandas 3,
  `save_next_zero` crashing on a fresh folder) is replaced by exact in-run de-duplication of identical PON vectors.
- `cleanpileup`: `-s N` was compared as a string (ignored with ≥ 10 samples); `N` in indel bases now counts as forward
  strand; the chromosome name is no longer used as a printf format.
- `pile2count`: the dot is counted with a regex literal, so gawk now gives the same counts as mawk.
- `cleanpileup` and `pile2count` split their input on tabs only (`FS = "\t"`). The legacy tools used mawk's default
  field splitting (runs of blanks), which collapsed empty fields and shifted the columns after them. Real `samtools
  mpileup` output has no empty fields (`0 * *` for no coverage) and no blanks in reference names, so it is unaffected.

### Pileup and PONmatrix (port behaviour, fixed after review)
- Chromosomes are the BED's (or `--chrom`'s, de-duplicated) that every BAM header has, in BAM header order. A BED
  contig missing from a BAM header is skipped with a warning (before, `-t 1` skipped it silently and `-t N` failed);
  if no chromosome is left, the run fails (exit 2) instead of writing an empty file. So `pileup -t 1` and `-t N` give
  the same rows in the same order also with a lexically sorted vendor BED (myCNV's `filterBed.mawk` relies on a
  defined order), and one PON BAM without a BED contig no longer aborts `pon-matrix`.
- Every pileup and PONmatrix file has its header line, also when no position is covered (before: a 0-byte gzip that
  `score` rejected; an empty chromosome list made the threaded pileup read stdin).
- BED files may be gzipped and blank-delimited (as `samtools -l` accepts).
- The threaded pileup writes gzipped per-chromosome parts (concatenated as a multi-member gzip) into a scratch folder
  next to the output (`--tmpdir` overrides), not ~15 GB of plain text per WES pair into the system temp. On the first
  failing chromosome the pending ones are cancelled.
- Repeated `--chrom` values count once (before: duplicated pileup rows, and a `pon-matrix` race on the same file).
- `pon-matrix`: relative BAM-list entries resolve against the folder the list is seen in, without following a symlink
  (Nextflow stages the list as a symlink next to the BAMs). `PON_list.txt` is written via a temp file + rename, so a
  staged input symlink of that name is replaced instead of rewriting its target.

### Scoring, inputs and outputs (port behaviour, fixed after review)
- **VCF alleles with right-hand context.** The shared trailing bases are trimmed before the shared leading ones (as
  `bcftools norm` does). In a multi-allelic record the shorter alleles keep the longest REF's context (Mutect2 writes
  `ATT AT,ATTT`); before, such indels were shifted to the right and scored EB 0 (artifact-like, also in the VCF
  annotation), and SNVs/insertions with context got NA. Indels inside repeats still need left-normalisation upstream
  (`bcftools norm -f <fasta>`).
- Fit keys are flat int tuples (about a fifth of the memory of a tuple of pairs).
- `--vcf-annotation`: sorted by numeric POS also for `vcf2table` tables (before, text order, which tabix rejects),
  and a `<file>.hdr` with `##INFO=<ID=EB,Number=A,Type=Float,...>` is written for `bcftools annotate -h` (with
  `Number=1` bcftools keeps only the first allele's EB of a multi-allelic record).
- No variants (empty VCF, header-only table, `--chrom` without a match) gives a header-only output and exit 0 (before:
  KeyError, exit 1).
- A `--ponmatrix` without `matrix/`, or whose `matrix/` has no file for any variant chromosome, fails with exit 2
  (before: exit 0 with EB NA everywhere); `--allow-missing-chrom` restores the NA output for the second case.
- Python API: `EB` is float64 (NaN = no data) instead of an object column with `pd.NA`; VCF `Start`/`End`/`vcf_pos`
  are `Int64`.
- The synthetic test VCF writes the multi-allelic sites as one record with two ALTs (and declares its INFO field).
