# Changelog

## 1.0.0.dev0 — branch `kumodev` (2026-11-22)

Port of `main` (HEAD `ef67010`) for kumo-variant-calling-nf. Values were checked against the unmodified `main` code
(py 3.9, pandas 1.5.2, scipy 1.9.3) on synthetic data (`tests/data/legacy/`).

### Structure
- Installable package `src/ebscore/` with the console script `ebscore` (`pileup`, `pon-matrix`, `ab-cache`, `score`,
  `vcf2table`). No config dict needed; every input is a path or option.
- mawk tools are package data, start with `#!/usr/bin/env bash` and always run through `bash`. Under `#!/bin/sh` (dash
  on Debian/Ubuntu) the bash-only option parsing silently failed, so every option fell back to its default.
- Current Python (≥ 3.10), pandas 2/3, numpy 2, scipy 1.18; unpinned `environment.yml`; Dockerfile with named
  `build`/`test`/`runtime` stages.
- Removed: `EBFilter.sh`, Mac envs, `.vscode`, notebook checkpoints, `tumor2matrix.mawk`, `csv2bed.mawk`, the zero
  cache (`zerocache.py`), the Snakemake-oriented `ebrun`/`file2matrix` entry points.
- The bundle FASTA with `-r <chrom>` replaces the per-chromosome split FASTAs.
- Splitting (plan D16): `pon-matrix` and `ab-cache` work per chromosome (and in row chunks); `pileup` and `score` run
  on the whole sample by default, `--chrom` restricts them.

### Unchanged (verified)
- `cleanpileup` and `pile2count` output is byte-identical to the legacy tools on the legacy fixtures and on synthetic
  pileups with indels and start markers carrying pileup-like MAPQ characters (`^$`, `^+`, `^-`).
- Beta-binomial likelihood, penalty term, start values, bounds, upper-tail p-value, Fisher combination (with
  `--fisher legacy`), EB rounding and the 0/60 caps.
- With the same AB cache, EB equals the legacy `run_ebscore` output for 33 of the 35 legacy rows. The other 2 are the
  sorted-merge bug below.
- With `--fit legacy` and scipy 1.9.3, all 102 test fits reproduce the legacy AB values exactly.

### Fixed (deliberate behaviour changes)
- **Lost PON data after a PON gap.** `tumor2matrix.mawk` consumed and dropped the PON line that follows a variant
  position missing from the PONmatrix, so the next variant got no PON data and EB 0 (in the test set a real SNV went
  from EB 3.036 to 0). Positions are now looked up with a hash filter (`filterpos.mawk`).
- **Variants silently dropped.** The legacy chain kept only tumor pileup lines with an SNV mismatch (`$7~/[ATCGatcg]/`),
  so pure-indel sites vanished. It collapsed multi-allelic positions to the last allele and lost rows on the inner
  merge. Now every input variant is in the output (20 of 55 test variants were missing before).
- **On-the-fly AB fitting crashed** (`tumor_matrix.py:80` unpacked 3 values from `zero2AB`, which returns 2). It is now
  the default path, and it is the only one for PONs without an AB cache (e.g. `somvar_HAEv7_NovaSeq_v1`).
- **Fit not reproducible across scipy versions.** The numerical-gradient L-BFGS-B stopped early on flat likelihoods
  (all-zero alt vectors). With scipy 1.18 only 19 of 102 test fits matched scipy 1.9.3 (relative differences up to
  0.6). The default `--fit exact` uses the analytic (digamma) gradient and tight tolerances: never a worse objective
  than legacy (up to 0.108 better), equal across scipy versions to the 5th decimal. EB changes by ≤ 0.004 on the test
  set. `--fit legacy` keeps the old optimiser settings.
- **EB precision ceiling.** `1 - chi2.cdf` cannot go below ~2.2e-16, so no EB fell between ~16 and 60. The default
  `--fisher sf` uses the survival function. Below ~16 both give the same values (relative difference < 1e-9), so
  thresholds in that range are unaffected (plan D-Q4).
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
- AB cache coordinates: `Chr Start Ref` are always written (the `61cba31` cache lacked them).
