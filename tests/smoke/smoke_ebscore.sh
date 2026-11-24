#!/usr/bin/env bash
# Smoke test of the ebscore runtime image on fixture tier A (run INSIDE the container).
#   podman run --rm --platform linux/amd64 -v "$PWD/tests/data/plumbing:/data:ro" -v "$PWD/tests/smoke:/smoke:ro" \
#       localhost/ebscore:<version> bash /smoke/smoke_ebscore.sh /data
# Writes to /tmp/smoke inside the container; compares with /data/expected/ebscore. Exit 0 = all PASS.
set -uo pipefail

D=${1:-/data}
S=$(cd "$(dirname "$0")" && pwd)
E=$D/expected/ebscore
O=/tmp/smoke/ebscore
mkdir -p "$O"
fail=0
step() { echo "== $*"; }
check() { if "$@"; then return 0; else fail=1; return 1; fi; }

step "environment"
echo "user=$(id -un) uid=$(id -u)  (must not be root)"
[ "$(id -u)" != "0" ] || { echo "FAIL runs as root"; fail=1; }
ebscore --version
samtools --version | head -1
mawk -W version 2>&1 | head -1
python -c "import pandas, numpy, scipy; print('pandas', pandas.__version__, 'numpy', numpy.__version__, 'scipy', scipy.__version__)"
ps -o rss= -p $$ >/dev/null && echo "ps OK" || { echo "FAIL ps (procps) missing"; fail=1; }

step "pon-matrix"
check ebscore pon-matrix --bam-list "$D/pon/PON_list.txt" --fasta "$D/ref/chr21.fa.gz" \
    --bed "$D/bed/testkit.padded.bed" -o "$O/PONmatrix" --quiet
check python "$S/compare_tables.py" "$E/PONmatrix/matrix/chr21.pon.gz" "$O/PONmatrix/matrix/chr21.pon.gz"
check cmp "$E/PONmatrix/PON_list.txt" "$O/PONmatrix/PON_list.txt"

step "pileup (1 thread and 2 threads must agree)"
check ebscore pileup --normal "$D/bam/normal.bam" --tumor "$D/bam/tumor.bam" --fasta "$D/ref/chr21.fa.gz" \
    --bed "$D/bed/testkit.padded.bed" -o "$O/pair.pileup.gz" --quiet
check zcmp "$E/pair.pileup.gz" "$O/pair.pileup.gz" && echo "PASS pair.pileup.gz identical to expected"
check ebscore pileup --normal "$D/bam/normal.bam" --tumor "$D/bam/tumor.bam" --fasta "$D/ref/chr21.fa.gz" \
    --bed "$D/bed/testkit.padded.bed" -o "$O/pair.t2.pileup.gz" -t 2 --quiet
check zcmp "$O/pair.pileup.gz" "$O/pair.t2.pileup.gz" && echo "PASS 2-thread pileup identical"

step "score (pileup route, VCF input, annotation side output)"
check ebscore score --variants "$D/variants/variants.vcf.gz" --ponmatrix "$O/PONmatrix" --pileup "$O/pair.pileup.gz" \
    -o "$O/variants.EB.tsv" --vcf-annotation "$O/variants.EB.annot.tsv" --keep-matrix --quiet
check python "$S/compare_tables.py" "$E/variants.EB.tsv" "$O/variants.EB.tsv" --col-atol EB=0.002
check python "$S/compare_tables.py" "$E/variants.EB.annot.tsv" "$O/variants.EB.annot.tsv" --col-atol EB=0.002
check cmp "$E/variants.EB.annot.tsv.hdr" "$O/variants.EB.annot.tsv.hdr" && echo "PASS variants.EB.annot.tsv.hdr identical"

step "score (tumor BAM route must equal the pileup route)"
check ebscore score --variants "$D/variants/variants.tsv" --ponmatrix "$O/PONmatrix" --tumor-bam "$D/bam/tumor.bam" \
    --fasta "$D/ref/chr21.fa.gz" -o "$O/variants.EB.frombam.tsv" --quiet
cut -f1-5,10-12 "$O/variants.EB.tsv" > "$O/eb_pileup_cols.tsv"
check python "$S/compare_tables.py" "$O/eb_pileup_cols.tsv" "$O/variants.EB.frombam.tsv"

step "error handling (missing input must exit non-zero)"
if ebscore score --variants /nonexistent.tsv --ponmatrix "$O/PONmatrix" --pileup "$O/pair.pileup.gz" -o "$O/x.tsv" --quiet 2>/dev/null; then
    echo "FAIL missing input exited 0"; fail=1
else
    echo "PASS missing input exits non-zero"
fi

echo
if [ "$fail" = 0 ]; then echo "ALL PASS (ebscore)"; else echo "SOME CHECKS FAILED (ebscore)"; fi
exit $fail
