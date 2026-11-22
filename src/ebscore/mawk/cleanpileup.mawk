#!/usr/bin/env bash
# cleanpileup v1.2.0 (ebscore kumodev)
# Cleans `samtools mpileup` text output (3 + 3*N columns, no extra output columns).
# - removes read-start markers (^ + mapping-quality char) and read-end markers ($)
# - replaces each "<base>[+-]<len><seq>" with I/i (insertion) or D/d (deletion);
#   upper case = indel bases on the forward strand
# - drops the quality columns and prints a header
#
# USAGE: samtools mpileup ... | bash cleanpileup.mawk [-d|--output_depth] [-s|--use-samples N]
#   -d  also print the depth column of each sample (Chr Start Ref Depth1 Read1 Depth2 Read2 ...)
#   -s  only print the first N samples
#
# Changes vs v1.1.0 (see CHANGELOG.md): bash shebang (the option parser is bash-only),
# numeric -s comparison, printf with an explicit format, N counts as a forward-strand indel base.
set -euo pipefail

outputDepth=0
useSamples=0
while (( "$#" )); do
    [[ $1 == --*=* ]] && set -- "${1%%=*}" "${1#*=}" "${@:2}"
    case "$1" in
        -d|--output_depth|--output-depth)
            outputDepth=1
            shift
            ;;
        -s|--use-samples)
            if [[ -n "${2:-}" && ${2:0:1} != "-" ]]; then
                useSamples=$2
                shift 2
            else
                echo "<cleanpileup> Error: sample count argument is missing [-s|--use-samples N]" >&2
                exit 1
            fi
            ;;
        -*)
            echo "<cleanpileup> Error: Unsupported flag $1" >&2
            exit 1
            ;;
        *)
            echo "<cleanpileup> Error: unexpected argument $1 (input is read from stdin)" >&2
            exit 1
            ;;
    esac
done

mawk -v printDepth="$outputDepth" -v useSamples="$useSamples" '
BEGIN {
    FS = "\t";
    printDepth = printDepth + 0;
    useSamples = useSamples + 0;
}
NR == 1 {
    samples = (NF - 3) / 3;
    if (useSamples > 0 && useSamples < samples) {
        samples = useSamples;
    }
    printf("Chr\tStart\tRef");
    for (s = 0; s++ < samples;) {
        if (printDepth) {
            printf("\tDepth%s\tRead%s", s, s);
        } else {
            printf("\tRead%s", s);
        }
    }
    printf("\n");
}
{
    printf("%s\t%s\t%s", $1, $2, $3);
    for (i = 0; i++ < samples;) {
        col = (i * 3) + 2;
        read = $col;
        # remove read-start (^ + MAPQ char) and read-end ($) markers
        gsub(/\^[^\t]|\$/, "", read);
        # replace <base><+|-><len><seq> with I/i or D/d
        while (match(read, /[+-][0-9]+/)) {
            pre = substr(read, 1, RSTART - 2);
            indel = substr(read, RSTART, 1);
            l = substr(read, RSTART + 1, RLENGTH - 1) + 0;
            indel_bases = substr(read, RSTART + RLENGTH, l);
            post = substr(read, RSTART + RLENGTH + l);
            fwd = match(indel_bases, /[ACGTN]/);
            if (indel == "-") {
                base = fwd ? "D" : "d";
            } else {
                base = fwd ? "I" : "i";
            }
            read = pre base post;
        }
        if (printDepth) {
            printf("\t%s\t%s", $(col - 1), read);
        } else {
            printf("\t%s", read);
        }
    }
    printf("\n");
}'
