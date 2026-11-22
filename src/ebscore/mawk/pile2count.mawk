#!/usr/bin/env bash
# pile2count v1.1.0 (ebscore kumodev)
# Per-sample, per-strand counts of A G C T I D plus depth from cleanpileup output.
#
# USAGE: samtools mpileup ... | bash cleanpileup.mawk | bash pile2count.mawk
# INPUT (with header, no depth columns):
#   Chr  Start  Ref  Read1  Read2 ...
# OUTPUT (the PONmatrix format):
#   Chr  Start  Ref  A  G  C  T  I  D  Depth
#   each cell: <sample1|sample2|...> for the + strand, "=", the same for the - strand
#   Depth = letters + "." (+ strand) / letters + "," (- strand); "*", "#", "N", "<", ">" are not counted.
#
# Changes vs v1.0 (see CHANGELOG.md): bash shebang; the dot is counted with a regex literal, so the
# output no longer depends on mawk's handling of "\." in strings (gawk now gives the same result).
set -euo pipefail

mawk '
BEGIN {
    FS = "\t";
    BASECOLS = 3;
    SEP = "=";
    letterCount = split("Aa-Gg-Cc-Tt-Ii-Dd", LETTERS, "-");
}
NR == 1 {
    sampleCount = NF - BASECOLS;
    printf("Chr\tStart\tRef");
    for (l = 0; l++ < letterCount;) {
        printf("\t%s", substr(LETTERS[l], 1, 1));
    }
    printf("\tDepth\n");
    next;
}
{
    for (r = 0; r++ < sampleCount;) {
        read = $(r + BASECOLS);
        for (strand = 0; strand++ < 2;) {
            for (l = 0; l++ < letterCount;) {
                lettersFound = gsub(substr(LETTERS[l], strand, 1), "", read);
                COUNT[strand "-" r "-" l] = lettersFound;
                DEPTH[strand "-" r] += lettersFound;
            }
        }
        DEPTH["1-" r] += gsub(/\./, "", read);
        DEPTH["2-" r] += gsub(/,/, "", read);
    }
    printf("%s\t%s\t%s\t", $1, $2, $3);
    for (l = 0; l++ < letterCount;) {
        for (strand = 0; strand++ < 2;) {
            printf("%s", COUNT[strand "-1-" l]);
            for (r = 1; r++ < sampleCount;) {
                printf("|%s", COUNT[strand "-" r "-" l]);
            }
            if (strand == 1) printf(SEP);
        }
        printf("\t");
    }
    for (strand = 0; strand++ < 2;) {
        printf("%s", DEPTH[strand "-1"]);
        for (r = 1; r++ < sampleCount;) {
            printf("|%s", DEPTH[strand "-" r]);
        }
        if (strand == 1) printf(SEP);
    }
    delete DEPTH;
    printf("\n");
}'
