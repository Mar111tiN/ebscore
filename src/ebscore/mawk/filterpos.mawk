#!/usr/bin/env bash
# filterpos v1.0.0 (ebscore kumodev)
# Keeps the header line and every line whose Chr + Start (columns 1 and 2) are listed in a positions file.
# Replaces the sorted merge of tumor2matrix.mawk: works on any number of chromosomes, any order, and never
# drops a line because a neighbouring position is missing.
#
# USAGE: <tab-separated stream with a "Chr" header> | bash filterpos.mawk positions.tsv
#   positions.tsv: two tab-separated columns Chr, Pos; no header
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "<filterpos> Error: USAGE: ... | bash filterpos.mawk positions.tsv" >&2
    exit 1
fi

mawk -v posfile="$1" '
BEGIN {
    FS = OFS = "\t";
    while ((getline line < posfile) > 0) {
        split(line, f, "\t");
        POS[f[1] SUBSEP f[2]] = 1;
    }
    close(posfile);
}
NR == 1 && $1 == "Chr" { print; next; }
(($1 SUBSEP $2) in POS)
'
