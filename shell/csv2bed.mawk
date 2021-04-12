#!/bin/sh

# converts a bed file with Chr  Start   End into bed file
# moves coords of deletion down for catching the deletion in mpileup format

mawk '
BEGIN {
    OFS="\t";
}
$1 == "'$2'" {
    start = $2 - 1;
    if ($5 == "-") {
        start = start -1;
    }
    print($1,start,start+1);
}' $1
