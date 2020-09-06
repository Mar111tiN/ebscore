#!/bin/sh

echo $1 | cat - $2 | mawk  '
BEGIN {
    FS="/";
    pon="'$3'";
}
# print tumor bam and store basename (before _A) in pattern
NR==1 {
    string=$NF;
    split(string,pat,"_");
    pattern="/" pat[1];
    print $0;
    next;
}
$0 !~ pattern {
    printf("%s/%s\n",pon,$0);
}'