#!/bin/sh

# mawk tool for extraction of base counts for clean pileup data
# 
# USAGE:
# samtools mpileup | cleanpileup | pile2count

# INPUT:
#  expects  a header from clean pileup 
# Chr   Start   Ref     Read1   Read2   Read3   Read4....
# output:
# Chr   Start   Ref     A G C T I D Depth
# Chr1  123124   A      sample|counts|in|+|strand-sample|counts|in|-|strand     depth|in|+|strand-depth|in|-|strand

mawk '
NR==1 {
    # get the ref-base
    base = $3;
    BASECOLS = 3; # how many columns left of the pileup reads
    sampleCount = NF - BASECOLS;
    ###### QUERY ##############
    # get the letters to look for
    letterCount = split("Aa-Gg-Cc-Tt-Ii-Dd",LETTERS,"-");
    ###### HEADER ################

    printf("Chr\tStart\tRef");
    # loop through the letter patterns and print first string (upper case)
    for (l = 0; l++ < letterCount;) {
        printf("\t%s", substr(LETTERS[l],1,1));
    }
    printf("\tDepth\n");
    # skip header
    next;
}
{   ######### LINES #############   
    # loop through the reads
    for (r=0;r++<sampleCount;) {

        read = $(r + BASECOLS);
        # print(reaLen);
        ######### --> COUNTS ############
        # start with +strand then -strand
        for (strand=0;strand++<2;) {
            # loop through the letters
            for (l=0;l++<letterCount;) {
                lettersFound = gsub(substr(LETTERS[l],strand,1), "", read);
                # COUNT[1-1-1] ..+strand-read1-letter1(A)
                COUNT[strand "-" r "-" l] = lettersFound;
                DEPTH[strand "-" r] += lettersFound;
            }
        }
        # fill up the depth with "." count for +strand aLen - +depth for -strand
        DEPTH["1-" r] += gsub("\.","",read);
        DEPTH["2-" r] += gsub(",","",read);
    }
    ######### OUTPUT #############
    # Basic output
    printf("%s\t%s\t%s\t",$1,$2,$3);
    # loop through the letters
    for (l = 0; l++<letterCount;) {
        for (strand=0; strand++<2;){
            printf("%s", COUNT[strand "-1-"l]);
            for (r=1;r++<sampleCount;){
                printf("|%s", COUNT[strand "-" r "-"l]);
            }
            if (strand==1) printf("-")
        }
        printf("\t")
    }
    # DEPTH OUTPUT
    for (strand=0; strand++<2;){
        printf("%s", DEPTH[strand "-1"]);
        for (r=1;r++<sampleCount;){
            printf("|%s", DEPTH[strand "-" r]);
        }
        if (strand==1) printf("-")
    }
    delete DEPTH
    printf("\n");
}'
