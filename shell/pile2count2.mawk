#!/bin/sh
mawk '
NR==1 {
    # get the ref-base
    base = $3;
    BASECOLS = 3; # how many columns left of the pileup reads
    poncount = NF - BASECOLS;
    ###### QUERY ##############
    # get the letters to look for
    letter_len = split("Aa-Gg-Cc-Tt-Ii-Dd",LETTERS,"-");
    ###### HEADER ################

    printf("Chr\tStart\tRef");
    # loop through the letter patterns and print first string (upper case)
    for (l = 0; l++ < letter_len;) {
        printf("\t%s", substr(LETTERS[l],1,1));
    }
    printf("\tDepth\n");
}
{   ######### LINES #############   
    # loop through the reads
    for (r=0;r++<poncount;) {

        read = $(r + BASECOLS);
        read_len=length(read);
        # print(read, read_len);
        ######### --> COUNTS ############
        # start with +strand then -strand
        for (strand=0;strand++<2;) {
            # loop through the letters
            for (l = 0; l++ < letter_len;) {
                letters_found = gsub(substr(LETTERS[l],strand,1), "", read);
                # COUNT[1-1-1] ..+strand-read1-letter1(A)
                COUNT[strand "-" r "-" l] = letters_found;
                DEPTH[strand "-" r] += letters_found;
            }
        }
        # fill up the depth with "." count for +strand and read_len - +depth for -strand
        DEPTH["1-" r] += gsub("\.","",read);
        DEPTH["2-" r] += gsub(",","",read);
    }
    ######### OUTPUT #############
    # Basic output
    printf("%s\t%s\t%s\t",$1,$2,$3);
    # loop through the letters
    for (l = 0; l++<letter_len;) {
        for (strand=0; strand++<2;){
            printf("%s", COUNT[strand "-1-"l]);
            for (r=1;r++<poncount;){
                printf("|%s", COUNT[strand "-" r "-"l]);
            }
            if (strand==1) printf("-")
        }
        printf("\t")
    }
    # DEPTH OUTPUT
    for (strand=0; strand++<2;){
        printf("%s", DEPTH[strand "-1"]);
        for (r=1;r++<poncount;){
            printf("|%s", DEPTH[strand "-" r]);
        }
        if (strand==1) printf("-")
    }
    delete COUNT
    delete DEPTH
    printf("\n");
}'
