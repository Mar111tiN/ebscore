#!/bin/sh

# takes a position as argument and moves the data from that position to the tumor column
mawk '
# KEEP HEADER
NR == 1 {
    SEP="=";
    BASECOLS = 5;
    for (col=0;col++<BASECOLS;) printf("%s\t",$col);
    printf("Tumor:Alt%sDepth\tPON:ALT%sDepth\n", SEP,SEP)
    # get samplepos as bash arg with default 1 for first position
    samplepos='"${1:-1}"';
}
NR > 1 {
    # go through fields
    for (col=0;col++<BASECOLS;) printf("%s\t",$col);
    # get the data
    data = $(BASECOLS+1);
    split(data, DATAAT, SEP);
    # scan through Alt and Depth (AT=1|2)
    for (AT=0;AT++<2;) {
        split(DATAAT[AT],DATASTRAND,"-");
        # scan through pos|neg Strand (ST=1|2)
        for (ST=0;ST++<2;) {
            ponCount = split(DATASTRAND[ST], DATA, "|")
            # counter for ordered DATA array accessible via for loop
            count=0;
            for (i=0; i++<ponCount;) {
                if (i == samplepos) {
                    TUMOR[AT "-" ST] = DATA[i];
                } else {
                    count++;
                    DDATA[AT "-" ST "-" count] = DATA[i];
                }
            }            
        }
    }

    ###### OUTPUT ################
    # output tumor
    printf("%s-%s%s%s-%s\t",TUMOR["1-1"],TUMOR["1-2"],SEP,TUMOR["2-1"],TUMOR["2-2"])
    for (AT=0;AT++<2;) {
        for (ST=0;ST++<2;) {
            printf("%s",DDATA[AT "-" ST "-1"])
            for (i=1; i++<ponCount-1;) {
               printf("|%s",DDATA[AT "-" ST "-" i]) 
            }
            if (ST==1) printf("-");           
        }
        if (AT==1) printf(SEP);
    }
    printf("\n");
}'