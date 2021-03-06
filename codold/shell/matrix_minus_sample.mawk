#!/bin/sh

# takes a position as argument and moves the data from that position to the tumor column
mawk '
# KEEP HEADER
NR == 1 {
    print $0;
    # get samplepos as bash arg with default 1 for first position
    samplepos='"${1:-1}"';
}
NR > 1 {
    # go through fields
    printf("%s\t%s", $1,$2);
    for (i = 2; i++ < NF;) {
        # split data by "|"
        split($i,array,"|");
        count=0;
        for (j = 0; j++ < length(array);) {
            # get the position argument
            if (j != samplepos) {
                count++;
                newarray[count]=array[j];
            } else {
                target=array[j];
            }
        }
        #output
        printf("\t%s", target);
        for (r=0; r++ < length(newarray);) {
            printf("|%s", newarray[r]);
        }
    }
    printf("\n");
}'