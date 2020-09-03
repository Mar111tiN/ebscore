#!/bin/sh

# takes a pon_list as argument and returns the column numbers containing base info and read info
# usage in pipeline
# samtools mpileup | cut -f $(pon2cols ponlist.txt)
# returns Chr   Start   End     Read1   Read2   Read3....

mawk '
END {
    printf("%s,%s,%s",1,2,3);
    for (i=0; i++<NR;) {
        printf(",%s",2+3*i)
    }
}' $1