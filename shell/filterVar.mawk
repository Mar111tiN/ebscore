#!/bin/sh

# takes a pile2count or ABcache file with Chr Start Ref A-datacols C-data G-data T-data I-data D-data
# and returns the columns matching the coords in mut.csv reduced to the specific Alt data
# mpileup | cut .. | cleanpileup | pile2count2 | filterVar <mut.csv> chrom
# mut.csv must have columns: Chr Start End Ref Alt
echo "STOP" > stop.file;
cat $1 stop.file - | mawk '
BEGIN {
  # INIT
  # read mutation file is active
  readMut = 1;
  readData = 0;
  chrom = "'$2'";
  step = 1; # the position counter 
}
# read the mut_file
readMut {
  if ($0 !~ "STOP") {
    # filter for the right chromosome
    if ($1 !~ chrom) next;
    # # store the current pos
    # if (step == 1) {
    #   current = $2;
    # }
    # # store ordered coords in POS array
    # POS[step++] = $2;
    # store coords and respective Alt data in POSALT
    # POSALT[1242352] = "A"
    POSALT[$2] = $5;
  }
  else {
    readMut = 0;
    writeHeader = 1;
  }
  next;
}
# HEADER
writeHeader {
  writeHeader = 0;
  readData = 1;
  # check if there is a header
  if ($0 ~ "Chr") {
    # get the col numbers for the vars from the header into COL array
    for (col=0; col++<NF;) {
      COL[$col] = col;
    }
    printf("Chr\tStart\tRef\tAlt\tData\tDepth\n");
    next;
  } else {
    print("<filterVar> No header detected") > "/dev/stderr";
  }
}
readData {
  pos=$2;
  if (pos in POSALT) {
    printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,POSALT[pos],$(COL[POSALT[pos]]), $NF);
  }
}'