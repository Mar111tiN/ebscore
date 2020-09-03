#!/bin/sh

# takes a pile2count or ABcache file with Chr Start Ref A-datacols C-data G-data T-data I-data D-data
# and returns the columns matching the coords in mut.csv reduced to the specific Alt data
# mpileup | cut .. | cleanpileup | pile2count2 | filterVar <mut.csv> chrom [showDepth=0|1]
# showDepth=0 -> do not print the last column
# showDepth=1 -> default: show last column (usually containing depth info)
# mut.csv must have columns: Chr Start End Ref Alt


echo "STOP" > stop.file;
cat $1 stop.file - | mawk '
BEGIN {
  # INIT
  # set the separator between Alt and Depth
  SEP = "=";
  # read mutation file is active
  readMut = 1;
  readData = 0;
  chrom = "'$2'";
  showDepth = "'${3:-1}'";
  step = 1; # the position counter 
}
# read the mut_file
readMut {
  if ($0 !~ "STOP") {
    # filter for the right chromosome
    if ($1 !~ chrom) next;
    # adjust coords for deletion
    ($5 == "-") ? pos=$2-1 : pos=$2;
    # store the current pos
    if (step == 1) {
      currentPos = pos;
    }
    # store ordered coords in POS array
    POS[step++] = pos;
    # store coords and respective Ref and Alt data in POSALT
    # POSALT[1242352] = "A"
    POSEND[pos] = $3;
    POSREF[pos] = $4;
    POSALT[pos] = $5;
  }
  else {
    # save last position in mutfile
    lastPos = POS[step-1];
    # reset step
    step=1;
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
    printf("Chr\tStart\tEnd\tRef\tAlt\tPONAlt");
    if (showDepth==1) {
      printf("%sDepth\n", SEP);
    } else {
      printf("\n");
    }
    next;
  } else {
    print("<filterVar> No header detected") > "/dev/stderr";
  }
}
readData {
  pos=$2;
  if (pos < currentPos) next;
  # print data at matching positions 
  if (pos in POSALT) {
    # get the right column depending on the ALT column
    # get the Ref and Alt from the arrays
    ref = POSREF[pos];
    alt = POSALT[pos];
    # get the right data column
    # for deletion, start has to be increased again
    if (alt == "-") {
      base = "D";
      pos++;
    } else {
      if (alt == "-") {
        base = "D"
      } else {
        base = alt;
      }
    }
    
    data = $(COL[base]);
    printf("%s\t%s\t%s\t%s\t%s\t%s",$1,pos,POSEND[pos],ref,alt,data);
    if (showDepth==1) {
      printf("%s%s\n",SEP,$NF);
    } else {
      printf("\n");
    }
    # bump currentPos to next mut position
    currentPos = POS[step++];
    next;
  }
  # stop if end is reached
  if (pos > lastPos) exit;
}'