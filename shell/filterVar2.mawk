#!/bin/sh

# filterEB v1.2
# takes a pile2count or ABcache file with Chr Start Ref A-datacols C-data G-data T-data I-data D-data
# and returns the columns matching the coords in mut.csv reduced to the specific Alt data
# mpileup | cut .. | cleanpileup | pile2count2 | filterVar <mut.csv> chrom
# mut.csv must have columns: Chr Start End Ref Alt

###################################################
####### PARAMS ####################################

####### ARGPARSE ##################
PARAMS=""
while (( "$#" )); do
    # allow for equal sign in long-format options
    [[ $1 == --*=* ]] && set -- "${1%%=*}" "${1#*=}" "${@:2}"
    case "$1" in
        # pileup input
        -p|-t|--pileup)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            tumorPileup=$2
            shift 2
        else
            echo "<filterEB> Error: tumor pileup file is missing\n[-p|-t|--pileup <path_to_tumor_pileup>]" >&2;
            exit 1;
        fi
        ;;
        # PON input
        -P|--pon-matrix)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            PONfile=$2;
            shift 2;
        else
            echo "<filterEB> Error: PON file is missing\n[-P|--pon-matrix <path_to_pon-matrix>]" >&2
            exit 1;
        fi
        ;;
        # pileup input
        -x|--exclude|--pon-exclude)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            PONexclude=$2
            shift 2
        else
            echo "<filterEB> Error: position for PON exclude is missing\n[-x|--exclude|--pon-exclude <INT>]" >&2
            exit 1
        fi
        ;;
        -*|--*=) # unsupported flags
        echo "<filterEB> Error: Unsupported flag $1" >&2
        exit 1
        ;;
        *) # preserve positional arguments
        PARAMS="$PARAMS $1"
        shift
        ;;
    esac
done
# # set positional arguments in their proper place
eval set -- "$PARAMS"

####### DEFAULT ARGS ####################
tumorPileup=${useExonicCoords-0};
filterChrom=${filterChrom-""};
mutFile=$1;

PONfile=${PONfile-"none"};
PONexclude=${PONexclude-0};

#############################################
###### START ################################
cat $mutFile - | mawk '
NR == 1 { # @HEADER of mutFile
    ### ARGS #################
    # set the separator between Alt and Depth
    SEP = "=";
    PONfile = "'$PONfile'";
    PONexclude='$PONexclude';
    chrom = "'$2'";
    # set STATE variables
    readMut = 1;
    readData = 0
    step = 1; # the position counter

    # open PON getline stream (Gstream) and skip first line
    if (PONfile != "none") {
        getline < PONfile;
    }

    # if mutfile contains header, skip to next line
    # else, start processing
    if ($0 ~ "Chr\tStart\t") next;
}
# read the mut_file
readMut { # read until the header of the next file
    if ($0 !~ "Chr\tStart\t") { #@ not header of stream data
        # filter for the right chromosome
        if ($1 != chrom) next;

        #### DEBUG ########
        # print($1, $2, $3, $4, readMut);
        ###################

        # adjust coords for deletion
        # for deletions, we need the pileup BEFORE the base
        ($5 == "-") ? pos=$2-1 : pos=$2;
        # store the current pos as start of search
        if (step == 1) {
            # STATE 
            currentPOS = pos;
        }
        # store ordered coords in POS array
        POS[step++] = pos;
        # store coords and respective Alt data in POSALT
        # get REF from 
        # POSALT[1242352] = "A"
        POSEND[pos] = $3;
        POSREF[pos] = $4;
        POSALT[pos] = $5;
        next;
    } else {
        # save last position in mutfile
        # in order to stop searching
        lastPos = POS[step-1];
        # reset step
        step=1;
        readMut = 0;
        writeHeader = 1;
    }
}
# HEADER
writeHeader { #@stream header
    writeHeader = 0;
    readData = 1;
    # check if there is a header
    if ($0 ~ "Chr") {
        # get the col numbers for the vars from the header into COL array
        # COL["A"] = 4 etc...
        for (col=0; col++<NF;) {
        COL[$col] = col;
        }
        printf("Chr\tStart\tEnd\tRef\tAlt\tTumor:Alt=Depth\tPON:Alt=Depth\n");
        ########
        printf("stdinBase\tstdinDepth\tPONData\tPONdepth\n");
        ########
        next;
    } else {
        print("<filterVar> No header detected") > "/dev/stderr";
    }
}
readData { #@ stream data
    pos=$2;
    # data is between last and currentPOS
    if (pos < currentPOS) next;
    # data for currentPOS is missing and moved past currentPOS
    if (pos > currentPOS) {
        # moved past the last POS
        if (pos > lastPos) exit;
        # go to next position in POS ARRAY
        currentPOS = POS[++step];
        next;
    }
    # @found stream data matching currentPOS:
    # print data at matching positions 
    # get the right stream column depending on POSALT
    # get the Ref and Alt from the arrays
    ref = POSREF[pos];
    alt = POSALT[pos];

    # modify altBase for Indels
    if (ref == "-") {
        altBase = "I";
    } else {
        if (alt == "-") {
            # for deletion, start has to be increased
            pos++;
            altBase = "D";
        } else {
            altBase = alt;
        }
    } 

    ### ########### BASE OUTPUT  ################
    printf("%s\t%s\t%s\t%s\t%s\t",$1,pos,POSEND[pos],ref,alt);
    
    # store the streamData and Depth
    streamData = $(COL[altBase]);
    # get the Depth from last column
    streamDepth = $NF;

    ############ get the PON cache data via Gstream#############
    if (PONfile != "none") {
        while ((getline < PONfile) > 0) { # Gstream

            if ($2 == currentPOS) { # found position in Gstream 
                # print the tumor data local stream
                printf("%s%s%s\t", streamData, SEP, streamDepth);
                PONdata = $(COL[altBase]);
                PONdepth = $NF; 
                if (PONexclude) {  #  
                    # split the PONData and PONdepth into SData array and retrieve the tumor data
                    split(PONdata, pdata, "-");
                    split(PONdepth, pdepth, "-");
                    # assume same structure of data and depth
                    for (strand in pdata) {
                        ponCount = split(pdata[strand], PDATASTRAND, "|");
                        split(pdepth[strand], PDEPTHSTRAND, "|");
                        i=1;
                        for (pon=0; pon++< ponCount;) {
                            # transfer the entire data via SDATASTRAND into SDATA
                            # same thing for depth
                            if (pon != PONexclude) {
                                PDATA[strand "-" i] = PDATASTRAND[pon];
                                PDEPTH[strand "-" i] = PDEPTHSTRAND[pon];
                                i++
                            }
                        }
                    }

                    ####### DEBUG ###########
                    # for (p in PDEPTH) {
                    #     print(p, PDEPTH[p]);
                    # }
                    ####### DEBUG ###########

                    ##### OUTPUT ################
                    ponData="";
                    ponDepth="";
                    for (strand=0; strand++<2;) {
                        # start with 1 because array[1]  is already printed
                        for (pon=0;pon++<ponCount-1;) {
                            ponData = ponData PDATA[strand "-" pon];
                            ponDepth = ponDepth PDEPTH[strand "-" pon];
                            if (pon != ponCount-1) {
                                ponData = ponData "|";
                                ponDepth = ponDepth "|";
                            } else {
                                if (strand == 1) {
                                    ponData = ponData "-";
                                    ponDepth = ponDepth "-";                                   
                                }
                            }
                        }
                    }
                    printf("%s%s%s\n", ponData, SEP, ponDepth);
                } else {
                    printf("%s%s%s\n", PONdata, SEP, PONdepth);
                }
                break;
            }
            # if the data is not in the PON file, write "NAN"
            if ($2 > currentPOS) {
                PONdata = "NAN";
                PONdepth = "NAN";  
                break;              
            };
        }
    } else {
    # PON and tumor combined
        # split the streamData into SData array and retrieve the tumor data
        split(streamData, sdata, "-");
        split(streamDepth, sdepth, "-");
        # assume same structure of data and depth
        for (strand in sdata) {
            ponCount = split(sdata[strand], SDATASTRAND, "|");
            split(sdepth[strand], SDEPTHSTRAND, "|");
            for (pon in SDATASTRAND) {
                # transfer the entire data via SDATASTRAND into SDATA
                # same thing for depth
                SDATA[strand "-" pon] = SDATASTRAND[pon];
                SDEPTH[strand "-" pon] = SDEPTHSTRAND[pon];
            }
        }
        ####### DEGUG #######
        # for (s in SDATA) {
        #     print(s, SDATA[s]);
        # }
        ####### DEGUG #######

        # print the tumor data as first elements of all arrays
        printf("%s-%s%s%s-%s\t", SDATA["1-1"], SDATA["2-1"], SEP, SDEPTH["1-1"], SDEPTH["2-1"]);

        # output the rest
        ponData="";
        ponDepth="";
        for (strand=0; strand++<2;) {
            # start with 1 because array[1]  is already printed
            for (pon=1;pon++<ponCount;) {
                ponData = ponData SDATA[strand "-" pon];
                ponDepth = ponDepth SDEPTH[strand "-" pon];
                if (pon != ponCount) {
                    ponData = ponData "|";
                    ponDepth = ponDepth "|";
                } else {
                    if (strand == 1) {
                        ponData = ponData "-";
                        ponDepth = ponDepth "-";
                    }
                }
            }

        }
        printf("%s%s%s\n", ponData, SEP, ponDepth);
    }

    ############# get the pileup data from pileup file ###


    printf("%s\t%s\t%s\t%s\n",streamData, streamDepth, PONdata, PONdepth);
    if (currentPOS == lastPos) exit;
    # bump currentPos to next mut position
    currentPOS = POS[++step];
    next;
    # stop if end is reached
}'