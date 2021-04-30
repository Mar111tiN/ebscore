#!/bin/sh

# filterEB v1.3

# powerhorse for EBraw data output:
# collects tumor (and optionally PONbam pileup) and combines it with PONcache and ABcache data
# USAGE: 
# pileup data comes from tumor bam only or a tumor_PONlist with tumor and PONfiles
# samtools mpileup -f $HG38 -q 20 -Q 25 -l $BED -r chr? [-b tumor_PONlist|tumor_bam] | 
# cleanpileup | pile2count | tumor2matrix [options] mutfile chrom
####### OPTIONS ####################################         
# [     -P | --pon-matrix]          =stream|<path_to_PONcache>|none     path to (opt. gzipped) PONcache_file, default is from input stream  ]
# [     -x | --pon-exclude          <INT=0>                             removes PONdata at INT position in case of tumor_bam matching PON   ]
# [     -A | --ABcache              =none|<path_to_PONcache>            add AB params from (opt. gzipped) ABcache to output stream          ]
# [     -m | --include-missing      FLAG                                if mut_positions missing in tumor should be printed                 ]
###################################################

####### ARGPARSE ##################
PARAMS=""
while (( "$#" )); do
    # allow for equal sign in long-format options
    [[ $1 == --*=* ]] && set -- "${1%%=*}" "${1#*=}" "${@:2}"
    case "$1" in
        -m|--include-missing)
        printMissing=1;
        shift
        ;;
        # PON input (can come from stream in case of combined pileup or pon_matrix or not at all [option None])
        -P|--pon-matrix)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            PONfile=$2;
            shift 2;
        else
            echo "<filterEB> Error: PON file is missing\n[-P|--pon-matrix <path_to_PONcache matrix>| None | stream=default]" >&2
            exit 1;
        fi
        ;;
        # exclude sample from PON at given position
        -x|--exclude|--pon-exclude|--PON-exclude)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            PONexclude=$2
            shift 2
        else
            echo "<filterEB> Error: position for PON exclude is missing\n[-x|--exclude|--pon-exclude <INT>]" >&2
            exit 1
        fi
        ;;
        # use ABcache
        -A|-a|--ABcache)
        if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
            ABfile=$2
            shift 2
        else
            echo "<filterEB> Error: ABcache file is missing\n[-p|-t|--pileup <path_to_ABcache file>]" >&2;
            exit 1;
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







#############################################
###### START line80 ################################
# read the mutation file followed by the stream data (pileup of tumor/tumor+PON)
cat $1 - | mawk '
NR == 1 { # @HEADER of mutFile

    ### ARGS/GLOBALS #################
    # set the separator between Alt and Depth
    SEP = "=";
    chrom = "'$2'";
    printMissing = '${printMissing-0}';
    ### PONmatrix #####################
    # PON matrix comes either from stream or from file
    PONfile = "'${PONfile-"stream"}'";
    PONfile = (tolower(PONfile) == "none") ? "none" : PONfile;
    PONfile = (tolower(PONfile) == "stream") ? "stream" : PONfile;
    PONexclude='${PONexclude-0}';
    # getting PONfile
    if (PONfile == "none") {
        printf("<tumor2matrix> PON output is disabled.\n") > "/dev/stderr";
    } else {
        if (PONfile != "stream") {
            if (PONfile ~ /.gz$/) {
                PONcmd = "cat " PONfile " 2>/dev/null | gunzip ";
            } else {
                PONcmd = "cat " PONfile " 2>/dev/null ";
            } 
            # open PON getline stream (Gstream) and skip first line
            # getline == 0 if file not found --> PON > "none
            if ((PONcmd | getline) == 0) {
                printf("<tumor2matrix> PONcache file %s not found!\n", PONfile) > "/dev/stderr";
                PONfile = "none";
            } else {                   
                printf("<tumor2matrix> Using PONcache %s.\n", PONfile) > "/dev/stderr";
                getPON=1; # stateVariable to activate getline for PONfile
            }
        }
    }

    ### ABcache #####################
    ABfile = "'${ABfile-"none"}'";
    ABfile = (tolower(ABfile) == "none") ? "none" : ABfile;
    if (ABfile != "none") {
        if (PONexclude) {
            # if normal matching to sample is in normal, ABcache cannot be used!
            printf("<tumor2matrix> ABcache cannot be used if -x is set!\n", ABfile) > "/dev/stderr";
            ABfile = "none";
        } else {
            if (ABfile ~ /.gz$/) {
                ABcmd = "cat " ABfile " 2>/dev/null | gunzip ";
            } else {
                ABcmd = "cat " ABfile " 2>/dev/null ";
            } 
            # open AB getline stream (Gstream) and skip first line
            # getline == 0 if file not found --> AB > "none
            if ((ABcmd | getline) == 0) {
                printf("<tumor2matrix> ABcache file %s not found!\n", ABfile) > "/dev/stderr";
                ABfile = "none";
            } else {                   
                printf("<tumor2matrix> Using ABcache %s.\n", ABfile) > "/dev/stderr";
            }
        }
    }

    # set STATE variables
    readMut = 1;
    readData = 0
    step = 1; # the position counter

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
        # PRINT HEADER
        printf("Chr\tStart\tEnd\tRef\tAlt\tTumor");
        if (PONfile != "none") printf("\tPON+\tPON-");
        if (ABfile != "none") printf("\tAB");
        printf("\n");
        ########
        # printf("Chr\tStart\tEnd\tRef\tAlt\tTumor:Alt=Depth\tPON:Alt=Depth\n");
        # printf("stdinBase\tstdinDepth\tPONData\tPONdepth\n");
        ########
        next;
    } else {
        print("<tumor2matrix> No header detected") > "/dev/stderr";
    }
}
readData { #@ stream data
    pos=$2;
    # data is between last and currentPOS
    if (pos < currentPOS) next;

    ####### MISSING DATA
    while (pos > currentPOS) { # streamdata for currentPOS is missing
        # get the Ref and Alt from missing position
        if (printMissing) {
            ref = POSREF[currentPOS];
            alt = POSALT[currentPOS];
            printf("%s\t%s\t%s\t%s\t%s\n",$1,currentPOS,POSEND[currentPOS],ref,alt);
        }
        # moved past the last POS
        if (pos > lastPos) exit;
        # go to next position in POS ARRAY
        currentPOS = POS[++step];
        # data is again between last and currentPOS
        if (pos < currentPOS) next;
    }
    # get the Ref and Alt from the arrays
    ref = POSREF[pos];
    alt = POSALT[pos];
    end = POSEND[pos];
    # @found stream data matching currentPOS:
    # print data at matching positions 
    # get the right stream column depending on POSALT
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
    # print only after pos has been raised for deletions
    printf("%s\t%s\t%s\t%s\t%s\t",$1,pos,end,ref,alt);

    ####### STREAMDATA
    # store the streamData and Depth
    streamData = $(COL[altBase]);
    # get the Depth from last column
    streamDepth = $NF;

    ####### TUMOR AND PON ##############
    if (PONfile == "none") {
        printf("%s%s%s", streamData, SEP, streamDepth);
    } else {
        if (getPON == 1) { ############ PONcache from Gstream #############
            while ((PONcmd | getline) > 0) { # Gstream
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

                        for (strand=0; strand++<2;) {
                            ponData="";
                            ponDepth="";
                            # start with 1 because array[1]  is already printed
                            for (pon=0;pon++<ponCount-1;) {
                                ponData = ponData PDATA[strand "-" pon];
                                ponDepth = ponDepth PDEPTH[strand "-" pon];
                                if (pon != ponCount-1) {
                                    ponData = ponData "|";
                                    ponDepth = ponDepth "|";
                                } else {
                                    STRAND[strand] = ponData SEP ponDepth;
                                }
                            }
                        }
                        printf("%s\t%s", STRAND[1], STRAND[2]);
                    } else {
                        split(PONdata, A, "-");
                        split(PONdepth, D, "-");
                        printf("%s%s%s\t%s%s%s", A[1], SEP, D[1], A[2], SEP, D[2]);
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
        } else { # PON and tumor are combined in STREAM
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
            for (strand=0; strand++<2;) {
                ponData="";
                ponDepth="";
                # start with 1 because array[1]  is already printed
                for (pon=1;pon++<ponCount;) {
                    ponData = ponData SDATA[strand "-" pon];
                    ponDepth = ponDepth SDEPTH[strand "-" pon];
                    if (pon != ponCount) {
                        ponData = ponData "|";
                        ponDepth = ponDepth "|";
                    } else {
                        STRAND[strand] = ponData SEP ponDepth;
                    }
                }

            }
            printf("%s\t%s", STRAND[1], STRAND[2]);
        }
    }

    ####### ABcache ##############
    if (ABfile != "none") {
        while ((ABcmd | getline) > 0) { # check for end of file
            if ($2 == currentPOS) { # found position in Gstream 
                ABdata = $(COL[altBase]);
                printf("\t%s", ABdata);
                break;
            }
            # if the data is not in the PON file, write "NAN"
            if ($2 > currentPOS) {
                PONdata = "NAN";
                PONdepth = "NAN";  
                break;              
            };
        }
    }

    ########## FINAL ############
    printf("\n");
    # stop if end is reached
    if (currentPOS == lastPos) exit;
    # bump currentPos to next mut position
    currentPOS = POS[++step];

    ###### DEBUG ##############
    # printf("%s\t%s\t%s\t%s\n",streamData, streamDepth, PONdata, PONdepth);
    ###### DEBUG ##############
}'