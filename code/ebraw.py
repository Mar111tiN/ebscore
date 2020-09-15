import pandas as pd
from io import StringIO
import os
import sys
from subprocess import PIPE, run
from ebutils import get_pon
from script_utils import show_output

# here come all the functions involved in converting raw data (bam or pileup) to dataframes


def bam2matrix(bam_file, mut_file, chrom, pon_list, EBconfig):
    '''
    converts the bam_file and the Pon-bams into a matrix for EB computations
    '''
    # unwrap the shell tools
    cleanpileup = EBconfig["cleanpileup"]
    csv2bed = EBconfig["csv2bed"]
    pon2cols = EBconfig["pon2cols"]
    pile2count = EBconfig["pile2count"]
    filterVar = EBconfig["filterVar"]
    pon2tumor = EBconfig["pon2tumor"]
    gsplit = EBconfig["genome_split"]
    pon_path = EBconfig["pon_path"]
    q = EBconfig["MAPQ"]
    Q = EBconfig["Q"]

    # create the bed file
    bed_file = f"{os.path.splitext(mut_file)[0]}_{chrom}.bed"
    run(f"{csv2bed} {mut_file} > {bed_file}", check=True, shell=True)

    # create the pon_list
    pon_list = os.path.join(EBconfig['pon_path'], pon_list)
    # make filename for pon_bam
    pon_bam = f"{os.path.splitext(mut_file)[0]}_{chrom}.pon"

    pon_df, pos, removed_ponbam = get_pon(
        bam_file, pon_list, pon_path=pon_path, prepend_bam=True)
    # output if bam is included in pon_list
    if pos:
        show_output(
            f"{removed_ponbam} matches target bam {bam_file} and is removed from pon_list. If not desired, change names in PoN", color="warning")

    pon_df.to_csv(pon_bam, index=False, header=None)

    pileup_cmd = f"samtools mpileup -Q {Q} -q {q} -l {bed_file} -f {gsplit}/{chrom}.fa -b {pon_bam} -r {chrom}"
    process_cmd = f"cut -f $({pon2cols} {pon_list}) | {cleanpileup} | {pile2count} | {filterVar} {mut_file} {chrom} | {pon2tumor} 1"
    cmd = f"{pileup_cmd} | {process_cmd}"
    matrix_df = pd.read_csv(StringIO(
        run(cmd, stdout=PIPE, check=True, shell=True).stdout.decode('utf-8')), sep='\t')
    return matrix_df
