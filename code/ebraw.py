import pandas as pd
from io import StringIO
import os
import sys
from subprocess import PIPE, run


# here come all the functions involved in converting raw data (bam or pileup) to dataframes

def bam2matrix(bam_file, mut_file, chrom, pon_list, EBconfig):
    '''
    converts the bam_file and the Pon-bams into a matrix for EB computations
    '''
    #unwrap the shell tools
    cleanpileup = EBconfig["cleanpileup"]
    makeponlist = EBconfig["makeponlist"]
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
    pon_bam = f"{os.path.splitext(mut_file)[0]}_{chrom}.pon"
    run(f"{makeponlist} {bam_file} {pon_list} {pon_path} > {pon_bam}", check=True, shell=True)
    
    
    pileup_cmd = f"samtools mpileup -Q {Q} -q {q} -l {bed_file} -f {gsplit}/{chrom}.fa -b {pon_bam} -r {chrom}"
    process_cmd = f"cut -f $({pon2cols} {pon_list}) | {cleanpileup} | {pile2count} | {filterVar} {mut_file} {chrom} | {pon2tumor} 1"
    cmd = f"{pileup_cmd} | {process_cmd}"
    matrix_df = pd.read_csv(StringIO(
    run(cmd, stdout=PIPE, check=True, shell=True).stdout.decode('utf-8')), sep='\t')
    return matrix_df