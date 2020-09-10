from io import StringIO
import os
import sys
from subprocess import PIPE, run
from ebcore import fit_bb, get_count_df, bb_value, fisher_combination


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


def matrix2AB(row, config=None):
    count_df = get_count_df(row['PON:ALT=Depth'], config["count_dict"])
    return fit_bb(count_df, config['fit_pen'])


######### EB2AB ######################

def AB2EBscore(row, EBconfig):
    '''
    takes a df containing an AB column of shape "A+|A--B+|B-""
    '''
    # retrieve the data from the row
    target_s, AB_list = retrieveABdata(row, EBconfig)
    
    # get the p_values for each strand
    p_values = {}
    p_values['+'] = bb_pvalue(get_obs_df(target_s, ['depth+', 'alt+']), params_split[:2])
    p_values['-'] = bb_pvalue(get_obs_df(target_s, ['depth-', 'alt-']), params_split[2:])
    # combine p_values with fisher combination
    EB_p = fisher_combination(p_values)
    if EB_p < 1e-60:
        return 60
    if EB_p > 1.0 - 1e-10:
        return 0
    return -round(math.log10(EB_p), 3)