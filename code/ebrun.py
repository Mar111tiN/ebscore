import os
import pandas as pd

from file2matrix import tumor2matrix
from tumor_matrix import tumor_matrix2AB_multi
from AB2EB import AB2EB_multi
from script_utils_EB import show_output


def run_ebscore(
    mut_file,
    tumor_bam="",
    pileup_file="",
    output_file="",
    pon_list="",
    chrom="",
    config=dict(
        debug=False,
        threads=8,
        use_cache=True,
        mawk_path="./shell",  # path relative to
        temp_dir="temp",
        pon_path=".",
        zero_path="./zero",
        genome_split="",
        MAPQ=20,  # params for pileup
        Q=25,  # params for pileup
        fit_pen=0.5,  # fitting penalty for beta_binomial parameter finding
        chunksize=5000,  # size of zero splits for zerocaching
        ZDFactor=13,  # how much complexity remains after flattening the tumor-zero lines
        min_zt=1000,  # minimum number of Tzero lines to bother with zero cache
    ),
):
    """
    master function for eb_computation
    """

    # ############## LOAD DATA ###############################
    tumor_file = tumor_bam if tumor_bam else pileup_file
    debug = config["debug"]

    show_output(f"Computing EBscore for chrom {chrom} on target {tumor_file}")

    # ############## bam2matrix #######################
    show_output("Retrieving tumor data and PONmatrix", time=False)
    matrix_df = tumor2matrix(
        mut_file,
        bam=tumor_bam,
        pileup=pileup_file,
        chrom=chrom,
        pon_list=pon_list,
        config=config,
    )
    if debug:
        matrix_file = f"{os.path.splitext(output_file)[0]}.EBmatrix"
        matrix_df.to_csv(matrix_file, sep="\t", index=False)

    # check if matrix_df is empty
    if matrix_df.empty:
        show_output("Something went wrong - created empty matrix file!")
        if debug:
            show_output(
                f"Created empty file {matrix_file}", color="warning", time=False
            )
        return

    # ########## ABparams ############################
    # # check for ABparams and compute if necessary
    if "AB" in matrix_df.columns:
        show_output(
            f"PONmatrix and AB computation on chrom {chrom} of {tumor_file} finished. Used ABcache!",
            color="success",
        )
        AB_df = matrix_df
    else:
        show_output(
            f"PONmatrix computation on chrom {chrom} of {tumor_file} finished. Starting AB computation..",
            color="success",
        )
        AB_df = tumor_matrix2AB_multi(matrix_df, config=config)
        show_output("Computing ABparams finished.", color="success")
        if debug:
            AB_file = f"{os.path.splitext(output_file)[0]}.AB"
            AB_df.to_csv(AB_file, sep="\t", index=False)

    # ######### AB2EB ################################
    show_output("Computing EBscore from ABparams", time=False)
    EB_df = AB2EB_multi(AB_df, config=config).drop(["AB", "Tumor"], axis=1)

    # ######### WRITE TO FILE ##############################################
    if output_file:
        EB_df.to_csv(output_file, sep="\t", index=False)
        show_output(
            f"Created EBscore for chrom {chrom} of {tumor_file} and written to {output_file}",
            color="success",
        )
    else:
        show_output(
            f"Computation of EBscore for chrom {chrom} of {tumor_file} finished!",
            color="success",
        )

    return EB_df
