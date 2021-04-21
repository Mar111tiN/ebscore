import os
from multiprocessing import Pool
from functools import partial
import math
import numpy as np
import pandas as pd

from ebcore import fit_bb
from ebutils import get_pon_df
from script_utils import show_output, run_cmd
from matrix2AB import matrix2AB


def PON2matrix(pon_list, chrom, config={}):
    """
    generates matrix file from pon_list and writes it to pon_path/matrix/<chrom>.pon.gz
    """

    # PARAMS
    # mawk tool unwrapper
    def mawk(tool):
        return os.path.join(config["mawk_path"], f"{tool}.mawk")

    gsplit = config["genome_split"]
    pon_path = config["pon_path"]
    q = config["MAPQ"]
    Q = config["Q"]
    bed = config["bed_file"]

    matrix_path = os.path.join(pon_path, "matrix")
    if not os.path.isdir(matrix_path):
        os.mkdir(matrix_path)

    # check the temp_folder with default pon_path/temp
    temp_folder = config.get("temp_dir", os.path.join(pon_path, "temp"))
    if not os.path.isdir(temp_folder):
        os.mkdir(temp_folder)
    pon_list_full = os.path.join(config["temp_dir"], f"full_{pon_list}")

    # create the pon_list with full path
    pon_list = os.path.join(config["pon_path"], pon_list)

    if not os.path.isfile(pon_list_full):
        get_pon_df(pon_list, pon_path).to_csv(
            pon_list_full, sep="\t", index=False, header=False
        )

    pileup_cmd = f"samtools mpileup -Q {Q} -q {q} -l {bed} -f {gsplit}/{chrom}.fa -b {pon_list_full} -r {chrom}"

    pon_matrix_file = os.path.join(matrix_path, f"{chrom}.pon")
    cmd = f"{pileup_cmd} | {mawk('cleanpileup')} | {mawk('pile2count')} | tee {pon_matrix_file} |  gzip  > {pon_matrix_file}.gz"
    run_cmd(cmd)
    return pon_matrix_file


def stackPONmatrix(PON_df):
    """
    converts the PON matrix into stacked version for easy computation
    all string operations are performed here
    """
    # stack the letters
    df = (
        PON_df.set_index(["Chr", "Start", "Ref", "Depth"])
        .stack()
        .reset_index()
        .rename({"level_4": "Alt", 0: "tumor"}, axis=1)
    )
    # split left and right
    df[["DL", "DR"]] = df["Depth"].str.split("-", expand=True)
    df[["TL", "TR"]] = df["tumor"].str.split("-", expand=True)
    df["L"] = df["TL"] + "=" + df["DL"]
    df["R"] = df["TR"] + "=" + df["DR"]
    # stack left and right
    df = (
        df.loc[:, ["Chr", "Start", "Ref", "Alt", "L", "R"]]
        .set_index(["Chr", "Start", "Ref", "Alt"])
        .stack()
        .reset_index()
        .rename({"level_4": "strand", 0: "PON"}, axis=1)
    )
    df[["T", "D"]] = df["PON"].str.split("=", expand=True)
    return df.loc[:, ["Chr", "Start", "Ref", "Alt", "strand", "T", "D"]]


def unstack_PONAB(stack_df):
    """
    convert the stacked PON_df back into the normal PON_df for human readable output
    """

    un1 = (
        stack_df.loc[:, ["Chr", "Start", "Ref", "Alt", "strand", "AB"]]
        .set_index(["Chr", "Start", "Ref", "Alt", "strand"])
        .unstack("strand")["AB"]
    )
    un1["AB"] = un1["L"] + "=" + un1["R"]
    return (
        un1.loc[:, "AB"].unstack("Alt").reset_index(drop=False).sort_values(["Start"])
    )


def PONmatrix2AB_multi(
    pon_matrix_df,
    config={
        "fit_pen": 0.5,
        "threads": 8,
        "chunk_size": 50000,
        "zero_path": "zero",
        "zero_condense_factor": 13,
    },
):
    """
    converts a PONmatrix into a PONAB file
    """

    threads = config["threads"]

    stack_pool = Pool(threads)
    show_output(f"Starting AB conversion of PON matrix using {threads} cores")

    # stack the pon_matrix_df using multithreads if possible
    if threads > 1:
        stack_split = np.array_split(pon_matrix_df, threads)
        stacked_dfs = stack_pool.imap(stackPONmatrix, stack_split)
        stack_pool.close()
        stack_df = pd.concat(stacked_dfs).reset_index(drop=True)
    else:
        # only one core
        stack_df = stackPONmatrix(pon_matrix_df)
    show_output("PON matrix has been stacked")
    pon_len = len(stack_df.index)

    config["len"] = pon_len
    # retrieve the zero_string and pon_size from the panel of normals of first row
    config["pon_size"] = len(stack_df.loc[0, "D"].split("|"))
    # retrieve the zero_string from the panel of normals of first row
    # and store as state in config
    config["zero_string"] = "|".join(
        np.array([0] * len(stack_df.loc[0, "D"].split("|"))).astype(str)
    )

    AB_pool = Pool(threads)
    # minimal length of 200 lines
    # split_factor = min(math.ceil(len(pon_matrix_df.index) / 200), threads)
    split_factor = math.ceil(pon_len / config["chunk_size"])
    # split the matrix
    split = np.array_split(stack_df, split_factor)
    dfs = AB_pool.imap(partial(matrix2AB, config), split)
    AB_pool.close()
    # out_df contains AB params
    pon_AB_df = unstack_PONAB(pd.concat(dfs).reset_index(drop=True))
    show_output("PON matrix successfully converted!", color="success")
    return pon_AB_df
