import os
from multiprocessing import Pool
from functools import partial
import math
import numpy as np
import pandas as pd

from ebcore import fit_bb
from ebutils import get_pon_df, get_zero_df
from script_utils import show_output, run_cmd
from ebconvert import matrix2AB


def PON2matrix(pon_list, chrom, EBconfig={}):
    """
    generates matrix file from pon_list
    """

    # PARAMS
    # mawk tool unwrapper
    def mawk(tool):
        return os.path.join(EBconfig["mawk_path"], f"{tool}.mawk")

    gsplit = EBconfig["genome_split"]
    pon_path = EBconfig["pon_path"]
    q = EBconfig["MAPQ"]
    Q = EBconfig["Q"]
    bed = EBconfig["bed_file"]

    matrix_path = os.path.join(pon_path, "matrix")
    if not os.path.isdir(matrix_path):
        os.mkdir(matrix_path)

    # check the temp_folder with default pon_path/temp
    temp_folder = EBconfig.get("temp_dir", os.path.join(pon_path, "temp"))
    if not os.path.isdir(temp_folder):
        os.mkdir(temp_folder)
    pon_list_full = os.path.join(EBconfig["temp_dir"], f"full_{pon_list}")

    # create the pon_list with full path
    pon_list = os.path.join(EBconfig["pon_path"], pon_list)

    if not os.path.isfile(pon_list_full):
        get_pon_df(pon_list, pon_path).to_csv(
            pon_list_full, sep="\t", index=False, header=False
        )

    pileup_cmd = f"samtools mpileup -Q {Q} -q {q} -l {bed} -f {gsplit}/{chrom}.fa -b {pon_list_full} -r {chrom}"

    pon_matrix_file = os.path.join(matrix_path, f"{chrom}.pon")
    cmd = f"{pileup_cmd} | {mawk('cleanpileup')} | {mawk('pile2count')} | tee {pon_matrix_file} |  gzip  > {pon_matrix_file}.gz"
    run_cmd(cmd)
    return pon_matrix_file


def tidyPONmatrix(PON_df):
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
    un1 = (
        stack_df.loc[:, ["Chr", "Start", "Ref", "Alt", "strand", "AB"]]
        .set_index(["Chr", "Start", "Ref", "Alt", "strand"])
        .unstack("strand")["AB"]
    )
    un1["AB"] = un1["L"] + "=" + un1["R"]
    return (
        un1.loc[:, "AB"]
        .unstack("Alt")
        .reset_index(drop=False)
        .sort_values(
            [
                "Chr",
                "Start",
            ]
        )
    )


def PONmatrix2AB_multi(
    pon_matrix_df,
    config={
        "fit_pen": 0.5,
        "threads": 8,
        "chunk_size": 50000,
        "zero_file": "zero.csv",
    },
):
    """
    converts a PONmatrix into a PONAB file
    """

    threads = config["threads"]

    pool = Pool(threads)

    # tidy the pon_matrix_df
    stack_df = tidyPONmatrix(pon_matrix_df)
    show_output(f"Stacked")

    pon_len = len(stack_df.index)

    config["len"] = pon_len
    # minimal length of 200 lines
    # split_factor = min(math.ceil(len(pon_matrix_df.index) / 200), threads)
    split_factor = math.ceil(pon_len / config["chunk_size"])
    # split the matrix
    split = np.array_split(stack_df, split_factor)
    dfs = pool.imap(partial(matrix2AB, config), split)
    pool.close()
    # out_df contains AB params
    pon_AB_df = unstack_PONAB(pd.concat(dfs))

    return pon_AB_df