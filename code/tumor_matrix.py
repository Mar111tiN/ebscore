from multiprocessing import Pool
from functools import partial
import math
import numpy as np
import pandas as pd

from script_utils_EB import show_output, run_cmd
from zerocache import flatten_df, zero2AB, update_zero_file, get_zerostring
from matrix2AB import matrix2AB

# separator between strands
STRANDSEP = "="
# separator between Alt and Depth
ADSEP = "<"

# ######## compute AB directly from tumor-matrix file
def stack_tumor_matrix(config, df):
    """
    create a stacked matrix_df for direct matrix2AB computation
    """

    df = (
        df.set_index(["Chr", "Start", "End", "Ref", "Alt", "Tumor"])
        .stack()
        .reset_index()
        .rename({"level_6": "strand", 0: "PON"}, axis=1)
    )
    df[["T", "D"]] = df["PON"].str.split(ADSEP, expand=True)

    df.loc[:, ["Chr", "Start", "End", "Ref", "Alt", "Tumor", "strand", "T", "D"]]
    df = flatten_df(df, ZDfactor=config["ZDfactor"])
    return df


def unstack_AB(stackAB_df):
    """
    converts the AB_matrix_df into unstacked version for further processing
    """
    unstack_df = (
        stackAB_df.loc[
            :, ["Chr", "Start", "End", "Ref", "Alt", "Tumor", "strand", "AB"]
        ]
        .set_index(["Chr", "Start", "End", "Ref", "Alt", "Tumor", "strand"])
        .unstack("strand")["AB"]
    )
    unstack_df["AB"] = unstack_df["PON+"] + STRANDSEP + unstack_df["PON-"]
    unstack_df = unstack_df.reset_index(drop=False).sort_values(["Start", "End"])
    return unstack_df.loc[:, ["Chr", "Start", "End", "Ref", "Alt", "Tumor", "AB"]]


def tumor_matrix2AB_multi(
    tumor_matrix_df,
    config={
        "fit_pen": 0.5,
        "threads": 8,
        "zero_path": "zero",
        "ZDfactor": 13,
        "min_zt": 1000,
        "chunksize": 20000
    },
):
    """
    converts a matrix_df into a AB file using multicore
    """

    threads = config["threads"]
    # stack the tumor_matrix from generator using multiprocessing
    stack_pool = Pool(threads)
    show_output(f"Stacking tumor matrix using {threads} cores.")
    stack_split = np.array_split(tumor_matrix_df, threads)
    dfs = stack_pool.imap(partial(stack_tumor_matrix, config), stack_split)
    stack_df = pd.concat(dfs)
    stack_pool.close()
    show_output("tumor matrix has been stacked")

    # retrieve the zerostring and ponsize from the panel of normals of first row
    # store in configs
    config["zerostring"], config["ponsize"] = get_zerostring(stack_df)
    # ######### reading AB from zerocache into AB_df
    stack_df, AB_df, useZero = zero2AB(stack_df, config)

    # case all AB have been retrieved from zero2AB
    if stack_df.empty:
        AB_df = unstack_AB(AB_df).reset_index(drop=True).sort_values(["Start"])
    else:

        # compute the AB_df for the stack_df
        AB_pool = Pool(threads)
        show_output(f"Starting AB conversion of matrix using {threads} cores")

        # minimal length of 200 lines
        # split_factor = min(math.ceil(len(pon_matrix_df.index) / 200), threads)
        split_factor = max(1, min(int(len(stack_df.index) / 200), threads))
        # split the matrix
        split = np.array_split(stack_df, split_factor)
        del stack_df
        AB_dfs = AB_pool.imap(partial(matrix2AB, config), split)
        AB_pool.close()

        # collect AB_dfs + get all the new zeros and write to file

        new_AB_df = pd.concat(AB_dfs)
        # check, if zero has been used (or zero_df was smaller then min_zt)
        if useZero:
            _ = update_zero_file(new_AB_df, config=config)

        AB_df = unstack_AB(
            pd.concat([AB_df, new_AB_df]).reset_index(drop=True)
        ).sort_values(["Start"])

    show_output("Matrix successfully converted!", color="success")
    # bring back the PONmatrix
    AB_df = AB_df.merge(tumor_matrix_df)

    # AB_df.loc[:, 'PON'] = AB_df['PON+'].str.split("=").str[0] + "-" + AB_df['PON-'].str.split("=").str[0] + "=" + AB_df['PON+'].str.split("=").str[1] + "-" + AB_df['PON-'].str.split("=").str[1]
    return AB_df.loc[:, ["Chr", "Start", "End", "Ref", "Alt", "Tumor", "PON+", "PON-", "AB"]]
