import pandas as pd
import numpy as np
import math

from multiprocessing import Pool
from functools import partial

from script_utils_EB import show_output
from ebcore import bb_pvalue, fisher_combination


def get_obs_array(t_pair):
    """
    turns an array [5,9] into array of observation pairs
       [[9, 5],
        [9, 6],
        [9, 7],
        [9, 8],
        [9, 9]]
    """
    return np.array([[t_pair[1], s] for s in range(t_pair[0], t_pair[1] + 1)])


######### EB2AB ######################
def AB2EB_row(row):
    """
    takes a df containing an AB column of shape "A+|A=B+|B-""
    """

    # get AB params from AB string
    AB_params = np.array([p.split("|") for p in row["AB"].split("-")]).astype(float)

    # get tumor matrix from Tumor string
    t_matrix = np.transpose([p.split("-") for p in row["Tumor"].split("=")]).astype(int)

    # convert tumor matrix in observation arrays
    obs_arrays = [get_obs_array(t_pair) for t_pair in t_matrix]

    # get the p_values for each strand
    p_values = [
        bb_pvalue(obs_array, AB) for (obs_array, AB) in zip(obs_arrays, AB_params)
    ]

    # combine p_values with fisher combination
    EB_p = fisher_combination(p_values)
    if EB_p < 1e-60:
        return 60
    if EB_p > 1.0 - 1e-10:
        return 0
    return -round(math.log10(EB_p), 3)


def AB2EB(AB_df):
    """
    single-threaded converter of AB params and Tumor to EBscore
    """
    show_output(
        f"Computing EBscore for {len(AB_df.index)} lines!",
        multi=True,
    )
    AB_df["EB"] = AB_df.apply(AB2EB_row, axis=1)
    show_output(
        f"EB computation finished",
        multi=True,
        color="success",
    )
    return AB_df


def AB2EB_multi(AB_df, config={"threads": 8}):
    """
    multi-threaded converter of AB params and Tumor to EBscore
    """

    threads = config["threads"]

    pool = Pool(threads)

    AB_len = len(AB_df.index)
    show_output(f"Starting EBscore computation of {AB_len} using {threads} cores")

    # split the matrix
    split = np.array_split(AB_df, threads)
    dfs = pool.imap(AB2EB, split)
    pool.close()
    # out_df contains AB params
    EB_df = pd.concat(dfs).reset_index(drop=True)
    show_output(f"EBscore finished!", color="success")
    return EB_df