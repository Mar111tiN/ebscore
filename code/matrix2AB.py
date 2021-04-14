import pandas as pd
import numpy as np
import math
import os
from multiprocessing import Pool
from functools import partial

from script_utils import show_output

from zerocache import extract_zero_df, load_zero_df, get_next_zero, flatten_zeros
from ebcore import fit_bb

### matrix to AB
def matrix2AB_row(row, pen=0.5):
    """
    takes a row of T and D, passes the matrix to fit_bb
    returns the AB string
    """

    d = row["D"]
    t = row["T"]
    count_matrix = np.transpose([d.split("|"), t.split("|")]).astype(int)
    return fit_bb(count_matrix)


##############   convert PONmatrix to PON AB
def matrix2AB(config, matrix_df):
    """
    config = {
        "fit_pen": 0.5,
        "threads": 8,
        "chunk_size": 1000,
        "zero_file": os.path.join(pon_path, "zero.csv"),
        "zero_condense_factor": 13 # how much complexity remains after flattening the tumor-zero lines
    }
    """

    # finished_line = matrix_df.iloc[0].name + config["chunk_size"]
    # percent = round(finished_line / config["len"] * 100, 1)

    zero_string = config["zero_string"]
    show_output(
        f"Computing {len(matrix_df.index)} lines!",
        multi=True,
    )
    # sort the depths at tumor == zero
    matrix_df.loc[matrix_df["T"] == zero_string, "D"] = flatten_zeros(
        matrix_df.loc[matrix_df["T"] == zero_string, "D"],
        zero_condense_factor=config["zero_condense_factor"],
    )

    # load the zero_df
    # check if zero_folder exists
    zero_path = config["zero_path"]
    if not os.path.isdir(zero_path):
        os.makedirs(zero_path)

    zero_df = load_zero_df(zero_path)
    if zero_df is None:
        use_zero = False
    else:
        # merge the AB params from the zero_df
        matrix_df = matrix_df.merge(zero_df, how="left")
        # store the joint data
        use_zero = True
        AB_df = matrix_df.query("AB == AB")
        # for all empty AB, split and compute
        matrix_df = matrix_df.query("AB != AB")

        show_output(
            f"{len(AB_df.index)} matching lines found in zero cache! Computing remaining {len(matrix_df.index)} lines.",
            multi=True,
        )

    matrix_df.loc[:, "AB"] = matrix_df.apply(
        matrix2AB_row, pen=config["fit_pen"], axis=1
    )

    # get all the new zeros and merge with zero_df and write to file

    # here comes the zero_file update
    updated_zero_df = extract_zero_df(matrix_df, zero_string=zero_string)

    # reload the current zero_df (could be updated in between)
    old_zero_df = load_zero_df(zero_path)

    if old_zero_df is None:
        # count the number of lines from previous zero_df
        old_lines = 0
    else:
        old_lines = len(old_zero_df.index)
        updated_zero_df = pd.concat([updated_zero_df, old_zero_df]).drop_duplicates("D")

    if len(zero_df.index) > old_lines:
        # write to new zero
        zero_file = get_next_zero(zero_path)
        zero_df.to_csv(zero_file, sep="\t", index=False)
        show_output(
            f"Saving updated zero cache to {os.path.basename(zero_file)}", multi=True
        )
    show_output(
        f"AB computation finished",
        multi=True,
        color="success",
    )
    return pd.concat([matrix_df, AB_df]) if use_zero else matrix_df


######### compute AB directly from tumor-matrix file
def stack_matrix(df):
    """
    create a stacked matrix_df for direct matrix2AB computation
    """

    df = (
        df.set_index(["Chr", "Start", "End", "Ref", "Alt", "Tumor"])
        .stack()
        .reset_index()
        .rename({"level_6": "strand", 0: "PON"}, axis=1)
    )
    df[["T", "D"]] = df["PON"].str.split("=", expand=True)
    return df.loc[:, ["Chr", "Start", "End", "Ref", "Alt", "Tumor", "strand", "T", "D"]]


def unstackAB(stackAB_df):
    """
    converts the AB_matrix_df into unstacked version for further processing
    """
    un1 = (
        stackAB_df.loc[
            :, ["Chr", "Start", "End", "Ref", "Alt", "Tumor", "strand", "AB"]
        ]
        .set_index(["Chr", "Start", "End", "Ref", "Alt", "Tumor", "strand"])
        .unstack("strand")["AB"]
    )
    un1["AB"] = un1["PON+"] + "=" + un1["PON-"]
    un2 = un1.reset_index(drop=False).sort_values(["Start", "End"])
    return un2.loc[:, ["Chr", "Start", "End", "Ref", "Alt", "Tumor", "AB"]]


def matrix2AB_multi(
    matrix_df,
    config={
        "fit_pen": 0.5,
        "threads": 8,
        "chunk_size": 500,
        "zero_path": "zero",
        "zero_condense_factor": 13,
    },
):
    """
    converts a matrix_df into a AB file using multicore
    """

    threads = config["threads"]

    pool = Pool(threads)
    show_output(f"Starting AB conversion of matrix using {threads} cores")

    # tidy the pon_matrix_df
    stack_df = stack_matrix(matrix_df)
    show_output(f"matrix has been stacked")

    pon_len = len(stack_df.index)

    config["len"] = pon_len

    # retrieve the zero_string from the panel of normals of first row
    config["zero_string"] = "|".join(
        np.array([0] * len(stack_df.loc[0, "D"].split("|"))).astype(str)
    )
    # minimal length of 200 lines
    # split_factor = min(math.ceil(len(pon_matrix_df.index) / 200), threads)
    split_factor = math.ceil(pon_len / config["chunk_size"])
    # split the matrix
    split = np.array_split(stack_df, split_factor)
    dfs = pool.imap(partial(matrix2AB, config), split)
    pool.close()
    # out_df contains AB params
    AB_df = unstackAB(pd.concat(dfs).reset_index(drop=True))
    show_output(f"matrix successfully converted!", color="success")
    return AB_df
