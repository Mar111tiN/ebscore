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


def stackPONmatrix(config, PON_df):
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
    df[["DL", "DR"]] = df["Depth"].str.split(STRANDSEP, expand=True)
    df[["TL", "TR"]] = df["tumor"].str.split(STRANDSEP, expand=True)
    df["L"] = df["TL"] + ADSEP + df["DL"]
    df["R"] = df["TR"] + ADSEP + df["DR"]
    # stack left and right
    df = (
        df.loc[:, ["Chr", "Start", "Ref", "Alt", "L", "R"]]
        .set_index(["Chr", "Start", "Ref", "Alt"])
        .stack()
        .reset_index()
        .rename({"level_4": "strand", 0: "PON"}, axis=1)
    )

    df[["T", "D"]] = df["PON"].str.split(ADSEP, expand=True)

    # reduce cols
    df = df.loc[:, ["Chr", "Start", "Ref", "Alt", "strand", "T", "D"]]

    df = flatten_df(df, ZDfactor=config["ZDfactor"])

    return df


def unstack_PONAB(stack_df):
    """
    convert the stacked PON_df back into the normal PON_df for human readable output
    """

    unstack_df = (
        stack_df.loc[:, ["Chr", "Start", "Ref", "Alt", "strand", "AB"]]
        .set_index(["Chr", "Start", "Ref", "Alt", "strand"])
        .unstack("strand")["AB"]
    )
    unstack_df["AB"] = unstack_df["L"] + STRANDSEP + unstack_df["R"]
    return (
        unstack_df.loc[:, "AB"]
        .unstack("Alt")
        .reset_index(drop=False)
        .sort_values(["Start"])
    )


def PONmatrix2AB_multi(
    pon_matrix_gen,
    config=dict(fit_pen=0.5, threads=8, zero_path="./zero", ZDfactor=13, min_zt=2000),
):
    """
    converts a PONmatrix (as generator to save memory) into a PONAB file
    """

    threads = config["threads"]
    # stack the PONmatrix from generator using multiprocessing
    stack_pool = Pool(threads)
    show_output(f"Stacking PON matrix using {threads} cores.")
    dfs = stack_pool.imap(partial(stackPONmatrix, config), pon_matrix_gen)
    stack_df = pd.concat(dfs)
    stack_pool.close()
    show_output("PON matrix has been stacked")
    # print(stack_df[:100])
    # retrieve the zerostring and ponsize from the panel of normals of first row
    # store in configs
    config["zerostring"], config["ponsize"] = get_zerostring(stack_df)

    # ######### reading AB from zerocache into AB_df
    stack_df, AB_df = zero2AB(stack_df, config)

    # case all AB have been retrieved from zero2AB
    if stack_df.empty:
        show_output("PON matrix successfully converted!", color="success")
        return unstack_PONAB(AB_df).reset_index(drop=True).sort_values(["Start"])

    # compute the AB_df for the stack_df
    AB_pool = Pool(threads)
    show_output(f"Starting AB conversion of matrix using {threads} cores")

    # minimal length of 200 lines
    split_factor = max(1, min(int(len(stack_df.index) / 200), threads))

    # split the matrix
    split = np.array_split(stack_df, split_factor)
    del stack_df
    AB_dfs = AB_pool.imap(partial(matrix2AB, config), split)
    AB_pool.close()

    # collect AB_dfs + get all the new zeros and write to file

    new_AB_df = pd.concat(AB_dfs)

    # check, if zero has been used (or zero_df was smaller then min_zt)
    _ = update_zero_file(new_AB_df, config=config)

    AB_df = unstack_PONAB(
        pd.concat([AB_df, new_AB_df]).reset_index(drop=True)
    ).sort_values(["Start"]).loc[:, list("AGCTID")]
    show_output("PON matrix successfully converted!", color="success")
    return AB_df
