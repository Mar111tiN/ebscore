import pandas as pd
import os
import re
import numpy as np
from script_utils_EB import show_output


def load_zero_list(zero_path, ponsize=10):
    """
    load the sorted zero_list for that ponsize and returns list and counter list
    """

    # get all the relevant zero files for that ponsize

    zero_files = [
        file
        for file in os.listdir(zero_path)
        if os.path.isfile(os.path.join(zero_path, file)) and f"zero{ponsize}." in file
    ]
    # get the numbers of the files
    # zero files are named <PONPATH>/zero/zero<ponsize>.INT.gz
    file_counter = sorted([int(f.split(".")[1]) for f in zero_files])
    # get the file list back from load_file_counter
    zero_file_list = [
        os.path.join(zero_path, f"zero{ponsize}.{i}.gz") for i in file_counter
    ]

    return zero_file_list, file_counter


def save_next_zero(zero_df, zero_path, ponsize=10):
    """
    saves the zero_df into the next free spot
    """
    # get all the zero files
    _, file_counter = load_zero_list(zero_path, ponsize)
    # get the next zero_file
    if len(file_counter):
        next_file_count = file_counter[-1] + 1
    else:
        next_file_count = 0

    next_zero_file = os.path.join(zero_path, f"zero{ponsize}.{next_file_count}.gz")

    # create the zero_path if not existing
    if not os.path.isdir(zero_path):
        os.makedirs(zero_path)

    zero_df.to_csv(next_zero_file, sep="\t", index=False, compression="gzip")
    show_output(
        f"Saving updated zero cache ({len(zero_df.index)} lines) to {os.path.basename(next_zero_file)}",
        multi=False,
    )


def zero2AB(stack_df, config):
    """
    go through zero files and get matching ABs from zero cache for tumor == zerostring
    """

    # split off the zeroT_df
    zeroT_df = stack_df.loc[stack_df["T"] == config["zerostring"], :].sort_values("D")
    # if zeroT_df is too small, just return stack_df and empty AB_df
    if len(zeroT_df) < config["min_zt"]:
        useZero = False
        return stack_df, pd.DataFrame(), useZero

    stack_df = stack_df.loc[stack_df["T"] != config["zerostring"], :]

    # load the zerofiles and merge incrementally to save memory
    AB_dfs = []
    for zfile in load_zero_list(config["zero_path"], config["ponsize"])[0]:
        try:
            zgen = pd.read_csv(
                zfile, sep="\t", compression="gzip", chunksize=config["chunksize"]
            )
            show_output(f"Loading {zfile}")
            for zdf in zgen:
                merge_df = zeroT_df.merge(zdf, how="left")
                AB_df = merge_df.query("AB == AB")
                AB_dfs.append(AB_df)
                zeroT_df = merge_df.query("AB != AB").drop("AB", axis=1)
        except Exception as e:
            show_output(f"{zfile} could not be loaded, {e}", color="warning")

    AB_df = pd.concat(AB_dfs)

    # load remaining
    stack_df = pd.concat([stack_df, zeroT_df]).reset_index(drop=True)
    show_output(
        f"{len(AB_df.index)} matching lines found in zero cache! Computing remaining {len(stack_df.index)} lines.",
        multi=False,
    )
    useZero = True
    return stack_df, AB_df, useZero


def update_zero_file(AB_df, config={}):
    """
    extracts the zero_df from the computed AB_df
    saves new zero file if exists
    """

    # extracts the zero_df from the computed AB_df
    new_zero_df = AB_df.loc[
        AB_df["T"] == config["zerostring"], ["D", "AB"]
    ].drop_duplicates("D")

    if not new_zero_df.empty:
        # write to new zero
        save_next_zero(
            new_zero_df, zero_path=config["zero_path"], ponsize=config["ponsize"]
        )
    return AB_df


def get_zerostring(df):
    """
    returns the zerostring and pon_length from first row of stacked df
    """

    zerostring = re.sub(r"[0-9]+", "0", df.iloc[0]["D"])
    pon_len = int((len(zerostring) + 1) / 2)
    return zerostring, pon_len


def flatten_zeros(col, ZDfactor=13):
    return col.str.split("|").apply(
        lambda li: "|".join(
            np.round(
                np.exp(
                    np.round(
                        (np.log(np.sort(np.array(li).astype(int))) * ZDfactor),
                        0,
                    )
                    / ZDfactor
                ),
                0,
            )
            .astype(int)
            .astype(str)
        )
    )


def flatten_df(df, ZDfactor=13):
    """
    flatten a stacked df using ZDfactor
    """

    # get the zerostring
    # reduce zeros with condense_factor
    zerostring = get_zerostring(df)

    # sort the depths at tumor == zero and reduce zero_complexity via flatten_zero
    df.loc[df["T"] == zerostring, "D"] = flatten_zeros(
        df.loc[df["T"] == zerostring, "D"], ZDfactor
    )
    return df


def collapse_zeros(zero_path, ponsize=10, reflat=False, ZDfactor=13):
    """
    reduces all zero_files to zero.0.csv
    """
    # get all the zero files

    zero_files, _ = load_zero_list(zero_path, ponsize)
    zero_df = pd.DataFrame()
    if len(zero_files) == 0:
        show_output(f"No zerofiles found in {zero_path}!", color="warning")
        return
    if len(zero_files) == 1:
        if not reflat:
            show_output(
                f"Only one file found in {zero_path}! No need to collapse!",
                color="success",
            )
            return
    show_output(
        f"Collapsing all {len(zero_files)} zero files in {zero_path} for PONsize {ponsize} into zero{ponsize}.0.gz"
    )

    zdfs = []
    for file in zero_files:
        show_output(f"Cleaning up zero file {file}", time=False)
        try:
            zdf = pd.read_csv(file, sep="\t", compression="gzip")
            # concat to zdf_all and drop duplicates
            zdfs.append(zdf)
            # reduce complexity using condense factor
        except:
            show_output(f"{file} could not be loaded", color="warning", time=False)

    show_output("Collapsing all zeros into one file")
    zero_df = pd.concat(zdfs).drop_duplicates("D").sort_values("D")

    # reflat if selected using ZDfactor
    if reflat:
        # reapply the flatten procedure
        show_output(f"Flatten zero file with condense_factor {ZDfactor}")
        zero_df.loc[:, "D"] = flatten_zeros(zero_df["D"], ZDfactor=ZDfactor)
        # remove the created duplicates
        zero_df = zero_df.drop_duplicates("D").sort_values("D")
    zero0_file = os.path.join(zero_path, f"zero{ponsize}.0.gz")
    zero_df.to_csv(
        zero0_file,
        sep="\t",
        index=False,
        compression="gzip",
    )
    show_output(
        f"Written collapsed zerofile ({len(zero_df.index)} lines) to {zero_path}/{zero0_file}",
        color="success",
    )
    for file in zero_files:
        if file == zero0_file:
            continue
        os.remove(file)
    return zero_df
