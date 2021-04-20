import pandas as pd
import os
import numpy as np
from script_utils import show_output


def extract_zero_df(stack_df, zero_string="0|0|0|0|0"):
    """
    takes the computed stack_df and filters out the tumor==0 data
    """
    zero_df = stack_df.loc[stack_df["T"] == zero_string, ["D", "AB"]].drop_duplicates(
        "D"
    )
    return zero_df


def load_zero_list(zero_path, pon_size=10):
    """
    load the sorted zero_list for that pon_size and returns list and counter list
    """

    # get all the relevant zero files for that pon_size

    zero_files = [
        file
        for file in os.listdir(zero_path)
        if os.path.isfile(os.path.join(zero_path, file)) and f"zero{pon_size}." in file
    ]
    # get the numbers of the files
    # zero files are named <PONPATH>/zero/zero<pon_size>.INT.gz
    file_counter = sorted([int(f.split(".")[1]) for f in zero_files])
    # get the file list back from load_file_counter
    zero_file_list = [
        os.path.join(zero_path, f"zero{pon_size}.{i}.gz") for i in file_counter
    ]

    return zero_file_list, file_counter


def load_zero_df(zero_path, pon_size=10):
    """
    looks into zero_path and gets the previous zero_files into one df
    zero files are named <zero_path>/zero/zero<pon_size>.INT.gz
    """

    # try getting the two files before the last 2
    zero_file_list, _ = load_zero_list(zero_path, pon_size)

    # empty list
    if not len(zero_file_list):
        return None
    # if more than 2 files get the 4th (if any) and (3rd) latest file else try at least zeroXX.0.gz
    if len(zero_file_list) > 2:
        load_file_list = zero_file_list[:-2]
    else:
        load_file_list = [os.path.join(zero_path, f"zero{pon_size}.0.gz")]

    zero_dfs = []

    ### DEBUG
    print(
        f"Try loading zero files {' '.join([os.path.basename(f) for f in load_file_list])}"
    )
    ###

    for f in load_file_list:
        try:
            zdf = pd.read_csv(f, sep="\t", compression="gzip")
            zero_dfs.append(zdf)
            show_output(f"Loading zero cache file {os.path.basename(f)}", multi=True)
        except:
            show_output(f"{f} could not be loaded", color="warning", time=False)

    # if nothing could be found
    if not len(zero_dfs):
        return None
    zero_df = pd.concat(zero_dfs).drop_duplicates("D")
    # show_output(f"{len(zero_df.index)} zero cache lines loaded", multi=True)

    return zero_df


def update_zero_file(AB_df, config={}):
    """
    extracts the zero_df from the computed AB_df
    reloads the zero_df (could be updated)
    compares new zero_df with reloaded zero_df and if any changes:
    saves new zero file
    """

    # extracts the zero_df from the computed AB_df
    updated_zero_df = extract_zero_df(AB_df, zero_string=config["zero_string"])

    zero_path = os.path.join(config["pon_path"], "zero")

    pon_size = config["pon_size"]
    # reload the current zero_df state(could be updated in between)
    latest_df = load_zero_df(zero_path)

    if latest_df is None:
        # count the number of lines from previous zero_df
        old_lines = 0
    else:
        old_lines = len(latest_df.index)
        new_zero_df = pd.concat([updated_zero_df, latest_df]).drop_duplicates(
            "D", keep=False
        )

    if len(new_zero_df.index) > old_lines:
        # write to new zero
        save_next_zero(new_zero_df, zero_path, pon_size=pon_size)


def save_next_zero(zero_df, zero_path, pon_size=10):
    """
    returns the next new zero file to write to
    """
    # get all the zero files
    zero_files, file_counter = load_zero_list(zero_path, pon_size)
    # get the next zero_file
    if len(file_counter):
        next_file_count = file_counter[-1] + 1
    else:
        next_file_count = 0
    next_zero_file = os.path.join(zero_path, f"zero{pon_size}.{next_file_count}.gz")

    zero_df.to_csv(next_zero_file, sep="\t", index=False, compression="gzip")
    show_output(
        f"Saving updated zero cache to {os.path.basename(next_zero_file)}", multi=True
    )


def flatten_zeros(col, zero_condense_factor=13):
    return col.str.split("|").apply(
        lambda li: "|".join(
            np.round(
                np.exp(
                    np.round(
                        (
                            np.log(np.sort(np.array(li).astype(int)))
                            * zero_condense_factor
                        ),
                        0,
                    )
                    / zero_condense_factor
                ),
                0,
            )
            .astype(int)
            .astype(str)
        )
    )


def collapse_zeros(zero_path, pon_size=10, zero_condense_factor=13):
    """
    reduces all zero_files to zero.0.csv
    """
    # get all the zero files

    zero_files, _ = load_zero_list(zero_path, pon_size)
    zero_df = pd.DataFrame()
    if len(zero_files):
        show_output(
            f"Collapsing all {len(zero_files)} zero files in {zero_path} for PONsize {pon_size} into zero{pon_size}.0.gz"
        )
    else:
        show_output(f"No zerofiles found in {zero_path}!")
        return
    for file in zero_files:
        show_output(f"Cleaning up zero file {file}", time=False)
        try:
            zdf = pd.read_csv(file, sep="\t", compression="gzip")
            # concat to zdf_all and drop duplicates
            zero_df = pd.concat([zero_df, zdf]).drop_duplicates("D")
            # reduce complexity using condense factor
        except:
            show_output(f"{file} could not be loaded", color="warning", time=False)

    zero_df = zero_df.sort_values("D")
    # zero_df.loc[:, "D"] = flatten_zeros(
    #     zero_df["D"], zero_condense_factor=zero_condense_factor
    # )
    zero0_file = os.path.join(zero_path, f"zero{pon_size}.0.gz")
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
