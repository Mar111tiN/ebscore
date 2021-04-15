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


def load_zero_df(zero_path):
    """
    looks into zero_path and gets the previous zero_files into one df
    """

    # get all the zero files
    zero_files = [
        file
        for file in os.listdir(zero_path)
        if os.path.isfile(os.path.join(zero_path, file))
    ]
    # get the numbers of the files
    file_counter = sorted(
        [int(f.replace("zero.", "").replace(".csv", "")) for f in zero_files]
    )
    # try getting the two files before the last 2
    load_file_counter = file_counter[-4:-2]
    if not len(load_file_counter):
        # try getting the second last file
        load_file_counter = file_counter[-2:-1]
        if not len(load_file_counter):
            # if nothing is there, load the zero.0.csv
            load_file_counter = [0]
    # get the file list
    load_files = [os.path.join(zero_path, f"zero.{i}.csv") for i in load_file_counter]
    zero_dfs = []
    for f in load_files:
        try:
            zdf = pd.read_csv(f, sep="\t")
            zero_dfs.append(zdf)
            # show_output(f"Loading zero cache file {os.path.basename(f)}", multi=True)
        except:
            show_output(f"{f} could not be loaded", color="warning", time=False)
            # show_output(f"Loading zero cache file {os.path.basename(f)}", multi=True)
            
    # if nothing could be found
    if len(zero_dfs):
        zero_df = pd.concat(zero_dfs).drop_duplicates("D")
        # show_output(f"{len(zero_df.index)} zero cache lines loaded", multi=True)
    else:
        return None
    return zero_df


def get_next_zero(zero_path):
    """
    returns the next new zero file to write to
    """
    # get all the zero files
    zero_files = [
        file
        for file in os.listdir(zero_path)
        if os.path.isfile(os.path.join(zero_path, file))
    ]
    # get the numbers of the files
    file_counter = sorted(
        [int(f.replace("zero.", "").replace(".csv", "")) for f in zero_files]
    )
    if len(file_counter):
        next_file_count = file_counter[-1] + 1
    else:
        next_file_count = 0
    return os.path.join(zero_path, f"zero.{next_file_count}.csv")


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


def collapse_zeros(zero_path, zero_condense_factor=13):
    """
    reduces all zero_files to zero.0.csv
    """
    # get all the zero files
    show_output(f"Collapsing all zero files in {zero_path} into zero.0.csv")
    zero_files = [
        os.path.join(zero_path, file)
        for file in os.listdir(zero_path)
        if os.path.isfile(os.path.join(zero_path, file)) and "zero" in file
    ]
    zdfs = []
    if len(zero_files) == 0:
        show_output(f"No zerofiles found in {zero_path}!")
        return
    for file in zero_files:
        show_output(f"Cleaning up zero file {file}", time=False)
        try:
            zdf = pd.read_csv(file, sep="\t")
            zdfs.append(zdf)
        except:
            show_output(f"{file} could not be loaded", color="warning", time=False)
        os.remove(file)
    # concat and drop duplicates
    zero_df = pd.concat(zdfs).drop_duplicates("D")
    # reduce complexity using condense factor
    zero_df.loc[:, "D"] = flatten_zeros(
        zero_df["D"], zero_condense_factor=zero_condense_factor
    )

    zero_df = zero_df.drop_duplicates("D").sort_values("D")

    zero_df.to_csv(os.path.join(zero_path, "zero.0.csv"), sep="\t", index=False)
    show_output(f"Written collapsed zerofile ({len(zero_df.index)} lines) to {zero_path}/zero.0.csv")
    return zero_df