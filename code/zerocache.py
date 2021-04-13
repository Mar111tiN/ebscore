import pandas as pd
import os


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
        if os.path.isfile(f):
            zdf = pd.read_csv(f, sep="\t")
            # show_output(f"Loading zero cache file {os.path.basename(f)}", multi=True)
            zero_dfs.append(zdf)
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


def collapse_zeros(zero_path):
    """
    reduces all zero_files to zero.0.csv
    """
    # get all the zero files
    zero_files = [
        os.path.join(zero_path, file)
        for file in os.listdir(zero_path)
        if os.path.isfile(os.path.join(zero_path, file))
    ]
    zdfs = []
    for file in zero_files:
        print(file)
        zdfs.append(pd.read_csv(file, sep="\t"))
        os.remove(file)
    zero_df = pd.concat(zdfs).drop_duplicates("D").sort_values("D")

    zero_df.to_csv(os.path.join(zero_path, "zero.0.csv"), sep="\t", index=False)
    return zero_df