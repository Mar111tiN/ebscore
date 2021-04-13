import pandas as pd
import os


def extract_zero_df(stack_df, zero_string="0|0|0|0|0"):
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
    print(file_counter)
    load_file_counter = file_counter[-4:-2]
    if not len(load_file_counter):
        load_file_counter = file_counter[-2:-1]
        if not len(load_file_counter):
            # if nothing is there, load the zero.0.csv
            load_file_counter = [0]
    load_files = [os.path.join(zero_path, f"zero.{i}.csv") for i in load_file_counter]
    zero_dfs = []
    for f in load_files:
        if os.path.isfile(f):
            print(f)
            zero_df = pd.read_csv(f, sep="\t")
            zero_dfs.append(zero_df)
    if len(zero_dfs):
        zero_df = pd.concat(zero_dfs).drop_duplicates("D")
    else:
        return None
    return zero_df


def get_next_zero(zero_path):
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
    if len(load_file_counter):
        next_file_count = file_counter[-1]
    else:
        next_file_count = 0
    return os.path.join(zero_path, f"zero.{next_file_count}.csv")


def get_pon_df(pon_list, pon_path):
    """
    complete the pon_list to a full_path pon_list
    """
    if pon_path and not pon_path.endswith("/"):
        pon_path += "/"
    # load the pon_list into df
    pon_df = pd.read_csv(pon_list, header=None, names=["bam"])
    # add the pon_path to all relative paths (not starting with /)
    pon_df.loc[~pon_df["bam"].str.startswith("/"), "bam"] = (
        f"{pon_path}" + pon_df["bam"]
    )
    return pon_df


def get_pon(file, pon_list, pon_path="", prepend_bam=False):
    """
    checks whether bam or pileup file is contained in the pon_list
    returns reduced pon_list and the matching 1-based position in the pon_list
    """
    # get the stripped sample name
    sample_base = os.path.basename(file).split(".")[0].split("_")[0].lstrip("0")
    # pon path is not "", add a backslash

    pon_df = get_pon_df(pon_list, pon_path)
    # extract the stripped sample into sample column
    pon_df["basename"] = pon_df["bam"].str.split("/").str[-1]
    pon_df["sample"] = (
        pon_df["basename"]
        .str.replace(".bam", "")
        .str.split("_", expand=True)[0]
        .str.lstrip("0")
    )
    # DEBUG
    # print(sample_base, pon_df['sample'])

    # add bam to beginning if prepend_bam == True
    if prepend_bam:
        pon_df = pd.concat([pd.DataFrame([{"bam": file}]), pon_df]).reset_index()

    if pon_df.query("sample == @sample_base").empty:
        return pon_df.loc[:, ["bam"]], 0, ""
    # bam is in pon
    reduced_pon_df = pon_df.query("sample != @sample_base").loc[:, ["bam"]]
    pos = pon_df.query("sample == @sample_base").index[0] + 1 - int(prepend_bam)
    removed_pon = pon_df.query("sample == @sample_base").iloc[0]["basename"]
    return reduced_pon_df, pos, removed_pon


######### AB2EB ######################


def retrieveABdata(row):
    # retrieve the data from the row
    # target_s:
    # turn string "0-5=12-42" into target_s:
    # pd.Series
    # alt+       0
    # alt-       5
    # depth+    12
    # depth-    42
    target = row["Tumor"]
    count_dict = {0: "alt+", 1: "alt-", 2: "depth+", 3: "depth-"}
    target_split = [s for ad in target.split("=") for s in ad.split("-")]
    target_s = pd.Series({count_dict[i]: v for i, v in enumerate(target_split)}).astype(
        int
    )
    # params:
    # turn string A+|B+-A-|B- into AB dict {'+':[A+, B+], '-':[A-, B-]}
    params = row["ABparams"]
    AB_list = [float(ab) for s in params.split("-") for ab in s.split("|")]
    AB_dict = {"+": AB_list[:2], "-": AB_list[2:]}
    return target_s, AB_dict


def get_obs_df(target_s, cols):
    """
    turn the target_s into obs_df
    """

    # cols is either ['depth+', 'alt+'] or ['depth-', 'alt-']
    alt_type = cols[1]
    # creates an observation df for each observation from depth-alt to depth-depth
    n_minus_k = target_s[cols[0]] - target_s[alt_type]
    # obs_df is instantiated from target_s dict with n_minus_k + 1 rows
    obs_df = pd.DataFrame(target_s[cols].to_dict(), index=range(n_minus_k + 1))
    # alt column is incremented using index
    obs_df[alt_type] = obs_df[alt_type] + obs_df.index
    return obs_df
