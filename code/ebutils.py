import pandas as pd
import os


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


def get_count_df(pon_string):
    """
    converts the base-wise read coverage to a matrix
    """

    count_dict = {0: "alt+", 1: "alt-", 2: "depth+", 3: "depth-"}

    data_split = [s for ad in pon_string.split("=") for s in ad.split("-")]
    count_df = pd.DataFrame.from_dict(
        {count_dict[i]: v.split("|") for i, v in enumerate(data_split)}
    )
    return count_df.astype(int)


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
