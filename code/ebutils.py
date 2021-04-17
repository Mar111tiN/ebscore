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
    pon_df.loc[:, "basename"] = pon_df["bam"].str.split("/").str[-1]
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
