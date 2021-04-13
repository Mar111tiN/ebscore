import pandas as pd
import numpy as np
import math
import os
from multiprocessing import Pool
from functools import partial

from script_utils import show_output
from ebcore import bb_pvalue, fisher_combination


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
    params = row["AB"]
    AB_list = [float(ab) for s in params.split("=") for ab in s.split("|")]
    AB_dict = {"+": AB_list[:2], "-": AB_list[2:]}
    return target_s, AB_dict


######### EB2AB ######################
def AB2EBscore(row):
    """
    takes a df containing an AB column of shape "A+|A--B+|B-""
    """
    # retrieve the data from the row
    target_s, AB_dict = retrieveABdata(row)

    # get the p_values for each strand
    p_values = {}
    p_values["+"] = bb_pvalue(get_obs_df(target_s, ["depth+", "alt+"]), AB_dict["+"])
    p_values["-"] = bb_pvalue(get_obs_df(target_s, ["depth-", "alt-"]), AB_dict["-"])
    # combine p_values with fisher combination
    EB_p = fisher_combination(p_values)
    if EB_p < 1e-60:
        return 60
    if EB_p > 1.0 - 1e-10:
        return 0
    return -round(math.log10(EB_p), 3)