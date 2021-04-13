import pandas as pd
import numpy as np
import math
import os

from script_utils import show_output
from ebutils import retrieveABdata, get_obs_df, get_zero_df
from ebcore import fit_bb, bb_pvalue, fisher_combination

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
        "zero_file": os.path.join(pon_path, "zero.csv")
    }
    """
    # finished_line = matrix_df.iloc[0].name + config["chunk_size"]
    # percent = round(finished_line / config["len"] * 100, 1)

    # sort the depths at tumor zero
    matrix_df.loc[matrix_df["T"] == "0|0|0|0|0", "D"] = (
        matrix_df.loc[matrix_df["T"] == "0|0|0|0|0", "D"]
        .str.split("|")
        .apply(lambda li: "|".join(np.sort(np.array(li).astype(int)).astype(str)))
    )

    # load the zero_df
    zero_file = config["zero_file"]
    if os.path.isfile(zero_file):
        zero_df = pd.read_csv(zero_file, sep="\t")
        while zero_df.empty:
            zero_df = pd.read_csv(zero_file, sep="\t")
        # merge the AB params from the zero_df
        matrix_df = matrix_df.merge(zero_df, how="left")

        # store the joint data
        AB_df = matrix_df.query("AB == AB")

        # for all empty AB, split and compute
        matrix_df = matrix_df.query("AB != AB")

        show_output(
            f"{len(AB_df.index)} lines found! Computing remaining {len(matrix_df.index)} lines.",
            multi=True,
        )

    matrix_df.loc[:, "AB"] = matrix_df.apply(
        matrix2AB_row, pen=config["fit_pen"], axis=1
    )

    # get all the new zeros and merge with zero_df and write to file

    # here comes the zero_file update
    zero_df = get_zero_df(matrix_df, zero_string="0|0|0|0|0")
    # reload the current zero_df (could be updated in between)

    if os.path.isfile(zero_file):
        old_zero_df = pd.read_csv(zero_file, sep="\t")
        while old_zero_df.empty:
            old_zero_df = pd.read_csv(zero_file, sep="\t")
        pd.concat([zero_df, old_zero_df]).drop_duplicates("D").to_csv(
            zero_file, sep="\t", index=False
        )

    return pd.concat([matrix_df, AB_df])


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
