import pandas as pd
import numpy as np
import math
import os
from multiprocessing import Pool
from functools import partial

from script_utils_EB import show_output

from zerocache import update_zero_file, flatten_df
from ebcore import fit_bb

# ## matrix to AB
def matrix2AB_row(row, pen=0.5, lines=1000):
    """
    takes a row of T and D, passes the matrix to fit_bb
    returns the AB string
    """
    steps = round(lines / 10)
    if (row.name % steps == 0) and row.name:
        show_output(f"{round(row.name/lines * 100)}% finished", multi=True)
    d = row["D"]
    t = row["T"]
    count_matrix = np.transpose([d.split("|"), t.split("|")]).astype(int)
    return fit_bb(count_matrix)


# #############   convert PONmatrix to PON AB
def matrix2AB(config, matrix_df):
    """
    config = {
        "fit_pen": 0.5,
        "threads": 8,
        "AB_chunksize": 1000,
        "ZDfactor": 13 # how much complexity remains after flattening the tumor-zero lines
    }
    """

    # reset for line counting
    matrix_df = matrix_df.reset_index(drop=True)
    lines = len(matrix_df.index)
    show_output(
        f"Computing {lines} lines!",
        multi=True,
    )
    matrix_df.loc[:, "AB"] = matrix_df.apply(
        matrix2AB_row, pen=config["fit_pen"], lines=lines, axis=1
    )

    show_output(
        f"AB computation finished",
        multi=True,
        color="success",
    )
    return matrix_df
