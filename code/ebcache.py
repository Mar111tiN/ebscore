import os
from multiprocessing import Pool
from functools import partial
import math
import numpy as np
import pandas as pd

from ebcore import fit_bb
from ebutils import get_count_df, get_pon_df
from script_utils import show_output, run_cmd


def PON2matrix(pon_list, chrom, EBconfig={}):
    """
    generates matrix file from pon_list
    """

    # PARAMS
    # mawk tool unwrapper
    def mawk(tool):
        return os.path.join(EBconfig["mawk_path"], f"{tool}.mawk")

    gsplit = EBconfig["genome_split"]
    pon_path = EBconfig["pon_path"]
    q = EBconfig["MAPQ"]
    Q = EBconfig["Q"]
    bed = EBconfig["bed_file"]

    matrix_path = os.path.join(pon_path, "matrix")
    if not os.path.isdir(matrix_path):
        os.mkdir(matrix_path)

    # check the temp_folder with default pon_path/temp
    temp_folder = EBconfig.get("temp_dir", os.path.join(pon_path, "temp"))
    if not os.path.isdir(temp_folder):
        os.mkdir(temp_folder)
    pon_list_full = os.path.join(EBconfig["temp_dir"], f"full_{pon_list}")

    # create the pon_list with full path
    pon_list = os.path.join(EBconfig["pon_path"], pon_list)

    if not os.path.isfile(pon_list_full):
        get_pon_df(pon_list, pon_path).to_csv(
            pon_list_full, sep="\t", index=False, header=False
        )

    pileup_cmd = f"samtools mpileup -Q {Q} -q {q} -l {bed} -f {gsplit}/{chrom}.fa -b {pon_list_full} -r {chrom}"

    pon_matrix_file = os.path.join(matrix_path, f"{chrom}.pon")
    cmd = f"{pileup_cmd} | {mawk('cleanpileup')} | {mawk('pile2count')} | tee {pon_matrix_file} |  gzip  > {pon_matrix_file}.gz"
    run_cmd(cmd)
    return pon_matrix_file


##############   convert PONmatrix to PON AB


def PONmatrix2AB_row(row, pen=0.5):
    for base in "AGCTID":
        row[base] = fit_bb(
            get_count_df(row[base] + "=" + row["Depth"]),
            pen,
        )
    return row.loc[["Chr", "Start", "Ref"] + list("AGCTID")]


def PONmatrix2AB(config, matrix_df):
    finished_line = matrix_df.iloc[0].name + config["chunk_size"]
    percent = round(finished_line / config["len"] * 100, 1)
    df = matrix_df.apply(PONmatrix2AB_row, pen=config["fit_pen"], axis=1)
    show_output(
        f"{finished_line} of {config['len']} lines ({percent}%) finished!",
        multi=True,
        time=True,
    )
    return df


def PONmatrix2AB_multi(
    pon_matrix_df, config={"fit_pen": 0.5, "threads": 8, "chunk_size": 50000}
):
    """
    converts a PONmatrix into a PONAB file
    """

    threads = config["threads"]

    pool = Pool(threads)

    pon_len = len(pon_matrix_df.index)
    config["len"] = pon_len
    # minimal length of 200 lines
    # split_factor = min(math.ceil(len(pon_matrix_df.index) / 200), threads)
    split_factor = math.ceil(pon_len / config["chunk_size"])
    # split the matrix
    split = np.array_split(pon_matrix_df, split_factor)
    dfs = pool.imap(partial(PONmatrix2AB, config), split)
    pool.close()
    # out_df contains EB_score
    out_df = pd.concat(dfs).sort_values(["Chr", "Start"])
    return out_df