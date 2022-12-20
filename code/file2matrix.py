import pandas as pd
import os

from ebutils import get_pon, get_pon_df
from script_utils_EB import show_output, cmd2df, run_cmd

# here come all the functions involved in converting raw data (bam or pileup) to dataframes


def tumor2matrix(
    mut_file, bam="", pileup="", cleanpileup="", pon_list="", chrom="", config={}
):
    """
    converts the bam_file and the Pon-bams into a matrix for EB computations
    does all the work behind the curtain:
    checks for PONcache file
    checks for ABcache file
    checks if PON is included in PON file
    """

    # PARAMS
    # mawk tool unwrapper
    def mawk(tool):
        return os.path.join(config["mawk_path"], f"{tool}.mawk")

    gsplit = config["genome_split"]
    pon_path = config["pon_path"]
    q = config["MAPQ"]
    Q = config["Q"]
    base_name = os.path.splitext(os.path.basename(mut_file))[0]
    use_cache = config["use_cache"]
    # check the temp_folder with default pon_path/temp
    temp_folder = config.get("temp_dir", os.path.join(pon_path, "temp"))
    os.makedirs(temp_folder, exist_ok=True)
        

    if bam:
        tumor_file = bam
        tumor_type = "bam"
        # create the bed file needed for pileup from bam-file
        # could be pythonized
        bed_file = os.path.join(
            temp_folder,
            f"{base_name}_{chrom}.bed",
        )
        run_cmd(f"{mawk('csv2bed')} {mut_file} {chrom} > {bed_file}")

        tumor_cmd = f"samtools mpileup -Q {Q} -q {q} -l {bed_file} -f {gsplit}/{chrom}.fa -r {chrom}"
    elif pileup:
        # for pileup, use the pileup file and only use first read (if mpileup comes from tumor-normal)
        tumor_file = pileup
        tumor_type = "pileup"
        if tumor_file.endswith(".gz"):
            tumor_cmd = "gunzip < "
        else:
            tumor_cmd = "cat "
        tumor_cmd += f"{tumor_file} | {mawk('cleanpileup')} -s 2"
    elif cleanpileup:
        tumor_file = cleanpileup
        tumor_type = "pileup"
        if tumor_file.endswith(".gz"):
            tumor_cmd = "gunzip < "
        else:
            tumor_cmd = "cat "
        # mawk snippet to convert
        mawk_cmd = 'BEGIN{OFS="\\t"}$7~/[ATCGatcg]/{print($1,$2,$3,$7)}'
        tumor_cmd += f"{tumor_file} | mawk '{mawk_cmd}'"
    else:
        show_output(
            "Neither bam nor pileup file is given! Either is required!",
            color="warning",
            multi=False,
        )

    # PON file
    # find out whether bam/pileup is contained in pon_list
    # create the pon_list
    pon_list = os.path.join(config["pon_path"], pon_list)

    # if normal matching the tumor is contained in PON file
    # normal_in_PON_pos > 0
    pon_df, normal_in_PON_pos, removed_ponbam = get_pon(
        tumor_file, pon_list, pon_path=pon_path, prepend_bam=True
    )

    if normal_in_PON_pos:
        show_output(
            f"{removed_ponbam} from PON list matches {tumor_file} and is removed from pon_list. If not desired, change names",
            color="warning",
        )

    # build the tumor2matrix command before other commands
    tumor2matrix_cmd = f"{mawk('tumor2matrix')} {mut_file} {chrom}"

    # #### use PON cache
    if use_cache:
        # check if PONcache matrix file is there
        matrix_file = os.path.join(config["pon_path"], f"matrix/{chrom}.pon.gz")
        if os.path.isfile(matrix_file):
            tumor2matrix_cmd += f" -P {matrix_file}"
        elif os.path.isfile(matrix_file := matrix_file.replace(".gz", "")):
            tumor2matrix_cmd += f" -P {matrix_file}"
        else:  # did not find the matrix file
            show_output(
                f"PON matrix {matrix_file} not found. Falling back to PON bam_files",
                color="warning",
                multi=False,
            )
            use_cache = False
        # check if ABcache file is there
        AB_file = os.path.join(config["pon_path"], f"AB/{chrom}.AB.gz")
        # check existence of matrix file
        if os.path.isfile(AB_file):
            tumor2matrix_cmd += f" -A {AB_file}"
        elif os.path.isfile(AB_file := AB_file.replace(".gz", "")):
            tumor2matrix_cmd += f" -A {AB_file}"

    if (
        use_cache
    ):  # re-check for bam files in case PONmatrix was not found (fallback to pileup of bam+PON)
        if tumor_type == "bam":
            # just add the bam_file as
            tumor_cmd += f" {tumor_file} | {mawk('cleanpileup')}"
        # add sample exclusion
        if normal_in_PON_pos:
            tumor2matrix_cmd += f" -x {normal_in_PON_pos}"

    else:  # no caching
        if tumor_type == "bam":
            # make filename for pon_bam
            pon_bam = os.path.join(temp_folder, f"{base_name}.pon")
            pon_df.to_csv(pon_bam, index=False, header=None)
            tumor_cmd += f" -b {pon_bam} | {mawk('cleanpileup')}"
        else:
            show_output(
                "If PON cache is not used, tumor pileup cannot be used!",
                color="warning",
                multi=False,
            )
            return None
    # print(pon_df)

    cmd = f"{tumor_cmd} | {mawk('pile2count')} | {tumor2matrix_cmd}"

    matrix_df = cmd2df(cmd)
    return matrix_df


def PON2matrix(pon_list, chrom, config={}):
    """
    generates matrix file from pon_list and writes it to pon_path/matrix/<chrom>.pon.gz
    """

    # PARAMS
    # mawk tool unwrapper
    def mawk(tool):
        return os.path.join(config["mawk_path"], f"{tool}.mawk")

    gsplit = config["genome_split"]
    pon_path = config["pon_path"]
    q = config["MAPQ"]
    Q = config["Q"]
    bed = config["bed_file"]

    matrix_path = os.path.join(pon_path, "matrix")
    os.makedirs(matrix_path, exist_ok=True)
        

    # check the temp_folder with default pon_path/temp
    temp_folder = config.get("temp_dir", os.path.join(pon_path, "temp"))
    os.makedirs(temp_folder, exist_ok=True)
    pon_list_full = os.path.join(config["temp_dir"], f"full_{pon_list}")

    # create the pon_list with full path
    pon_list = os.path.join(config["pon_path"], pon_list)

    get_pon_df(pon_list, pon_path).to_csv(
        pon_list_full, sep="\t", index=False, header=False
    )

    pileup_cmd = f"samtools mpileup -Q {Q} -q {q} -l {bed} -f {gsplit}/{chrom}.fa -b {pon_list_full} -r {chrom}"

    pon_matrix_file = os.path.join(matrix_path, f"{chrom}.pon")
    cmd = f"{pileup_cmd} | {mawk('cleanpileup')} | {mawk('pile2count')} | gzip  > {pon_matrix_file}.gz"
    run_cmd(cmd, show=True)
    return pon_matrix_file
