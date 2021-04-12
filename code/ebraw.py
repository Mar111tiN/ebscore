import pandas as pd
from io import StringIO
import os
import sys
from subprocess import PIPE, run
from ebutils import get_pon
from script_utils import show_output, show_command, run_cmd

# here come all the functions involved in converting raw data (bam or pileup) to dataframes


def tumor2matrix(mut_file, bam="", pileup="", pon_list="", chrom="", EBconfig={}):
    """
    converts the bam_file and the Pon-bams into a matrix for EB computations
    """

    # PARAMS
    # mawk tool unwrapper
    def mawk(tool):
        return os.path.join(EBconfig["mawk_path"], f"{tool}.mawk")

    gsplit = EBconfig["genome_split"]
    pon_path = EBconfig["pon_path"]
    q = EBconfig["MAPQ"]
    Q = EBconfig["Q"]
    base_name = os.path.splitext(os.path.basename(mut_file))[0]
    use_cache = EBconfig["use_cache"]
    # check the temp_folder with default pon_path/temp
    temp_folder = EBconfig.get("temp_dir", os.path.join(pon_path, "temp"))
    if not os.path.isdir(temp_folder):
        os.mkdir(temp_folder)

    if bam:
        tumor_file = bam
        tumor_type = "bam"
        # create the bed file needed for pileup from bam-file
        # could be pythonized
        bed_file = os.path.join(
            temp_folder,
            f"{base_name}_{chrom}.bed",
        )
        run(
            f"{mawk('csv2bed')} {mut_file} {chrom} > {bed_file}", check=True, shell=True
        )
        # run(f"{mawk('csv2bed')} {mut_file} > {bed_file}", check=True, shell=True)
        tumor_cmd = f"samtools mpileup -Q {Q} -q {q} -l {bed_file} -f {gsplit}/{chrom}.fa -r {chrom}"
    else:
        if pileup:
            tumor_file = pileup
            tumor_type = "pileup"
            if tumor_file.endswith(".gz"):
                tumor_cmd = f"gunzip < {tumor_file} | {mawk('cleanpileup')} -s 1"
            else:
                tumor_cmd = f"cat {tumor_file} | {mawk('cleanpileup')} -s 1"
        else:
            show_output(
                "Neither bam nor pileup file is given! Either is required!",
                color="warning",
                multi=False,
            )
    # find out whether bam/pileup is contained in pon_list
    # create the pon_list
    pon_list = os.path.join(EBconfig["pon_path"], pon_list)

    pon_df, pos, removed_ponbam = get_pon(
        tumor_file, pon_list, pon_path=pon_path, prepend_bam=True
    )

    if pos:
        show_output(
            f"{removed_ponbam} from PON list matches {tumor_file} and is removed from pon_list. If not desired, change names",
            color="warning",
        )

    # build the filter command before other commands
    filter_cmd = f"{mawk('tumor2matrix')} {mut_file} {chrom}"
    ##### use PON cache
    if use_cache:
        # check if gzip matrix file is there
        matrix_file = os.path.join(EBconfig["pon_path"], f"matrix/{chrom}.pon")
        # check existence of matrix file
        # if os.path.isfile(f"{matrix_file}.gz"):
        #     tumor_cmd = (
        #         f"mkfifo ponzip; gunzip < {matrix_file}.gz > ponzip; {tumor_cmd}"
        #     )
        #     filter_cmd += f" -P ponzip"
        # el
        if os.path.isfile(matrix_file):
            filter_cmd += f" -P {matrix_file}"

        else:  # did not find the matrix file
            show_output(
                f"PON matrix {matrix_file} not found. Falling back to PON bam_files",
                color="warning",
                multi=False,
            )
            use_cache = False

    if use_cache:
        if tumor_type == "bam":
            # just add the bam_file as
            tumor_cmd += f" {tumor_file} | {mawk('cleanpileup')}"
        # add sample exclusion
        if pos:
            filter_cmd += f" -x {pos}"
    else:
        if tumor_type == "bam":
            # make filename for pon_bam
            pon_bam = os.path.join(temp_folder, f"{base_name}_{chrom}.pon")
            pon_df.to_csv(pon_bam, index=False, header=None)
            tumor_cmd += f" -b {pon_bam} | {mawk('cleanpileup')}"
        else:
            show_command(
                "If PON cache is not used, tumor pileup cannot be used!",
                color="warning",
                multi=False,
            )
            return None
    # print(pon_df)

    cmd = f"{tumor_cmd} | {mawk('pile2count')} | {filter_cmd}"
    show_command(cmd, multi=False)
    matrix_df = pd.read_csv(
        StringIO(run(cmd, stdout=PIPE, check=True, shell=True).stdout.decode("utf-8")),
        sep="\t",
    )
    return matrix_df
