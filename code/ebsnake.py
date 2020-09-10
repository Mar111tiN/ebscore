import os
from ebrun import run_eb


def main(s):
    '''
    wrapped into function lest module be called by each subprocess
    '''

    w = s.wildcards
    p = s.params
    static_path = s.config['paths']['mystatic']
    EBconfig = s.config['EBFilter']
    pon_list = os.path.join(static_path, EBconfig['pon_list'])

    # load configs into EBconfig
    EBconfig = {
        "cleanpileup": "../shell/cleanpileup.mawk",
        "makeponlist": "../shell/makeponlist.sh",
        "csv2bed":"../shell/csv2bed.mawk",
        "pon2cols": "../shell/pon2cols.mawk",
        "pile2count": "../shell/pile2count2.mawk",
        "filterVar": "../shell/filterVar.mawk",
        "pon2tumor": "../shell/pon2tumor.mawk",
        "pon_path": pon_path,
        "genome_split": "/Users/mahtin/Dropbox/Icke/Work/static/genome/gatk/hg38/split",
        "MAPQ": 20,
        "Q": 25,
        "fit_pen": 0.5,
        "count_dict": {0:"alt+", 1:"alt-", 2:"depth+", 3:"depth-"}
    }


    run_eb(
        table=s.input.table,
        tumor_bam=s.input.tumor_bam,
        output=s.output,
        pon_list=pon_list,
        chrom=w.chrom,
        log=s.log,
        threads=s.threads,
        EBparams=EBconfig['params'],
        full_output=EBconfig['full_pon_output'],
        # import the scripts
        cleanpileup=p.cleanpileup,
        csv2bed=p.csv2bed,
        pon2cols=p.pon2cols,
        pile2count=p.pile2count,
        matrix2EBinput=p.matrix2EBinput,
        makeponlist=p.makeponlist
    )


if __name__ == "__main__":
    main(snakemake)
