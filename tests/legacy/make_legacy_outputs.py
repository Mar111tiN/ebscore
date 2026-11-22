"""Produce the legacy "before" outputs from the unmodified ebscore main code (HEAD ef67010).

Run in the somvar-legacy env (py 3.9, pandas 1.5.2) against a main worktree:
    git -C <ebscore> worktree add --detach /tmp/ebscore-main main
    conda run -n somvar-legacy python tests/legacy/make_legacy_outputs.py /tmp/ebscore-main

Writes tests/data/legacy/:
- pon.matrix.tsv, pair.cpile.tsv   legacy cleanpileup | pile2count of the synthetic PON, cleanpileup -d of the pair
- ponmatrix/matrix/<chrom>.pon.gz  legacy PONmatrix per chromosome (+ PON_list.txt)
- ponmatrix/AB/<chrom>.AB.gz       legacy AB cache (ebcache.PONmatrix2AB_multi, zero cache disabled)
- eb/<chrom>.EB                    legacy run_ebscore (cleanpileup mode, PON matrix + AB cache)
- fits.tsv                         legacy ebcore.fit_bb on the PON count vectors at the variant positions
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

MAIN = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(MAIN / "code"))

from ebcache import PONmatrix2AB_multi  # noqa: E402
from ebcore import fit_bb  # noqa: E402
from ebrun import run_ebscore  # noqa: E402

HERE = Path(__file__).resolve().parent
SYN = HERE.parent / "data" / "synthetic"
OUT = HERE.parent / "data" / "legacy"
SHELL = MAIN / "shell"
CHROMS = ["chr7", "chr12"]
N_PON = 6


def sh(cmd):
    subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True)


def main():
    pon_dir = OUT / "ponmatrix"
    (pon_dir / "matrix").mkdir(parents=True, exist_ok=True)
    (pon_dir / "AB").mkdir(parents=True, exist_ok=True)
    (OUT / "eb").mkdir(parents=True, exist_ok=True)
    (pon_dir / "PON_list.txt").write_text("".join(f"N{i:02d}.bam\n" for i in range(1, N_PON + 1)))

    # mawk tools (legacy copies, always via bash)
    sh(f"bash {SHELL}/cleanpileup.mawk < {SYN}/pon.mpileup | bash {SHELL}/pile2count.mawk > {OUT}/pon.matrix.tsv")
    sh(f"bash {SHELL}/cleanpileup.mawk -d < {SYN}/pair.mpileup > {OUT}/pair.cpile.tsv")
    for chrom in CHROMS:
        sh(
            f"mawk -F'\\t' '$1 == \"{chrom}\"' {SYN}/pon.mpileup | bash {SHELL}/cleanpileup.mawk "
            f"| bash {SHELL}/pile2count.mawk | gzip -n > {pon_dir}/matrix/{chrom}.pon.gz"
        )
    pair_gz = OUT / "pair.cpile.tsv.gz"
    sh(f"gzip -n -c {OUT}/pair.cpile.tsv > {pair_gz}")

    with tempfile.TemporaryDirectory() as tmp:
        # The legacy Python runs the tools by path; their #!/bin/sh is dash here (Ubuntu), which silently breaks the
        # bash-only option parsing (AGENTS.md trap). The HPC's /bin/sh was bash, so run copies with a bash shebang.
        bash_shell = Path(tmp) / "shell"
        bash_shell.mkdir()
        for tool in SHELL.glob("*.mawk"):
            lines = tool.read_text().split("\n")
            lines[0] = "#!/bin/bash"
            (bash_shell / tool.name).write_text("\n".join(lines))
            (bash_shell / tool.name).chmod(0o755)
        config = dict(
            debug=False, threads=2, use_cache=True, mawk_path=str(bash_shell), temp_dir=tmp, pon_path=str(pon_dir),
            zero_path=os.path.join(tmp, "zero"), genome_split="", MAPQ=20, Q=25, fit_pen=0.5, chunksize=5000,
            ZDfactor=13, min_zt=10**9,
        )
        # legacy bug: zerocache.save_next_zero lists zero_path before creating it (zerocache.py:17 vs :46)
        os.makedirs(config["zero_path"], exist_ok=True)
        # AB cache
        for chrom in CHROMS:
            gen = pd.read_csv(pon_dir / "matrix" / f"{chrom}.pon.gz", sep="\t", chunksize=1000)
            ab = PONmatrix2AB_multi(gen, dict(config))
            ab.to_csv(pon_dir / "AB" / f"{chrom}.AB.gz", sep="\t", index=False, compression="gzip")

        # EB (legacy cannot read MNVs: drop them from its input)
        var = pd.read_csv(SYN / "variants.tsv", sep="\t")
        var = var.query("scenario != 'mnv'").loc[:, ["Chr", "Start", "End", "Ref", "Alt"]]
        mut_file = Path(tmp) / "variants.csv"
        var.sort_values(["Chr", "Start"]).to_csv(mut_file, sep="\t", index=False)
        for chrom in CHROMS:
            run_ebscore(
                str(mut_file), cleanpileup=str(pair_gz), output_file=str(OUT / "eb" / f"{chrom}.EB"),
                pon_list="PON_list.txt", chrom=chrom, config=dict(config),
            )

    # direct fits on the PON vectors at the variant positions (alt letter, both strands)
    matrix = pd.concat(
        [pd.read_csv(pon_dir / "matrix" / f"{c}.pon.gz", sep="\t", dtype=str) for c in CHROMS]
    ).set_index(["Chr", "Start"])
    rows = []
    var = pd.read_csv(SYN / "variants.tsv", sep="\t", dtype=str)
    for _, v in var.iterrows():
        start = int(v["Start"])
        if v["Ref"] == "-":
            letter, pos = "I", start
        elif v["Alt"] == "-":
            letter, pos = "D", start - 1
        elif len(v["Alt"]) == 1:
            letter, pos = v["Alt"], start
        else:
            continue
        key = (v["Chr"], str(pos))
        if key not in matrix.index:
            continue
        alts, deps = matrix.loc[key, letter].split("="), matrix.loc[key, "Depth"].split("=")
        for strand in (0, 1):
            d, a = deps[strand], alts[strand]
            cm = np.transpose([d.split("|"), a.split("|")]).astype(int)
            rows.append((d, a, fit_bb(cm)))
    # edge cases: all-zero alts, one PON sample, zero depth in some samples
    for d, a in [("30|40|50|60|70|80", "0|0|0|0|0|0"), ("100", "3"), ("0|0|25|0|0|10", "0|0|1|0|0|0")]:
        cm = np.transpose([d.split("|"), a.split("|")]).astype(int)
        rows.append((d, a, fit_bb(cm)))
    pd.DataFrame(rows, columns=["depths", "alts", "AB"]).drop_duplicates().to_csv(
        OUT / "fits.tsv", sep="\t", index=False
    )
    print("legacy outputs written to", OUT)


if __name__ == "__main__":
    main()
