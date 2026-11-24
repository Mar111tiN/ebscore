"""CLI contract (plan C1): the subcommands and flags that the kumo-variant-calling-nf modules call must not drift.

tests/data/cli_contract.json is the snapshot (made by tests/make_cli_contract.py from ebscore.cli.build_parser()).
"""

import argparse
import json
import shlex

import pytest

from ebscore.cli import build_parser
from make_cli_contract import CONTRACT, affected_modules, build_contract, contract_diff, render


def test_cli_contract_unchanged():
    expected = json.loads(CONTRACT.read_text())
    actual = build_contract()
    changes = contract_diff(expected, actual)
    assert actual == expected, "\n".join([
        f"the ebscore CLI no longer matches its contract {CONTRACT}:",
        *(f"  {c}" for c in changes),
        "If the change is deliberate: rerun `python tests/make_cli_contract.py` and commit the JSON with the change,",
        "and review the kumo-variant-calling-nf modules that call the changed flags (design section 3 process table):",
        *(f"  {m}" for m in affected_modules(expected, actual)),
    ])


def test_contract_file_is_canonical():
    """The JSON is the script's own output (one line per option), so a rerun gives a minimal diff."""
    text = CONTRACT.read_text()
    assert text == render(json.loads(text)), "rerun `python tests/make_cli_contract.py` instead of editing the JSON"


# The command lines of design section 3 (process table, 2026-11-23), with placeholders for the Nextflow variables.
# The snapshot catches every change; these show directly that the module calls still parse.
MODULE_CALLS = {
    "EBSCORE_PONMATRIX": "pon-matrix --bam-list bam_list.txt --fasta genome.fa --bed kit.padded.bed -o ponmatrix "
                         "--chrom chr21 --min-mapq 20 --min-baseq 25 -t 2",
    "EBSCORE_PILEUP paired": "pileup --normal N.bam --tumor T.bam --fasta genome.fa --bed kit.padded.bed "
                             "-o P1.pileup.gz --tmpdir tmp --min-mapq 20 --min-baseq 25 -t 8",
    "EBSCORE_PILEUP tumor-only": "pileup --tumor T.bam --fasta genome.fa --bed kit.padded.bed -o T1.pileup.gz "
                                 "--tmpdir tmp --min-mapq 20 --min-baseq 25 -t 8",
    "EBSCORE_UNION": "vcf2table --vcf P1.mutect2.vcf.gz -o P1.mutect2.tsv",
    "EBSCORE_PONSUBSET run": "pon-subset --ponmatrix ponmatrix -o run --variants unions/P1.union.tsv "
                             "--variants unions/P2.union.tsv -t 8",
    "EBSCORE_PONSUBSET pair": "pon-subset --ponmatrix run -o pairs/P1 --variants unions/P1.union.tsv -t 8",
    "EBSCORE_SCORE pileup": "score --variants P1.union.tsv --ponmatrix subset --pileup P1.pileup.gz -o P1.EB.tsv "
                            "--vcf-annotation P1.EB.annot.tsv --exclude-sample P1_N.bam --min-mapq 20 --min-baseq 25 "
                            "--tmpdir tmp -t 8",
    "EBSCORE_SCORE tumor-bam": "score --variants T1.union.tsv --ponmatrix subset --tumor-bam T.bam --fasta genome.fa "
                               "-o T1.EB.tsv --vcf-annotation T1.EB.annot.tsv --min-mapq 20 --min-baseq 25 "
                               "--tmpdir tmp -t 8",
}


@pytest.mark.parametrize("module", MODULE_CALLS)
def test_module_calls_parse(module):
    argv = shlex.split(MODULE_CALLS[module])
    args = build_parser().parse_args(argv)
    assert args.command == argv[0]


# ---------------------------------------------------------------- the generator itself, on a toy parser
def toy_parser(count_default=0):
    p = argparse.ArgumentParser(prog="toy")
    sub = p.add_subparsers(dest="command", required=True)
    s = sub.add_parser("run")
    s.add_argument("infile")
    s.add_argument("-n", "--count", action="count", default=count_default)
    s.add_argument("--mode", action="store_const", const="fast")
    s.add_argument("--tag", action="extend", nargs="+", default=["a"])
    s.add_argument("--level", type=float, choices=[0.5, 1.0], default=0.5)
    s.add_argument("--clock", default=object())  # a computed default: not part of the contract
    g = s.add_mutually_exclusive_group()
    g.add_argument("--x", action="store_true")
    g.add_argument("-y", action="store_false")
    return p


def test_contract_of_toy_parser():
    c = build_contract(toy_parser())
    assert c["command_required"] is True and list(c["commands"]) == ["run"]
    run = c["commands"]["run"]
    assert run["option_strings"] == ["--clock", "--count", "--help", "--level", "--mode", "--tag", "--x", "-h", "-n",
                                     "-y"]
    assert run["positionals"] == [{"dest": "infile", "flags": [], "required": True, "action": "store", "nargs": None,
                                   "type": None, "choices": None, "default": None}]
    o = run["options"]
    assert {k: v["action"] for k, v in o.items()} == {
        "--help": "help", "--count": "count", "--mode": "store_const", "--tag": "extend", "--level": "store",
        "--clock": "store", "--x": "store_true", "-y": "store_false"}
    assert o["--count"]["flags"] == ["--count", "-n"] and o["--count"]["default"] == 0
    assert o["--mode"]["const"] == "fast"
    assert o["--tag"]["nargs"] == "+" and o["--tag"]["default"] == ["a"]
    assert o["--level"]["type"] == "float" and o["--level"]["choices"] == [0.5, 1.0]
    assert "default" not in o["--clock"]
    assert run["mutually_exclusive"] == [{"required": False, "options": ["--x", "-y"]}]
    assert json.loads(render(c)) == c


def test_contract_diff_names_the_change():
    old, new = build_contract(toy_parser()), build_contract(toy_parser(count_default=3))
    assert contract_diff(old, old) == []
    assert contract_diff(old, new) == ["ebscore run --count default: 0 -> 3"]
    del new["commands"]["run"]["options"]["--mode"]
    new["commands"]["run"]["option_strings"].remove("--mode")
    new["commands"]["run"]["options"]["--new"] = {"action": "store_true"}
    new["commands"]["run"]["option_strings"].append("--new")
    assert contract_diff(old, new) == [
        "ebscore run option_strings: added --new; removed --mode",
        "ebscore run --count default: 0 -> 3",
        "ebscore run --mode: removed",
        'ebscore run --new: added {"action": "store_true"}',
    ]
    assert affected_modules(old, new) == ["`run`: no module calls it yet"]
