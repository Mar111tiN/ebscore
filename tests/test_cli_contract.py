"""CLI contract (plan C1): the subcommands and flags that the kumo-variant-calling-nf modules call must not drift.

tests/data/cli_contract.json is the snapshot (made by tests/make_cli_contract.py from ebscore.cli.build_parser()).
"""

import argparse
import json
import shlex
import shutil
import types

import pytest

import ebscore
from ebscore.cli import build_parser
from make_cli_contract import (CONTRACT, SRC_PKG, TOP_LEVEL_MODULES, VERSION_MASK, affected_modules, build_contract,
                               contract_diff, render, source_note, stale_install)


def mismatch_message(expected: dict, actual: dict, module=ebscore) -> str:
    lines = [f"the ebscore CLI no longer matches its contract {CONTRACT}:",
             *(f"  {c}" for c in contract_diff(expected, actual) or ["(see the rendered JSON)"]),
             *source_note(module)]
    if not stale_install(module):
        lines += [
            "If the change is deliberate: rerun `python tests/make_cli_contract.py` and commit the JSON with the change,",
            "and review the kumo-variant-calling-nf modules that call the changed flags (design section 3 process table):",
            *(f"  {m}" for m in affected_modules(expected, actual)),
        ]
    return "\n".join(lines)


def test_cli_contract_unchanged():
    expected = json.loads(CONTRACT.read_text())
    actual = build_contract()
    # as canonical JSON text, not `==`: in Python 1 == 1.0 == True, in the snapshot they differ
    assert render(actual) == render(expected), mismatch_message(expected, actual)


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


def exact_flags(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """`parser` and all its subcommand parsers without prefix matching (argparse's allow_abbrev).

    The real CLI accepts `--vcf-annotation` for a flag renamed to `--vcf-annotation-file`. A module call must use the
    exact flags, else a stale flag keeps parsing until a later option with the same prefix makes it ambiguous.
    """
    todo = [parser]
    while todo:
        p = todo.pop()
        p.allow_abbrev = False
        todo += [s for a in p._actions if isinstance(a, argparse._SubParsersAction) for s in a.choices.values()]
    return parser


@pytest.mark.parametrize("module", MODULE_CALLS)
def test_module_calls_parse(module):
    argv = shlex.split(MODULE_CALLS[module])
    args = exact_flags(build_parser()).parse_args(argv)
    assert args.command == argv[0]


def test_version_output_matches_contract(capsys):
    """versions.yml parses `ebscore --version`: the printed line is the snapshot's template with the version."""
    template = json.loads(CONTRACT.read_text())["options"]["--version"]["version"]
    assert template == f"ebscore {VERSION_MASK}"
    with pytest.raises(SystemExit) as exit_:
        build_parser().parse_args(["--version"])
    assert exit_.value.code == 0
    assert capsys.readouterr().out == template.replace(VERSION_MASK, ebscore.__version__) + "\n"


# ---------------------------------------------------------------- the generator itself, on a toy parser
def toy_parser(count_default=0, version="toy 1.2", tag_flag="--tag"):
    p = argparse.ArgumentParser(prog="toy")
    p.add_argument("--version", action="version", version=version)
    sub = p.add_subparsers(dest="command", required=True)
    s = sub.add_parser("run")
    s.add_argument("infile")
    s.add_argument("-n", "--count", action="count", default=count_default)
    s.add_argument("--mode", action="store_const", const="fast")
    s.add_argument(tag_flag, dest="tag", action="extend", nargs="+", default=["a"])
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
    assert c["options"]["--version"]["version"] == "toy 1.2"  # no version given: nothing masked
    assert json.loads(render(c)) == c


def test_version_template_masks_the_version():
    """A changed --version format is a contract change (versions.yml parses it), a version bump is not."""
    v12 = build_contract(toy_parser(version="toy 1.2"), version="1.2")
    assert v12["options"]["--version"]["version"] == f"toy {VERSION_MASK}"
    assert contract_diff(v12, build_contract(toy_parser(version="toy 1.3"), version="1.3")) == []
    changed = build_contract(toy_parser(version="toy v1.2"), version="1.2")
    assert contract_diff(v12, changed) == [f'ebscore --version version: "toy {VERSION_MASK}" -> "toy v{VERSION_MASK}"']
    assert affected_modules(v12, changed) == [TOP_LEVEL_MODULES]


@pytest.mark.parametrize("new_default, shown", [(0.0, "0 -> 0.0"), (False, "0 -> false")])
def test_contract_diff_is_type_exact(new_default, shown):
    """0 == 0.0 == False in Python; in the snapshot (and for a Nextflow module) they are different defaults."""
    old, new = build_contract(toy_parser()), build_contract(toy_parser(count_default=new_default))
    assert old == new  # what a plain `==` comparison would have let through
    assert render(old) != render(new)
    assert contract_diff(old, new) == [f"ebscore run --count default: {shown}"]
    assert affected_modules(old, new) == ["`run`: no module calls it yet"]


def test_exact_flags_reject_a_renamed_flag():
    """A flag renamed to a longer one (--tag -> --tag-name) must fail a module call that still uses the old flag."""
    argv = ["run", "in.txt", "--tag", "b"]
    renamed = toy_parser(tag_flag="--tag-name")
    assert renamed.parse_args(argv).tag == ["a", "b"]  # argparse's prefix matching hides the rename
    with pytest.raises(SystemExit):
        exact_flags(toy_parser(tag_flag="--tag-name")).parse_args(argv)
    assert exact_flags(toy_parser()).parse_args(argv).tag == ["a", "b"]


def test_stale_install_is_named(tmp_path):
    """A test that imports another, different ebscore than src/ says so instead of 'rerun the script'."""
    copy = tmp_path / "site-packages" / "ebscore"
    shutil.copytree(SRC_PKG, copy, ignore=shutil.ignore_patterns("__pycache__"))
    installed = types.SimpleNamespace(__file__=str(copy / "__init__.py"), __version__="0.9")
    checkout = types.SimpleNamespace(__file__=str(SRC_PKG / "__init__.py"), __version__="1.0")
    assert not stale_install(checkout)
    assert not stale_install(installed)  # a non-editable install of this checkout (CI, Docker test stage)
    (copy / "cli.py").write_text((copy / "cli.py").read_text() + "\n# older\n")
    assert stale_install(installed)
    expected = build_contract()
    message = mismatch_message(expected, expected, module=installed)
    assert f"(ebscore 0.9 from {copy.resolve()})" in message
    assert "stale install" in message and "PYTHONPATH=src" in message
    assert "rerun `python tests/make_cli_contract.py`" not in message
    message = mismatch_message(expected, expected, module=checkout)
    assert f"(ebscore 1.0 from {SRC_PKG.resolve()})" in message
    assert "rerun `python tests/make_cli_contract.py`" in message and "stale" not in message


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
