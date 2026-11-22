"""mawk tools: byte parity with the legacy tools and the documented fixes."""

import subprocess

import pytest

from conftest import DATA, LEGACY, SYN, sh
from ebscore.tools import mawk_tool, tool


def run_tool(cmd, stdin_text):
    return subprocess.run(["bash", "-o", "pipefail", "-c", cmd], input=stdin_text, capture_output=True,
                          text=True, check=True).stdout


def test_fixture_cleanpileup():
    out = sh(f"{tool('cleanpileup')} < {DATA / 'mawk/test.pile'}").stdout
    assert out == (DATA / "mawk/test.cpile").read_text()


def test_fixture_pile2count():
    out = sh(f"{tool('pile2count')} < {DATA / 'mawk/test.cpile'}").stdout
    assert out == (DATA / "mawk/test.bamcount").read_text()


def test_synthetic_pon_matrix_matches_legacy():
    out = sh(f"{tool('cleanpileup')} < {SYN / 'pon.mpileup'} | {tool('pile2count')}").stdout
    assert out == (LEGACY / "pon.matrix.tsv").read_text()


def test_synthetic_paired_cleanpileup_matches_legacy():
    out = sh(f"{tool('cleanpileup', '-d')} < {SYN / 'pair.mpileup'}").stdout
    assert out == (LEGACY / "pair.cpile.tsv").read_text()


def test_cleanpileup_n_in_indel_is_forward_strand():
    # legacy classified +2NN as a reverse-strand insertion (i); N is an upper-case (forward) base
    line = "chr1\t10\tA\t3\t.+2NN,-1t.\tIII\n"
    out = run_tool(tool("cleanpileup"), line).splitlines()[1]
    assert out.split("\t")[3] == "Id."


def test_cleanpileup_use_samples_is_numeric():
    # legacy compared -s as a string, so "-s 2" was ignored with >= 10 samples ("2" > "10")
    fields = "\t".join(["1\t.\tI"] * 11)
    out = run_tool(tool("cleanpileup", "-s", "2"), f"chr1\t5\tA\t{fields}\n").splitlines()
    assert out[0] == "Chr\tStart\tRef\tRead1\tRead2"


def test_cleanpileup_start_markers_with_special_mapq_chars():
    # ^$ ^+ ^- ^5 are read starts with MAPQ chars that look like pileup syntax
    line = "chr1\t10\tA\t5\t^$.^+,^-A^5c.$\tIIIII\n"
    out = run_tool(tool("cleanpileup"), line).splitlines()[1]
    assert out.split("\t")[3] == ".,Ac."


def test_pile2count_counts_dots_like_mawk_in_gawk_too():
    text = "Chr\tStart\tRef\tRead1\nchr1\t1\tA\t..,,AaGg.\n"
    expected = "Chr\tStart\tRef\tA\tG\tC\tT\tI\tD\tDepth\nchr1\t1\tA\t1=1\t1=1\t0=0\t0=0\t0=0\t0=0\t5=4\n"
    assert run_tool(tool("pile2count"), text) == expected
    gawk = subprocess.run(["bash", "-c", "command -v gawk"], capture_output=True, text=True).stdout.strip()
    if gawk:
        script = open(mawk_tool("pile2count")).read().replace("mawk '", "gawk '", 1)
        assert run_tool(script, text) == expected


def test_filterpos_keeps_header_and_listed_positions(tmp_path):
    pos = tmp_path / "pos.tsv"
    pos.write_text("chr2\t7\nchr1\t3\n")
    text = "Chr\tStart\tX\nchr1\t2\ta\nchr1\t3\tb\nchr2\t3\tc\nchr2\t7\td\n"
    assert run_tool(tool("filterpos", pos), text) == "Chr\tStart\tX\nchr1\t3\tb\nchr2\t7\td\n"


@pytest.mark.parametrize("name", ["cleanpileup", "pile2count", "filterpos"])
def test_tools_have_bash_shebang(name):
    assert open(mawk_tool(name)).readline().strip() == "#!/usr/bin/env bash"
