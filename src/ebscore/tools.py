"""Shell helpers: locate the packaged mawk tools and run bash pipelines."""

from __future__ import annotations

import gzip
import io
import logging
import os
import shlex
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from importlib.resources import files
from pathlib import Path

import pandas as pd

log = logging.getLogger("ebscore")


def mawk_tool(name: str) -> str:
    """Path of a packaged mawk tool (``src/ebscore/mawk/<name>.mawk``)."""
    path = files("ebscore") / "mawk" / f"{name}.mawk"
    if not path.is_file():
        raise FileNotFoundError(f"mawk tool {name!r} not found at {path}")
    return str(path)


def tool(name: str, *args) -> str:
    """Shell snippet that runs a packaged mawk tool with bash (its option parser is bash-only)."""
    return " ".join(["bash", shlex.quote(mawk_tool(name)), *(shlex.quote(str(a)) for a in args)])


def q(path) -> str:
    return shlex.quote(str(path))


def cat_cmd(path) -> str:
    """Shell snippet that writes a plain or gzipped file to stdout."""
    return f"gzip -dc {q(path)}" if str(path).endswith(".gz") else f"cat {q(path)}"


def find_gz(path_base: Path):
    """``x.gz`` if it exists, else ``x`` without .gz, else None (PONmatrix files may be plain)."""
    path_base = Path(path_base)
    if path_base.is_file():
        return path_base
    plain = path_base.with_suffix("") if path_base.suffix == ".gz" else None
    return plain if plain and plain.is_file() else None


def write_text_atomic(path, text: str) -> None:
    """Write a file via a temp file + rename, so that an existing symlink at ``path`` is replaced, not followed."""
    path = Path(path)
    tmp = path.with_name(f".{path.name}.part")
    tmp.write_text(text)
    os.replace(tmp, path)


def run(cmd: str) -> None:
    """Run a shell pipeline with ``bash -o pipefail``; raise on any failing stage."""
    log.debug("$ %s", cmd)
    subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True)


def cmd2df(cmd: str) -> pd.DataFrame:
    """Run a shell pipeline and read its tab-separated stdout (with header) as strings.

    Returns an empty DataFrame when the pipeline prints nothing (e.g. no covered position).
    """
    log.debug("$ %s", cmd)
    out = subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True, stdout=subprocess.PIPE).stdout
    if not out.strip():
        return pd.DataFrame()
    return pd.read_csv(io.BytesIO(out), sep="\t", dtype=str, keep_default_na=False)


def open_text(path):
    """Open a plain or gzipped (detected by its magic bytes, not the name) text file for reading."""
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path)


def first_lines(path, n: int = 2) -> list[str]:
    """Up to ``n`` first lines of a plain or gzipped text file (without newlines)."""
    out = []
    with open_text(path) as fh:
        for line in fh:
            out.append(line.rstrip("\n"))
            if len(out) == n:
                break
    return out


def first_line(path) -> str:
    """First line of a plain or gzipped text file ("" if empty)."""
    lines = first_lines(path, 1)
    return lines[0] if lines else ""


def bed_chroms(bed) -> list[str]:
    """Chromosomes of a BED file (plain or gzipped, tab- or space-delimited) in order of first appearance.

    Header, track and browser lines are skipped.
    """
    seen: dict[str, None] = {}
    with open_text(bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            seen.setdefault(line.split()[0], None)
    return list(seen)


def bam_contigs(bam) -> list[str]:
    """Reference sequence names (@SQ SN) of a BAM/CRAM header, in header order."""
    header = subprocess.run(["samtools", "view", "-H", str(bam)], check=True, stdout=subprocess.PIPE, text=True).stdout
    return [
        field[3:]
        for line in header.splitlines()
        if line.startswith("@SQ\t")
        for field in line.split("\t")[1:]
        if field.startswith("SN:")
    ]


def unique(values) -> list:
    """Values in first-appearance order without duplicates (e.g. repeated --chrom)."""
    return list(dict.fromkeys(values or []))


def run_parallel(fn, items, threads: int) -> list:
    """``[fn(x) for x in items]`` on a thread pool, results in input order.

    On the first failure the items that have not started yet are cancelled, and the error is raised once the running
    ones are done.
    """
    items = list(items)
    if threads <= 1 or len(items) <= 1:
        return [fn(x) for x in items]
    with ThreadPoolExecutor(max_workers=threads) as pool:
        futures = [pool.submit(fn, x) for x in items]
        try:
            for fut in as_completed(futures):
                fut.result()
        except BaseException:
            for fut in futures:
                fut.cancel()
            raise
        return [fut.result() for fut in futures]


def ensure_dir(path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def setup_logging(verbose: bool = False, quiet: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.WARNING if quiet else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s %(name)s %(levelname)s: %(message)s", datefmt="%H:%M:%S")
