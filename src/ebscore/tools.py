"""Shell helpers: locate the packaged mawk tools and run bash pipelines."""

from __future__ import annotations

import io
import logging
import shlex
import subprocess
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


def first_line(path) -> str:
    """First line of a plain or gzipped text file ("" if empty)."""
    import gzip

    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        return fh.readline().rstrip("\n")


def bed_chroms(bed) -> list[str]:
    """Chromosomes of a BED file in order of first appearance (header/track lines skipped)."""
    seen: dict[str, None] = {}
    with open(bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            seen.setdefault(line.split("\t", 1)[0], None)
    return list(seen)


def ensure_dir(path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def setup_logging(verbose: bool = False, quiet: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.WARNING if quiet else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s %(name)s %(levelname)s: %(message)s", datefmt="%H:%M:%S")
