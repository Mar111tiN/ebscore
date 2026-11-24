# syntax=docker/dockerfile:1
# ebscore: EB score of somatic variants against a pileup PON (PONmatrix). CLI only: `ebscore <command> ...`.
# Build from the repo root:
#   docker build --target test -t ebscore:test .        # runs pytest, not shipped
#   docker build -t ebscore:<version> .                 # runtime image (last stage)

# ---- build: solve the env once at its final prefix, install the package ----
FROM mambaorg/micromamba:2.10.0-debian13-slim AS build
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
 && micromamba clean -a -y
COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml README.md /tmp/ebscore/
COPY --chown=$MAMBA_USER:$MAMBA_USER src /tmp/ebscore/src
# pip is needed only for this install. It is removed afterwards: its vendored libraries (pip/_vendor/vendor.txt:
# msgpack, urllib3, setuptools/pkg_resources) carried 5 fixed HIGH CVEs in the first CI Trivy gate (2026-11-23), and
# nothing calls pip at runtime. The lock file is written after the removal, so it lists what ships.
RUN /opt/conda/bin/python -m pip install --no-cache-dir --no-deps --no-build-isolation /tmp/ebscore \
 && /opt/conda/bin/ebscore --version \
 && micromamba remove -y -n base --force pip \
 && micromamba env export -n base > /opt/conda/env.lock.yml \
 && micromamba clean -a -y \
 && /opt/conda/bin/ebscore --version

# ---- test: docker build --target test . ; never shipped ----
FROM build AS test
COPY --chown=$MAMBA_USER:$MAMBA_USER tests /tmp/ebscore/tests
WORKDIR /tmp/ebscore
RUN /opt/conda/bin/python -m pytest -q -p no:cacheprovider tests

# ---- slim: the tested env minus what nothing runs (plan 3b, 2026-11-24; about half the size) ----
# docker/slim_env.py: force-removes the listed packages that are installed (pytest and its helpers, tk), deletes test
# suites, headers, static libs, docs, translations and __pycache__, and strips the ELF files (--strip-unneeded keeps the
# dynamic symbols). Checked on a copy of the env: full test suite and tier-A smoke unchanged.
FROM build AS slim
USER root
RUN apt-get update \
 && apt-get install -y --no-install-recommends binutils \
 && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER docker/slim_env.py /tmp/slim_env.py
RUN /opt/conda/bin/python -I /tmp/slim_env.py /opt/conda --remove pytest pygments packaging pluggy colorama \
        exceptiongroup typing_extensions tomli iniconfig tk pip \
 && micromamba env export -n base > /opt/conda/env.lock.yml \
 && micromamba clean -a -y \
 && /opt/conda/bin/ebscore --version

# ---- runtime: the shipped image (must stay the last stage) ----
FROM debian:13.7-slim AS runtime
RUN apt-get update \
 && apt-get install -y --no-install-recommends bash procps \
 && rm -rf /var/lib/apt/lists/*
COPY --from=slim /opt/conda /opt/conda
# no __pycache__ in the env (slim stage) and /opt/conda is not writable for the task user: never try to write .pyc
ENV PATH=/opt/conda/bin:$PATH \
    PYTHONDONTWRITEBYTECODE=1
RUN useradd -m ebscore_user
USER ebscore_user
WORKDIR /home/ebscore_user
