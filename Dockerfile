# syntax=docker/dockerfile:1
# ebscore: EB score of somatic variants against a pileup PON (PONmatrix). CLI only: `ebscore <command> ...`.
# Build from the repo root:
#   docker build --target test -t ebscore:test .        # runs pytest, not shipped
#   docker build -t ebscore:<version> .                 # runtime image (last stage)

# ---- build: solve the env once at its final prefix, install the package ----
FROM mambaorg/micromamba:2.10.0-debian13-slim AS build
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
 && micromamba env export -n base > /opt/conda/env.lock.yml \
 && micromamba clean -a -y
COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml README.md /tmp/ebscore/
COPY --chown=$MAMBA_USER:$MAMBA_USER src /tmp/ebscore/src
RUN /opt/conda/bin/python -m pip install --no-cache-dir --no-deps --no-build-isolation /tmp/ebscore \
 && /opt/conda/bin/ebscore --version

# ---- test: docker build --target test . ; never shipped ----
FROM build AS test
RUN micromamba install -y -n base -c conda-forge pytest && micromamba clean -a -y
COPY --chown=$MAMBA_USER:$MAMBA_USER tests /tmp/ebscore/tests
WORKDIR /tmp/ebscore
RUN /opt/conda/bin/python -m pytest -q -p no:cacheprovider tests

# ---- runtime: the shipped image (must stay the last stage) ----
FROM debian:13.7-slim AS runtime
RUN apt-get update \
 && apt-get install -y --no-install-recommends bash procps \
 && rm -rf /var/lib/apt/lists/*
COPY --from=build /opt/conda /opt/conda
ENV PATH=/opt/conda/bin:$PATH
RUN useradd -m ebscore_user
USER ebscore_user
WORKDIR /home/ebscore_user
