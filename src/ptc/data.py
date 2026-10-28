"""Datasets of the R package PTC, exported by scripts/export_from_R.R."""

import csv
import gzip
import lzma
from functools import lru_cache
from importlib import resources

import numpy as np
import pandas as pd


def _path(name):
    return resources.files("ptc") / "data" / name


@lru_cache(maxsize=1)
def tscan():
    """TargetScan 7.0 predicted interactions (columns miRNA, mRNA; miRBase v21)."""
    return pd.read_csv(_path("TScan.csv.gz"), dtype=str, keep_default_na=False)


@lru_cache(maxsize=1)
def context_scores():
    """TargetScan 7.0 conserved site context++ scores
    (columns Gene.Symbol, miRNA, context.score)."""
    ctx = pd.read_csv(_path("TS7.0_context_scores.csv.gz"), dtype=str,
                      keep_default_na=False)
    # Python's float() reads each value exactly as R does, with any pandas version
    ctx["context.score"] = np.array([float(x) for x in ctx["context.score"]])
    return ctx


@lru_cache(maxsize=1)
def ground_truth():
    """Experimentally confirmed interactions from miRTarBase 6.1, TarBase 7.0
    and miRWalk 2.0 (GroundT in R), as {mRNA: set of miRNAs}."""
    out = {}
    with gzip.open(_path("GroundT.csv.gz"), "rt", newline="") as f:
        reader = csv.reader(f)
        next(reader)
        for mr, mir in reader:
            out.setdefault(mr, set()).add(mir)
    return out


def _matrix(name):
    """An expression matrix saved exactly (raw 64-bit values, column by column)."""
    with gzip.open(_path(name + ".names.txt.gz"), "rt") as f:
        rows, cols = f.read().rstrip("\n").split("\n")
    rows, cols = rows.split("\t"), cols.split("\t")
    with lzma.open(_path(name + ".f64.xz")) as f:
        values = np.frombuffer(f.read(), dtype="<f8")
    return pd.DataFrame(values.reshape(len(cols), len(rows)).T, index=rows, columns=cols)


def load_tcga_brca():
    """Matched miRNA-mRNA bulk data from TCGA BRCA (TCGA_BRCAdata in R):
    518 miRNAs and 17403 mRNAs in 503 samples. Returns (miRNAs, mRNAs)."""
    return _matrix("TCGA_BRCA_miRNAs"), _matrix("TCGA_BRCA_mRNAs")


def load_single_cell():
    """Matched single-cell miRNA-mRNA data from GSE114071 (Wang et al., 2019),
    prepared as in inst/PTC.MainApplication.R of the R package (19 half-cells).
    Returns (miRNAs, mRNAs); samples are matched by row order."""
    return _matrix("SC_miRNAs"), _matrix("SC_mRNAs")


def load_wanderlust_time():
    """Pseudotime of the single-cell data obtained with Wanderlust (Bendall et
    al., 2014; k=4, l=2, ng=2, snn=0), one value per cell in the order of
    load_single_cell() (wandTime.csv in the R package)."""
    with open(_path("wandTime.csv")) as f:
        return np.array([float(x) for x in f.read().strip().split(",")])
