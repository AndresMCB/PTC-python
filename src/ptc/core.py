"""PTC: Pseudo-temporal Causality, a Python version of the R package PTC 1.1.2.

Each function follows the R function of the same name, and returns the same
results (see tests/). Expression data are pandas DataFrames with samples in
rows and genes in columns, as matrices in the R package.
"""

import math
import warnings
from dataclasses import dataclass, field
from itertools import combinations

import numpy as np
import pandas as pd

from . import data as ptc_data
from .decoupled_test import decoupled_test
from .mirbase import mirv21
from .rrandom import RRandom


@dataclass
class PTCResult:
    """Results of ptc(), as the list returned by PTC() in R.

    index: for each mRNA with at least one parent, the positions of the
        parents in its set of plausible parents (1-based, as in R).
    names: for each of those mRNAs, the names of the parent miRNAs.
    genes: names of those mRNAs.
    summary: inferred miRNA-mRNA interactions ranked by TargetScan context++
        score, with columns rank, miR, mR and Score.
    failed: mRNAs for which the invariance test could not be computed, with
        the error. As in R, they are in genes and names (with no parent name)
        but not in summary.
    """
    index: dict
    names: dict
    genes: list
    summary: pd.DataFrame
    failed: dict = field(default_factory=dict)


def ptc(miRNAs, mRNAs, VIM, nmiR=30, nmR=1500, ngrid=2, alpha=0.02,
        complements=True, explore_all=True, silent=True, TScan=None):
    """Estimates the causal miRNA parents of a set of nmR mRNAs (PTC() in R).

    miRNAs, mRNAs: DataFrames with matched samples in rows (same order) and
        miRNAs / mRNAs in columns.
    VIM: expression of the marker gene (or any pseudotime), one value per
        sample, used to order the samples.
    nmiR, nmR: number of miRNAs and mRNAs selected by largest MAD.
    ngrid: number of segments of the pseudotime series (environments).
    alpha: significance level of the invariance test.
    complements: compare each segment against its complement (True) or all
        segments pairwise (False).
    explore_all: return the union of all invariant sets (True) or the first
        invariant set found (False).
    silent: if False, print each evaluated and accepted set.
    TScan: DataFrame with columns miRNA and mRNA (miRBase v21 names); if
        None, TargetScan 7.0 is used.
    """
    mRNAs = mRNAs.copy()
    # for compatibility, mRNAs are transformed to uppercase
    mRNAs.columns = [str(c).upper() for c in mRNAs.columns]

    seq = ptime([miRNAs, mRNAs], VIM)
    sel = gene_sel(seq, nmiR=nmiR, nmR=nmR, TScan=TScan)
    d, columns = sel["d"], sel["columns"]

    outcome, failed = {}, {}
    rng = RRandom(1)
    for gene in sel["mRs"]:
        parents = sel["PParents"].get(gene)
        if not parents:
            continue
        if not silent:
            print(f" Current mRNA  = {{ {gene} }}")
        rng.set_seed(1)  # set.seed(1) before each gene, as in R
        try:
            outcome[gene] = test_invariance(
                d[:, columns[gene]], d[:, [columns[p] for p in parents]],
                ngrid=ngrid, alpha=alpha, explore_all=explore_all,
                silent=silent, complements=complements, rng=rng)
        except (np.linalg.LinAlgError, ValueError) as e:
            failed[gene] = str(e)
            warnings.warn(f"The invariance test failed for {gene}; it is reported "
                          f"without parents. Error: {e}", stacklevel=2)

    res = extract_parents(outcome, failed, sel["PParents"])
    summary = inter_list_to_matrix(res["names"])
    summary = rank_by_context(ptc_data.context_scores(), summary)
    return PTCResult(index=res["index"], names=res["names"], genes=res["genes"],
                     summary=summary, failed=failed)


def ptime(matched_data, VIM):
    """Orders the samples of [miRNAs, mRNAs] by VIM (PTC.ptime in R)."""
    order = np.argsort(np.asarray(VIM, dtype=float), kind="stable")
    return [m.iloc[order] for m in matched_data]


def gene_sel(seq_data, nmiR=30, nmR=1500, TScan=None):
    """Selects predictor candidates and targets (PTC.GeneSel in R).

    Returns a dict with d (numeric matrix, miRNAs then mRNAs in columns),
    columns (name -> column of d, first match as in R), miRs (miRBase v21
    names), mRs and PParents (plausible parents of each mRNA).
    """
    if not (isinstance(seq_data, (list, tuple)) and len(seq_data) == 2):
        raise ValueError("seq_data must be a list with 2 elements. A matrix with miRNAs "
                         "gene expression (1st element), and a matrix with mRNAs gene "
                         "expression (2nd element)")
    sel = get_data_by_mad(seq_data, nmiR, nmR)
    # change selected miRNA names to miRBase v21
    names = list(sel["names"])
    sel["miRs"] = mirv21(names[:nmiR])
    names[:nmiR] = sel["miRs"]
    sel["names"] = names
    columns = {}
    for i, n in enumerate(names):
        columns.setdefault(n, i)
    sel["columns"] = columns
    sel["PParents"] = find_pp(TScan, sel["miRs"], sel["mRs"])
    return sel


def get_data_by_mad(seq_data, nmiR, nmR):
    """The nmiR miRNAs and nmR mRNAs with largest MAD (getDatabyMAD in R)."""
    miRs, mRNAs = seq_data
    mir_idx = _select_by_mad(miRs.to_numpy(dtype=float), nmiR)
    mr_idx = _select_by_mad(mRNAs.to_numpy(dtype=float), nmR)
    mr_names = [str(mRNAs.columns[i]) for i in mr_idx]
    # remove duplicated mRNA names (keep the first)
    seen, keep = set(), []
    for k, n in enumerate(mr_names):
        if n not in seen:
            seen.add(n)
            keep.append(k)
    mr_idx = [mr_idx[k] for k in keep]
    mir_names = [str(miRs.columns[i]) for i in mir_idx]
    mr_names = [mr_names[k] for k in keep]
    d = np.hstack([miRs.to_numpy(dtype=float)[:, mir_idx],
                   mRNAs.to_numpy(dtype=float)[:, mr_idx]])
    return {"d": d, "names": mir_names + mr_names, "miRs": mir_names, "mRs": mr_names}


def _select_by_mad(x, value):
    """Columns of x with the largest MAD (FSbyMAD, cut.type = "topk")."""
    mads = np.array([_mad(x[:, j]) for j in range(x.shape[1])])
    if value > len(mads):
        value = len(mads)
        print("Warning: the feature selection number is beyond the original feature numnber")
    order = np.argsort(-mads, kind="stable")  # decreasing, ties in original order
    return list(order[:value])


def _median(x):
    """R's median(): the mean of the two middle values is computed as R's
    mean(), with extended precision."""
    s = np.sort(x)
    n = len(s)
    half = (n + 1) // 2
    if n % 2 == 1:
        return float(s[half - 1])
    return _r_mean(s[half - 1:half + 1])


def _mad(x):
    return 1.4826 * _median(np.abs(x - _median(x)))


def _r_mean(x):
    """R's mean(): sum in extended precision, then one correction step."""
    x = np.asarray(x, dtype=np.longdouble)
    n = len(x)
    s = x.sum() / n
    if math.isfinite(float(s)):
        s = s + (x - s).sum() / n
    return float(s)


def find_pp(TScan, miRs, mRs):
    """Plausible parents of each mRNA: the miRNAs that TargetScan predicts
    can bind it (PTC.findPP in R). Returns {mRNA: [miRNAs]}."""
    if TScan is None:
        TScan = ptc_data.tscan()
    targets = {}
    for mir, mr in zip(TScan.iloc[:, 0].astype(str), TScan.iloc[:, 1].astype(str)):
        targets.setdefault(mr, set()).add(mir)
    pparents = {}
    for mr in mRs:
        hits = [j for j in miRs if j in targets.get(mr, ())]
        pparents[mr] = hits if hits else None
    return pparents


def test_invariance(Y, X, ngrid=2, alpha=0.02, explore_all=True, silent=True,
                    complements=True, rng=None):
    """Indexes (1-based) of the parents of a target gene: the union of the
    invariant sets of predictors (PTC.TestInvariance in R).

    Y: expression of the target (n values); X: n x p expression of its
    plausible parents. rng: an RRandom generator (R's set.seed(1) if None).
    """
    if rng is None:
        rng = RRandom(1)
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float).reshape(-1, 1)
    n, p = X.shape
    # environments: ngrid consecutive segments of the pseudotime series
    grid = _r_seq_int(0, n, ngrid + 1)

    s_union = []
    aux = p
    rejected = True
    finish = False
    while not finish:
        sets = list(combinations(range(1, p + 1), aux))
        # remove all subsets of the current union so they are not evaluated again
        if s_union:
            known = {frozenset(c) for c in combinations(s_union, aux)}
            if not silent:
                print(f" (PTC) nSets to REMOVE in current level = {{ {len(known)} }}")
            sets = [s for s in sets if frozenset(s) not in known]
        if not silent:
            print(f" (PTC) nSets to explore in current level = {{ {len(sets)} }}")
        if not sets:
            raise ValueError("subscript out of bounds")  # as R's Sets[,1] with no sets
        for s in sets:
            if not silent:
                print(f" Current evaluated set  = {{ {', '.join(map(str, s))} }}")
            p_value = decoupled_test(X, Y, s, grid, complements, B=100, rng=rng)
            if p_value > alpha:
                rejected = False
                s_union += [k for k in s if k not in s_union]
                if not silent:
                    print(f" Accepted set  = {{ {', '.join(map(str, s))} }}")
                if not explore_all:
                    finish = True
                    break
            if len(s_union) == p:
                finish = True
                break
        aux -= 1
        if aux == 0:
            finish = True
    return [] if rejected else s_union


def _r_seq_int(start, end, length_out):
    """as.integer(seq(start, end, length.out = length_out)) as in R."""
    by = (end - start) / (length_out - 1)
    values = [start + i * by for i in range(length_out)]
    values[-1] = end
    return [int(v) for v in values]


def extract_parents(outcome, failed, predictors):
    """All mRNAs with at least one causal parent (Extract.Parents in R)."""
    index, names, genes = {}, {}, []
    for gene in predictors:
        if gene in failed:
            # as in R, a failed test is kept with no parent name (NA)
            index[gene], names[gene] = failed[gene], [None]
        elif outcome.get(gene):
            index[gene] = outcome[gene]
            names[gene] = [predictors[gene][k - 1] for k in outcome[gene]]
        else:
            continue
        genes.append(gene)
    return {"index": index, "names": names, "genes": genes}


def inter_list_to_matrix(inter_list):
    """miRNA-mRNA interactions as a two-column table (InterList.toMatrix)."""
    rows = [(mir, mr) for mr, mirs in inter_list.items()
            if mirs is not None and mirs != [None] for mir in mirs]
    return pd.DataFrame(rows, columns=["miR", "mR"])


def inter_matrix_to_list(inter_matrix):
    """{mRNA: [miRNAs]} from a two-column table (InterMatrix.toList)."""
    out = {}
    for mir, mr in inter_matrix.dropna().itertuples(index=False):
        out.setdefault(mr, []).append(mir)
    return out


def rank_by_context(context, relationships):
    """Ranks interactions by the mean TargetScan context++ score of each
    pair, from the most negative; pairs without score (NaN) are last
    (PTC.RankByContext in R). Score is text with 15 significant digits, as
    in the R results."""
    if relationships.shape[1] != 2:
        raise ValueError("relationshipsM has to be a matrix with 2 columns [miR,mR]")
    rel = relationships.copy()
    rel.columns = ["miR", "mR"]
    if len(rel) == 0:
        return pd.DataFrame(columns=["rank", "miR", "mR", "Score"])
    # only the context-score rows of the genes being ranked are needed
    ctx = context[context["Gene.Symbol"].isin(set(rel["mR"]))]
    mirs = mirv21(ctx["miRNA"].astype(str).tolist())
    scores = {}
    for mir, gene, score in zip(mirs, ctx["Gene.Symbol"], ctx["context.score"]):
        scores.setdefault((mir, gene), []).append(score)
    mean_scores = [_r_mean(scores[k]) if k in scores else float("nan")
                   for k in zip(rel["miR"], rel["mR"])]
    text = [_r_format(s) for s in mean_scores]
    rel["Score"] = text
    order = np.argsort(np.array([float(t) for t in text]), kind="stable")  # NaN last
    rel = rel.iloc[order].reset_index(drop=True)
    rel.insert(0, "rank", [str(i) for i in range(1, len(rel) + 1)])
    return rel


def _r_format(x):
    """A number as R's as.character(): 15 significant digits."""
    if math.isnan(x):
        return "NaN"
    return f"{x:.15g}"


def confirmed_from_list(inter_list, ground_truth=None):
    """The inferred interactions that are experimentally confirmed
    (Confirmed.fromList in R). Returns {mRNA: [miRNAs]}."""
    if ground_truth is None:
        ground_truth = ptc_data.ground_truth()
    out = {}
    for mr, mirs in inter_list.items():
        if mirs is None or mirs == [None]:
            continue
        hits = [m for m in mirs if m in ground_truth.get(mr, ())]
        if hits:
            out[mr] = hits
    return out


def confirmed_from_matrix(inter_matrix, ground_truth=None):
    """Rows of a miRNA-mRNA table that are experimentally confirmed
    (Confirmed.fromMatrix in R)."""
    if ground_truth is None:
        ground_truth = ptc_data.ground_truth()
    m = inter_matrix.iloc[:, :2]
    rows = []
    for mir in pd.unique(m.iloc[:, 0]):
        for mr in m.iloc[:, 1][m.iloc[:, 0] == mir]:
            if mir in ground_truth.get(mr, ()):
                rows.append((mir, mr))
    return pd.DataFrame(rows, columns=["miR", "mR"])
