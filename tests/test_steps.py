"""The building blocks give the same values as in R."""

import numpy as np
import pandas as pd

import ptc
from ptc.core import _r_seq_int
from ptc.decoupled_test import decoupled_test
from ptc.rrandom import RRandom
from conftest import EXPECTED, read_expected


def test_random_numbers_are_r_ones():
    expected = [float.fromhex(x) for x in (EXPECTED / "rnorm.txt").read_text().split()]
    got = np.concatenate([RRandom(1).rnorm(5000), RRandom(12345).rnorm(1000)])
    assert got.tolist() == expected


def test_mirbase_v21_names():
    for batch in ["tcga", "sc", "odd"]:
        expected = read_expected(f"mirv21_{batch}.csv")
        assert ptc.mirv21(expected["name"]) == expected["v21"].tolist()


def test_decoupled_test_p_values(tcga):
    miRNAs, mRNAs = tcga
    sel = ptc.gene_sel(ptc.ptime([miRNAs, mRNAs], mRNAs["VIM"]), 30, 1500)
    d, col = sel["d"], sel["columns"]
    expected = pd.read_csv(EXPECTED / "tcga_pvalues.csv")
    got = []
    for r in expected.itertuples():
        X = d[:, [col[p] for p in sel["PParents"][r.gene]]]
        Y = d[:, [col[r.gene]]]
        S = [int(s) for s in r.S.split()]
        grid = _r_seq_int(0, X.shape[0], r.ngrid + 1)
        got.append(decoupled_test(X, Y, S, grid, r.complements, B=100, rng=RRandom(r.seed)))
    # p-values are multiples of 1/101; the CSV keeps 15 significant digits
    np.testing.assert_allclose(got, expected["p"], rtol=1e-12, atol=0)


def test_ranking_by_context_score():
    rel = pd.DataFrame({"miR": ["a", "b", "c", "d"], "mR": ["A", "B", "C", "D"]})
    ctx = pd.DataFrame({"Gene.Symbol": ["A", "B", "C"], "miRNA": ["a", "b", "c"],
                        "context.score": [-0.05, -0.9, 0.2]})
    ranked = ptc.rank_by_context(ctx, rel)
    assert ranked["miR"].tolist() == ["b", "a", "c", "d"]
    assert ranked["Score"].tolist() == ["-0.9", "-0.05", "0.2", "NaN"]


def test_no_interactions_gives_empty_summary():
    empty = ptc.inter_list_to_matrix({})
    ranked = ptc.rank_by_context(ptc.context_scores(), empty)
    assert len(ranked) == 0
    assert ranked.columns.tolist() == ["rank", "miR", "mR", "Score"]
