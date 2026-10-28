"""PTC() on the datasets of the paper and with other parameter settings,
compared with the outputs of the R package (tests/expected)."""

import pytest

import ptc
from conftest import parents, read_expected


@pytest.mark.filterwarnings("ignore:The invariance test failed")
def test_single_cell(single_cell):
    miRNAs, mRNAs = single_cell
    res = ptc.ptc(miRNAs, mRNAs, VIM=mRNAs["VIM"])
    assert res.summary.equals(read_expected("sc_summary.csv"))
    assert parents(res) == list(read_expected("sc_parents.csv").itertuples(index=False, name=None))
    # as in the paper, the test cannot be computed for these genes
    assert sorted(res.failed) == ["BZW1", "EIF4G2", "SUB1"]


def test_tcga_brca(tcga):
    miRNAs, mRNAs = tcga
    res = ptc.ptc(miRNAs, mRNAs, VIM=mRNAs["VIM"])
    assert res.summary.equals(read_expected("tcga_summary.csv"))
    assert parents(res) == list(read_expected("tcga_parents.csv").itertuples(index=False, name=None))
    assert res.failed == {}


VARIANTS = {
    "first_set": dict(nmiR=15, nmR=150, explore_all=False),
    "pairwise": dict(nmiR=15, nmR=150, complements=False, ngrid=3, alpha=0.05),
    "own_tscan": dict(nmiR=15, nmR=300, TScan=ptc.tscan().iloc[::2]),
    "ngrid4": dict(nmiR=20, nmR=200, ngrid=4, alpha=0.01),
}


@pytest.mark.parametrize("name", VARIANTS)
def test_parameters(tcga, name):
    miRNAs, mRNAs = tcga
    res = ptc.ptc(miRNAs, mRNAs, VIM=mRNAs["VIM"], **VARIANTS[name])
    assert res.summary.equals(read_expected(f"tcga_{name}_summary.csv"))
    expected = read_expected(f"tcga_{name}_parents.csv")
    assert parents(res) == list(expected.itertuples(index=False, name=None))
