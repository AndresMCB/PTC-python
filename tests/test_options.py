"""Other parameter settings give the same results as the R package."""

import pytest

import ptc
from conftest import parents, read_expected

VARIANTS = {
    "first_set": dict(nmiR=15, nmR=150, explore_all=False),
    "pairwise": dict(nmiR=15, nmR=150, complements=False, ngrid=3, alpha=0.05),
    "own_tscan": dict(nmiR=15, nmR=300, TScan=ptc.tscan().iloc[::2]),
    "ngrid4": dict(nmiR=20, nmR=200, ngrid=4, alpha=0.01),
}


@pytest.mark.parametrize("name", VARIANTS)
def test_option(tcga, name):
    miRNAs, mRNAs = tcga
    res = ptc.ptc(miRNAs, mRNAs, VIM=mRNAs["VIM"], **VARIANTS[name])
    assert res.summary.equals(read_expected(f"tcga_{name}_summary.csv"))
    expected = read_expected(f"tcga_{name}_parents.csv")
    assert parents(res) == list(expected.itertuples(index=False, name=None))
