"""PTC with the settings of the paper gives the same results as the R package."""

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
