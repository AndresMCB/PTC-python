"""The example script reproduces the experiments of the paper (R results of
inst/PTC.MainApplication.R)."""

import runpy
from pathlib import Path

import numpy as np
import pytest

import ptc
from conftest import EXPECTED, parents, read_expected

SCRIPT = Path(__file__).parent.parent / "examples" / "paper_experiments.py"


def test_wanderlust_time_is_read_as_in_r():
    expected = [float.fromhex(x) for x in (EXPECTED / "wandTime.txt").read_text().split()]
    assert ptc.load_wanderlust_time().tolist() == expected


@pytest.mark.filterwarnings("ignore:The invariance test failed")
def test_example_script_gives_r_results():
    script = runpy.run_path(str(SCRIPT))
    results = script["run_experiments"]()
    for key, name in [("test1", "tcga"), ("test2", "sc"), ("test3", "sc_wanderlust")]:
        assert results[key].summary.equals(read_expected(f"{name}_summary.csv"))
        expected = read_expected(f"{name}_parents.csv")
        assert parents(results[key]) == list(expected.itertuples(index=False, name=None))
    # numbers reported by inst/PTC.MainApplication.R in R
    counts = [script["confirmed"](results[k]) for k in ["test1", "test2", "test3"]]
    assert counts == [(30, 25), (294, 121), (293, 122)]
