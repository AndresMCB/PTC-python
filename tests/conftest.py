from pathlib import Path

import pandas as pd
import pytest

import ptc

EXPECTED = Path(__file__).parent / "expected"


def read_expected(name):
    """An expected result exported from R by scripts/PTC.ExportData.R."""
    return pd.read_csv(EXPECTED / name, dtype=str, keep_default_na=False)


def parents(result):
    """(mRNA, miRNA) pairs of result.names; R writes a missing parent as NA."""
    return [(g, m if m is not None else "NA") for g, ms in result.names.items() for m in ms]


@pytest.fixture(scope="session")
def tcga():
    return ptc.load_tcga_brca()


@pytest.fixture(scope="session")
def single_cell():
    return ptc.load_single_cell()
