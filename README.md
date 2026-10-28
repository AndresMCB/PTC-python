# PTC for Python

PTC-python is a Python implementation of the R package [PTC](https://github.com/AndresMCB/PTC) (version 1.1.2), which finds causal miRNA-mRNA relationships during biological processes. Both implementations produce identical outputs: interactions, ranks and scores (see [Validation](#validation)).

Detailed information about the method can be found in:

A Pseudo-Temporal Causality Approach to Identifying miRNA-mRNA Interactions During Biological Processes
Andres M. Cifuentes-Bernal, Vu VH Pham, Xiaomei Li, Lin Liu, Jiuyong Li, Thuc Duy Le
Bioinformatics, 2021; doi: https://doi.org/10.1093/bioinformatics/btaa899

## Installation
PTC needs Python 3.9 or higher and depends only on numpy (1.26 or higher) and pandas (1.5 or higher). It was tested with Python 3.11 to 3.13.

```bash
pip install git+https://github.com/AndresMCB/ptc-python
```

## Usage
Expression data are pandas DataFrames with samples in rows and genes in columns. The miRNA and mRNA tables must contain the same samples, in the same order.

```python
import ptc

# breast cancer bulk data from TCGA, as in the paper
miRNAs, mRNAs = ptc.load_tcga_brca()

# samples are ordered by the expression of the EMT marker VIM (pseudotime)
result = ptc.ptc(miRNAs, mRNAs, VIM=mRNAs["VIM"])

result.summary      # inferred interactions ranked by TargetScan context++ score
result.names        # parent miRNAs of each mRNA

# experimentally confirmed interactions (miRTarBase, TarBase, miRWalk)
ptc.confirmed_from_list(result.names)
```

To use your own data, read the expression tables with pandas, for example:

```python
import pandas as pd
miRNAs = pd.read_csv("mirna.csv", index_col=0)
mRNAs = pd.read_csv("mrna.csv", index_col=0)
result = ptc.ptc(miRNAs, mRNAs, VIM=mRNAs["VIM"])
```

Parameters and default values follow the R package: `nmiR=30`, `nmR=1500`, `ngrid=2`, `alpha=0.02`, `complements=True`, `explore_all=True`, `silent=True` and `TScan=None` (TargetScan 7.0). `TScan` can be a DataFrame with your own miRNA-mRNA predictions (columns miRNA and mRNA, miRBase v21 names).

On the TCGA BRCA data, `ptc.ptc()` takes a few seconds.

## Experiments of the paper
[examples/PTC.MainApplication.py](examples/PTC.MainApplication.py) runs the experiments of the paper, as [inst/PTC.MainApplication.R](https://github.com/AndresMCB/PTC/blob/master/inst/PTC.MainApplication.R) in the R package:

```bash
python examples/PTC.MainApplication.py
```

| Experiment | Interactions | Experimentally confirmed | Confirmed in the top 200 |
|---|---|---|---|
| test1: TCGA BRCA, VIM pseudotime | 292 | 30 | 25 |
| test2: single-cell, VIM pseudotime | 567 | 294 | 121 |
| test3: single-cell, Wanderlust pseudotime | 565 | 293 | 122 |

The R script reports the same values. In the single-cell data, the invariance test cannot be computed for BZW1, EIF4G2 and SUB1 (singular design matrix); they are listed in `result.failed`.

## Functions
| Python | R |
|---|---|
| `ptc` | `PTC` |
| `ptime` | `PTC.ptime` |
| `gene_sel` | `PTC.GeneSel` |
| `get_data_by_mad` | `getDatabyMAD` |
| `find_pp` | `PTC.findPP` |
| `test_invariance` | `PTC.TestInvariance` |
| `extract_parents` | `Extract.Parents` |
| `inter_list_to_matrix`, `inter_matrix_to_list` | `InterList.toMatrix`, `InterMatrix.toList` |
| `rank_by_context` | `PTC.RankByContext` |
| `mirv21` | `PTC.miRv21` |
| `confirmed_from_list`, `confirmed_from_matrix` | `Confirmed.fromList`, `Confirmed.fromMatrix` |
| `load_tcga_brca`, `load_single_cell`, `load_wanderlust_time` | `TCGA_BRCAdata`, `SC_miRNAsdata` and `SC_mRNAsdata`, `wandTime.csv` |
| `tscan`, `context_scores`, `ground_truth` | `TScan`, `TS7.0_Conserved_Site_Context_Scores`, `GroundT` |

## Validation
PTC-python was validated against the R package. To obtain identical outputs, it:
* uses R's random numbers: `set.seed(1)` and `rnorm()` (Mersenne-Twister and inversion) are reproduced in `ptc/rrandom.py`;
* converts miRNA names to miRBase v21 with the same table and rules as miRBaseConverter;
* computes medians and means as R does (extended precision).

The data in `src/ptc/data` and the expected results in `tests/expected` were exported from the R package with `scripts/PTC.ExportData.R`. To run the tests:

```bash
pip install "ptc-mirna[test] @ git+https://github.com/AndresMCB/ptc-python"
pytest
```

The results were checked on Linux (x86-64). On Windows, numpy's extended precision is not available, so in rare cases the last digit of a median could differ from R.

## License
GPL-3, as the R package. The decoupled invariance test is adapted from the R package seqICP by Niklas Pfister and Jonas Peters (GPL-3).
