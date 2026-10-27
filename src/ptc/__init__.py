"""PTC: Pseudo-temporal Causality approach to identify miRNA-mRNA interactions.

Python version of the R package PTC (https://github.com/AndresMCB/PTC).
"""

from .core import (PTCResult, confirmed_from_list, confirmed_from_matrix,
                   extract_parents, find_pp, gene_sel, get_data_by_mad,
                   inter_list_to_matrix, inter_matrix_to_list, ptc, ptime,
                   rank_by_context, test_invariance)
from .data import (context_scores, ground_truth, load_single_cell,
                   load_tcga_brca, tscan)
from .mirbase import mirv21

__version__ = "1.1.2"
