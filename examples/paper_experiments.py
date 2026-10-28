"""Experiments of the paper, as inst/PTC.MainApplication.R in the R package.

Cifuentes-Bernal AM, Pham VVH, Li X, Liu L, Li J, Le TD. A pseudo-temporal
causality approach to identifying miRNA-mRNA interactions during biological
processes. Bioinformatics, 2021. https://doi.org/10.1093/bioinformatics/btaa899

Run with:  python examples/paper_experiments.py
"""

import ptc


def run_experiments():
    results = {}

    ###################### Test with bulk data ######################
    miRNAs, mRNAs = ptc.load_tcga_brca()
    test1 = ptc.ptc(miRNAs, mRNAs, VIM=mRNAs["VIM"])
    results["test1"] = test1

    ###################### Single-cell data ######################
    # Data of GSE114071 (Wang et al., 2019), prepared as in
    # inst/PTC.MainApplication.R: matched half-cells only, miRNA expression
    # back from log2 and amplified, and mRNAs with mean expression > 0.1
    # that are expressed in 20% of the cells. Cells are matched by row order.
    SC_miRNAs, SC_mRNAs = ptc.load_single_cell()
    test2 = ptc.ptc(SC_miRNAs, SC_mRNAs, VIM=SC_mRNAs["VIM"])
    results["test2"] = test2

    #### Additional test by using pseudotime from Wanderlust ####
    test3 = ptc.ptc(SC_miRNAs, SC_mRNAs, VIM=ptc.load_wanderlust_time())
    results["test3"] = test3

    return results


def confirmed(result):
    """Experimentally confirmed interactions: all of them, and among the
    200 interactions with the best TargetScan context++ score."""
    all_confirmed = ptc.confirmed_from_list(result.names)
    top200 = ptc.confirmed_from_matrix(result.summary.iloc[:200][["miR", "mR"]])
    return sum(len(m) for m in all_confirmed.values()), len(top200)


if __name__ == "__main__":
    import warnings
    # three single-cell genes cannot be tested (singular design matrix),
    # as in the paper; PTC reports them in result.failed
    warnings.filterwarnings("ignore", message="The invariance test failed")

    names = {"test1": "TCGA BRCA, VIM pseudotime",
             "test2": "single-cell, VIM pseudotime",
             "test3": "single-cell, Wanderlust pseudotime"}
    for key, result in run_experiments().items():
        n_conf, n_top = confirmed(result)
        print(f"{key} ({names[key]}): {len(result.summary)} interactions, "
              f"{n_conf} confirmed, {n_top} confirmed in the top 200")
