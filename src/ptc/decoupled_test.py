"""Decoupled invariance test used by test_invariance().

Python version of R/decoupledTest.R in the R package PTC, which is the part
of seqICP::seqICP.s() (seqICP 1.1, by Niklas Pfister and Jonas Peters, GPL-3)
that PTC uses: test = "decoupled", model = "iid", link = sum and a parametric
bootstrap. The random draws are R's (see rrandom.py), so p-values match
those computed in R.

Pfister, N., P. Buhlmann and J. Peters (2019). Invariant Causal Prediction for
Sequential Data. Journal of the American Statistical Association.
https://doi.org/10.1080/01621459.2018.1491403
"""

from itertools import combinations
from math import comb

import numpy as np

_EPS = np.finfo(float).eps


def _solve(a, b):
    """solve(a, b) with R's check: R stops when the reciprocal condition
    number of a is below the machine precision."""
    rcond = 1.0 / np.linalg.cond(a, 1)
    if not rcond >= _EPS:
        raise np.linalg.LinAlgError(
            f"system is computationally singular: reciprocal condition number = {rcond:g}")
    return np.linalg.solve(a, b)


def decoupled_test(X, Y, S, grid, complements, B=100, rng=None):
    """Bonferroni-adjusted p-value of H0: "the predictors S are invariant".

    X: n x p matrix of candidate predictors, Y: n x 1 response,
    S: columns of X to test (1-based, as in R), grid: change points that
    define the environments, complements: compare each environment against
    its complement (True) or against all non-overlapping environments,
    B: bootstrap samples, rng: an RRandom generator.
    """
    n = X.shape[0]
    # iid model: intercept plus the predictors in S
    X = np.column_stack([np.ones(n), X[:, [s - 1 for s in S]]])
    d = X.shape[1]

    ###
    # Environments defined by the grid, and the matrices Z to fit each of them
    ###
    grid = list(grid)
    L = len(grid)
    if L <= 1:
        raise ValueError("not enough grid points given, specify correct grid")
    if grid[0] == 1:
        grid[0] = 0
    L_tot = comb(L, 2)
    cuts = list(combinations(range(1, L + 1), 2))
    env = []  # 0-based rows of each environment
    if complements:
        comp = []
        for c1, c2 in cuts:
            cp1, cp2 = grid[c1 - 1], grid[c2 - 1]
            env.append(np.arange(cp1, cp2))
            if cp1 == 0 and cp2 != n:
                comp.append(np.arange(cp2, n))
            elif cp2 == n and cp1 != 0:
                comp.append(np.arange(0, cp1))
            elif cp1 == 0 and cp2 == n:
                comp.append(np.arange(0, n))
            else:
                comp.append(np.concatenate([np.arange(0, cp1), np.arange(cp2, n)]))
        env += comp
        fullenv = L - 1
        idx = [i for i in range(1, L_tot + 1) if i != fullenv]
        pairwise = [(i, i + L_tot) for i in idx]
        L_tot = 2 * L_tot
    else:
        for c1, c2 in cuts:
            env.append(np.arange(grid[c1 - 1], grid[c2 - 1]))
        first, second = [], []
        for i in range(1, L_tot + 1):
            for k in range(1, L_tot + 1):
                if cuts[k - 1][0] >= cuts[i - 1][1]:
                    first.append(i)
                    second.append(k)
        pairwise = list(zip(first + second, second + first))
    X_env = [X[e] for e in env]
    Z = [_solve(xe.T @ xe, xe.T) for xe in X_env]
    pairs = np.array(pairwise) - 1

    ###
    # Test statistics: differences in coefficients (row 1) and in residual
    # variances (row 2) between pairs of environments. R holds one column
    # per residual vector (observed or bootstrap).
    ###
    # (as in R, a segment with no residual degrees of freedom gives an
    # infinite or undefined variance, without warning)
    @np.errstate(divide="ignore", invalid="ignore")
    def test_stat(R):
        gamma = [Z[i] @ R[env[i]] for i in range(L_tot)]
        var = np.array([((R[env[i]] - X_env[i] @ gamma[i]) ** 2).sum(axis=0) / (len(env[i]) - d)
                        for i in range(L_tot)])
        T1 = sum(((gamma[i] - gamma[j]) ** 2).sum(axis=0) for i, j in pairs)
        T2 = sum(np.abs(var[i] / var[j] - 1) for i, j in pairs)
        return np.vstack([T1, T2])

    # Scaled residuals of the pooled linear fit; M = I - X (X'X)^-1 X'
    XtX_Xt = _solve(X.T @ X, X.T)
    M = np.eye(n) - X @ XtX_Xt
    R = M @ Y
    R = R / np.sqrt((R ** 2).sum())
    observed = test_stat(R)

    ###
    # Parametric bootstrap of the scaled residuals for the null distribution
    ###
    e = rng.rnorm(n * B).reshape(B, n).T  # column b = rnorm(n) of bootstrap b
    tmp = e - X @ (XtX_Xt @ e)
    boot = test_stat(tmp / np.sqrt((tmp ** 2).sum(axis=0)))
    p_value = ((boot >= observed).sum(axis=1) + 1) / (B + 1)

    # Bonferroni adjustment for the two statistics
    return min(2 * p_value.min(), 1.0)
