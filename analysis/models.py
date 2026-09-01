"""
models.py
=========
Genomic-prediction models for multi-environment GxE benchmarking (pure numpy).

Models (all BLUP variants, with known/oracle variance components so the methods
can be compared on equal footing):
  - gblup_single           : single-environment GBLUP (RR-BLUP equivalent)
  - gblup_pooled           : across-environment GBLUP (line main effect)
  - mt_gblup               : multi-trait GBLUP (environments = correlated traits)
  - mt_predict_new_env     : CV1 prediction of an untested environment
  - mt_predict_new_lines   : CV2 prediction of untested lines
  - rn_covariance          : reduced-rank reaction-norm covariance (rank-2)
"""
import numpy as np
from scipy.linalg import solve

from utils import grm_vanraden


# ----------------------------------------------------------------- single env
def gblup_single(G_train, y_train, G_cross, lam):
    """Single-env GBLUP. y_train is a centered phenotype vector (n_train,).
    Returns predictions for the rows of G_cross (n_test x n_train)."""
    G_train = np.asarray(G_train, dtype=float)
    n = G_train.shape[0]
    K = G_train + lam * np.eye(n)
    alpha = solve(K, y_train, assume_a="pos")
    return G_cross @ alpha


# ----------------------------------------------------------------- pooled
def gblup_pooled(Y_tr, X, lam_pool):
    """Across-environment GBLUP: a single line main effect shared across
    environments. Y_tr is n x r (env-centered phenotypes in the r training
    environments). Returns the n-vector of predicted line main effects."""
    Y_tr = np.asarray(Y_tr, dtype=float)
    G = grm_vanraden(X, coding="01")
    n = G.shape[0]
    r = Y_tr.shape[1]
    # BLUP of the line main effect from the r within-line means
    ybar = Y_tr.mean(axis=1)
    K = G + (lam_pool / r) * np.eye(n)   # residual var of the mean = var_res / r
    return solve(K, ybar, assume_a="pos")


# ----------------------------------------------------------------- multi-trait
def mt_gblup(Y, G, Sig_g, Sig_e, eig=None):
    """Exact multi-trait GBLUP for complete n x e data with covariance
    Cov(vec U) = Sig_g (x) G and Cov(vec E) = Sig_e (x) I.

    Uses the eigen-decomposition of G to reduce the Kronecker system to n
    independent e x e solves. Returns U_hat (n x e) = BLUP of breeding values.
    Pass `eig=(w, V)` (precomputed np.linalg.eigh(G)) to avoid recomputation.
    """
    Y = np.asarray(Y, dtype=float)
    n, e = Y.shape
    # center by column means (fixed env effects)
    Yc = Y - Y.mean(axis=0, keepdims=True)

    if eig is None:
        w, V = np.linalg.eigh(G)          # G = V diag(w) V'
    else:
        w, V = eig
    YV = V.T @ Yc                         # rotate to eigen-space (n x e)

    U_hat_V = np.empty_like(YV)
    for l in range(n):
        A_l = w[l] * Sig_g + Sig_e        # e x e
        try:
            rhs = solve(A_l, YV[l], assume_a="pos")
        except np.linalg.LinAlgError:
            # estimated covariance may be singular (moment estimator / rank-reduced
            # reaction-norm); generalized inverse gives the legitimate BLUP solution
            rhs = np.linalg.pinv(A_l) @ YV[l]
        U_hat_V[l] = w[l] * (Sig_g @ rhs)

    return V @ U_hat_V


def mt_predict_new_env(U_hat_train, Sig_g, target_env, train_envs):
    """CV1: predict breeding values in an untested environment.

    U_hat_train : n x r BLUP in the training environments
    Sig_g       : full e x e genetic covariance
    Returns n-vector of predicted breeding values in `target_env`.
    """
    train_envs = list(train_envs)
    S_tt = Sig_g[np.ix_(train_envs, train_envs)]        # r x r
    S_tj = Sig_g[np.ix_(train_envs, [target_env])]      # r x 1
    try:
        w = solve(S_tt, S_tj, assume_a="pos")           # r x 1 (full rank)
    except np.linalg.LinAlgError:
        w = np.linalg.pinv(S_tt) @ S_tj                 # reduced-rank (reaction norm)
    return U_hat_train @ w.ravel()


def mt_predict_new_lines(G_cross, G_train, U_hat_train):
    """CV2: predict breeding values of untested lines from their genomic
    relationship to the training lines."""
    return G_cross @ solve(G_train, U_hat_train, assume_a="pos")


# ----------------------------------------------------------------- factor analytic
def fa_covariance(Sig_g, k):
    """Factor-analytic rank-k covariance (Burgueno et al. 2012 style).

    Returns  Sig_fa = Lambda Lambda' + Psi, where Lambda is the e x k loading
    matrix from the top-k eigenvalues of Sig_g and Psi is the diagonal matrix of
    residual (environment-specific) variances, so the diagonal of Sig_fa equals
    the diagonal of Sig_g. For k = e this reproduces Sig_g exactly, so the FA
    models are nested within the full multi-trait model.
    """
    Sig_g = np.asarray(Sig_g, dtype=float)
    e = Sig_g.shape[0]
    S = (Sig_g + Sig_g.T) / 2.0
    w, V = np.linalg.eigh(S)
    w = np.clip(w, 0.0, None)
    order = np.argsort(w)[::-1][:k]
    Lam = V[:, order] * np.sqrt(w[order])          # e x k loadings
    resid = np.diag(S) - np.sum(Lam ** 2, axis=1)
    resid = np.clip(resid, 0.0, None)
    return Lam @ Lam.T + np.diag(resid)


# ----------------------------------------------------------------- reaction norm
def rn_covariance(Sig_g, z):
    """Reduced-rank (rank <= 2) reaction-norm covariance on environmental
    covariate z:  Sig_g_rr = P_z Sig_g P_z, where P_z projects onto span([1, z]).

    Implements a linear reaction norm g_ij = g0_i + g1_i * z_j, with the
    covariance structure obtained by projecting the oracle Sig_g.
    """
    z = np.asarray(z, dtype=float)
    e = z.shape[0]
    Z = np.column_stack([np.ones(e), z])
    P = Z @ solve(Z.T @ Z, Z.T, assume_a="pos")
    return P @ Sig_g @ P
