"""
Plant-breeding quantitative-genetics utilities (pure Python / numpy).

All methods are standard and citable:
  - VanRaden (2008) genomic relationship matrix (method 1)
  - GBLUP (equivalent to RR-BLUP under the corresponding GRM; Habier et al. 2007, 2013)
  - ANOVA-based variance components for two-way (line x environment) tables
    (Comstock & Moll 1963; standard GxE practice)
  - Multivariate Haseman-Elston (1972) moment estimator for GxE variance components
"""
import numpy as np
import pandas as pd
from scipy.linalg import solve


# ---------------------------------------------------------------- GRM
def grm_vanraden(X, coding="01"):
    """VanRaden (2008) method-1 G matrix from marker matrix X (n x m).

    coding = "01"  -> markers are 0/1 (e.g. dominant DArT); allele freq p = colmean
                      marker variance = p(1-p)  ->  denom = sum p(1-p)
    coding = "012" -> markers are 0/1/2 dosage;      allele freq p = colmean/2
                      marker variance = 2p(1-p) ->  denom = 2 sum p(1-p)
    """
    X = np.asarray(X, dtype=float)
    p = X.mean(axis=0)
    if coding == "012":
        p = p / 2.0
    Z = X - p  # centering by allele frequency (method 1)
    denom = np.sum(p * (1.0 - p))
    if coding == "012":
        denom *= 2.0
    G = Z @ Z.T / denom
    return G


def standardize_markers(X, coding="01"):
    """Center markers by sample allele frequency and scale column-wise by sqrt(p(1-p)).

    This is the RR-BLUP design matrix. Note: column-wise scaling corresponds to
    the VanRaden *method-2* GRM (Zs @ Zs.T); the method-1 GRM uses a single
    global scaling (sum of p(1-p)) instead. `gblup_cv`/`grm_vanraden` use method 1.
    """
    X = np.asarray(X, dtype=float)
    p = X.mean(axis=0)
    if coding == "012":
        p = p / 2.0
    Z = X - p
    sd = np.sqrt(p * (1.0 - p))  # note: for 0/1 coding, scale uses sqrt(p(1-p))
    sd[sd < 1e-12] = 1.0
    return Z / sd


# ---------------------------------------------------------------- RR-BLUP / GBLUP
def rrblup_fit(Xtr, ytr, lam):
    """Ridge regression BLUP on standardized marker matrix. Returns intercept + beta."""
    X = np.asarray(Xtr, dtype=float)
    y = np.asarray(ytr, dtype=float).ravel()
    n, m = X.shape
    Xc = np.column_stack([np.ones(n), X])
    XtX = Xc.T @ Xc
    XtX[np.arange(m + 1), np.arange(m + 1)] += lam  # do not shrink intercept
    beta = solve(XtX, Xc.T @ y, assume_a="pos")
    return beta


def rrblup_predict(Xte, beta):
    X = np.asarray(Xte, dtype=float)
    Xc = np.column_stack([np.ones(X.shape[0]), X])
    return Xc @ beta


def rrblup_cv(X, y, k=5, seed=42, lam_grid=None):
    """k-fold cross-validated RR-BLUP; lambda tuned by inner k-fold on the training fold."""
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float).ravel()
    n = X.shape[0]
    rng = np.random.default_rng(seed)
    idx = rng.permutation(n)
    folds = np.array_split(idx, k)
    if lam_grid is None:
        lam_grid = np.logspace(-2, 4, 13)
    pred = np.full(n, np.nan)
    best_lams = []
    for tr in folds:
        va_mask = np.ones(n, bool)
        va_mask[tr] = False
        # inner CV for lambda on the training set only
        inner = np.array_split(rng.permutation(np.flatnonzero(va_mask)), 5)
        scores = []
        for lam in lam_grid:
            accs = []
            for iv in inner:
                itr = np.ones(n, bool)
                itr[iv] = False
                itr[tr] = False
                beta = rrblup_fit(X[itr], y[itr], lam)
                accs.append(np.corrcoef(rrblup_predict(X[iv], beta), y[iv])[0, 1])
            scores.append(np.mean(accs))
        lam_best = lam_grid[int(np.argmax(scores))]
        best_lams.append(lam_best)
        beta = rrblup_fit(X[va_mask], y[va_mask], lam_best)
        pred[tr] = rrblup_predict(X[tr], beta)
    corr = np.corrcoef(pred, y)[0, 1]
    return {"corr": corr, "pred": pred, "obs": y, "lam_grid": lam_grid, "best_lams": best_lams}


# ---------------------------------------------------------------- GBLUP CV
def gblup_cv(G, y, k=5, seed=42, lam_grid=None, folds=None):
    """k-fold cross-validated GBLUP on a precomputed GRM G, with lambda tuned by
    an inner k-fold on each training set (method-1 GBLUP, consistent with the
    simulation / validation scripts). Returns dict with 'corr', 'pred', 'obs'."""
    G = np.asarray(G, dtype=float)
    y = np.asarray(y, dtype=float).ravel()
    n = G.shape[0]
    rng = np.random.default_rng(seed)
    if folds is None:
        idx = rng.permutation(n)
        folds = np.array_split(idx, k)
    if lam_grid is None:
        lam_grid = np.logspace(-2, 4, 13)
    pred = np.full(n, np.nan)
    best_lams = []
    for tr in folds:
        tr = np.asarray(tr)
        va = np.ones(n, bool)
        va[tr] = False
        va_idx = np.flatnonzero(va)
        inner = np.array_split(rng.permutation(va_idx), 5)
        scores = []
        for lam in lam_grid:
            accs = []
            for iv in inner:
                itr = np.ones(n, bool)
                itr[iv] = False
                itr[tr] = False
                itr_idx = np.flatnonzero(itr)
                Gtr = G[np.ix_(itr_idx, itr_idx)]
                Gcr = G[np.ix_(iv, itr_idx)]
                ytr = y[itr_idx]
                K = Gtr + lam * np.eye(itr_idx.size)
                alpha = solve(K, ytr, assume_a="pos")
                accs.append(np.corrcoef(Gcr @ alpha, y[iv])[0, 1])
            scores.append(np.mean(accs))
        lam_best = lam_grid[int(np.argmax(scores))]
        best_lams.append(lam_best)
        Gtr = G[np.ix_(va_idx, va_idx)]
        Gcr = G[np.ix_(tr, va_idx)]
        alpha = solve(Gtr + lam_best * np.eye(va_idx.size), y[va_idx], assume_a="pos")
        pred[tr] = Gcr @ alpha
    corr = np.corrcoef(pred, y)[0, 1]
    return {"corr": corr, "pred": pred, "obs": y, "best_lams": best_lams}


# ---------------------------------------------------------------- variance components
def estimate_components_haseman_elston(Y, G):
    """Multivariate Haseman-Elston (GCTA-style) moment estimator.

    Estimates the genetic covariance matrix Sig_g (e x e) and the residual
    covariance matrix Sig_e (diagonal) of a multi-environment trait vector from
    an n x e phenotype matrix Y and an n x n relationship matrix G, under
      Cov(y_ij, y_i'k) = Sig_g[j, k] * G[i, i'] + Sig_e[j, k] * 1[i=i'].

    Uses off-diagonal entries only (unbiased by residual), i.e. classical
    Haseman-Elston regression of the empirical cross-product on G.
    """
    Y = np.asarray(Y, dtype=float)
    G = np.asarray(G, dtype=float)
    n, e = Y.shape

    # off-diagonal mask
    tri = np.triu_indices(n, k=1)
    a = G[tri]                        # off-diagonal relationship values

    Yc = Y - Y.mean(axis=0, keepdims=True)
    Sig_g = np.empty((e, e))
    Sig_e = np.empty((e, e))
    g_bar_diag = float(np.mean(np.diag(G)))

    for j in range(e):
        for k in range(j, e):
            sj = Yc[:, j]
            sk = Yc[:, k]
            # cross-products of two different individuals, weighted by relatedness
            s_off = sj[tri[0]] * sk[tri[1]]
            if j != k:
                s_off = 0.5 * (s_off + sj[tri[1]] * sk[tri[0]])   # symmetrize
            sg = float(np.sum(s_off * a) / np.sum(a * a))
            Sig_g[j, k] = sg
            Sig_g[k, j] = sg
            if j == k:
                se = float(np.var(Y[:, j]) - sg * g_bar_diag)
                Sig_e[j, k] = max(se, 1e-6)
            else:
                Sig_e[j, k] = 0.0
                Sig_e[k, j] = 0.0

    return Sig_g, Sig_e


def make_psd(S, jitter=1e-4):
    """Return the positive semi-definite 'bending' of a symmetric matrix S by
    clipping negative eigenvalues to `jitter`. Standard practice for moment-
    estimated covariance matrices (which need not be PSD in finite samples)."""
    S = np.asarray(S, dtype=float)
    S = (S + S.T) / 2.0
    w, V = np.linalg.eigh(S)
    w = np.clip(w, jitter, None)
    return V @ np.diag(w) @ V.T


def two_way_anova_components(Y):
    """Variance components from a two-way (line x env) table with r=1 per cell.

    Returns dict with sigma2_g, sigma2_res (=sigma2_ge + sigma2_e, confounded at r=1),
    sigma2_env, and broad-sense heritability on entry means H2.
    """
    Y = np.asarray(Y, dtype=float)
    g, e = Y.shape
    grand = Y.mean()
    SS_env = g * np.sum((Y.mean(axis=0) - grand) ** 2)
    SS_g = e * np.sum((Y.mean(axis=1) - grand) ** 2)
    SS_ge = np.sum((Y - Y.mean(axis=1, keepdims=True) - Y.mean(axis=0, keepdims=True) + grand) ** 2)
    df_env = e - 1
    df_g = g - 1
    df_ge = (g - 1) * (e - 1)
    MS_env = SS_env / df_env
    MS_g = SS_g / df_g
    MS_ge = SS_ge / df_ge
    sigma2_res = MS_ge
    sigma2_g = (MS_g - MS_ge) / e
    sigma2_g = max(sigma2_g, 0.0)
    sigma2_env = (MS_env - MS_ge) / g
    sigma2_env = max(sigma2_env, 0.0)
    H2 = sigma2_g / (sigma2_g + sigma2_res / e)
    return {
        "SS_env": SS_env, "SS_g": SS_g, "SS_ge": SS_ge,
        "MS_env": MS_env, "MS_g": MS_g, "MS_ge": MS_ge,
        "sigma2_g": sigma2_g, "sigma2_res": sigma2_res,
        "sigma2_env": sigma2_env, "H2": H2,
    }
