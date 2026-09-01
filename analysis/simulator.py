"""
simulator.py
============
Forward population + GxE simulator (pure numpy), calibrated to the CIMMYT
wheat panel used in Phase 1.

Population structure (standard breeding-genetics forward simulation):
  1. nf diploid founders with biallelic markers in linkage equilibrium;
     marker allele frequencies resampled from the real wheat panel (0/1 DArT-like).
  2. ngen generations of random mating with recombination (Poisson crossovers
     per chromosome) -> drift builds linkage disequilibrium and relatedness.
  3. n doubled-haploid (DH) inbred lines are sampled from the final generation,
     so lines are RELATED and markers are coded 0/1 (homozygous), like wheat.

Genetic value (infinitesimal model, all markers causal):
      c_i  = sum_j (X_ij - p_j) a_j ,  a_j ~ N(0, 1), standardized to Var(c)=1

GxE (compound-symmetry factor model, genetic correlation rG):
      g_ij = sqrt(rG) * c_i + sqrt(1 - rG) * s_ij
      y_ij = mu_j + g_ij + e_ij,   Var(e) = (1 - h2)/h2

Oracle variance matrices are returned for the BLUP benchmark:
      Sig_g = (1-rG) I + rG J   (e x e),   Sig_e = sigma2_e I.
"""
import numpy as np


# ------------------------------------------------------------------ genome
def assign_chromosomes(m, C=21):
    """Split marker indices into C ordered chromosomes."""
    base, rem = divmod(m, C)
    chroms, idx = [], 0
    for c in range(C):
        size = base + (1 if c < rem else 0)
        chroms.append(np.arange(idx, idx + size))
        idx += size
    return chroms


def gamete(H_parent, chroms, rng, lam=1.0):
    """One meiosis: a gamete (m-vector of 0/1) from a diploid parent H_parent (2 x m),
    with Poisson(lam) crossovers per chromosome."""
    m = H_parent.shape[1]
    gam = np.empty(m, dtype=float)
    for ci in chroms:
        n = len(ci)
        strand = rng.integers(0, 2)
        if n == 1:
            gam[ci[0]] = H_parent[strand, ci[0]]
            continue
        k = rng.poisson(lam)
        breaks = (np.sort(rng.choice(n - 1, size=min(k, n - 1), replace=False) + 1)
                  if k > 0 else np.array([], dtype=int))
        start = 0
        for b in np.append(breaks, n):
            gam[ci[start:b]] = H_parent[strand, ci[start:b]]
            strand = 1 - strand
            start = b
    return gam


def forward_population(nf, m, p, ngen, chroms, rng, lam=1.0):
    """Founder haplotypes -> ngen random mating -> final diploid population (nf x 2 x m)."""
    H = rng.random((nf, 2, m)) < p[None, None, :]
    H = H.astype(float)
    for _ in range(ngen):
        mom = rng.integers(0, nf, size=nf)
        dad = rng.integers(0, nf, size=nf)
        newH = np.empty_like(H)
        for i in range(nf):
            newH[i, 0] = gamete(H[mom[i]], chroms, rng, lam)
            newH[i, 1] = gamete(H[dad[i]], chroms, rng, lam)
        H = newH
    return H


# ------------------------------------------------------------------ main
def simulate(n=599, m=1279, e=4, h2=0.4, rg=0.5, sigma2_env=1.0,
             nf=100, ngen=15, lam=1.0, C=21, maf_empirical=None, seed=0):
    """Simulate one GxE dataset with related inbred lines.

    Returns dict with X (n x m, 0/1), c (n main genetic values),
    G_true (n x e breeding values), Y (n x e phenotypes), mu, sigma2_e,
    Sig_g (e x e) and Sig_e (e x e) oracle covariances, plus marker freqs p.
    """
    rng = np.random.default_rng(seed)

    if maf_empirical is not None:
        maf = np.asarray(maf_empirical, dtype=float)
        p = rng.choice(maf, size=m, replace=True)
    else:
        p = rng.beta(0.3, 0.7, size=m)
    p = np.clip(p, 0.05, 0.95)

    chroms = assign_chromosomes(m, C)
    H = forward_population(nf, m, p, ngen, chroms, rng, lam)

    # sample n DH lines: one random gamete from a random diploid
    X = np.empty((n, m))
    for i in range(n):
        X[i] = gamete(H[rng.integers(0, nf)], chroms, rng, lam)

    # infinitesimal genetic value (all markers causal), centered by the SAMPLE
    # allele frequencies p_hat of the simulated lines so Cov(c) == VanRaden G
    # (grm_vanraden also centers by the sample frequencies).
    p_hat = X.mean(axis=0)
    Z = X - p_hat[None, :]
    denom = float(np.sum(p_hat * (1.0 - p_hat)))
    sd = np.sqrt(1.0 / denom) if denom > 0 else 1.0
    a = rng.normal(0.0, sd, size=m)
    c = Z @ a                                     # Cov(c) = VanRaden G exactly

    # GxE factor model with heritable interaction (standard factor/CS structure):
    #   g_ij = sqrt(rg) c_i + sqrt(1-rg) t_ij ,  t_j = Z b_j (b_j ~ N(0, I/denom))
    # so Cov(g_ij, g_i'k) = Sig_g[j,k] * G[i,i']  with  Sig_g = (1-rg) I + rg J.
    B = rng.normal(0.0, sd, size=(m, e))          # env-specific marker effects
    T = Z @ B                                     # n x e interaction (Cov = G per env)
    G_true = np.sqrt(rg) * c[:, None] + np.sqrt(1.0 - rg) * T

    mu = rng.normal(0.0, np.sqrt(sigma2_env), size=e)
    sigma2_e = (1.0 - h2) / h2
    E = rng.normal(0.0, np.sqrt(sigma2_e), size=(n, e))
    Y = mu[None, :] + G_true + E

    Sig_g = (1.0 - rg) * np.eye(e) + rg * np.ones((e, e))
    Sig_e = sigma2_e * np.eye(e)

    return {
        "X": X, "c": c, "a": a,
        "G_true": G_true, "Y": Y, "mu": mu,
        "sigma2_e": sigma2_e, "Sig_g": Sig_g, "Sig_e": Sig_e,
        "p": p, "n": n, "m": m, "e": e, "h2": h2, "rg": rg,
    }


# ------------------------------------------------------------------ structured
def simulate_structured(n=599, m=1279, e=4, h2=0.4, kappa=0.5,
                        z=None, nf=100, ngen=15, lam=1.0, C=21,
                        maf_empirical=None, seed=0):
    """Simulate STRUCTURED (reaction-norm / covariate-driven) GxE data.

    The true generative model is a linear reaction norm
        g_ij = c_i + b_i * z_j
    with marker-based intercept (c) and slope (b) breeding values, both having
    covariance proportional to G and being mutually independent. z is an
    observable environmental covariate (it also drives the environment mean, so
    it can be read off from the environment-mean phenotype). The genetic
    covariance is therefore
        Sig_g = sigma2_c * J + sigma2_b * z z'        (rank 2)
    with kappa = sigma2_b / sigma2_c controlling GxE strength, scaled so the
    average genetic variance across environments equals 1 (nominal h2 applies
    on average). Returns the same dict layout as `simulate`, plus z, sigma2_c,
    sigma2_b.
    """
    rng = np.random.default_rng(seed)

    if maf_empirical is not None:
        maf = np.asarray(maf_empirical, dtype=float)
        p = rng.choice(maf, size=m, replace=True)
    else:
        p = rng.beta(0.3, 0.7, size=m)
    p = np.clip(p, 0.05, 0.95)

    if z is None:
        z = np.array([-1.5, -0.5, 0.5, 1.5])[:e]
    z = np.asarray(z, dtype=float)

    chroms = assign_chromosomes(m, C)
    H = forward_population(nf, m, p, ngen, chroms, rng, lam)
    X = np.empty((n, m))
    for i in range(n):
        X[i] = gamete(H[rng.integers(0, nf)], chroms, rng, lam)

    p_hat = X.mean(axis=0)
    Z = X - p_hat[None, :]
    denom = float(np.sum(p_hat * (1.0 - p_hat)))
    s = np.sqrt(1.0 / denom) if denom > 0 else 1.0

    # reaction-norm variance components (average genetic variance = 1)
    z2bar = float(np.mean(z ** 2))
    sigma2_c = 1.0 / (1.0 + kappa * z2bar)
    sigma2_b = kappa * sigma2_c

    a_c = rng.normal(0.0, s * np.sqrt(sigma2_c), size=m)
    a_b = rng.normal(0.0, s * np.sqrt(sigma2_b), size=m)
    c = Z @ a_c                             # Cov(c) = sigma2_c * G
    b = Z @ a_b                             # Cov(b) = sigma2_b * G
    G_true = c[:, None] + b[:, None] * z[None, :]

    # environment mean follows the covariate z (so z is observable via env means)
    mu = z.copy()
    sigma2_e = (1.0 - h2) / h2
    E = rng.normal(0.0, np.sqrt(sigma2_e), size=(n, e))
    Y = mu[None, :] + G_true + E

    Sig_g = sigma2_c * np.ones((e, e)) + sigma2_b * np.outer(z, z)
    Sig_e = sigma2_e * np.eye(e)

    return {
        "X": X, "c": c, "b": b, "z": z,
        "G_true": G_true, "Y": Y, "mu": mu,
        "sigma2_e": sigma2_e, "Sig_g": Sig_g, "Sig_e": Sig_e,
        "p": p, "n": n, "m": m, "e": e, "h2": h2, "kappa": kappa,
        "sigma2_c": sigma2_c, "sigma2_b": sigma2_b,
    }
