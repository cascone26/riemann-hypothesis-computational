#!/usr/bin/env python3
"""
may_criterion.py — PRISM III item #1: May's ecosystem-stability structural diagnostic.

May (1972): a random n×n interaction matrix with connectivity C, variance σ²
is stable iff σ√(nC) < 1. Real food webs escape this by NOT being random:
  - weak-interaction dominant (heavy-tailed |a_ij| distribution)
  - modular (block structure: strong within clusters, weak between)
  - few strong channels buffering many weak links

Question: is τ_N's off-diagonal coupling "food-web-structured" in the same way?
If positivity survives because τ_N is food-web-structured (not a generic argument),
that's a new structural constraint invisible to spectral tools.

Diagnostic battery:
  A. May's parameter σ√(DC) for τ_N's off-diagonal block
  B. Interaction-strength distribution: heavy-tailed vs uniform?
  C. Modularity: block structure of |a_ij| matrix (hierarchical clustering proxy)
  D. Strong-interaction fraction: what fraction carry >50% of total weight?
  E. Comparison to a random matrix of the same size/variance
"""

import math
import numpy as np
import scipy.linalg
from scipy.special import digamma as _digamma
import json
import time

TWO_PI = 2*math.pi
LOGPI  = math.log(math.pi)

def omega_vec(r):
    return np.real(_digamma(0.25+0.5j*r)) - LOGPI

def build_sieve(lsq):
    sieve = np.ones(lsq+1, dtype=bool)
    sieve[0] = sieve[1] = False
    for i in range(2, int(lsq**0.5)+1):
        if sieve[i]: sieve[i*i::i] = False
    lkpp, wts = [], []
    for p in range(2, lsq+1):
        if not sieve[p]: continue
        pk = p
        while pk <= lsq:
            lkpp.append(math.log(pk))
            wts.append(math.log(p)/math.sqrt(pk))
            pk *= p
    return np.array(lkpp), np.array(wts)

def build_tau(N, lsq=14, G=8000, K_BATCH=2000):
    L = 2*math.log(math.sqrt(lsq))
    ns = np.arange(-N, N+1)
    D = len(ns)
    R = TWO_PI*N/L + 60.0
    grid = np.linspace(-R, R, G)
    h = grid[1]-grid[0]
    w = np.full(G, h); w[0]=w[-1]=0.5*h
    Om = omega_vec(grid)
    c = TWO_PI*ns[:,None]/L
    d = c - grid[None,:]
    sm = np.abs(d)<1e-12
    E = np.where(sm, np.sqrt(L)+0j,
                 (1/np.sqrt(L))*np.exp(0.5j*d*L)*2*np.sin(0.5*d*L)/d)
    arch = (E*(w*Om) @ E.conj().T).real / TWO_PI
    sh = np.sinh(L/4)**2
    ni = ns[:,None]; nj = ns[None,:]
    num = L**2 - 16*np.pi**2*ni*nj
    den = (L**2+16*np.pi**2*ni**2)*(L**2+16*np.pi**2*nj**2)
    B = 32*L*sh*num/den
    lkpp, wts_p = build_sieve(lsq)
    K = len(lkpp)
    S = np.zeros(D); S_diag = np.zeros(D)
    for k0 in range(0, K, K_BATCH):
        k1 = min(k0+K_BATCH, K)
        ph = TWO_PI*np.outer(ns, lkpp[k0:k1])/L
        S += np.sin(ph) @ wts_p[k0:k1]
        S_diag += np.cos(ph) @ (2*wts_p[k0:k1]*(1-lkpp[k0:k1]/L))
    Si = S[:,None]; Sj = S[None,:]
    dni = ns[:,None] - ns[None,:]
    with np.errstate(divide='ignore', invalid='ignore'):
        P = np.where(dni != 0, (Sj-Si)/(math.pi*dni), 0.0)
    np.fill_diagonal(P, S_diag)
    tau = arch+B-P
    return (tau+tau.T)/2

def may_analysis(M, label=""):
    """Run May's criterion + food-web diagnostics on matrix M."""
    D = M.shape[0]
    # Extract off-diagonal block
    off_mask = ~np.eye(D, dtype=bool)
    off_vals = M[off_mask]  # D²-D values

    # May's parameter components
    sigma = off_vals.std()       # interaction strength std
    abs_off = np.abs(off_vals)
    # Connectivity C = fraction of off-diag entries "non-negligible"
    noise_floor = abs_off.max() * 0.01
    C = (abs_off > noise_floor).mean()
    may_param = sigma * math.sqrt(D * C)

    # Heavy-tailed test: Gini coefficient of |a_ij| distribution
    abs_sorted = np.sort(abs_off)
    n = len(abs_sorted)
    gini = (2*np.sum((np.arange(1,n+1))*abs_sorted) / (n*abs_sorted.sum()) - (n+1)/n)

    # Strong-interaction fraction (top 10% carry what % of total weight?)
    thresh_90 = np.percentile(abs_off, 90)
    frac_weight_top10 = abs_off[abs_off > thresh_90].sum() / abs_off.sum()

    # Modularity proxy: variance across rows vs within rows
    # If modular: some rows have large variance (they're in a "hub"), others small
    row_stds = np.array([np.abs(M[i, off_mask[i]]).std() for i in range(D)])
    row_std_cv = row_stds.std() / (row_stds.mean() + 1e-12)  # CV of row stds

    # Weak-interaction dominant: fraction with |a_ij| < median
    median_off = np.median(abs_off)
    frac_weak = (abs_off < median_off).mean()  # always ~0.5 by definition, use quartile
    frac_weak_q1 = (abs_off < np.percentile(abs_off, 25)).mean()

    # Off-diagonal decay with distance: |a_{i,j}| vs |i-j|
    ni = np.arange(D)[:,None]; nj = np.arange(D)[None,:]
    dist = np.abs(ni - nj)[off_mask]
    # Bin by distance
    max_dist = D - 1
    dist_bins = {}
    for d in range(1, min(D, 15)):
        mask_d = (dist == d)
        if mask_d.sum() > 0:
            dist_bins[d] = float(abs_off[mask_d].mean())

    # Power-law fit to decay: |a_{ij}| ~ dist^{-alpha}
    dists_arr = np.array(sorted(dist_bins.keys()))
    means_arr = np.array([dist_bins[d] for d in dists_arr])
    if len(dists_arr) >= 3:
        logd = np.log(dists_arr)
        logm = np.log(np.maximum(means_arr, 1e-20))
        alpha, logA = np.polyfit(logd, logm, 1)
    else:
        alpha = 0.0

    # May stability check: sigma*sqrt(DC) < 1?
    may_stable = may_param < 1.0

    print(f"\n  {label} (D={D}):")
    print(f"    σ_off = {sigma:.4e},  C = {C:.4f},  σ√(DC) = {may_param:.4f}  "
          f"({'STABLE' if may_stable else 'UNSTABLE'} under May's criterion)")
    print(f"    Gini coefficient of |a_ij|: {gini:.4f}  "
          f"({'heavy-tailed' if gini > 0.5 else 'uniform-like'})")
    print(f"    Top-10% entries carry {frac_weight_top10*100:.1f}% of total off-diag weight  "
          f"({'concentrated (food-web-like)' if frac_weight_top10 > 0.5 else 'diffuse'})")
    print(f"    Row-std CV = {row_std_cv:.4f}  "
          f"({'modular (high CV)' if row_std_cv > 0.3 else 'homogeneous (low CV)'})")
    print(f"    Off-diag decay exponent α = {alpha:.4f}  "
          f"({'power-law decay (structured)' if abs(alpha) > 0.5 else 'slow decay (long-range)'})")

    # Distance-decay table
    print(f"    Off-diag mean by distance:")
    for d in list(dists_arr[:8]):
        print(f"      |i-j|={int(d):2d}: mean|a_ij| = {dist_bins[int(d)]:.4e}")

    return {
        "D": D, "sigma_off": float(sigma), "C": float(C),
        "may_param": float(may_param), "may_stable": bool(may_stable),
        "gini": float(gini), "frac_weight_top10": float(frac_weight_top10),
        "row_std_cv": float(row_std_cv), "decay_exponent": float(alpha),
        "dist_bins": {str(k): v for k, v in dist_bins.items()}
    }

def foodweb_vs_random(D, sigma, rng_seed=42):
    """Generate a random Gaussian matrix of same size/variance for comparison."""
    rng = np.random.default_rng(rng_seed)
    R = rng.standard_normal((D, D)) * sigma
    R = (R + R.T) / 2
    np.fill_diagonal(R, 0.0)
    return R

if __name__ == "__main__":
    print("May's Criterion — Food-Web Structural Diagnostic on τ_N")
    print("==========================================================")
    print()
    print("May (1972): σ√(nC) < 1 for stability. Real food webs escape by being")
    print("non-random: weak-interaction dominant, modular, few strong channels.")
    print("Question: is τ_N's off-diagonal structure food-web-structured?")

    results = {}

    # Test at multiple N and lsq values
    for (N, lsq) in [(12, 14), (16, 14), (20, 14), (16, 50), (20, 50)]:
        t0 = time.time()
        tau = build_tau(N, lsq)
        r = may_analysis(tau, f"τ_N (lsq={lsq}, N={N})")
        results[f"lsq{lsq}_N{N}"] = r

        # Compare to random matrix of same variance
        sigma_tau = r["sigma_off"]
        R = foodweb_vs_random(2*N+1, sigma_tau)
        r_rand = may_analysis(R, f"RANDOM baseline (same σ={sigma_tau:.3e})")
        results[f"random_D{2*N+1}"] = r_rand
        print(f"  ({time.time()-t0:.1f}s)")

    print()
    print("=== SUMMARY: τ_N vs Random Matrix ===")
    print()

    # Load gram data for larger N comparison
    print("Gram matrix (Nyman-Beurling operator, bvec data) — larger N:")
    for n_pts in [300, 1000]:
        try:
            gram = np.load(f"results/gram_{n_pts}.npy")
            # Extract off-diagonal only
            r_gram = may_analysis(gram, f"Gram (N={n_pts})")
            results[f"gram_{n_pts}"] = r_gram
        except FileNotFoundError:
            print(f"  gram_{n_pts}.npy not found, skip")

    print()
    print("=== KEY FINDINGS ===")
    tau_cases = {k: v for k, v in results.items() if k.startswith("lsq")}
    rand_cases = {k: v for k, v in results.items() if k.startswith("random")}

    print(f"  {'Case':>20}  {'may_param':>10}  {'Gini':>6}  {'top10%':>8}  {'CV_row':>7}  {'α_decay':>8}")
    print(f"  {'─'*70}")
    for k, r in {**tau_cases, **rand_cases}.items():
        print(f"  {k:>20}  {r['may_param']:>10.4f}  {r['gini']:>6.4f}  "
              f"{r['frac_weight_top10']*100:>7.1f}%  {r['row_std_cv']:>7.4f}  {r['decay_exponent']:>8.4f}")

    print()
    tau_vals = list(tau_cases.values())
    rand_vals = list(rand_cases.values())
    if tau_vals and rand_vals:
        avg_gini_tau = np.mean([r["gini"] for r in tau_vals])
        avg_gini_rand = np.mean([r["gini"] for r in rand_vals])
        avg_may_tau = np.mean([r["may_param"] for r in tau_vals])
        avg_may_rand = np.mean([r["may_param"] for r in rand_vals])
        avg_top10_tau = np.mean([r["frac_weight_top10"] for r in tau_vals])
        avg_top10_rand = np.mean([r["frac_weight_top10"] for r in rand_vals])

        print(f"  τ_N mean Gini = {avg_gini_tau:.4f}  vs  Random mean Gini = {avg_gini_rand:.4f}")
        print(f"  τ_N mean May param = {avg_may_tau:.4f}  vs  Random = {avg_may_rand:.4f}")
        print(f"  τ_N top-10% weight = {avg_top10_tau*100:.1f}%  vs  Random = {avg_top10_rand*100:.1f}%")
        print()

        if avg_gini_tau > avg_gini_rand * 1.2:
            print("  FINDING: τ_N is MORE heavy-tailed than a random matrix (food-web-like)")
        else:
            print("  FINDING: τ_N is NOT more heavy-tailed than random — not food-web-structured")

        if all(r["may_stable"] for r in tau_vals):
            print("  All τ_N cases: STABLE under May's criterion (σ√(DC) < 1)")
        else:
            unstable = [k for k, r in tau_cases.items() if not r["may_stable"]]
            print(f"  UNSTABLE cases: {unstable}")

        print()
        print("  Physical interpretation:")
        print("  May's stability ↔ prime-comb coupling is 'sparse + weak-dominated'")
        print("  If Gini(τ_N) >> Gini(random): positivity is STRUCTURAL, not generic")
        print("  If Gini(τ_N) ≈ Gini(random): May's lens gives no additional insight")

    out = "results/may_criterion.json"
    with open(out, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"\nSaved → {out}")
