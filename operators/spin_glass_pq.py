#!/usr/bin/env python3
"""
spin_glass_pq.py — PRISM EXTREME item #5: Spin-glass overlap distribution P(q).

Treats tau_N's Rayleigh-quotient landscape as a spin-glass Hamiltonian.
Computes the overlap distribution P(q) between independent replicas.

P(q) diagnostic:
  Two delta peaks at {0, 1} → single basin (monofractal), landscape is convex
  Continuous P(q) → replica-symmetry broken (RSB), rugged multi-basin landscape
  One peak at q=1 → degenerate ground state (multi-eigenvector near-zero manifold)

Method:
1. Build tau_N at two lsq values: near-resonant (lsq≈50, positive spectrum) and
   exact resonant from spectral_gap_analysis (N=32, 22 near-zero modes).
2. For each: sample K random unit vectors {x_i}, compute Rayleigh quotients R(x_i).
3. Compute pairwise overlaps q_{ij} = (x_i · x_j)^2.
4. For the near-zero subspace: project onto eigenvectors with |λ|<threshold,
   compute overlaps within that subspace.
5. Bin P(q) and analyze: width, peaks, shape.
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

def build_tau(N, lsq, G=8000, K_BATCH=3000):
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

def compute_pq(tau, K_samples=2000, near_zero_thresh=None, rng_seed=42):
    """
    Compute P(q) overlap distribution.
    Returns: overlaps array, Rayleigh quotients array, near-zero overlap array.
    """
    rng = np.random.default_rng(rng_seed)
    D = tau.shape[0]
    evals, evecs = scipy.linalg.eigh(tau)
    min_ev = evals[0]

    # Shift tau to be PSD (add |min_ev|+epsilon)
    shift = max(0, -min_ev + 1e-10)
    tau_shifted = tau + shift * np.eye(D)
    evals_shifted = evals + shift

    # Sample random unit vectors
    X = rng.standard_normal((K_samples, D))
    X = X / np.linalg.norm(X, axis=1, keepdims=True)

    # Rayleigh quotients
    R = np.einsum('ki,ij,kj->k', X, tau_shifted, X)  # x^T tau_s x

    # Pairwise overlaps q = (x_i · x_j)^2
    # Subsample to 500 for efficiency
    idx = rng.choice(K_samples, min(500, K_samples), replace=False)
    X_sub = X[idx]
    dots = X_sub @ X_sub.T  # (500, 500)
    overlaps = dots**2
    # Upper triangle only
    tri_idx = np.triu_indices(len(idx), k=1)
    q_vals = overlaps[tri_idx]

    # Near-zero subspace overlaps
    if near_zero_thresh is not None:
        nz_mask = evals_shifted < near_zero_thresh
    else:
        # Adaptive: find spectral gap
        gaps = np.diff(evals_shifted)
        gap_idx = np.argmax(gaps[:D//2])  # largest gap in lower half
        nz_mask = np.zeros(D, dtype=bool)
        nz_mask[:gap_idx+1] = True

    n_nz = nz_mask.sum()
    if n_nz > 0 and n_nz < D:
        V_nz = evecs[:, nz_mask]  # D × n_nz
        # Project X_sub onto near-zero subspace
        X_proj = X_sub @ V_nz  # (500, n_nz)
        # Overlaps in projected space: q_nz = ||proj(x_i)||^2 * ||proj(x_j)||^2 * (x_i·x_j in subspace)^2
        dots_nz = X_proj @ X_proj.T
        q_nz = (dots_nz / (np.linalg.norm(X_proj, axis=1, keepdims=True) * np.linalg.norm(X_proj, axis=1)))**2 if n_nz > 0 else None
        q_nz_vals = q_nz[tri_idx] if q_nz is not None else np.array([])
    else:
        q_nz_vals = np.array([])
        n_nz = 0

    return q_vals, R, q_nz_vals, n_nz, evals_shifted, evecs

def analyze_pq(q_vals, label):
    """Print P(q) statistics."""
    if len(q_vals) == 0:
        print(f"  {label}: no overlaps")
        return
    bins = np.linspace(0, 1, 21)
    hist, _ = np.histogram(q_vals, bins=bins)
    hist = hist / hist.sum()

    # Find peaks
    peak_at_0 = hist[0]   # q ∈ [0, 0.05]
    peak_at_1 = hist[-1]  # q ∈ [0.95, 1]
    mid_mass = hist[4:-4].sum()  # q ∈ [0.2, 0.8]

    print(f"  {label}: {len(q_vals)} overlap pairs")
    print(f"    q_mean={q_vals.mean():.4f}, q_std={q_vals.std():.4f}")
    print(f"    P(q∈[0,0.05])={peak_at_0:.3f}, P(q∈[0.2,0.8])={mid_mass:.3f}, P(q∈[0.95,1])={peak_at_1:.3f}")
    if peak_at_0 > 0.5:
        interpretation = "RANDOM (most vectors are orthogonal — high-D concentration)"
    elif peak_at_1 > 0.5:
        interpretation = "ALIGNED (most vectors overlap — degenerate ground state)"
    elif mid_mass > 0.3:
        interpretation = "RUGGED (continuous P(q) — RSB, spin-glass-like)"
    else:
        interpretation = "BIMODAL (two-peak structure — simple landscape)"
    print(f"    Interpretation: {interpretation}")

if __name__ == "__main__":
    print("Spin-Glass P(q) Overlap Distribution — CvS Operator")
    print("======================================================")
    print()

    results = {}

    # Case 1: lsq=50 (near-resonant j≈8.8), multiple N values
    print("=== Case 1: lsq=50, N=20 and N=40 ===")
    for N, lsq in [(20, 50), (40, 50)]:
        t0 = time.time()
        tau = build_tau(N, lsq)
        D = 2*N+1
        evals = scipy.linalg.eigvalsh(tau)
        print(f"\nN={N} (D={D}): min_ev={evals[0]:.4e}, max_ev={evals[-1]:.4e}  ({time.time()-t0:.1f}s)")

        q_vals, R_vals, q_nz, n_nz, evs_s, vecs = compute_pq(tau, K_samples=1000)

        print(f"  Near-zero subspace: {n_nz} modes (gap-detected)")
        print(f"  Rayleigh quotient R: mean={R_vals.mean():.4e}, std={R_vals.std():.4e}, min={R_vals.min():.4e}")
        analyze_pq(q_vals, "Full overlap P(q)")
        if len(q_nz) > 0:
            analyze_pq(q_nz, "Near-zero subspace P(q)")

        results[f"lsq50_N{N}"] = {
            "D": D, "min_ev": float(evals[0]), "n_nz": int(n_nz),
            "q_mean": float(q_vals.mean()), "q_std": float(q_vals.std()),
            "R_mean": float(R_vals.mean()), "R_std": float(R_vals.std())
        }

    # Case 2: lsq=14, multiple N, shifted to PSD
    print("\n=== Case 2: lsq=14 (has negative eigenvalues → shifted to PSD) ===")
    for N in [8, 12, 16, 20]:
        t0 = time.time()
        tau = build_tau(N, 14)
        D = 2*N+1
        evals = scipy.linalg.eigvalsh(tau)
        shift = max(0, -evals[0] + 1e-8)
        print(f"\nN={N} (D={D}): min_ev={evals[0]:.4e}, shift={shift:.4e}  ({time.time()-t0:.1f}s)")

        q_vals, R_vals, q_nz, n_nz, evs_s, vecs = compute_pq(tau, K_samples=1000)
        print(f"  Near-zero (post-shift): {n_nz} modes")
        analyze_pq(q_vals, "Full overlap P(q)")
        if len(q_nz) > 0:
            analyze_pq(q_nz, "Near-zero P(q)")

        results[f"lsq14_N{N}"] = {
            "D": D, "min_ev": float(evals[0]), "shift": float(shift), "n_nz": int(n_nz),
            "q_mean": float(q_vals.mean()), "q_std": float(q_vals.std())
        }

    # Case 3: spectral_gap_analysis data (N=32, resonant)
    print("\n=== Case 3: Resonant N=32 from spectral_gap_analysis.json ===")
    sg = json.load(open("/Users/scones/projects/riemann/results/spectral_gap_analysis.json"))
    if sg:
        entry = sg[0]  # N=32
        evs_res = np.array(entry['smallest_30'])
        n_nz_res = entry.get('n_near_zero_1e4', 22)
        print(f"N={entry['N']}: {n_nz_res} near-zero eigenvalues (|λ|<1e-4)")
        print(f"  Spectral gap ratio: {entry['spectral_gap_ratio']:.3e}")
        print(f"  First_above_1e-4: {entry['first_above_1e4']:.4e}")
        print(f"  Near-zero eigenvalues: {evs_res[:n_nz_res]}")

        # The near-zero space has n_nz_res ≈ 22 modes.
        # For P(q) we need eigenvectors, not just eigenvalues.
        # The overlap structure in a D-dimensional space with k near-zero modes:
        # Random vector projected onto k-dim subspace: E[q] = k/D (by Marchenko-Pastur)
        D_res = 2*entry['N']+1
        k = n_nz_res
        E_q_random = k / D_res  # expected overlap for random vectors
        print(f"\n  D={D_res}, k={k} near-zero modes")
        print(f"  Expected random overlap in near-zero space: E[q] = k/D = {E_q_random:.4f}")
        print(f"  If P(q) peaks at q>{E_q_random:.2f}: near-zero modes are coherently aligned")
        print(f"  If P(q) ~ δ(q-k/D): no preferential alignment (Haar-random within subspace)")

    print()
    print("=== SUMMARY ===")
    print()
    print("P(q) structure across cases:")
    for key, r in results.items():
        print(f"  {key}: D={r['D']}, n_nz={r['n_nz']}, q_mean={r['q_mean']:.4f}, q_std={r['q_std']:.4f}")
    print()
    print("Key finding:")
    # For high-dimensional random vectors, E[q] = 1/D (concentration of measure)
    # If q_mean >> 1/D: there IS preferential alignment (ground state manifold structure)
    for key, r in results.items():
        D = r['D']
        E_q_random = 1.0 / D
        ratio = r['q_mean'] / E_q_random
        print(f"  {key}: q_mean/E[random] = {ratio:.2f}×  "
              f"({'ALIGNED (>2×)' if ratio > 2 else 'random-like (≈1×)'})")

    out = "/Users/scones/projects/riemann/results/spin_glass_pq.json"
    with open(out, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"\nSaved → {out}")
