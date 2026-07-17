#!/usr/bin/env python3
"""
tda_eigenvectors.py — PRISM EXTREME item #3: Persistent homology / TDA.

Computes persistent homology on the near-zero eigenvectors of tau_N across
multiple N values. Looks for topological invariants (Betti numbers,
persistence-diagram features) invisible to scalar statistics.

Two TDA passes:
A. EIGENVALUE SPECTRUM TDA: treat {λ_i} sorted eigenvalues as a 1D point cloud.
   Sublevel-set persistence captures the "staircase" structure of the spectrum.

B. EIGENVECTOR COMPONENT TDA: treat the near-zero eigenvectors (at different N)
   as point clouds in R^k. Vietoris-Rips persistent homology reveals topology
   of the eigenvector shape manifold.

C. GRAM MATRIX TDA: use gram_N as a precomputed distance matrix for TDA.
   gram_N[i,j] = inner product between eigenvector components at Lambda_sq=i,j.
   Convert to distance via D = sqrt(gram[i,i] + gram[j,j] - 2*gram[i,j]).
"""

import numpy as np
import json
import time
import math
import scipy.linalg
from scipy.special import digamma as _digamma
from ripser import ripser
from persim import plot_diagrams

TWO_PI = 2*math.pi
LOGPI  = math.log(math.pi)
LAMBDA_SQ = 14
L = 2 * math.log(math.sqrt(LAMBDA_SQ))

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

def build_tau(N, G=6000):
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

    lkpp, wts_p = build_sieve(LAMBDA_SQ)
    K = len(lkpp)
    S = np.zeros(D); S_diag = np.zeros(D)
    for k0 in range(0, K, 2000):
        k1 = min(k0+2000, K)
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

def persistence_summary(dgm, dim_name):
    """Print summary of a persistence diagram."""
    if len(dgm) == 0:
        print(f"    {dim_name}: no features")
        return []
    finite = dgm[np.isfinite(dgm[:,1])]
    lifetimes = (finite[:,1] - finite[:,0]) if len(finite) > 0 else np.array([])
    inf_count = np.sum(~np.isfinite(dgm[:,1]))
    max_lt = f"{lifetimes.max():.4e}" if len(lifetimes) > 0 else "n/a"
    mean_lt = f"{lifetimes.mean():.4e}" if len(lifetimes) > 0 else "n/a"
    print(f"    {dim_name}: {len(dgm)} features ({inf_count} infinite), "
          f"max lifetime={max_lt}, mean lifetime={mean_lt}")
    if len(lifetimes) > 0:
        top = np.argsort(lifetimes)[-min(5, len(lifetimes)):][::-1]
        print(f"      Top features by lifetime:")
        for idx in top:
            b, d = finite[idx]
            print(f"        born={b:.4e}, died={d:.4e}, lifetime={d-b:.4e}")
    return lifetimes.tolist()

if __name__ == "__main__":
    print("Persistent Homology / TDA — CvS Operator")
    print("==========================================")
    print()

    # === PASS A: Eigenvalue spectrum TDA ===
    print("=== PASS A: Sublevel-set persistence on eigenvalue spectra ===")
    print("(1D point cloud = sorted eigenvalues; captures spectral staircase)")
    print()

    evsA_all = {}
    for N in [4, 6, 8, 10, 12, 16, 20, 24]:
        t0 = time.time()
        tau = build_tau(N)
        evs = scipy.linalg.eigvalsh(tau)
        evsA_all[N] = evs
        print(f"N={N:>3} (D={2*N+1:>3}): eigenvalues computed ({time.time()-t0:.1f}s)")

    print()
    results_A = {}
    for N, evs in evsA_all.items():
        # Treat eigenvalues as a 1D point cloud: Vietoris-Rips on {λ_i}
        # H0 lifetime = gap between consecutive eigenvalues (merge times)
        pts = evs.reshape(-1, 1)
        dgms = ripser(pts, maxdim=0)['dgms']
        H0 = dgms[0]
        finite_H0 = H0[np.isfinite(H0[:,1])]
        lifetimes = (finite_H0[:,1] - finite_H0[:,0]).tolist() if len(finite_H0) > 0 else []
        # Largest H0 lifetime = largest consecutive eigenvalue gap (topological spectral gap)
        spectral_gap = max(lifetimes) if lifetimes else 0
        # Find WHERE in the spectrum the largest gap is
        gaps = np.diff(np.sort(evs))
        gap_pos_idx = np.argmax(gaps)
        gap_position = (evs[gap_pos_idx] + evs[gap_pos_idx+1]) / 2  # midpoint of largest gap
        # Negative eigenvalues (arch+B modes not cancelled by P)
        n_neg = np.sum(evs < -1e-10)
        print(f"N={N:>3}: #neg_ev={n_neg}, largest_gap={gaps.max():.4e} at λ≈{gap_position:.4f}, "
              f"spectral_gap_TDA={spectral_gap:.4e}")
        results_A[N] = {
            "H0_count": len(H0),
            "spectral_gap_topology": spectral_gap,
            "largest_gap": float(gaps.max()),
            "gap_position": float(gap_position),
            "n_negative_ev": int(n_neg),
            "lifetimes_top5": sorted(lifetimes, reverse=True)[:5]
        }

    print()
    print("=== PASS B: Eigenvector point cloud TDA (near-zero eigenvectors) ===")
    print("(Each near-zero eigenvector is a point in R^D; TDA on the manifold of shapes)")
    print()

    # Collect near-zero eigenvectors across N values
    # Stack them as rows in a point cloud
    near_zero_vecs = []
    vec_labels = []
    for N, evs in evsA_all.items():
        tau = build_tau(N)
        evals, evecs = scipy.linalg.eigh(tau)
        D = 2*N+1
        # Near-zero: |eigenvalue| < 1e-3
        near_zero_idx = np.where(np.abs(evals) < 1e-3)[0]
        if len(near_zero_idx) > 0:
            # Pad to the largest D (2*24+1=49) for uniform embedding
            for idx in near_zero_idx:
                v = evecs[:, idx]
                # Center to D=49 by zero-padding
                v_padded = np.zeros(49)
                offset = (49 - D) // 2
                v_padded[offset:offset+D] = v
                near_zero_vecs.append(v_padded)
                vec_labels.append(N)

    if len(near_zero_vecs) > 1:
        X = np.array(near_zero_vecs)
        print(f"Point cloud: {X.shape[0]} near-zero eigenvectors in R^49")
        print(f"(From N ∈ {sorted(set(vec_labels))})")
        print()

        # Persistent homology on eigenvector manifold
        dgms_B = ripser(X, maxdim=1)['dgms']
        print(f"  Persistent homology:")
        lifetimes_B0 = persistence_summary(dgms_B[0], "H0 (connected components)")
        lifetimes_B1 = persistence_summary(dgms_B[1], "H1 (loops/cycles)")

        # H0: number of connected components at death = number of clusters
        B0_finite = dgms_B[0][np.isfinite(dgms_B[0][:,1])]
        if len(B0_finite) > 0:
            largest_gap = (B0_finite[:,1] - B0_finite[:,0]).max()
            print(f"  Largest H0 gap: {largest_gap:.4e} — "
                  f"{'strong clustering' if largest_gap > 0.1 else 'continuous manifold'}")
    else:
        print("Not enough near-zero eigenvectors for manifold TDA")
        dgms_B = None

    print()
    print("=== PASS C: Gram matrix TDA (existing bvec/gram data) ===")
    results_dir = "/Users/scones/projects/riemann/results/"
    gram_results = {}
    for n_pts in [300, 1000, 2000]:
        try:
            bvec = np.load(f"{results_dir}bvec_{n_pts}.npy")
            gram = np.load(f"{results_dir}gram_{n_pts}.npy")

            # gram is n_pts × n_pts. Interpret as inner product matrix.
            # Distance: D[i,j] = sqrt(g[i,i] + g[j,j] - 2*g[i,j])
            diag = np.diag(gram)
            D_mat = np.sqrt(np.maximum(diag[:,None] + diag[None,:] - 2*gram, 0))

            # For TDA with large matrices, subsample
            max_pts = 200
            if n_pts > max_pts:
                idx = np.linspace(0, n_pts-1, max_pts).astype(int)
                D_sub = D_mat[np.ix_(idx, idx)]
            else:
                D_sub = D_mat

            t0 = time.time()
            dgms_C = ripser(D_sub, metric="precomputed", maxdim=1)['dgms']
            dt = time.time() - t0
            print(f"\nGram N={n_pts} (subsample {D_sub.shape[0]}×{D_sub.shape[0]}, {dt:.1f}s):")
            lC0 = persistence_summary(dgms_C[0], "  H0")
            lC1 = persistence_summary(dgms_C[1], "  H1")
            gram_results[n_pts] = {
                "H0_top5": sorted(lC0, reverse=True)[:5] if lC0 else [],
                "H1_top5": sorted(lC1, reverse=True)[:5] if lC1 else []
            }
        except FileNotFoundError:
            print(f"Gram N={n_pts}: file not found, skip")

    print()
    print("=== SUMMARY ===")
    print("Key TDA findings:")
    print(f"  Pass A (spectral): Topological spectral gap vs N")
    for N, r in results_A.items():
        print(f"    N={N:>3}: #neg_ev={r['n_negative_ev']}, largest_gap={r['largest_gap']:.4e} "
              f"at λ≈{r['gap_position']:.4f}")

    # Save results
    save_data = {
        "LAMBDA_SQ": LAMBDA_SQ,
        "pass_A": {str(k): v for k, v in results_A.items()},
        "pass_C_gram": {str(k): v for k, v in gram_results.items()}
    }
    out = "/Users/scones/projects/riemann/results/tda_eigenvectors.json"
    with open(out, "w") as fh:
        json.dump(save_data, fh, indent=2)
    print(f"\nSaved → {out}")
