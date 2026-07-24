#!/usr/bin/env python3
"""
§199.4: SCHUR COMPLEMENT TAIL BOUND
-----------------------------------
Prove p=2-only Weil positivity beyond the 4×4 (k=1,3,5,7) by Schur complement.
Build even matrix on modes k=1,3,5,7,...,up to K, partition into:
  - Q_44: the 4×4 block (k=1,3,5,7)
  - Q_rr: the tail block (k≥9)
  - Q_4r: the coupling
Then Schur(Q_4r, Q_rr) = Q_44 - Q_4r * Q_rr^{-1} * Q_4r^T remains PD iff:
  (a) Q_rr is PD (eigenvalues bounded below)
  (b) coupling norm ||Q_4r|| is bounded
  (c) min_eig(Q_44) > ||Q_4r||²/min_eig(Q_rr)

Report: Schur-corrected minimum eigenvalue vs original 4×4 minimum,
        and whether the correction closes rigorously.
"""
import sys, math
import mpmath as mp
sys.path.insert(0, '/Users/scones/projects/riemann/operators/dbn_dynamics')
from mpmath_weil_hiprec import weil_matrices

# Settings
mp.mp.dps = 40
K_max = 20                    # modes up to k=2*K_max-1 odd indices
M = 800
even_modes = [c for c in range(K_max) if (c+1) % 2 == 1]  # k=1,3,5,...

def block_extract(Q, Mass, idx):
    """Extract block of Q and Mass for mode indices idx."""
    n = len(idx)
    B = mp.matrix(n, n)
    M_block = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            B[i, j] = Q[idx[i], idx[j]]
            M_block[i, j] = Mass[idx[i], idx[j]]
    return B, M_block

def gen_eig_sym(Q, Mass):
    """Lowest eigenpair of Q c = λ Mass c."""
    L = mp.cholesky(Mass)
    Linv = L ** (-1)
    A = Linv * Q * Linv.T
    A = (A + A.T) / 2
    E, V = mp.eigsy(A)
    return E[0], (Linv.T * V)[:, 0]  # return lowest eigenvalue and corresponding eigenvector

def matrix_norm(M):
    """Spectral norm (largest singular value) of matrix M."""
    n = M.rows
    # Compute M^T M and find its largest eigenvalue
    MtM = M.T * M
    eigs, _ = mp.eigsy(MtM)
    return mp.sqrt(eigs[-1])

def schur_complement_analysis(a, K_list):
    """
    For each K in K_list:
      - Extract 4×4 block (k=1,3,5,7)
      - Extract tail block (k≥9)
      - Compute Schur complement: S = Q_44 - Q_4r * Q_rr^{-1} * Q_4r^T
      - Report min eigenvalues and bounds.
    """
    am = mp.mpf(str(a))
    Qf, Qa, Mass, Vv, xs, w = weil_matrices(am, K_max, M)

    print(f"\n{'='*90}")
    print(f"a = {a:.4f}  (p=2-only zone)")
    print(f"{'='*90}")
    print(f"{'K':>4} {'modes in tail':>18} {'λ_min(4×4)':>16} {'λ_min(tail)':>16} "
          f"{'||Q_4r||':>12} {'Schur λ_min':>16} {'bound closes':>14}")
    print(f"{'-'*110}")

    results = []

    for K in K_list:
        n_modes = 2 * K - 1  # highest odd index
        modes_to_use = [i for i in range(K_max) if i < K and (i+1) % 2 == 1]

        if len(modes_to_use) <= 4:
            print(f"{K:>4}  (skip: too few modes)")
            continue

        # Extract blocks
        idx_44 = [0, 1, 2, 3]  # k=1,3,5,7 (indices in modes_to_use)
        idx_tail = [i for i in range(4, len(modes_to_use))]  # k≥9

        if len(idx_tail) == 0:
            print(f"{K:>4}  (no tail modes)")
            continue

        Q_44, M_44 = block_extract(Qf, Mass, [modes_to_use[i] for i in idx_44])
        Q_rr, M_rr = block_extract(Qf, Mass, [modes_to_use[i] for i in idx_tail])

        # Coupling: rows from {1,3,5,7}, cols from tail
        Q_4r = mp.matrix(4, len(idx_tail))
        for i in range(4):
            for j in range(len(idx_tail)):
                Q_4r[i, j] = Qf[modes_to_use[idx_44[i]], modes_to_use[idx_tail[j]]]

        # Compute generalized eigenvalues for the 4×4 and tail
        lam_44, _ = gen_eig_sym(Q_44, M_44)
        lam_rr, _ = gen_eig_sym(Q_rr, M_rr)

        # Schur complement (in original inner product)
        # S = Q_44 - Q_4r M_rr^{-1} Q_4r^T (generalized inverse, via Cholesky)
        try:
            L_rr = mp.cholesky(M_rr)
            L_rr_inv = L_rr ** (-1)

            # Transform to make M_rr = I: Q_rr' = L_rr_inv Q_rr L_rr_inv^T, Q_4r' = Q_4r L_rr_inv^T
            Q_rr_white = L_rr_inv * Q_rr * L_rr_inv.T
            Q_4r_white = Q_4r * L_rr_inv.T

            # Schur in whitened space: S' = Q_44 - Q_4r' (Q_rr')^{-1} Q_4r'^T
            Q_rr_white_inv = Q_rr_white ** (-1)
            S_white = Q_44 - Q_4r_white * Q_rr_white_inv * Q_4r_white.T

            # Eigenvalues of S wrt M_44: same as of M_44^{-1/2} S M_44^{-1/2}
            L_44 = mp.cholesky(M_44)
            L_44_inv = L_44 ** (-1)
            S_white_44 = L_44_inv * S_white * L_44_inv.T
            E_S, _ = mp.eigsy(S_white_44)
            lam_S = E_S[0]

            # Bounds for closure
            norm_Q_4r = matrix_norm(Q_4r_white)  # spectral norm in the whitened basis
            coupling_contribution = norm_Q_4r ** 2 / lam_rr  # upper bound on Schur correction's effect
            bound_closes = (lam_44 - coupling_contribution > 0)

            tail_str = f"k={modes_to_use[idx_tail[0]]+1}..{modes_to_use[idx_tail[-1]]+1}"
            print(f"{K:>4}  {tail_str:>18} {mp.nstr(lam_44, 6):>16} {mp.nstr(lam_rr, 6):>16} "
                  f"{mp.nstr(norm_Q_4r, 4):>12} {mp.nstr(lam_S, 6):>16} "
                  f"{'YES' if bound_closes else 'NO':>14}")

            results.append({
                'K': K,
                'modes': modes_to_use,
                'lam_44': lam_44,
                'lam_rr': lam_rr,
                'norm_Q_4r': norm_Q_4r,
                'lam_S': lam_S,
                'coupling_contrib': coupling_contribution,
                'closes': bound_closes
            })
        except Exception as e:
            print(f"{K:>4}  (error: {str(e)[:40]})")

    return results

def main():
    a_list = [0.44, 0.48, 0.52, 0.549]
    K_range = list(range(5, 12))  # K=5..11 (modes up to k=1,3,5,7,9,11,13,15,17,19,21)

    all_results = {}

    for a in a_list:
        results = schur_complement_analysis(a, K_range)
        all_results[a] = results

    print(f"\n{'='*90}")
    print("SUMMARY: Schur complement tail-bound closure across the p=2-only zone")
    print(f"{'='*90}")

    for a in a_list:
        if a in all_results and all_results[a]:
            results = all_results[a]
            all_close = all(r['closes'] for r in results)
            print(f"\na={a:.4f}: {'CLOSES' if all_close else 'DOES NOT CLOSE (see details above)'}")
            for r in results:
                print(f"  K={r['K']:2d}: λ_min(S)={mp.nstr(r['lam_S'], 6)}, bound={mp.nstr(r['coupling_contrib'], 6)}, "
                      f"margin={mp.nstr(r['lam_44'] - r['coupling_contrib'], 6)}")

if __name__ == '__main__':
    main()
