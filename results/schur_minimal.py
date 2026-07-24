#!/usr/bin/env python3
"""
§199.4 MINIMAL SCHUR COMPLEMENT TEST
-------------------------------------
Quick test: for a representative a, build matrices incrementally from K=4 (4×4 block)
to K=6,7,8,... and track how min-eigenvalue changes with tail.

Uses same infrastructure as lowdim_coercivity.py but tracks absolute eigenvalues
across increasing K instead of just truncated dimensions.
"""
import sys, math
import mpmath as mp
sys.path.insert(0, '/Users/scones/projects/riemann/operators/dbn_dynamics')
from mpmath_weil_hiprec import weil_matrices

mp.mp.dps = 40
K_max = 12  # only go up to K=12 (modes up to k=23)
M = 800
even_modes = [c for c in range(K_max) if (c+1) % 2 == 1]

def gen_eig_lowest(Q, Mass):
    """Lowest eigenvalue of Q c = λ Mass c."""
    try:
        L = mp.cholesky(Mass)
        Linv = L ** (-1)
        A = Linv * Q * Linv.T
        A = (A + A.T) / 2
        E, _ = mp.eigsy(A)
        return E[0]
    except:
        return None

def block_extract(Q, Mass, idx):
    """Extract block for indices."""
    n = len(idx)
    B = mp.matrix(n, n)
    M_block = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            B[i, j] = Q[idx[i], idx[j]]
            M_block[i, j] = Mass[idx[i], idx[j]]
    return B, M_block

def minimal_test(a_list):
    """Build matrices incrementally and track eigenvalues."""
    for a in a_list:
        am = mp.mpf(str(a))
        print(f"\n{'='*90}")
        print(f"a = {a:.4f}  (K_max = {K_max}, modes up to k={2*K_max-1})")
        print(f"{'='*90}")
        print(f"{'K':>3} {'n_odd_modes':>8} {'modes incl.':>15} {'λ_min(full)':>18} {'diff from K-1':>18}")
        print(f"{'-'*90}")

        Qf, Qa, Mass, Vv, xs, w = weil_matrices(am, K_max, M)

        prev_lam = None
        for K in range(1, K_max + 1):
            modes_to_use = even_modes[:K]
            if len(modes_to_use) == 0:
                continue

            # Extract blocks
            Q_K, M_K = block_extract(Qf, Mass, modes_to_use)
            lam = gen_eig_lowest(Q_K, M_K)

            if lam is not None:
                if K <= 4:
                    kind = "4×4 block"
                else:
                    kind = f"tail K={K}"

                diff_str = "---"
                if prev_lam is not None:
                    diff = lam - prev_lam
                    diff_str = f"{mp.nstr(diff, 6)}"

                modes_str = f"k=1..{2*len(modes_to_use)-1}"
                print(f"{K:>3} {len(modes_to_use):>8} {modes_str:>15} {mp.nstr(lam, 10):>18} {diff_str:>18}")
                prev_lam = lam

        # Final assessment: is the 4×4 (K=4) minimum stayed positive with tail?
        Q_4x4, M_4x4 = block_extract(Qf, Mass, even_modes[:4])
        lam_4x4 = gen_eig_lowest(Q_4x4, M_4x4)

        Q_full, M_full = block_extract(Qf, Mass, even_modes)
        lam_full = gen_eig_lowest(Q_full, M_full)

        if lam_4x4 and lam_full:
            correction = lam_4x4 - lam_full
            print(f"\n  λ_min(4×4 alone):         {mp.nstr(lam_4x4, 6)}")
            print(f"  λ_min(full K={K_max}):       {mp.nstr(lam_full, 6)}")
            print(f"  Tail correction:          {mp.nstr(correction, 6)}")
            if lam_full > 0:
                print(f"  → FULL FORM STAYS POSITIVE ✓")
            else:
                print(f"  → WARNING: full form went negative!")

if __name__ == '__main__':
    a_list = [0.44, 0.48, 0.52, 0.549]
    minimal_test(a_list)
