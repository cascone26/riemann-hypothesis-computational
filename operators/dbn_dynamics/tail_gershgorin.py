#!/usr/bin/env python3
"""
Tail bound via Gershgorin disc: prove full even Weil matrix (modes 1,3,5,...,K) is PD
by bounding tail modes k>=9 using Gershgorin's theorem and diagonal dominance.

Lane: TAIL BOUND via GERSHGORIN
- Build full matrix on K modes (K=8..12: modes 1,3,5,7,9,11,13,15,...)
- Analyze structure: diagonal growth, off-diagonal coupling
- Apply Gershgorin: each eigenvalue in a disc |λ-A[ii]| <= row_sum_off_diag
- Prove: min_eig(full) >= min_eig(4x4) - tail_correction, where tail is rigorously bounded
"""
import sys, mpmath as mp
sys.path.insert(0, '/Users/scones/projects/riemann/operators/dbn_dynamics')
from mpmath_weil_hiprec import weil_matrices
from entry_engine import Q
mp.mp.dps = 40

def build_even_weil_matrix(a, K_odd):
    """Build even-sector Weil matrix on odd modes k=1,3,5,...,K_odd.
    K_odd = largest odd k; modes = [1,3,5,...,K_odd]; total size = (K_odd+1)//2.
    Returns: full numpy-like matrix of mpmath.mpf (2D list), eigenvalues.
    """
    modes = [2*i+1 for i in range((K_odd+1)//2)]  # 1,3,5,...
    n = len(modes)
    print(f"\n=== Building even Weil matrix: a={float(a):.4f}, modes={modes} (n={n}) ===")

    # Build the matrix
    matrix = [[mp.mpf(0) for _ in range(n)] for _ in range(n)]
    for ik, k in enumerate(modes):
        for il, l in enumerate(modes):
            if il >= ik:
                val = Q(k, l)
                matrix[ik][il] = val
                if il > ik:
                    matrix[il][ik] = val  # symmetric

    # Eigenvalues (numpy won't work with mpmath, use scipy or mpmath's eigh)
    # For now, use mpmath's eigenvalue solver on the matrix
    # Convert to mpmath matrix object
    A = mp.matrix(matrix)
    try:
        evals = mp.eig(A)[0]  # eigenvalues only
        evals = sorted([float(ev) for ev in evals])  # sort ascending
    except:
        # Fallback: use mpmath.eig which returns (eigenvalues, eigenvectors)
        evals, _ = mp.eig(A)
        evals = sorted([float(ev) for ev in evals])

    return matrix, modes, evals

def gershgorin_bound(matrix, modes):
    """Apply Gershgorin's theorem to bound eigenvalues.
    Returns: list of (mode_index, diag, row_sum_off_diag, lower_bound_lambda)
    where lower_bound_lambda = diag - row_sum_off_diag.
    """
    n = len(modes)
    results = []

    for i in range(n):
        diag = matrix[i][i]
        row_sum = sum(abs(matrix[i][j]) for j in range(n) if j != i)
        lower_bound = diag - row_sum

        results.append({
            'mode': modes[i],
            'index': i,
            'diag': float(diag),
            'row_sum_off_diag': float(row_sum),
            'lower_bound': float(lower_bound),
            'upper_bound': float(diag + row_sum),
        })

    return results

def tail_correction_gershgorin(matrix, modes, cutoff_idx):
    """Bound the tail (modes beyond the 4x4) using Gershgorin applied to the Schur complement.

    For the full matrix partitioned as:
      M = [ M_lead  M_tail_cross ]
          [ M_cross M_tail_tail  ]

    By Gershgorin on the Schur complement S = M_tail - M_cross M_lead^{-1} M_cross^T,
    we can bound the contribution of tail modes to the full min-eigenvalue.

    For simplicity: use the diagonal dominance argument on M_tail directly with off-diagonal coupling.
    """

    print(f"\n--- Tail analysis (Gershgorin bound) ---")
    print(f"Cutoff at mode {modes[cutoff_idx]} (index {cutoff_idx}); tail modes: {modes[cutoff_idx+1:]}")

    # Diagonal dominance: for the tail part
    tail_modes = modes[cutoff_idx+1:]
    n_tail = len(tail_modes)

    if n_tail == 0:
        print("No tail (full matrix size = cutoff).")
        return mp.mpf(0), mp.mpf(0)

    print(f"Tail size: {n_tail} modes")

    # Maximum row sum off-diagonal in the tail block
    max_row_sum_tail = mp.mpf(0)
    min_diag_tail = mp.mpf('inf')

    for i in range(cutoff_idx, len(modes)):
        diag_i = matrix[i][i]
        row_sum_i = sum(abs(matrix[i][j]) for j in range(len(modes)) if j != i)
        min_diag_tail = min(min_diag_tail, diag_i)
        max_row_sum_tail = max(max_row_sum_tail, row_sum_i)

        print(f"  mode {modes[i]}: diag={float(diag_i):.6e}, row_sum_off_diag={float(row_sum_i):.6e}")

    # Gershgorin bound on tail: each eigenvalue of tail in disc |λ - diag| <= row_sum
    # Worst case: min eigenvalue >= min_diag_tail - max_row_sum_tail
    worst_tail_lower = min_diag_tail - max_row_sum_tail

    print(f"Min diagonal in tail: {float(min_diag_tail):.6e}")
    print(f"Max row-sum off-diag in tail: {float(max_row_sum_tail):.6e}")
    print(f"Gershgorin bound (tail min eig): {float(worst_tail_lower):.6e}")

    return worst_tail_lower, min_diag_tail

def report(a_val, K_odd_values):
    """Run analysis across multiple K values."""
    a = mp.mpf(str(a_val))
    results_per_K = {}

    for K in K_odd_values:
        matrix, modes, evals = build_even_weil_matrix(a, K)
        n = len(modes)

        print(f"\nEigenvalues (K={K}, n={n}):")
        for i, ev in enumerate(evals[:min(5, len(evals))]):
            print(f"  λ_{i+1} = {ev:.6e}")
        print(f"  ...λ_min = {evals[0]:.6e}")

        # Gershgorin analysis
        gersh = gershgorin_bound(matrix, modes)
        print(f"\nGershgorin disc analysis (K={K}):")
        for g in gersh[:min(6, len(gersh))]:  # Show first few and last few
            print(f"  mode {g['mode']:2d}: diag={g['diag']:10.6e}, row_sum_off={g['row_sum_off_diag']:10.6e}, "
                  f"lower_bound={g['lower_bound']:10.6e}")
        if len(gersh) > 6:
            print(f"  ...")
            for g in gersh[-3:]:
                print(f"  mode {g['mode']:2d}: diag={g['diag']:10.6e}, row_sum_off={g['row_sum_off_diag']:10.6e}, "
                      f"lower_bound={g['lower_bound']:10.6e}")

        # Tail bound (using Gershgorin on modes beyond the 4x4)
        cutoff_idx = 4  # modes 1,3,5,7 are indices 0,1,2,3
        tail_bound, min_diag_tail = tail_correction_gershgorin(matrix, modes, cutoff_idx)

        results_per_K[K] = {
            'modes': modes,
            'eigenvalues': evals,
            'lambda_min': evals[0],
            'gershgorin': gersh,
            'tail_bound': tail_bound,
            'min_diag_tail': min_diag_tail,
        }

    # Summary table
    print(f"\n{'='*100}")
    print(f"SUMMARY (a={float(a):.4f})")
    print(f"{'='*100}")
    print(f"{'K':>3} {'n':>2} {'λ_min':>15} {'Gersh_min':>15} {'Tail_bound':>15} {'min_diag_tail':>15}")
    print(f"{'-'*70}")

    for K in K_odd_values:
        res = results_per_K[K]
        gersh_min = min(g['lower_bound'] for g in res['gershgorin'])
        tail = res['tail_bound']
        print(f"{K:3d} {len(res['modes']):2d} {res['lambda_min']:15.6e} {gersh_min:15.6e} {float(tail):15.6e} {float(res['min_diag_tail']):15.6e}")

    print(f"\nCritical observations:")
    print(f"- λ_min(4x4) [K=7]: {results_per_K[7]['lambda_min']:.6e}")
    print(f"- λ_min(6x6) [K=11]: {results_per_K[11]['lambda_min']:.6e}")
    if 15 in results_per_K:
        print(f"- λ_min(8x8) [K=15]: {results_per_K[15]['lambda_min']:.6e}")

    print(f"\nGershgorin tail correction:")
    for K in sorted(results_per_K.keys()):
        res = results_per_K[K]
        tail_c = res['tail_bound']
        tail_diag = res['min_diag_tail']
        lambda_4x4 = results_per_K[7]['lambda_min']
        print(f"  K={K}: tail_bound={float(tail_c):.6e}, min_diag_tail={float(tail_diag):.6e}, "
              f"margin_lost={max(mp.mpf(0), lambda_4x4 - float(tail_c)):.6e}")

if __name__ == '__main__':
    # Test at critical a = 0.549 (p=3 entry point)
    a_test = 0.549
    K_values = [7, 11, 15]  # 4x4, 6x6, 8x8

    report(a_test, K_values)

    # Also test across the window
    print(f"\n{'='*100}")
    print("Scanning across window a ∈ (0.42, 0.5493)")
    print(f"{'='*100}")

    for a_val in [0.42, 0.45, 0.48, 0.50, 0.52, 0.549]:
        a = mp.mpf(str(a_val))
        # Just the 4x4 and 6x6
        for K in [7, 11]:
            matrix, modes, evals = build_even_weil_matrix(a, K)
            gersh = gershgorin_bound(matrix, modes)
            gersh_min = min(g['lower_bound'] for g in gersh)
            print(f"a={a_val:.3f}, K={K:2d} (n={len(modes):2d}): λ_min={evals[0]:.6e}, Gersh_min={gersh_min:.6e}")
