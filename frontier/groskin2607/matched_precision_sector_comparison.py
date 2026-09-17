#!/usr/bin/env python3
"""
MATCHED-PRECISION SECTOR COMPARISON

At the EXACT (c,N) points margin_law.py used, compute:
1. The full matrix's smallest eigenvalue (via margin_law.py's method: build_arb_tau + mpmath eigsy)
2. The sector-split min eigenvalues (via same Arb precision)
3. Compare: does min(λ_even, λ_odd) equal λ_full?

This resolves whether the parity-sector decomposition is exact or if there's a
numerical/algorithmic discrepancy.
"""
import sys
import math
import numpy as np
from flint import arb, acb, arb_mat, ctx
import mpmath as mp

sys.path.insert(0, '/Users/scones/projects/riemann/frontier/groskin2607')
from arb_ldlt_certify import build_arb_tau

def spectrum_full(c, N, prec):
    """Compute full matrix's smallest eigenvalue at given Arb precision.

    Replicates margin_law.py's method exactly.
    """
    ctx.prec = prec
    A, DIM = build_arb_tau(c, N, prec)

    # Convert to mpmath, same as margin_law.py does
    mp.mp.dps = int(prec * 0.29)
    Amp = mp.matrix(DIM, DIM)
    for i in range(DIM):
        for j in range(DIM):
            Amp[i, j] = mp.mpf(A[i, j].str(int(prec*0.28), radius=False))

    # Eigenvalues via mpmath eigsy
    E = sorted(float(mp.log10(abs(e))) for e in mp.eigsy(Amp, eigvals_only=True))

    # Return smallest eigenvalue (in real space, not log10)
    return E[0], E, Amp  # (log10 of smallest, list of all log10 eigs, mpmath matrix)

def extract_parity_blocks_mpmath(A_mp, N):
    """Extract even/odd blocks from mpmath matrix."""
    dim = 2*N+1
    even_idx = [i for i in range(dim) if (i-N) % 2 == 0]
    odd_idx = [i for i in range(dim) if (i-N) % 2 != 0]

    A_even = mp.matrix(len(even_idx), len(even_idx))
    A_odd = mp.matrix(len(odd_idx), len(odd_idx))

    for ii, i in enumerate(even_idx):
        for jj, j in enumerate(even_idx):
            A_even[ii, jj] = A_mp[i, j]

    for ii, i in enumerate(odd_idx):
        for jj, j in enumerate(odd_idx):
            A_odd[ii, jj] = A_mp[i, j]

    return A_even, A_odd, even_idx, odd_idx

def spectrum_sectors(A_mp, N):
    """Compute sector min eigenvalues from mpmath matrix."""
    A_even, A_odd, even_idx, odd_idx = extract_parity_blocks_mpmath(A_mp, N)

    E_even_raw = sorted(float(mp.log10(abs(e))) for e in mp.eigsy(A_even, eigvals_only=True))
    E_odd_raw = sorted(float(mp.log10(abs(e))) for e in mp.eigsy(A_odd, eigvals_only=True))

    return E_even_raw[0], E_odd_raw[0], E_even_raw, E_odd_raw

def main():
    print("="*110)
    print("MATCHED-PRECISION SECTOR COMPARISON")
    print("="*110)
    print("\nComputing at the EXACT (c,N) points from margin_law.py with prec=1400 bits:\n")

    # These are margin_law.py's exact test points
    test_points = [
        (30, 4),
        (30, 6),
        (30, 8),
        (30, 10),
        (30, 12),
    ]

    prec = 1400  # Same as margin_law.py uses

    print(f"{'c':>3} {'N':>2} {'log10(λ_full)':>15} {'log10(λ_even)':>15} {'log10(λ_odd)':>15} "
          f"{'min(even,odd)':>15} {'Δ=full-min':>12} {'which wins':>10}")
    print("-"*110)

    for c, N in test_points:
        # Compute full matrix spectrum
        log10_full, E_full, A_mp = spectrum_full(c, N, prec)

        # Compute sector spectra
        log10_even, log10_odd, E_even, E_odd = spectrum_sectors(A_mp, N)

        # min of sectors
        log10_min_sectors = min(log10_even, log10_odd)
        which_wins = "even" if log10_even <= log10_odd else "odd"

        # Difference
        delta = log10_full - log10_min_sectors

        print(f"{c:>3} {N:>2} {log10_full:>15.2f} {log10_even:>15.2f} {log10_odd:>15.2f} "
              f"{log10_min_sectors:>15.2f} {delta:>12.2f} {which_wins:>10}")

    print("\n" + "="*110)
    print("INTERPRETATION")
    print("="*110)
    print("""
If Δ ≈ 0 (< 0.1): The decomposition is EXACT. min(λ_even, λ_odd) = λ_full within numerical error.
  → Puzzle is NOT a block-diagonalization issue
  → Need to explain why sector exponents (0.60-0.63) differ from full exponent (0.44)

If Δ is large (> 0.5): There's a DISCREPANCY between block-split and full computation
  → Possible: one of the constructions has a subtle bug or different indexing
  → Need to investigate which sector's minimum is being computed incorrectly

If Δ shows a PATTERN (e.g., alternates or grows with N): The sectors may be COUPLED
  → The eigenvector for the global minimum may NOT respect the even/odd parity exactly
  → This would explain why min(sector eigs) ≠ global min eig
""")

if __name__ == '__main__':
    main()
