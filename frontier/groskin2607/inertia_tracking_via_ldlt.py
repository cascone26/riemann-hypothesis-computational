#!/usr/bin/env python3
"""
Track inertia (number of negative/positive eigenvalues) as we homotopy τ(t) = (W02-WR) - t·Wp.

This is MORE ROBUST than tracking eigenvalues directly, since Sylvester's criterion
(counting LDL^T pivots) is numerically stable even when eigenvalues are tiny.

Key idea:
- At t=0: τ(0) = W02 - WR has negative eigenvalues (known from S2636: τ_even has 2, τ_odd has 4)
- At t=1: τ(1) = W02 - WR - Wp = τ is PSD (0 negative eigenvalues)
- In between: something happens that flips all the signs

By tracking WHERE the sign flips occur (which pivots in the LDL^T become zero or change sign),
we can understand the mechanism better than raw eigenvalues.
"""

import sys
import argparse
import numpy as np
from pathlib import Path

sys.path.insert(0, '/Users/scones/projects/riemann/frontier/groskin2607')

from flint import arb, ctx
import mpmath as mp


def build_arb_tau_homotopy(c, N, t_weight, prec):
    """
    Build τ(t) = (W02 - WR) - t·Wp using Arb high precision.
    """
    ctx.prec = prec

    from arb_ldlt_certify import arb_closed_forms, arb_J, arb_kappa, prime_powers_up_to

    S, CC, XC, L = arb_closed_forms(N, c, prec)
    PI = arb.pi()
    sp2 = 16 * PI * PI
    l2 = L * L
    pref02 = 32 * L * (L / 4).sinh() ** 2
    kappa = arb_kappa(L)
    J = arb_J(L)

    pdata = prime_powers_up_to(c)
    weights = [arb(p).log() * (arb(q) ** arb('-0.5')) for (q, p) in pdata]
    positions = [arb(q).log() for (q, p) in pdata]
    M = len(weights)

    def S_signed(nn):
        return S[nn] if nn >= 0 else -S[-nn]

    DIM = 2 * N + 1
    from flint import arb_mat
    A = arb_mat(DIM, DIM)

    t_weight_arb = arb(t_weight)

    for i in range(DIM):
        n = i - N
        for j in range(i, DIM):
            m = j - N

            num = l2 - sp2 * m * n
            den = (l2 + sp2 * m * m) * (l2 + sp2 * n * n)
            W02 = pref02 * num / den

            if n == m:
                WR = kappa + 2 * CC[abs(n)] + J - (2 / L) * XC[abs(n)]
            else:
                WR = (S_signed(m) - S_signed(n)) / (PI * (n - m))

            Wp = arb(0)
            for idx in range(M):
                y = positions[idx]
                if n == m:
                    q = 2 * (1 - y / L) * (2 * PI * n * y / L).cos()
                else:
                    q = ((2 * PI * m * y / L).sin() - (2 * PI * n * y / L).sin()) / (PI * (n - m))
                Wp += weights[idx] * q

            val = W02 - WR - t_weight_arb * Wp
            A[i, j] = val
            A[j, i] = val

    return A, DIM


def count_negative_eigenvalues_ldlt(A, DIM):
    """
    Count the number of negative eigenvalues using LDL^T factorization (Sylvester's criterion).
    Returns (n_pos, n_neg) or (None, None) if the factorization encounters an undetermined pivot.
    """
    d = [None] * DIM
    Lf = [[arb(0)] * DIM for _ in range(DIM)]
    n_pos = 0
    n_neg = 0

    for i in range(DIM):
        s = A[i, i]
        for k in range(i):
            s = s - Lf[i][k] * Lf[i][k] * d[k]
        d[i] = s

        if s > 0:
            n_pos += 1
        elif s < 0:
            n_neg += 1
        else:
            # Undetermined pivot (zero within the ball)
            return None, None

        for j in range(i + 1, DIM):
            t = A[j, i]
            for k in range(i):
                t = t - Lf[j][k] * Lf[i][k] * d[k]
            Lf[j][i] = t / d[i]

    return n_pos, n_neg


def track_inertia(c, N, num_steps=21, prec=500):
    """
    Track inertia (number of negative eigenvalues) as t goes from 0 to 1.

    Returns:
      dict with t_values, n_pos_array, n_neg_array
    """
    print(f"Tracking inertia for c={c}, N={N}, dim={2*N+1}")
    print(f"  precision: {prec} bits")
    print(f"  t ∈ [0,1], {num_steps} steps\n")

    t_values = np.linspace(0, 1, num_steps)
    n_pos_array = []
    n_neg_array = []
    crossover_t = None

    for step_idx, t in enumerate(t_values):
        A, DIM = build_arb_tau_homotopy(c, N, float(t), prec)
        n_pos, n_neg = count_negative_eigenvalues_ldlt(A, DIM)

        if n_pos is None:
            print(f"  step {step_idx+1:2d}/{num_steps}: t={t:.3f} — LDL^T factorization undetermined (zero pivot)")
            n_pos_array.append(None)
            n_neg_array.append(None)
        else:
            print(f"  step {step_idx+1:2d}/{num_steps}: t={t:.3f}  n_pos={n_pos:2d}  n_neg={n_neg:2d}")
            n_pos_array.append(n_pos)
            n_neg_array.append(n_neg)

            # Track first sign change (when n_neg reaches 0)
            if crossover_t is None and n_neg == 0:
                crossover_t = t

    return {
        't_values': t_values,
        'n_pos': n_pos_array,
        'n_neg': n_neg_array,
        'crossover_t': crossover_t,
        'metadata': {'c': c, 'N': N, 'num_steps': num_steps, 'prec': prec},
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--c', type=int, default=30, help='cutoff parameter')
    p.add_argument('--N', type=int, default=8, help='number of frequency modes')
    p.add_argument('--steps', type=int, default=21, help='number of t values to sample')
    p.add_argument('--prec', type=int, default=500, help='Arb precision in bits')
    args = p.parse_args()

    result = track_inertia(args.c, args.N, num_steps=args.steps, prec=args.prec)

    print("\n" + "="*100)
    print("INERTIA TRACKING SUMMARY")
    print("="*100)

    # Find where n_neg becomes 0
    if result['crossover_t'] is not None:
        print(f"\nFirst zero-negativity crossing: t ≈ {result['crossover_t']:.4f}")
    else:
        # Find where n_neg first decreases (even if doesn't reach 0)
        n_neg_array = result['n_neg']
        n_neg_clean = [x for x in n_neg_array if x is not None]
        if len(n_neg_clean) > 1:
            min_n_neg = min(n_neg_clean)
            print(f"\nMinimum n_neg reached: {min_n_neg} (should be 0 at t=1)")

    # Interpretation
    print(f"\nInterpretation:")
    print(f"  At t=0: {result['n_neg'][0]} negative eigenvalues (background matrix W02-WR)")
    print(f"  At t=1: {result['n_neg'][-1]} negative eigenvalues (full matrix W02-WR-Wp, should be 0 if PSD)")
    print(f"\nKey finding:")
    print(f"  The structure of how these negative eigenvalues 'disappear' as t→1")
    print(f"  tells us how the primes Wp accomplish the cancellation.")

    # Save data
    output_json = f'inertia_tracking_c{args.c}_N{args.N}.json'
    import json
    with open(output_json, 'w') as f:
        json.dump({
            'c': args.c,
            'N': args.N,
            't_values': result['t_values'].tolist(),
            'n_pos': result['n_pos'],
            'n_neg': result['n_neg'],
            'crossover_t': result['crossover_t'],
        }, f, indent=2)
    print(f"\nSaved to {output_json}")


if __name__ == '__main__':
    main()
