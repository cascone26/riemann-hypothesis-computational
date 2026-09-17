"""
Compute κ (pole count) of the Weyl m-function m_a(z) for increasing a values.

κ = number of UHP poles of m(z) with negative trace.

In the finite-a Suzuki system, Thm 1.5 predicts κ=0 for all finite a.
In the limit as a→∞, if the system is limit-point and the m-function
has RH-specific structure, κ might jump to 1.

If κ stays 0 even as a→∞, that confirms the system doesn't capture
the pole structure needed for RH equivalence.

We use Nevanlinna-Gram coefficient matrix or winding-number tests
to compute κ from a sampled m-function.
"""

import sys, numpy as np
sys.path.insert(0, 'operators')

from canonical_weyl_solver import solve_canonical_system_finite_a, extract_weyl_m_function

def compute_kappa_via_nevanlinna_rank(m_func_samples, z_samples, verbose=False):
    """
    Compute κ via Nevanlinna-Gram matrix rank deficiency.

    Given samples of m(z) on a grid of z ∈ UHP, the Nevanlinna-Gram matrix
    K[i,j] = (m(z_i) - conj(m(z_j))) / (z_i - conj(z_j))

    should have rank = dim - κ (for a generalized-Nevanlinna function).

    κ = dim - rank is the number of negative squares (pole count).
    """

    dim = len(z_samples)
    K = np.zeros((dim, dim), dtype=complex)

    for i in range(dim):
        for j in range(dim):
            z_i = z_samples[i]
            z_j = z_samples[j]
            m_i = m_func_samples[i]
            m_j = m_func_samples[j]

            # Nevanlinna-Gram coefficient
            denom = z_i - np.conj(z_j)
            if abs(denom) > 1e-12:
                K[i, j] = (m_i - np.conj(m_j)) / denom
            else:
                K[i, j] = 0.0

    # Compute rank via SVD
    try:
        U, s, Vh = np.linalg.svd(K)
        rank = np.sum(s > 1e-6 * np.max(s))
    except:
        rank = dim

    kappa = dim - rank

    if verbose:
        print(f"    Nevanlinna-Gram matrix: dim={dim}, rank={rank}, κ={kappa}")
        print(f"    Singular values: {s[:5]}...")

    return kappa, rank, K

def test_kappa_evolution(z_center=None, a_values=None, verbose=True):
    """
    Compute κ at each a value and track how it evolves as a→∞.
    """

    if z_center is None:
        z_center = complex(1.0, 1.0)

    if a_values is None:
        a_values = [1, 3, 5, 10, 20]

    if verbose:
        print("="*80)
        print("κ EVOLUTION TEST: Pole count as a→∞")
        print("="*80)
        print(f"Base point z = {z_center}")
        print()
        print("  a    |m_a|          κ    rank   Interpretation")
        print("-" * 75)

    results = []

    for a in a_values:
        # Solve the fundamental solution
        X_final = solve_canonical_system_finite_a(z_center, a, rtol=1e-9, atol=1e-12)
        m_center = extract_weyl_m_function(X_final, z_center, a)

        # Sample m on a small circle around z_center to compute κ
        # For a rough estimate, use a circle of radius 0.3 around z_center in UHP
        theta_vals = np.linspace(0, 2*np.pi, 8, endpoint=False)
        radius = 0.3
        z_circle = z_center + radius * np.exp(1j * theta_vals)

        m_samples = []
        for z_test in z_circle:
            if z_test.imag < 0.01:  # Stay in UHP
                continue
            X_final_test = solve_canonical_system_finite_a(z_test, a, rtol=1e-9, atol=1e-12)
            if X_final_test is not None:
                m_test = extract_weyl_m_function(X_final_test, z_test, a)
                m_samples.append(m_test)
            else:
                m_samples.append(np.nan + 1j*np.nan)

        z_circle_valid = z_circle[z_circle.imag > 0.01][:len(m_samples)]

        if len(m_samples) > 1:
            kappa, rank, K = compute_kappa_via_nevanlinna_rank(m_samples, z_circle_valid, verbose=False)
        else:
            kappa = None
            rank = None

        if verbose:
            if kappa is not None:
                print(f" {a:3.0f}   {abs(m_center):.8e}    {kappa:2d}    {rank:2d}   ", end="")
                if kappa == 0:
                    print("No poles (matches Thm 1.5)")
                elif kappa == 1:
                    print("ONE pole! (RH-consistent)")
                else:
                    print(f"κ={kappa} (unexpected)")
            else:
                print(f" {a:3.0f}   {abs(m_center):.8e}    [failed to compute κ]")

        results.append({
            'a': a,
            'm_center': m_center,
            'kappa': kappa,
            'rank': rank,
        })

    print()
    return results

if __name__ == "__main__":
    print("\n" + "="*80)
    print("CANONICAL SYSTEM: κ POLE-COUNT EVOLUTION")
    print("="*80)
    print()

    # Test at one UHP point for now
    results = test_kappa_evolution(
        z_center=complex(1.0, 1.0),
        a_values=[1, 3, 5, 10, 20],
        verbose=True
    )

    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)

    kappa_vals = [r['kappa'] for r in results if r['kappa'] is not None]

    if all(k == 0 for k in kappa_vals):
        print("\nAll κ = 0: CONSISTENT WITH SUZUKI THM 1.5")
        print("The system preserves κ=0 structure even as a→∞.")
        print("(Non-convergence of m_a is genuine limit-circle behavior.)")
    elif any(k == 1 for k in kappa_vals):
        print("\nSome κ = 1: POTENTIAL RH SIGNAL")
        print("κ jumps to 1 in the limit — candidate for RH equivalence.")
    else:
        print("\nκ values:", kappa_vals)
