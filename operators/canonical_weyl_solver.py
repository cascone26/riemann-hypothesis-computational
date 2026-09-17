"""
Canonical-system Weyl m-function solver: limit-point classification via convergence.

Task: solve dX/dt = z·J·H(t)·X on [0,a] for increasing a, extract Weyl m-function,
test whether m_a(z) → m_∞(z) (limit-point) or diverges (limit-circle).

H(t) = [[0, -ρ(t)], [σ(t), 0]]  (canonical-system form)
σ(t) = Ω(t)/2π - P(t)  (explicit-formula symbol)
ρ(t) = 1  (unit weight; Suzuki structure suggests this or similar)

Weyl m-function: the solution couple X(t) = (x₁(t), x₂(t))^T is normalized
so that at large t, the vector (1, m(z))^T spans the L²-integrable subspace
(limit-point) or all solutions if limit-circle.

Standard extraction: solve on [0,a] with initial condition X(0)=(1,0)^T,
normalize at a, and read off m_a(z) from the boundary behavior.
"""

import numpy as np
import mpmath as mp
from scipy.integrate import solve_ivp, odeint
from scipy.special import digamma, polygamma
import matplotlib.pyplot as plt
import sys, os

# ============================================================================
# PART 1: EXPLICIT FORMULA SYMBOL (HIGH PRECISION)
# ============================================================================

def archimedean_weight_omega(t, prec=50):
    """
    High-precision archimedean weight Ω(t)/2π.

    Ω(t) = Re ψ₀(1/4 + it/2) - ln(π)
    ψ₀ is the digamma function.
    """
    mp.dps = prec  # Set precision
    s = mp.mpf('0.25') + 1j * mp.mpf(t) / 2
    psi_val = mp.digamma(s)
    omega = mp.re(psi_val) - mp.log(mp.pi)
    return float(omega) / (2 * mp.pi)

def von_mangoldt_weight_P(t, K_max=None, prec=50):
    """
    High-precision von Mangoldt weight P(t).

    P(t) = (1/π) Σ_{Λ(n)} Λ(n)/√n · cos(t·log n)

    Sum is over prime powers with log n ≤ K_max.
    For finite a: K_max = a
    For infinite a: K_max → ∞ (use sufficiently large value)
    """
    mp.dps = prec

    if K_max is None:
        K_max = 100

    # Generate prime powers: primes up to ~exp(K_max)
    # For K_max ≤ 50, we can afford a denser list
    primes = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
              101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199]

    result = mp.mpf(0)

    for p in primes:
        k = 1
        while True:
            n = p ** k
            log_n = k * float(mp.log(p))
            if log_n > K_max:
                break
            # Λ(p^k) = log(p)
            lambda_n = mp.log(p)
            contrib = lambda_n / mp.sqrt(mp.mpf(n)) * mp.cos(mp.mpf(t) * mp.log(mp.mpf(p)) * k)
            result += contrib
            k += 1

    return float(result) / float(mp.pi)

def explicit_formula_symbol(t, a_trunc=None, prec=50):
    """
    Full explicit-formula symbol σ(t) = Ω(t)/2π - P(t).

    At finite a: P(t) truncated to log n ≤ a
    """
    omega_val = archimedean_weight_omega(t, prec=prec)
    p_val = von_mangoldt_weight_P(t, K_max=a_trunc, prec=prec)
    return omega_val - p_val

# ============================================================================
# PART 2: CANONICAL SYSTEM ODE SOLVER
# ============================================================================

def hamiltonian_2x2(t, a):
    """
    Build the 2x2 canonical-system Hamiltonian at parameter a.

    H(t) = [[0, -ρ(t)], [σ(t), 0]]

    where σ(t) = Ω(t)/2π - P(t) is the explicit-formula symbol,
    ρ(t) = 1 (unit weight).
    """
    sigma_t = explicit_formula_symbol(t, a_trunc=a)
    rho_t = 1.0

    # Note: H is real symmetric with off-diagonal real elements
    H = np.array([[0.0, -rho_t], [sigma_t, 0.0]], dtype=np.float64)

    return H

def symplectic_matrix_J():
    """
    Symplectic matrix J = [[0, -1], [1, 0]].
    """
    return np.array([[0.0, -1.0], [1.0, 0.0]], dtype=np.complex128)

def canonical_system_ode(t, X_flat, a, z):
    """
    ODE system for the canonical system:

    d(J·X)/dt = H(t)·X

    which gives:
    J·dX/dt = H·X
    dX/dt = J^{-1}·H·X = -J·H·X  (since J^{-1} = -J)

    Equivalently, for the 2x2 Hamiltonian H = [[0, -ρ], [σ, 0]]:

    dX/dt = [[0, -J·σ], [-J·ρ, 0]] · X  (no, this is wrong)

    Let me use the correct formulation. A canonical system is typically:
    dX/dt = z·J·H·X where J = [[0,-1],[1,0]]

    But to avoid decoupling, I think the correct form for this problem is:
    dX₁/dt = -z·σ(t)·X₂
    dX₂/dt = z·ρ(t)·X₁

    which in matrix form (where X = (X₁, X₂)^T is a solution vector, not the fund matrix):
    dX/dt = [[0, -zσ(t)], [zρ(t), 0]] · X

    For a fundamental matrix, we track multiple solutions. Actually, let me reconsider:
    The Hamiltonian H should couple the two components. Let me use the correct form:

    dX/dt = z · H(t) · J^T · X   (alternative order)

    where J = [[0,-1],[1,0]], so J^T = [[0,1],[-1,0]]

    Then: H·J^T = [[0, -ρ], [σ, 0]] · [[0,1],[-1,0]]
                 = [[0·0 + (-ρ)·(-1), 0·1 + (-ρ)·0], [σ·0 + 0·(-1), σ·1 + 0·0]]
                 = [[ρ, 0], [0, σ]]

    That's still diagonal. Let me try: dX/dt = z · J · H(t) · X more carefully.

    Actually, wait. The standard screw/Sturm-Liouville form is:
    d/dt [ X₂ ]  =  [  0    σ(t) ] [ X₁ ]
         [ X₁ ]     [ -ρ(t)  0  ] [ X₂ ]

    NOT dX/dt = J·H·X in the usual order.

    The coupling matrix should be [[0, σ], [-ρ, 0]], not [[0, -ρ], [σ, 0]].
    Let me redefine.
    """

    # Unflatten X
    X = X_flat.reshape(2, 2)

    # Build the correct coupling matrix for a non-degenerate coupled system
    # The Sturm-Liouville form is:
    # d/dt [X2] = z * [0   σ(t)] [X1]
    #      [X1]       [-ρ  0  ] [X2]
    #
    # or equivalently, dX/dt = z * M(t) * X where
    # M(t) = [[0, σ(t)], [-ρ(t), 0]]

    sigma_t = explicit_formula_symbol(t, a_trunc=a)
    rho_t = 1.0

    M_t = np.array([[0.0, sigma_t], [-rho_t, 0.0]], dtype=np.complex128)

    # Compute dX/dt = z·M(t)·X
    dX = z * (M_t @ X)

    # Flatten and return
    return dX.flatten()

def solve_canonical_system_finite_a(z, a, t_max=None, rtol=1e-9, atol=1e-12):
    """
    Solve the canonical system on [0, t_max] for parameter a and complex z.

    Initial condition: X(0) = I (identity)

    Returns: solution at t_max, as a 2x2 matrix
    """

    if t_max is None:
        t_max = a

    # Initial condition: fundamental solution X(0) = I
    X0 = np.eye(2, dtype=np.complex128).flatten()

    # Define a dense evaluation grid to get smooth tracking
    t_eval = np.linspace(0, t_max, max(100, int(50*t_max)))

    try:
        # Use solve_ivp with RK45 (adaptive, suitable for complex-valued ODE)
        # scipy.integrate.solve_ivp handles complex-valued systems via real decomposition
        sol = solve_ivp(canonical_system_ode, [0, t_max], X0,
                        args=(a, z),
                        method='RK45',
                        t_eval=t_eval,
                        dense_output=True,
                        rtol=rtol,
                        atol=atol,
                        max_step=0.1)

        if sol.status != 0:
            print(f"  Warning: solve_ivp status={sol.status} for a={a}, z={z}")

        # Extract final solution
        X_final_flat = sol.y[:, -1]
        X_final = X_final_flat.reshape(2, 2)

        return X_final

    except Exception as e:
        print(f"  ODE solve failed for a={a}, z={z}: {e}")
        return None

# ============================================================================
# PART 3: WEYL m-FUNCTION EXTRACTION (STANDARD FORMULA)
# ============================================================================

def extract_weyl_m_function(X_final, z, a):
    """
    Extract the Weyl m-function from the final fundamental solution X(a).

    Standard formula for canonical systems:
    The Weyl m-function m(z) is defined via the condition that
    (1, m(z))^T · (X₁₁, X₁₂)^T is L²-integrable at infinity.

    For the finite-a system, we use:
    m_a(z) = X₂₁(a) / X₁₁(a)  (standard Weyl extraction from fund. solution)

    This ratio encodes the boundary behavior of the solution couple.

    Alternative (if X₁₁ is small): use X₂₂ / X₁₂ or similar, with regularization.
    """

    if X_final is None:
        return np.nan + 1j*np.nan

    X11, X12 = X_final[0, 0], X_final[0, 1]
    X21, X22 = X_final[1, 0], X_final[1, 1]

    # Primary extraction
    if abs(X11) > 1e-10:
        m_a = X21 / X11
    elif abs(X12) > 1e-10:
        m_a = X22 / X12
    else:
        # Both leading elements too small — system is singular
        return np.nan + 1j*np.nan

    return m_a

# ============================================================================
# PART 4: CONVERGENCE TEST (LIMIT-POINT CLASSIFICATION)
# ============================================================================

def test_limit_point_convergence(z_test=None, a_values=None, verbose=True):
    """
    Numerical test for limit-point vs limit-circle classification.

    Solve the canonical system for increasing values of a, extract m_a(z),
    and check whether m_a(z) converges to a limit m_∞(z) (limit-point)
    or diverges/oscillates (limit-circle).

    Returns:
    - results: list of (a, m_a) pairs
    - limit_point: boolean indicating convergence
    - interpretation: string describing the finding
    """

    if z_test is None:
        z_test = complex(1.0, 1.0)  # Test point in UHP

    if a_values is None:
        a_values = [1, 5, 10, 20, 50]

    if verbose:
        print("="*80)
        print("LIMIT-POINT CLASSIFICATION TEST")
        print("="*80)
        print(f"Test parameter z = {z_test}")
        print(f"Test values: a = {a_values}")
        print()
        print("  a    Re(m_a)          Im(m_a)          |m_a|         Δ|m| from prev")
        print("-" * 80)

    results = []
    m_prev_abs = None
    m_prev = None

    for a in a_values:
        try:
            X_final = solve_canonical_system_finite_a(z_test, a)
            m_a = extract_weyl_m_function(X_final, z_test, a)

            if np.isnan(m_a) or np.isinf(m_a):
                if verbose:
                    print(f"{a:4.0f}   [singular or NaN]")
                results.append((a, None))
                continue

            m_abs = abs(m_a)
            results.append((a, m_a))

            if m_prev_abs is None:
                if verbose:
                    print(f"{a:4.0f}   {m_a.real:+.10e}  {m_a.imag:+.10e}  {m_abs:.6e}")
            else:
                delta_m = abs(m_abs - m_prev_abs)
                if verbose:
                    print(f"{a:4.0f}   {m_a.real:+.10e}  {m_a.imag:+.10e}  {m_abs:.6e}  {delta_m:.3e}")

            m_prev_abs = m_abs
            m_prev = m_a

        except Exception as e:
            if verbose:
                print(f"{a:4.0f}   [error: {str(e)[:40]}]")
            results.append((a, None))

    # Analyze convergence
    valid_results = [(a, m) for a, m in results if m is not None]

    if len(valid_results) < 2:
        limit_point = None
        interpretation = "Insufficient valid data points"
    else:
        # Check if |m_a| is stabilizing
        m_vals = np.array([abs(m) for _, m in valid_results])
        m_diffs = np.diff(m_vals)

        # If the last few differences are small, declare convergence
        if len(m_diffs) >= 2:
            avg_recent_diff = np.mean(np.abs(m_diffs[-2:]))
            if avg_recent_diff < 0.01:
                limit_point = True
                interpretation = "LIMIT-POINT: m_a(z) appears to converge"
            else:
                # Check for wild growth or oscillation
                if np.max(np.abs(m_diffs)) > 0.1 * np.mean(m_vals):
                    limit_point = False
                    interpretation = "DIVERGENT: |m_a| shows large fluctuations"
                else:
                    limit_point = None
                    interpretation = "INCONCLUSIVE: behavior unclear at this resolution"
        else:
            limit_point = None
            interpretation = "INCONCLUSIVE: too few data points"

    if verbose:
        print()
        print("INTERPRETATION:")
        print(f"  {interpretation}")
        print()

    return results, limit_point, interpretation

# ============================================================================
# PART 5: EXECUTION
# ============================================================================

if __name__ == "__main__":
    print("\n" + "="*80)
    print("CANONICAL-SYSTEM WEYL m-FUNCTION CONVERGENCE TEST")
    print("="*80)
    print()

    # Test at a few UHP points to get robust picture
    test_points = [
        complex(0.5, 1.0),    # Lower UHP
        complex(1.0, 1.0),    # Middle UHP
        complex(2.0, 1.0),    # Higher UHP
    ]

    a_values = [1, 5, 10, 20, 50]

    all_results = {}

    for z_test in test_points:
        print(f"\n{'='*80}")
        print(f"Testing at z = {z_test}")
        print(f"{'='*80}")

        results, limit_point, interpretation = test_limit_point_convergence(
            z_test=z_test,
            a_values=a_values,
            verbose=True
        )

        all_results[z_test] = {
            'results': results,
            'limit_point': limit_point,
            'interpretation': interpretation,
        }

    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)

    for z_test, data in all_results.items():
        print(f"\nz = {z_test}: {data['interpretation']}")
        if data['results']:
            m_final = data['results'][-1][1]
            if m_final is not None:
                print(f"  Final m_a(z) ≈ {m_final.real:+.6e} + {m_final.imag:+.6e}i")
