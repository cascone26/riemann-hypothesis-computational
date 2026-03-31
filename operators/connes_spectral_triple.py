"""
Implementation of Connes-Consani-Moscovici Zeta Spectral Triples (2025)
arXiv:2511.22755

Constructs D_log^(lambda,N) operators as rank-one perturbations of the
scaling operator, whose eigenvalues numerically match zeta zeros.

Parameters:
  lambda = sqrt(14), so L = log(14) ≈ 2.6391
  N = 30 (smaller than paper's 120 due to precision constraints)
  Primes: {2, 3, 5, 7, 11, 13} (all p <= lambda^2 = 14)

Uses mpmath for high-precision arithmetic (50 digits — lower than
paper's 200 but should show the pattern).
"""

import mpmath
import numpy as np
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
PLOTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "plots")
os.makedirs(RESULTS_DIR, exist_ok=True)

# Set precision
mpmath.mp.dps = 50

# Parameters
LAMBDA = mpmath.sqrt(14)
L = 2 * mpmath.log(LAMBDA)  # = log(14)
PRIMES = [2, 3, 5, 7, 11, 13]  # primes <= lambda^2 = 14
N_TRUNC = 30  # truncation (paper uses 120, we use 30 for speed)


def load_zeros(count=50):
    path = os.path.join(DATA_DIR, "zeros_100k.txt")
    zeros = []
    with open(path, "r") as f:
        for i, line in enumerate(f):
            if i >= count:
                break
            line = line.strip()
            if line:
                zeros.append(mpmath.mpf(line))
    return zeros


def von_mangoldt(k):
    """Von Mangoldt function: Lambda(k) = log(p) if k = p^m, else 0."""
    if k <= 1:
        return mpmath.mpf(0)
    for p in PRIMES:
        if p > k:
            break
        pk = p
        while pk <= k:
            if pk == k:
                return mpmath.log(p)
            pk *= p
    return mpmath.mpf(0)


def W_boundary(n, m):
    """
    Boundary evaluation term W_{0,2}(V_n, V_m).
    Formula: 32*L * sinh^2(L/4) * (L^2 - 16*pi^2*m*n) /
             [(L^2 + 16*pi^2*m^2) * (L^2 + 16*pi^2*n^2)]
    """
    pi = mpmath.pi
    num = 32 * L * mpmath.sinh(L/4)**2 * (L**2 - 16 * pi**2 * m * n)
    den = (L**2 + 16 * pi**2 * m**2) * (L**2 + 16 * pi**2 * n**2)
    if den == 0:
        return mpmath.mpf(0)
    return num / den


def W_prime_contribution(n, m):
    """
    Non-archimedean (prime) contribution:
    sum_{1 < k <= exp(L)} Lambda(k) * k^{-1/2} * q(U_n, U_m)(log k)

    For n != m: q(U_n, U_m)(y) = [sin(2*pi*m*y/L) - sin(2*pi*n*y/L)] / [pi*(n-m)]
    For n = m: q(U_n, U_n)(y) = 2*(1 - |y|/L) * cos(2*pi*n*y/L)
    """
    pi = mpmath.pi
    total = mpmath.mpf(0)
    k_max = int(float(mpmath.exp(L))) + 1

    for k in range(2, min(k_max + 1, 200)):  # cap at 200 for speed
        lam_k = von_mangoldt(k)
        if lam_k == 0:
            continue
        y = mpmath.log(k)
        weight = lam_k * k ** mpmath.mpf('-0.5')

        if n != m:
            q_val = (mpmath.sin(2*pi*m*y/L) - mpmath.sin(2*pi*n*y/L)) / (pi * (n - m))
        else:
            q_val = 2 * (1 - abs(y)/L) * mpmath.cos(2*pi*n*y/L)

        total += weight * q_val

    return total


def W_archimedean_approx(n, m):
    """
    Archimedean contribution W_R(V_n, V_m).
    This involves hypergeometric functions and digamma functions.
    For a first implementation, we use a simplified version.

    The key structure: W_R has both diagonal and off-diagonal parts.
    The diagonal part dominates and is related to log(|s_n|) where
    s_n = 2*pi*n/L are the unperturbed eigenvalues.

    Approximation: W_R(V_n, V_n) ≈ (1/2)*psi(1/4 + i*pi*n/L) + (1/2)*psi(1/4 - i*pi*n/L) + C
    where psi is the digamma function.
    """
    pi = mpmath.pi
    if n == 0 and m == 0:
        # Special case: involves Euler-Mascheroni constant
        return -mpmath.euler - mpmath.log(2)

    if n == m:
        s = pi * n / L
        if abs(s) < 0.01:
            return -mpmath.euler - mpmath.log(2)
        # Digamma approximation for the archimedean contribution
        val = mpmath.re(mpmath.digamma(mpmath.mpf('0.25') + 1j * s))
        return val + mpmath.log(2)

    # Off-diagonal: smaller, related to the overlap of shifted digamma functions
    # For simplicity, use a leading-order approximation
    if n - m != 0:
        return mpmath.mpf(0)  # Leading order: diagonal dominant

    return mpmath.mpf(0)


def build_weil_matrix():
    """Build the (2N+1) x (2N+1) Weil quadratic form matrix."""
    dim = 2 * N_TRUNC + 1
    print(f"Building {dim}x{dim} Weil matrix...")
    print(f"  lambda = sqrt(14), L = {float(L):.6f}")
    print(f"  Primes: {PRIMES}")

    tau = mpmath.matrix(dim, dim)

    for i in range(dim):
        n = i - N_TRUNC  # n ranges from -N to N
        for j in range(i, dim):
            m = j - N_TRUNC

            # Weil form = archimedean + boundary - primes
            w_arch = W_archimedean_approx(n, m)
            w_bound = W_boundary(n, m)
            w_prime = W_prime_contribution(n, m)

            val = w_arch + w_bound - w_prime
            tau[i, j] = val
            tau[j, i] = val  # symmetric

        if (i + 1) % 10 == 0:
            print(f"  Row {i+1}/{dim} done")

    return tau


def find_min_eigenvector(tau):
    """Find the minimum eigenvalue and eigenvector of tau."""
    dim = tau.rows
    print(f"Finding minimum eigenvector of {dim}x{dim} matrix...")

    # Convert to numpy for eigenvalue computation
    tau_np = np.array([[float(tau[i,j]) for j in range(dim)] for i in range(dim)])

    eigenvalues, eigenvectors = np.linalg.eigh(tau_np)

    # Minimum eigenvalue
    min_idx = 0
    epsilon_N = eigenvalues[min_idx]
    xi = eigenvectors[:, min_idx]

    print(f"  Minimum eigenvalue: {epsilon_N:.10f}")
    print(f"  First few eigenvector components: {xi[:5]}")

    # Normalize: delta_N(xi) = sum xi_j / sqrt(L) = 1
    delta_sum = np.sum(xi) / float(mpmath.sqrt(L))
    if abs(delta_sum) > 1e-10:
        xi = xi / delta_sum

    return epsilon_N, xi


def xi_hat(z, xi_vec):
    """
    Fourier transform of xi:
    xi_hat(z) = (2/sqrt(L)) * sin(z*L/2) * sum_{j=-N}^{N} xi_j / (z - 2*pi*j/L)
    """
    L_f = float(L)
    prefix = (2 / np.sqrt(L_f)) * np.sin(z * L_f / 2)
    total = 0.0
    for j_idx in range(len(xi_vec)):
        j = j_idx - N_TRUNC
        pole = 2 * np.pi * j / L_f
        if abs(z - pole) < 1e-10:
            # Near pole — use L'Hopital or skip
            continue
        total += xi_vec[j_idx] / (z - pole)
    return prefix * total


def find_spectral_zeros(xi_vec, zeros_target, search_range=100):
    """Find zeros of xi_hat(z) near known zeta zeros."""
    from scipy.optimize import brentq

    z_test = np.linspace(0.1, search_range, 10000)
    xi_vals = [xi_hat(z, xi_vec) for z in z_test]

    # Find sign changes
    found_zeros = []
    for i in range(len(z_test) - 1):
        if xi_vals[i] * xi_vals[i+1] < 0:
            try:
                z_root = brentq(lambda z: xi_hat(z, xi_vec), z_test[i], z_test[i+1])
                found_zeros.append(z_root)
            except:
                pass

    return found_zeros


def main():
    print("=" * 60)
    print("CONNES-CONSANI-MOSCOVICI ZETA SPECTRAL TRIPLE")
    print("arXiv:2511.22755 (2025)")
    print("=" * 60)

    zeros = load_zeros(30)

    # Step 1: Build Weil matrix
    tau = build_weil_matrix()

    # Step 2: Find minimum eigenvector
    epsilon_N, xi = find_min_eigenvector(tau)

    # Step 3: Find zeros of xi_hat
    print("\nFinding zeros of xi_hat(z)...")
    found = find_spectral_zeros(xi, zeros, search_range=80)
    print(f"Found {len(found)} zeros")

    # Step 4: Compare to zeta zeros
    print(f"\n{'#':>3} {'Found':>12} {'Zeta Zero':>12} {'Error':>12}")
    print("-" * 42)
    matched = 0
    comparisons = []
    for i, fz in enumerate(found[:20]):
        # Find nearest zeta zero
        diffs = [abs(fz - float(z)) for z in zeros]
        nearest_idx = np.argmin(diffs)
        nearest = float(zeros[nearest_idx])
        err = abs(fz - nearest)
        if err < 2.0:
            matched += 1
        comparisons.append({"found": fz, "zero": nearest, "error": err})
        print(f"{i+1:>3} {fz:>12.4f} {nearest:>12.4f} {err:>12.6f}")

    print(f"\nMatched (within 2.0): {matched}/{min(len(found), 20)}")
    print(f"\nNote: With N={N_TRUNC} and 50-digit precision, we expect")
    print(f"moderate accuracy. The paper uses N=120 and 200 digits.")

    # Plot
    z_range = np.linspace(0.1, 80, 5000)
    xi_vals = [xi_hat(z, xi) for z in z_range]

    fig, ax = plt.subplots(figsize=(14, 5))
    ax.plot(z_range, xi_vals, 'b-', linewidth=0.5)
    ax.axhline(y=0, color='black', linewidth=0.5)
    for z in zeros[:20]:
        ax.axvline(x=float(z), color='red', alpha=0.3, linewidth=0.5)
    ax.set_xlabel('z')
    ax.set_ylabel('xi_hat(z)')
    ax.set_title(f'Connes Spectral Triple: xi_hat(z) zeros vs zeta zeros (N={N_TRUNC}, lambda=sqrt(14))')
    ax.set_xlim(0, 80)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOTS_DIR, 'connes_spectral_triple.png'), dpi=150)
    plt.close()
    print(f"\nPlot saved.")

    # Save
    results = {
        "N": N_TRUNC,
        "lambda": float(LAMBDA),
        "L": float(L),
        "primes": PRIMES,
        "epsilon_N": float(epsilon_N),
        "found_zeros": [float(z) for z in found[:20]],
        "comparisons": comparisons[:20],
        "matched": matched,
    }
    with open(os.path.join(RESULTS_DIR, "connes_spectral_triple.json"), "w") as f:
        json.dump(results, f, indent=2)


if __name__ == "__main__":
    main()
