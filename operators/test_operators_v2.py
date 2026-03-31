"""
Operator Testing v2 — Better discretization approaches.

The v1 tests all failed because naive finite-difference discretization
of xp-type operators on bounded intervals doesn't work. The issue is
that these operators fundamentally need specific boundary conditions
and domain choices to produce discrete spectra related to zeta zeros.

This version tries:
1. Spectral methods instead of finite differences
2. The specific regularization schemes from the literature
3. Working in Fourier space where p is diagonal
4. The explicit construction from Sierra (2008) that reportedly
   reproduces low-lying zeros
"""

import numpy as np
from scipy import linalg
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")


def load_target_zeros(count=200):
    path = os.path.join(DATA_DIR, "zeros_100k.txt")
    zeros = []
    with open(path, "r") as f:
        for i, line in enumerate(f):
            if i >= count:
                break
            line = line.strip()
            if line:
                zeros.append(float(line))
    return np.array(zeros)


def compare(eigs, zeros, label=""):
    """Compare eigenvalues to zeros and print report."""
    eigs = np.sort(np.abs(eigs))
    eigs = eigs[eigs > 0.5]  # filter noise

    n = min(len(eigs), len(zeros), 20)
    if n == 0:
        print(f"  {label}: No valid eigenvalues")
        return None

    e = eigs[:n]
    z = zeros[:n]
    rel_err = np.abs(e - z) / z

    print(f"\n  {label}")
    print(f"  {'#':>3} {'Eigenvalue':>12} {'Zero':>12} {'Rel Err':>10}")
    print(f"  {'-'*40}")
    for i in range(min(n, 10)):
        print(f"  {i+1:>3} {e[i]:>12.4f} {z[i]:>12.4f} {rel_err[i]:>10.6f}")

    print(f"  Mean rel error (first {n}): {np.mean(rel_err):.6f}")
    return {"eigenvalues": e.tolist(), "mean_rel_error": float(np.mean(rel_err))}


# =============================================================================
# APPROACH 1: Berry-Keating via quantization condition
# =============================================================================

def berry_keating_quantization():
    """
    Berry & Keating's key insight: the classical Hamiltonian H = xp
    has trajectories xp = E (hyperbolas). Bohr-Sommerfeld quantization
    on a bounded phase space gives:

    The number of eigenvalues below E is approximately N(E),
    the zero counting function.

    This means the eigenvalues should BE the zeta zeros, by construction.
    But this only works if you choose the right phase space cutoff.

    Let's verify: does N(t_n) ≈ n for the known zeros?
    """
    print("=" * 60)
    print("APPROACH 1: Berry-Keating Quantization Condition")
    print("=" * 60)

    zeros = load_target_zeros(200)

    def N_smooth(t):
        """Smooth part of zero counting function."""
        return (t / (2 * np.pi)) * np.log(t / (2 * np.pi)) - t / (2 * np.pi) + 7 / 8

    print("\n  Verifying N(t_n) ≈ n for known zeros:")
    print(f"  {'n':>5} {'t_n':>12} {'N(t_n)':>12} {'Error':>10}")
    print(f"  {'-'*42}")

    for n in [1, 2, 5, 10, 20, 50, 100, 200]:
        if n <= len(zeros):
            t_n = zeros[n - 1]
            N_val = N_smooth(t_n)
            err = N_val - n
            print(f"  {n:>5} {t_n:>12.4f} {N_val:>12.4f} {err:>10.4f}")

    print("\n  The N(t_n) - n offset is the oscillatory part S(t).")
    print("  It averages to 0 but fluctuates. This IS the zero distribution.")


# =============================================================================
# APPROACH 2: Sierra's explicit Hamiltonian (2008)
# =============================================================================

def sierra_explicit():
    """
    Sierra (2008) proposed a specific quantum mechanical model where
    the phase shifts of scattering states match the zeta zeros.

    The model uses the Hamiltonian H = p^2/2m + V(x) where V(x) is
    constructed from the zero counting function.

    In the WKB approximation, the eigenvalues E_n satisfy:
    integral sqrt(E_n - V(x)) dx = (n + 1/2) * pi

    If V(x) is chosen so that this integral equals N(E_n),
    then E_n = t_n (the zeta zeros).

    The explicit construction: V(x) = inverse function of
    (1/pi) * integral_0^E sqrt(E - V) dx
    """
    print("\n" + "=" * 60)
    print("APPROACH 2: Sierra Explicit Construction")
    print("=" * 60)

    zeros = load_target_zeros(200)

    # The WKB potential that reproduces zeta zeros
    # We need: pi * N(E) = integral_0^L sqrt(E - V(x)) dx
    # For a linear potential V(x) = cx on [0, L]:
    # integral = (2/3c) * E^(3/2) when L >= E/c

    # Better: use the smooth N(E) to construct the potential numerically
    # N(E) ~ (E/2pi) * log(E/2pi) - E/2pi + 7/8

    # Eigenvalues from WKB: E_n where N(E_n) = n + 1/2
    # This is just inverting the counting function!

    from scipy.optimize import brentq

    def N_smooth(E):
        if E < 1:
            return 0
        return (E / (2 * np.pi)) * np.log(E / (2 * np.pi)) - E / (2 * np.pi) + 7 / 8

    print("\n  WKB eigenvalues (from inverting N(E) = n + 1/2):")
    wkb_eigs = []
    for n in range(200):
        target = n + 0.5
        try:
            E_n = brentq(lambda E: N_smooth(E) - target, 1, 1e6)
            wkb_eigs.append(E_n)
        except:
            break

    wkb_eigs = np.array(wkb_eigs)
    compare(wkb_eigs, zeros, "Sierra WKB (N(E_n) = n + 1/2)")

    # Also try N(E_n) = n (no half-integer shift)
    wkb_eigs2 = []
    for n in range(1, 201):
        target = n
        try:
            E_n = brentq(lambda E: N_smooth(E) - target, 1, 1e6)
            wkb_eigs2.append(E_n)
        except:
            break

    wkb_eigs2 = np.array(wkb_eigs2)
    compare(wkb_eigs2, zeros, "Sierra WKB (N(E_n) = n)")


# =============================================================================
# APPROACH 3: Direct matrix from zeta function
# =============================================================================

def montgomery_matrix():
    """
    Construct a matrix whose eigenvalues match zeta zeros using
    the explicit formula connecting primes and zeros.

    The GUE connection suggests: build a random matrix from the
    GUE ensemble and check if its eigenvalue statistics match.

    More specifically: use the fact that the pair correlation of
    zeros matches GUE to construct a matrix model.

    The simplest version: take the known zeros and construct the
    matrix that has them as eigenvalues (trivially), then see if
    this matrix has special structure.
    """
    print("\n" + "=" * 60)
    print("APPROACH 3: Reverse Engineering — Matrix from Zeros")
    print("=" * 60)

    zeros = load_target_zeros(100)
    N = len(zeros)

    # Trivially: diagonal matrix with zeros as eigenvalues
    # But let's look for structure by putting them through GUE format

    # GUE matrix: H = (A + A†)/2 where A has i.i.d. complex Gaussian entries
    # We want a matrix with eigenvalues = zeros[0..N-1]
    # Construct via H = U D U† where D = diag(zeros), U = random unitary

    # But the INTERESTING question: is there a STRUCTURED (not random) U
    # such that H has a recognizable form?

    # Let's analyze the structure of the required matrix in various bases
    D = np.diag(zeros)

    # Check: what does the matrix look like in the Fourier basis?
    F = np.fft.fft(np.eye(N)) / np.sqrt(N)
    H_fourier = F @ D @ F.conj().T

    # Check sparsity / decay of off-diagonal elements
    off_diag = np.abs(H_fourier) - np.diag(np.diag(np.abs(H_fourier)))
    max_off = np.max(off_diag)
    mean_off = np.mean(off_diag)

    print(f"\n  Matrix in Fourier basis:")
    print(f"    Max off-diagonal: {max_off:.4f}")
    print(f"    Mean off-diagonal: {mean_off:.4f}")
    print(f"    Diagonal range: {np.min(np.abs(np.diag(H_fourier))):.4f} to {np.max(np.abs(np.diag(H_fourier))):.4f}")

    # What about a basis from prime numbers?
    # Construct a unitary from the prime counting function
    primes = []
    for n in range(2, 600):
        if all(n % p != 0 for p in primes):
            primes.append(n)
        if len(primes) >= N:
            break

    # log(p_k) values
    log_primes = np.log(np.array(primes[:N], dtype=float))

    # Build a matrix from prime structure
    # M[j,k] = exp(2*pi*i * zeros[j] * log(primes[k]) / (2*pi))
    #         = primes[k]^(i * zeros[j])
    M = np.zeros((N, N), dtype=complex)
    for j in range(N):
        for k in range(N):
            M[j, k] = primes[k] ** (1j * zeros[j]) / np.sqrt(N)

    # Is M close to unitary?
    MMt = M @ M.conj().T
    unitarity_error = np.linalg.norm(MMt - np.eye(N)) / N
    print(f"\n  Prime-zero matrix p_k^(i*t_j):")
    print(f"    Unitarity error (||MM†-I||/N): {unitarity_error:.6f}")

    # The Selberg/Guinand explicit formula connects primes and zeros via:
    # sum_zeros f(t_n) = f-hat(0)*log(pi) - sum_primes sum_m log(p)*f-hat(m*log(p))/p^(m/2)
    # This is the mathematical structure that any operator must encode


# =============================================================================
# APPROACH 4: Spectral zeta function approach
# =============================================================================

def spectral_approach():
    """
    Instead of guessing operators, work backwards:

    Given eigenvalues {t_n}, what can we deduce about the operator?

    The spectral zeta function of operator H is:
    Z(s) = sum_n t_n^(-s)

    For the Riemann zeros, this is:
    Z(s) = sum_n (1/t_n^s)

    This is related to the Hadamard product of zeta.
    Compute it and look for structure.
    """
    print("\n" + "=" * 60)
    print("APPROACH 4: Spectral Zeta Function of the Zeros")
    print("=" * 60)

    zeros = load_target_zeros(10000)

    print("\n  Z(s) = sum_n t_n^(-s) for various s:")
    print(f"  {'s':>6} {'Z(s)':>16} {'converged':>10}")
    print(f"  {'-'*35}")

    for s in [2.0, 3.0, 4.0, 5.0, 1.5, 1.1]:
        Z = np.sum(zeros ** (-s))
        # Check convergence by comparing N=5000 vs N=10000
        Z_half = np.sum(zeros[:5000] ** (-s))
        rel_change = abs(Z - Z_half) / abs(Z)
        converged = "yes" if rel_change < 0.01 else "no"
        print(f"  {s:>6.1f} {Z:>16.10f} {converged:>10}")

    # The trace of H^(-s) relates to the heat kernel
    # Tr(e^{-tH}) = sum_n e^{-t*t_n}
    print("\n  Heat trace Tr(e^{-tH}) = sum_n e^{-t*t_n}:")
    for t_val in [0.001, 0.01, 0.05, 0.1, 0.5, 1.0]:
        heat = np.sum(np.exp(-t_val * zeros))
        print(f"    t={t_val:.3f}: {heat:.10f}")

    # The heat trace encodes ALL information about the operator's spectrum
    # Its short-time asymptotics tell us about the geometry of the underlying space
    # Tr(e^{-tH}) ~ a_0 * t^(-d/2) + a_1 * t^(-(d-1)/2) + ...
    # where d is the "dimension" of the space

    print("\n  Short-time heat kernel asymptotics:")
    ts = np.array([0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005])
    heats = np.array([np.sum(np.exp(-t * zeros[:50000])) for t in ts])

    # Fit log(heat) vs log(t) to get the exponent
    log_t = np.log(ts)
    log_h = np.log(heats)
    slope, intercept = np.polyfit(log_t, log_h, 1)
    print(f"    Fitted exponent: {slope:.4f}")
    print(f"    (d/2 = {-slope:.4f}, so effective dimension d = {-2*slope:.4f})")
    print(f"    For comparison: d=1 would give slope -0.5")


def main():
    berry_keating_quantization()
    sierra_explicit()
    montgomery_matrix()
    spectral_approach()

    print("\n" + "=" * 60)
    print("KEY TAKEAWAYS")
    print("=" * 60)
    print("""
  1. The WKB/quantization approach works: inverting N(E) = n
     gives eigenvalues that match zeros to the accuracy of the
     smooth counting function. The REAL problem is the oscillatory
     correction S(t).

  2. The naive discretization of xp-type operators fails completely.
     These operators need very specific domains and boundary conditions.

  3. The prime-zero matrix structure is non-trivial but not unitary —
     the explicit formula is the mathematical link.

  4. The heat kernel reveals an effective dimension of the underlying
     space. This constrains what kind of operator we're looking for.
    """)


if __name__ == "__main__":
    main()
