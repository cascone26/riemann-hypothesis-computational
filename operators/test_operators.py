"""
Hilbert-Polya Operator Testing Framework

Tests candidate self-adjoint operators against known Riemann zeta zeros.
For each candidate, computes eigenvalues numerically and compares to zeros.

The Hilbert-Polya conjecture: there exists an unbounded self-adjoint operator H
whose eigenvalues are the imaginary parts of the non-trivial zeta zeros.
If we find it, RH follows automatically (eigenvalues of self-adjoint operators are real).
"""

import numpy as np
from scipy import linalg, sparse
from scipy.sparse.linalg import eigsh
import math
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")
os.makedirs(RESULTS_DIR, exist_ok=True)


def load_target_zeros(count=1000):
    """Load zeta zeros as our target eigenvalues."""
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


# =============================================================================
# OPERATOR CANDIDATES
# =============================================================================

class OperatorCandidate:
    """Base class for Hilbert-Polya operator candidates."""

    name = "base"
    description = ""
    reference = ""

    def build_matrix(self, N, **params):
        """Build NxN matrix approximation of the operator."""
        raise NotImplementedError

    def compute_eigenvalues(self, N, n_eigs=100, **params):
        """Compute the first n_eigs eigenvalues."""
        H = self.build_matrix(N, **params)

        # Ensure Hermitian (self-adjoint)
        H = (H + H.conj().T) / 2

        if sparse.issparse(H):
            eigenvalues = eigsh(H, k=min(n_eigs, N - 2), which='SM',
                               return_eigenvectors=False)
        else:
            eigenvalues = np.linalg.eigvalsh(H)

        eigenvalues = np.sort(np.abs(eigenvalues))

        # Return only positive eigenvalues (the operator should have
        # eigenvalues at ±t_n due to symmetry)
        positive = eigenvalues[eigenvalues > 0.1]
        return positive[:n_eigs]

    def compare_to_zeros(self, zeros, N=2000, n_compare=50, **params):
        """Compare eigenvalues to known zeta zeros."""
        try:
            eigs = self.compute_eigenvalues(N, n_eigs=n_compare, **params)
        except Exception as e:
            return {"error": str(e), "name": self.name}

        n_match = min(len(eigs), len(zeros), n_compare)
        if n_match == 0:
            return {"error": "No eigenvalues computed", "name": self.name}

        eigs = eigs[:n_match]
        target = zeros[:n_match]

        abs_errors = np.abs(eigs - target)
        rel_errors = abs_errors / target

        # Find where it diverges
        diverge_idx = n_match
        for i in range(n_match):
            if rel_errors[i] > 0.5:  # 50% relative error = diverged
                diverge_idx = i
                break

        return {
            "name": self.name,
            "description": self.description,
            "reference": self.reference,
            "N_matrix": N,
            "n_eigenvalues": n_match,
            "eigenvalues": eigs.tolist(),
            "target_zeros": target.tolist(),
            "abs_errors": abs_errors.tolist(),
            "rel_errors": rel_errors.tolist(),
            "mean_abs_error": float(np.mean(abs_errors)),
            "mean_rel_error": float(np.mean(rel_errors)),
            "max_abs_error": float(np.max(abs_errors)),
            "diverges_at": diverge_idx,
            "first_10_comparison": [
                {"n": i + 1, "eigenvalue": float(eigs[i]),
                 "zero": float(target[i]),
                 "abs_err": float(abs_errors[i]),
                 "rel_err": float(rel_errors[i])}
                for i in range(min(10, n_match))
            ],
        }


class BerryKeatingXP(OperatorCandidate):
    """
    Berry-Keating Hamiltonian: H = xp (symmetrized)
    H = (1/2)(xp + px) = -i(x d/dx + 1/2)

    On L2(R+), this gives eigenvalues on a continuous spectrum.
    Needs regularization/boundary conditions to get discrete spectrum.

    Berry & Keating (1999): "H = xp and the Riemann Zeros"
    """

    name = "Berry-Keating xp"
    description = "H = (xp + px)/2 with hard-wall boundary conditions"
    reference = "Berry & Keating, 1999"

    def build_matrix(self, N, x_max=100.0, **params):
        dx = x_max / N
        # Grid points (avoid x=0)
        x = np.linspace(dx, x_max, N)

        # H = -i(x d/dx + 1/2) discretized
        # Using central differences for d/dx
        H = np.zeros((N, N), dtype=complex)
        for i in range(N):
            H[i, i] = -1j * 0.5  # the 1/2 term
            if i > 0:
                H[i, i - 1] = -1j * x[i] / (2 * dx) * (-1)
            if i < N - 1:
                H[i, i + 1] = -1j * x[i] / (2 * dx) * (1)

        return H


class BerryKeatingRegularized(OperatorCandidate):
    """
    Regularized Berry-Keating: H = (1/2)(xp + px) + V(x)
    where V(x) provides confinement to give discrete spectrum.

    Sierra & Rodriguez-Laguna (2011) studied V(x) = x^2/2 + l(l+1)/(2x^2)
    """

    name = "Berry-Keating + harmonic"
    description = "H = (xp+px)/2 + x^2/2 (harmonic regularization)"
    reference = "Sierra & Rodriguez-Laguna, 2011"

    def build_matrix(self, N, x_max=50.0, **params):
        dx = x_max / N
        x = np.linspace(dx, x_max, N)

        H = np.zeros((N, N), dtype=complex)
        for i in range(N):
            # xp part
            H[i, i] = -1j * 0.5
            if i > 0:
                H[i, i - 1] = 1j * x[i] / (2 * dx)
            if i < N - 1:
                H[i, i + 1] = -1j * x[i] / (2 * dx)

            # Harmonic potential
            H[i, i] += x[i] ** 2 / 2

        return H


class BBMOperator(OperatorCandidate):
    """
    Bender-Brody-Müller (2017) PT-symmetric operator.

    H = (1 - e^{-i p̂})(x̂ + ℓ) - i(1 - e^{-i p̂})
    (simplified form)

    They claimed this operator has eigenvalues at the zeta zeros.
    Criticism: the operator may not be well-defined on a suitable domain.

    Bender, Brody & Müller, "Hamiltonian for the zeros of the
    Riemann zeta function", PRL 2017
    """

    name = "Bender-Brody-Müller"
    description = "PT-symmetric H involving shift operator"
    reference = "Bender, Brody & Müller, PRL 2017"

    def build_matrix(self, N, **params):
        # BBM involves a shift operator e^{-ip}: f(x) -> f(x-1)
        # Discretize on [0, L] with N points
        L = float(N)  # domain [0, N]
        dx = L / N
        x = np.linspace(dx / 2, L - dx / 2, N)

        # Shift operator (shift by 1 in x-space)
        shift_pixels = int(1.0 / dx)
        S = np.zeros((N, N))
        for i in range(N):
            j = i - shift_pixels
            if 0 <= j < N:
                S[i, j] = 1.0

        I = np.eye(N)
        X = np.diag(x)

        # H = (I - S)(X + ℓI) - i(I - S)
        # Simplified: H = (I - S)(X + ℓI - iI)
        ell = 0.5  # free parameter
        H = (I - S) @ (X + (ell - 1j) * I)

        return H


class QuantumHarmonicXP(OperatorCandidate):
    """
    Connes-inspired: use the quantum harmonic oscillator framework.
    H = a†a formalism but modified to target zeta zeros.

    This tests: can we find eigenvalues near zeta zeros from a simple
    quantum mechanical system with modified potential?

    V(x) = (1/2)x^2 + alpha/x^2 (Calogero-type)
    Eigenvalues: E_n = 2n + 1 + sqrt(alpha + 1/4)
    """

    name = "Calogero-type"
    description = "H = -d^2/dx^2 + x^2/2 + alpha/x^2"
    reference = "Calogero, 1969; tested as RH candidate"

    def build_matrix(self, N, alpha=10.0, x_max=30.0, **params):
        dx = x_max / N
        x = np.linspace(dx, x_max, N)

        # Kinetic energy (second derivative, finite differences)
        T = np.zeros((N, N))
        for i in range(N):
            T[i, i] = 2.0 / dx ** 2
            if i > 0:
                T[i, i - 1] = -1.0 / dx ** 2
            if i < N - 1:
                T[i, i + 1] = -1.0 / dx ** 2

        # Potential
        V = np.diag(x ** 2 / 2 + alpha / x ** 2)

        H = T + V
        return H


class SierraTownsend(OperatorCandidate):
    """
    Sierra-Townsend (2008): Quantum mechanical model related to
    the Landau problem. They found a Hamiltonian whose spectrum
    approximates low-lying zeta zeros.

    H = -d^2/dx^2 + (pi*N(x))^2 where N(x) is the zero counting function.
    This is a WKB-inspired construction.
    """

    name = "Sierra-Townsend WKB"
    description = "H = -d^2/dx^2 + (pi*N(x))^2, WKB-inspired"
    reference = "Sierra & Townsend, 2008"

    def build_matrix(self, N, x_max=100.0, **params):
        dx = x_max / N
        x = np.linspace(dx, x_max, N)

        # Zero counting function approximation
        def N_approx(t):
            if t < 1:
                return 0
            return (t / (2 * np.pi)) * np.log(t / (2 * np.pi)) - t / (2 * np.pi) + 7 / 8

        # Kinetic energy
        T = np.zeros((N, N))
        for i in range(N):
            T[i, i] = 2.0 / dx ** 2
            if i > 0:
                T[i, i - 1] = -1.0 / dx ** 2
            if i < N - 1:
                T[i, i + 1] = -1.0 / dx ** 2

        # Potential from zero counting function
        V_vals = np.array([(np.pi * N_approx(xi)) ** 2 for xi in x])
        V = np.diag(V_vals)

        H = T + V
        return H


# =============================================================================
# PARAMETERIZED SEARCH: H = f(x)p + g(x)
# =============================================================================

class ParameterizedFG(OperatorCandidate):
    """
    General family: H = f(x)p + g(x)
    Symmetrized: H = (1/2)(f(x)p + pf(x)) + g(x)

    Berry-Keating is f(x) = x, g(x) = 0.
    We parameterize f and g as polynomials and search.
    """

    name = "Parameterized f(x)p + g(x)"
    description = "General family with polynomial f,g"
    reference = "Generalization of Berry-Keating"

    def __init__(self, f_coeffs=None, g_coeffs=None):
        if f_coeffs is None:
            f_coeffs = [0, 1]  # f(x) = x (Berry-Keating)
        if g_coeffs is None:
            g_coeffs = [0]  # g(x) = 0
        self.f_coeffs = f_coeffs
        self.g_coeffs = g_coeffs
        self.name = f"f(x)p+g(x) [f={f_coeffs}, g={g_coeffs}]"

    def _eval_poly(self, x, coeffs):
        result = np.zeros_like(x)
        for i, c in enumerate(coeffs):
            result += c * x ** i
        return result

    def build_matrix(self, N, x_max=50.0, **params):
        dx = x_max / N
        x = np.linspace(dx, x_max, N)

        f_vals = self._eval_poly(x, self.f_coeffs)
        g_vals = self._eval_poly(x, self.g_coeffs)

        H = np.zeros((N, N), dtype=complex)

        for i in range(N):
            # g(x) diagonal
            H[i, i] = g_vals[i]

            # Symmetrized f(x)p term: (1/2)(f(x)p + pf(x))
            # = f(x) p + (1/2) f'(x) (-i hbar)
            # Discretized p = -i d/dx
            if i > 0:
                H[i, i - 1] += -1j * f_vals[i] / (2 * dx) * (-1)
            if i < N - 1:
                H[i, i + 1] += -1j * f_vals[i] / (2 * dx) * (1)

            # f'(x) term for symmetrization
            f_prime = 0
            for k, c in enumerate(self.f_coeffs):
                if k > 0:
                    f_prime += k * c * x[i] ** (k - 1)
            H[i, i] += -1j * f_prime / 2

        return H


# =============================================================================
# TEST RUNNER
# =============================================================================

def run_all_tests():
    """Test all operator candidates against known zeros."""
    zeros = load_target_zeros(count=200)
    print(f"Loaded {len(zeros)} target zeta zeros")
    print(f"First zero: {zeros[0]:.6f}, Last: {zeros[-1]:.6f}")
    print("=" * 70)

    candidates = [
        BerryKeatingXP(),
        BerryKeatingRegularized(),
        BBMOperator(),
        QuantumHarmonicXP(),
        SierraTownsend(),
        ParameterizedFG([0, 1], [0]),         # Berry-Keating: f=x, g=0
        ParameterizedFG([0, 1], [0, 0, 0.1]), # f=x, g=0.1*x^2
        ParameterizedFG([0, 1, 0.01], [0]),    # f=x+0.01*x^2, g=0
    ]

    all_results = []

    for candidate in candidates:
        print(f"\nTesting: {candidate.name}")
        print(f"  {candidate.description}")
        print(f"  Ref: {candidate.reference}")
        print("-" * 50)

        # Try different matrix sizes
        for N in [500, 1000]:
            print(f"  N={N}...")
            result = candidate.compare_to_zeros(zeros, N=N, n_compare=30)

            if "error" in result:
                print(f"    ERROR: {result['error']}")
                continue

            print(f"    Eigenvalues computed: {result['n_eigenvalues']}")
            print(f"    Mean relative error: {result['mean_rel_error']:.6f}")
            print(f"    Diverges at eigenvalue #{result['diverges_at']}")

            # Show first few comparisons
            for comp in result["first_10_comparison"][:5]:
                print(f"    #{comp['n']}: eig={comp['eigenvalue']:.4f} "
                      f"vs zero={comp['zero']:.4f} "
                      f"(err={comp['rel_err']:.4f})")

            all_results.append(result)

    # Save results
    outfile = os.path.join(RESULTS_DIR, "operator_tests.json")
    with open(outfile, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nResults saved to {outfile}")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    for r in all_results:
        if "error" not in r:
            print(f"  {r['name']} (N={r['N_matrix']}): "
                  f"mean_rel_err={r['mean_rel_error']:.4f}, "
                  f"diverges_at={r['diverges_at']}")

    return all_results


if __name__ == "__main__":
    run_all_tests()
