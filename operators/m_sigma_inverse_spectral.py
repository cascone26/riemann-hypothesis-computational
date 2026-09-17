#!/usr/bin/env python3
"""§990 (CANDIDATE A) — Inverse-spectral reconstruction from M_σ's spectral data.

Approach: Given M_σ = arch − prime's eigendecomposition (eigenvalues λ_k, eigenvectors v_k),
treat its spectral measure as an indefinite measure (one negative atom at κ=1) and apply
Krein's indefinite moment problem to reconstruct a Jacobi matrix whose m-function inherits
the explicit-formula structure. Test whether the pole enters naturally.

Key difference from Suzuki (truncation): directly inverse-spectralize M_σ itself, not via
a characteristic function. If successful, the κ=1 negative square should appear as a pole
in the reconstructed m-function.

Status: REPORT (implementation + numerical test). This is a probe on the §981 frontier.
"""

import numpy as np
from scipy.special import digamma
from scipy.linalg import solve
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# ────────────────────────────────────────────────────────────────────────────────────────
# STEP 1: Compute M_σ's full eigendecomposition (extending kappa_pole_location_FIXED.py)
# ────────────────────────────────────────────────────────────────────────────────────────

def vonmangoldt(N):
    lam = np.zeros(N+1)
    isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True
            lp = np.log(p)
            pk = p
            while pk <= N:
                lam[pk] = lp
                pk *= p
    return lam

def chi4(n):
    n %= 4
    return {1: 1, 3: -1}.get(n, 0)

def chi8_2(n):
    n %= 8
    return {1: 1, 7: 1, 3: -1, 5: -1}.get(n, 0)

def chi8_m2(n):
    n %= 8
    return {1: 1, 3: 1, 5: -1, 7: -1}.get(n, 0)

def jacobi_odd(a, n):
    a %= n
    r = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                r = -r
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            r = -r
        a %= n
    return r if n == 1 else 0

def chi_odd(q):
    return lambda n: jacobi_odd(n % q, q)

def build_M_sigma(cen_lo=14, cen_hi=44, fourW=5.0, dr=0.004, frac=0.8, is_zeta=True):
    """Build M_σ = arch − prime, return both matrix and data for spectral analysis."""
    NMAX = 8000
    LAM = vonmangoldt(NMAX)
    NS = np.arange(2, NMAX + 1)
    LOGN = np.log(NS)
    LAMv = LAM[2:]

    # Default: zeta
    chi = None
    a_par = 0
    q = 1

    W = fourW / 4
    L = cen_hi - cen_lo
    N = max(4, round(frac * L * 2 * W / np.pi))

    RG = np.arange(-(cen_hi + 120), (cen_hi + 120), dr)
    OM = np.real(digamma((a_par + 0.5) / 2 + 0.5j * RG)) + np.log(q / np.pi)
    CEN = np.linspace(cen_lo, cen_hi, N)

    def fej(x):
        y = W * x
        o = np.ones_like(y)
        nz = np.abs(y) > 1e-12
        o[nz] = (np.sin(y[nz]) / y[nz]) ** 2
        return o

    P = np.array([fej(RG - c) + fej(RG + c) for c in CEN])
    G = (P @ P.T) * dr
    arch = (P * OM[None, :]) @ P.T * dr / (2 * np.pi)

    mask = LOGN < fourW
    ln = LOGN[mask]
    chiv = np.ones(np.sum(mask))  # zeta: all 1
    wts = LAMv[mask] * chiv / np.sqrt(NS[mask])

    cosu = np.cos(np.outer(ln, RG))
    prime = np.zeros((N, N))
    for j in range(N):
        gp = ((P[j][None, :] * P) @ cosu.T) * dr / (2 * np.pi)
        prime[j, :] = 2.0 * (gp * wts[None, :]).sum(axis=1)

    M_sigma = arch - 0.5 * (prime + prime.T)
    return M_sigma, arch, prime, CEN, G

def compute_eigendecomp(M_sigma):
    """Full eigendecomposition of M_σ, sorted descending by eigenvalue."""
    eigvals, eigvecs = np.linalg.eigh(M_sigma)
    # Sort descending
    idx = np.argsort(-eigvals)
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    return eigvals, eigvecs

# ────────────────────────────────────────────────────────────────────────────────────────
# STEP 2: Construct indefinite spectral measure and inverse-spectral Jacobi matrix
# ────────────────────────────────────────────────────────────────────────────────────────

def lanczos_indefinite(diag_vals, mu, K, ref_vec=None, full_orth=True):
    """
    Generalized Lanczos for an indefinite (signed) measure.

    Given a diagonal matrix D = diag(diag_vals) and a weight/spectral measure μ
    (can have negative components), construct Jacobi matrix coefficients (a_n, b_n).

    ref_vec: reference direction for spectral weights (default: first basis vector).
    Returns: (a_coeffs, b_coeffs) as arrays of length K.

    This is an INDEFINITE generalization of krein_string_zeta.py's positive Lanczos.
    For an indefinite measure, the resulting Jacobi matrix may have indefinite inertia.
    """
    n_dof = len(diag_vals)
    if ref_vec is None:
        ref_vec = np.zeros(n_dof)
        ref_vec[0] = 1.0

    # Initial vector: weighted by sqrt(abs(mu)), sign-corrected
    v = np.sqrt(np.abs(mu)) * ref_vec
    norm_v = np.linalg.norm(v)
    if norm_v < 1e-14:
        # Degenerate; use random
        v = np.random.randn(n_dof)
        v = v / np.linalg.norm(v)
    else:
        v = v / norm_v

    # Indefinite inner product: <u, w> = Σ_i sign(mu_i) * u_i * w_i
    def indef_inner(u, w):
        return np.sum(np.sign(np.maximum(np.abs(mu), 1e-15)) * u * w)

    # Start Lanczos
    V = [v]  # Orthogonal (in indefinite sense) basis
    a = []
    b = []

    w = diag_vals * v

    for k in range(K):
        # Diagonal entry
        ak = indef_inner(v, w)
        a.append(ak)

        # Residual
        w = w - ak * v
        if k > 0 and len(b) > 0:
            w = w - b[-1] * V[-2]

        # Full indefinite reorthogonalization
        if full_orth:
            for u in V:
                coeff = indef_inner(u, w)
                w = w - coeff * u
            for u in V:
                coeff = indef_inner(u, w)
                w = w - coeff * u

        # Off-diagonal coefficient
        bk = np.sqrt(np.abs(indef_inner(w, w)))
        b.append(bk)

        if bk < 1e-14:
            # Breakdown; pad with zeros
            for _ in range(K - k - 1):
                a.append(0.0)
                b.append(0.0)
            break

        vn = w / bk
        V.append(vn)
        v = vn
        w = diag_vals * v

    return np.array(a), np.array(b)

# ────────────────────────────────────────────────────────────────────────────────────────
# STEP 3: Weyl m-function from Jacobi matrix
# ────────────────────────────────────────────────────────────────────────────────────────

def jacobi_weyl_mfunction(a, b, z, backward=True):
    """
    Weyl m-function for a Jacobi matrix from its coefficients.

    Tridiagonal: J = diag(a) + diag(b,1) + diag(b,-1)
    Weyl m-function via continued fraction (backward recurrence).
    m(z) = -b_0²/(z - a_0 - b_1²/(z - a_1 - b_2²/(z - a_2 - ...)))

    For indefinite inertia, m can have poles and zeros off the real axis.
    """
    K = len(a)

    if backward:
        # Start from the end
        m_chain = 0.0  # m_K = 0
        for k in range(K - 1, -1, -1):
            denom = z - a[k]
            if k < K - 1:
                denom = denom - b[k] ** 2 / (m_chain + 1e-300j)
            m_chain = -b[k - 1] ** 2 / (denom + 1e-300j) if k > 0 else m_chain

    return m_chain

# More robust: forward recurrence in a different form
def jacobi_weyl_forward(a, b, z):
    """
    Compute m-function via a stable forward recurrence
    (the solution of the ODE form of the recurrence relation).
    """
    K = len(a)
    # P_n(z) / P_{n-1}(z) recurrence
    # Characteristic: P_n(z) = (z - a_n) P_{n-1}(z) - b_{n-1}² P_{n-2}(z)

    # Boundary: P_0 = 1, P_{-1} = 0
    # m_n = -P'_n(z) / P_n(z) = -(1/(z-a_n) + something)

    # Use Weyl continued fraction: m = -b_0² / (z - a_0 - b_1²/(z - a_1 - ...))
    # But compute stably via a vector recurrence

    # Standard: numerator/denominator approach
    # m_n(z) satisfies: m_n = -b_n² / (z - a_{n+1} - m_{n+1})
    # Start from m_K = 0 (or ∞), step backward

    m = 0.0
    for n in range(K - 1, -1, -1):
        if n < K - 1 and abs(b[n]) > 1e-15:
            m = -b[n] ** 2 / (z - a[n] - m + 1e-300j)
        else:
            m = 0.0
    return m

# ────────────────────────────────────────────────────────────────────────────────────────
# STEP 4: Nevanlinna kernel and κ count
# ────────────────────────────────────────────────────────────────────────────────────────

def nevanlinna_kernel_kappa(m_func, z_test, tol=1e-9):
    """
    Compute κ = #negative eigenvalues of the Nevanlinna reproducing kernel
    K(z, w) = (m(z) - conj(m(w))) / (z - conj(w))
    for a generalized-Nevanlinna function m.

    This is the indefinite-metric version used in §974-984.
    """
    Z = np.array(z_test, dtype=complex)
    n = len(Z)

    # Evaluate m at test points
    mv = np.array([m_func(z) for z in Z], dtype=complex)

    # Build Nevanlinna kernel
    K = np.empty((n, n), dtype=complex)
    for i in range(n):
        for j in range(n):
            d = Z[i] - np.conj(Z[j])
            if abs(d) < 1e-12:
                # Diagonal: K(z,z) = m'(z) / (2i)  (via L'Hopital)
                # Approximate via finite difference
                eps = 1e-6
                m_plus = m_func(Z[i] + eps)
                m_minus = m_func(Z[i] - eps)
                K[i, j] = (m_plus - m_minus) / (2.0j * eps)
            else:
                K[i, j] = (mv[i] - np.conj(mv[j])) / d

    # Symmetrize
    K = 0.5 * (K + K.conj().T)

    # Eigenvalues
    ev = np.linalg.eigvalsh(K)
    s = max(1.0, np.abs(ev).max())
    kappa = int(np.sum(ev < -tol * s))

    return kappa, ev

# ────────────────────────────────────────────────────────────────────────────────────────
# MAIN: Run the construction
# ────────────────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 100)
    print("§990-A: Inverse-Spectral M_σ Construction")
    print("=" * 100)

    # Build M_σ
    print("\n[1] Building M_σ = arch − prime (window [14,44], 4W=5.0)...", flush=True)
    M_sigma, arch, prime, CEN, G = build_M_sigma(cen_lo=14, cen_hi=44, fourW=5.0)
    print(f"    M_σ shape: {M_sigma.shape}")

    # Eigendecompose
    print("\n[2] Eigendecomposition of M_σ...", flush=True)
    eigvals, eigvecs = compute_eigendecomp(M_sigma)
    kappa_exact = int(np.sum(eigvals < -1e-9 * np.abs(eigvals).max()))
    print(f"    Eigenvalues (sorted desc): {eigvals[:5]} ... {eigvals[-2:]}")
    print(f"    κ(M_σ) = {kappa_exact} (exact count of negatives)")
    print(f"    Smallest negative eig: {eigvals[kappa_exact] if kappa_exact < len(eigvals) else 'none'}")

    # The pole direction: the principal direction (first eigenvector)
    # and the negative direction: the (κ)-th eigenvector (first negative)
    pole_approx = eigvecs[:, 0]  # Principal direction (largest eigenvalue)
    neg_eig = eigvecs[:, kappa_exact] if kappa_exact < len(eigvecs) else None

    print(f"\n    Principal direction (pole approx): norm = {np.linalg.norm(pole_approx):.4f}")
    if neg_eig is not None:
        alignment = np.abs(np.dot(pole_approx, neg_eig))
        print(f"    Alignment (pole, negative-eig): {alignment:.6f}")

    # Spectral measure
    print("\n[3] Building indefinite spectral measure from M_σ...", flush=True)

    # The spectral measure is μ = Σ_k λ_k δ_{λ_k}
    # weighted by spectral projection coefficients
    # Natural choice: use |⟨e_0, v_k⟩|² (first basis vector) as weight
    weights_ref = np.abs(eigvecs[0, :]) ** 2
    print(f"    Spectral weights from 1st coordinate: min={weights_ref.min():.3e}, max={weights_ref.max():.3e}")
    print(f"    Sum of weights: {weights_ref.sum():.6f}")

    # Construct Jacobi matrix from this spectral data
    K_jacobi = min(20, len(eigvals) - 1)
    print(f"\n[4] Lanczos on indefinite spectral measure (K={K_jacobi} terms)...", flush=True)

    # Indefinite Lanczos: diagonal values = eigenvalues, weights = spectral weights
    # Use eigenvectors as the basis (orthonormal in standard metric)
    mu_weighted = eigvals * weights_ref  # Signed weights

    a_jacob, b_jacob = lanczos_indefinite(eigvals, mu_weighted, K=K_jacobi,
                                          ref_vec=weights_ref, full_orth=True)

    print(f"    a_n (diagonal): {a_jacob[:5]} ...")
    print(f"    b_n (off-diag): {b_jacob[:5]} ...")
    print(f"    b_n² (squared): {b_jacob[:5]**2} ...")

    # Weyl m-function
    print(f"\n[5] Computing Weyl m-function from Jacobi matrix...", flush=True)

    def m_func(z):
        return jacobi_weyl_forward(a_jacob, b_jacob, z)

    # Test points: UHP
    z_test_uhp = [complex(x, y) for x in np.linspace(5, 50, 30) for y in [0.1, 0.5, 1.0]]
    m_vals_sample = [m_func(z) for z in z_test_uhp[:5]]
    print(f"    Sample m-values at UHP: {m_vals_sample[:3]}")

    # κ via Nevanlinna kernel
    print(f"\n[6] Computing κ via Nevanlinna kernel (corrected)...", flush=True)
    z_test = [complex(x, y) for x in np.linspace(10, 50, 15) for y in np.linspace(0.1, 2.0, 10)]
    kappa_neo, evals_neo = nevanlinna_kernel_kappa(m_func, z_test)
    print(f"    κ(constructed m-function) = {kappa_neo}")
    print(f"    Kernel eigs: min={evals_neo.min():.3e}, max={evals_neo.max():.3e}")

    # Comparison with original
    print(f"\n[RESULT]")
    print(f"  κ(M_σ direct) = {kappa_exact}")
    print(f"  κ(reconstructed m-function) = {kappa_neo}")
    print(f"  Match: {kappa_neo == kappa_exact}")

    if kappa_neo == kappa_exact:
        print(f"\n  ✓ SUCCESS: Inverse-spectral reconstruction REPRODUCES κ.")
        print(f"    The negative square in M_σ is faithfully encoded in the m-function.")
    else:
        print(f"\n  ✗ MISMATCH: κ values do not agree.")
        print(f"    This indicates either a category error (Nevanlinna sampling) or a")
        print(f"    structural issue with the indefinite Lanczos reconstruction.")

    print("\n" + "=" * 100)
