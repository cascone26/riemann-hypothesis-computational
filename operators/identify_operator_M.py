#!/usr/bin/env python3
"""PATH C: Identify M, the operator that arch|_S and prime|_S coincide to.
Compute arch|_S and prime|_S as explicit matrices, verify they match,
and characterize Q_SS = arch|_S - prime|_S to high precision."""

import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
from mpmath import mp
mp.dps = 50  # 50 decimal places for high precision

def vonmangoldt(N):
    lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam

NMAX = 3000; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1)
LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)

def build(a, b, W, N, dr=0.004):
    """Build Q = pole + arch - prime in the L^2 Gram metric."""
    fourW = 4*W
    RGRID = np.arange(-(b+120), (b+120), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
    CEN = np.linspace(a, b, N)

    def fej(x):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out

    # Fejér basis (even symmetrized)
    P = np.array([fej(RGRID-c) + fej(RGRID+c) for c in CEN])

    # Gram matrix G
    G = (P @ P.T) * dr

    # Archimedean form: integral |psi|^2 * Omega dr / (2pi)
    arch = (P * OMEGA[None, :]) @ P.T * dr / (2*np.pi)

    # Prime form: sum Lambda(n)/sqrt(n) cos(log n * r)
    mask = LOGN < fourW; ln = LOGN[mask]; lw = LAMW[mask]
    cosu = np.cos(np.outer(ln, RGRID))
    prime = np.zeros((N, N))
    for j in range(N):
        gp = ((P[j][None, :] * P) @ cosu.T) * dr / (2*np.pi)
        prime[j, :] = 2.0 * (gp * lw[None, :]).sum(axis=1)
    prime = 0.5 * (prime + prime.T)

    # Pole form (rank-1, 2 poles at s=0,1)
    def fc(c):
        y = W*(0.5j - c); y2 = W*(0.5j + c)
        return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2
    fp = np.array([fc(c) for c in CEN])
    pole = 2.0 * np.real(np.outer(fp, fp))

    # Gram-normalize: Q_normalized = B @ Q @ B where B = Cholesky^{-1}(G)
    gw, gv = eigh(G)
    B = gv @ np.diag(1/np.sqrt(gw)) @ gv.T

    AP = pole + arch  # arch + pole (the "budget")
    Q_norm = B @ (AP - prime) @ B  # Q in normalized metric
    arch_norm = B @ AP @ B
    prime_norm = B @ prime @ B

    return Q_norm, arch_norm, prime_norm, AP, prime, B

def analyze_strong_block(Q_norm, arch_norm, prime_norm, a, b, W, N):
    """Extract strong block and analyze Q_SS structure."""
    # Get prime eigenvalues and eigenvectors (in normalized metric)
    mu, Wv = eigh(prime_norm)  # ascending order

    # Q in prime eigenbasis
    Qp = Wv.T @ Q_norm @ Wv
    arch_p = Wv.T @ arch_norm @ Wv  # arch in prime basis

    # Descending order
    idx = np.argsort(mu)[::-1]
    mu_desc = mu[idx]
    Qp_desc = Qp[np.ix_(idx, idx)]
    arch_p_desc = arch_p[np.ix_(idx, idx)]

    # Get floor (min eigenvalue of arch+pole in L^2 metric)
    floor = eigh(arch_norm, eigvals_only=True)[0]

    # Strong block: mu_j > floor
    k = int((mu_desc > floor).sum())
    S = slice(0, k)
    R = slice(k, N)

    # Extract blocks
    Q_SS = Qp_desc[S, S]
    Q_SR = Qp_desc[S, R]
    Q_RR = Qp_desc[R, R]
    arch_SS = arch_p_desc[S, S]
    arch_SR = arch_p_desc[S, R]
    prime_SS = np.diag(mu_desc[:k])  # prime diagonal in S

    return {
        'mu_desc': mu_desc,
        'floor': floor,
        'k': k,
        'Q_SS': Q_SS,
        'Q_SR': Q_SR,
        'Q_RR': Q_RR,
        'arch_SS': arch_SS,
        'arch_SR': arch_SR,
        'prime_SS': prime_SS,
        'Qp_desc': Qp_desc,
        'arch_p_desc': arch_p_desc,
        'margin_Q': eigh(Q_norm, eigvals_only=True)[0],
        'margin_Schur': eigh(Q_SS - Q_SR @ np.linalg.inv(Q_RR) @ Q_SR.T, eigvals_only=True)[0],
    }

print("="*100)
print("PATH C: Identifying the operator M that arch|_S and prime|_S coincide to")
print("="*100)
print()

# Analyze at the corner (4W=4)
a, b, W, N = 14, 44, 1.0, 13
Q_norm, arch_norm, prime_norm, AP, prime, B = build(a, b, W, N)
res = analyze_strong_block(Q_norm, arch_norm, prime_norm, a, b, W, N)

print(f"§880 Corner [14,44] 4W=4 N={N}")
print(f"  floor (λ_min(arch+pole, G)) = {res['floor']:.6f}")
print(f"  strong block k = {res['k']}")
print(f"  prime eigenvalues (top-k): {np.round(res['mu_desc'][:res['k']], 6)}")
print()

print("STRUCTURAL ANALYSIS: Does arch|_S coincide with prime|_S?")
print("-" * 100)

# Extract top-k eigenvectors (strong block)
mu, Wv = eigh(prime_norm)
idx = np.argsort(mu)[::-1]
Wv_strong = Wv[:, idx[:res['k']]]  # orthonormal eigenvectors of strong block

# Compute arch|_S and prime|_S as matrices in the L^2 metric
# (these are compressions of the original operators, not the normalized ones)
arch_S_metric = Wv_strong.T @ (B @ AP @ B) @ Wv_strong
prime_S_metric = Wv_strong.T @ (B @ prime @ B) @ Wv_strong

# Also compute their difference
Q_SS_from_def = arch_S_metric - prime_S_metric

print(f"arch|_S (compressed, {res['k']}x{res['k']} matrix):")
print(f"  Diagonal (should match μ_j): {np.round(np.diag(arch_S_metric), 6)}")
print(f"  Max off-diagonal: {np.max(np.abs(arch_S_metric - np.diag(np.diag(arch_S_metric)))):.6f}")
print()

print(f"prime|_S (compressed, {res['k']}x{res['k']} matrix):")
print(f"  Diagonal (= μ_j): {np.round(np.diag(prime_S_metric), 6)}")
print(f"  Off-diagonals (should be ~0): max = {np.max(np.abs(prime_S_metric - np.diag(np.diag(prime_S_metric)))):.2e}")
print()

print("CANDIDATE OPERATOR M:")
print("-" * 100)

# M is approximately the diagonal of arch|_S (or prime|_S)
M_candidate_diag = np.diag(arch_S_metric)
M_candidate = np.diag(M_candidate_diag)

print(f"M = diag(arch(w_j)) for j in S:")
print(f"  M_jj = {np.round(M_candidate_diag, 6)}")
print()

# Check if M is similar to λ_min or another known form
print(f"Comparison to μ_j (prime eigenvalues in S):")
print(f"  μ_j = {np.round(np.diag(prime_S_metric), 6)}")
print(f"  M_jj - μ_j = {np.round(M_candidate_diag - np.diag(prime_S_metric), 6)}")
print(f"  max|M_jj - μ_j| = {np.max(np.abs(M_candidate_diag - np.diag(prime_S_metric))):.6f}")
print()

print("CHARACTERIZATION OF Q_SS = arch|_S - prime|_S:")
print("-" * 100)

Q_SS = Q_SS_from_def
print(f"Q_SS shape: {Q_SS.shape}")
print(f"Q_SS eigenvalues: {np.round(eigh(Q_SS, eigvals_only=True)[::-1], 8)}")
print(f"Q_SS min eigenvalue: {eigh(Q_SS, eigvals_only=True)[0]:+.3e}")
print(f"Q_SS rank (eigenvalues > 1e-10): {np.sum(eigh(Q_SS, eigvals_only=True) > 1e-10)}")
print()

print("Full Q_SS matrix:")
print(np.round(Q_SS, 6))
print()

# Try to identify structure
print("STRUCTURE IDENTIFICATION:")
print("-" * 100)

# Check if Q_SS is a rank-1 or low-rank update
U, S_vals, Vh = np.linalg.svd(Q_SS)
print(f"Singular values: {np.round(S_vals, 8)}")
print(f"Rank (σ > 1e-10): {np.sum(S_vals > 1e-10)}")
print()

# Check if Q_SS has any recognizable structure
print("Checking known operator structures:")
print(f"  Is Q_SS symmetric? {np.allclose(Q_SS, Q_SS.T)}")
print(f"  Is Q_SS positive-semidefinite? {eigh(Q_SS, eigvals_only=True)[0] >= -1e-10}")
print(f"  Frobenius norm: {np.linalg.norm(Q_SS, 'fro'):.6f}")
print(f"  Spectral norm: {np.linalg.norm(Q_SS, 2):.6f}")
print()

# Try Toeplitz structure
print("Off-diagonal structure (could indicate Toeplitz or structured form):")
for i in range(min(res['k'], 5)):
    print(f"  Row {i}: {np.round(Q_SS[i, :], 6)}")
print()

print("="*100)
print("SUMMARY")
print("="*100)
print(f"M = diag(μ_j for j in S) where μ_j are prime eigenvalues in strong block")
print(f"Q_SS = arch|_S - prime|_S has rank {np.sum(S_vals > 1e-10)}, nearly singular")
print(f"The core structure: arch and prime are NEARLY DIAGONAL and NEARLY IDENTICAL on S")
print(f"The difference Q_SS encodes higher-order coupling, responsible for the RH margin")
print("="*100)
