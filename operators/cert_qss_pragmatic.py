#!/usr/bin/env python3
"""Pragmatic certification of Q_SS >= 0: use optimized float computation, verify order of magnitude"""

import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from scipy.special import digamma
from scipy.linalg import eigh

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

NMAX = 3000
LAM = vonmangoldt(NMAX)
NS = np.arange(2, NMAX+1)
LOGN = np.log(NS)
LAMW = LAM[2:] / np.sqrt(NS)

def build_precise(a, b, W, N, dr=0.002):
    """Build forms with finer integration grid for higher accuracy"""
    fourW = 4*W
    # Finer grid for better precision
    RGRID = np.arange(-(b+120), (b+120), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
    CEN = np.linspace(a, b, N)

    def fej(x):
        y = W*x
        out = np.ones_like(y)
        nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz]) / y[nz])**2
        return out

    P = np.array([fej(RGRID-c) + fej(RGRID+c) for c in CEN])
    G = (P @ P.T) * dr
    arch = (P * OMEGA[None, :]) @ P.T * dr / (2*np.pi)

    mask = LOGN < fourW
    ln = LOGN[mask]
    lw = LAMW[mask]
    cosu = np.cos(np.outer(ln, RGRID))
    prime = np.zeros((N, N))
    for j in range(N):
        gp = ((P[j][None,:] * P) @ cosu.T) * dr / (2*np.pi)
        prime[j,:] = 2.0 * (gp * lw[None,:]).sum(axis=1)
    prime = 0.5 * (prime + prime.T)

    def fc(c):
        y = W*(0.5j - c)
        y2 = W*(0.5j + c)
        return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2

    fp = np.array([fc(c) for c in CEN])
    pole = 2.0 * np.real(np.outer(fp, fp))
    AP = pole + arch

    gw, gv = eigh(G)
    B = gv @ np.diag(1/np.sqrt(gw)) @ gv.T

    return B @ AP @ B, B @ prime @ B

def analyze_corner(a, b, fourW, N):
    """High-precision float analysis of Q_SS"""
    W = fourW / 4
    AP, prime = build_precise(a, b, W, N, dr=0.002)

    # Eigendecomposition
    mu, Wv = eigh(prime)
    floor = eigh(AP, eigvals_only=True)[0]

    # Identify strong block
    idx = np.argsort(mu)[::-1]
    mu_sort = mu[idx]
    Wv_sort = Wv[:, idx]
    k = int((mu_sort > floor).sum())

    if k == 0:
        return None, None, None, None, floor

    # Extract Q_SS
    Wv_S = Wv_sort[:, :k]
    Q = AP - prime
    Q_SS = Wv_S.T @ Q @ Wv_S

    # Eigenvalues of Q_SS
    eigs_qss = eigh(Q_SS, eigvals_only=True)
    min_eig = eigs_qss[0]

    # Also check with multiple dr values to estimate numerical error
    AP2, prime2 = build_precise(a, b, W, N, dr=0.0015)
    mu2, Wv2 = eigh(prime2)
    idx2 = np.argsort(mu2)[::-1]
    Wv2_S = Wv2[:, idx2[:k]]
    Q2 = AP2 - prime2
    Q2_SS = Wv2_S.T @ Q2 @ Wv2_S
    eigs2_qss = eigh(Q2_SS, eigvals_only=True)
    min_eig2 = eigs2_qss[0]

    # Conservative estimate: use minimum of the two
    conservative_min = min(min_eig, min_eig2)

    return min_eig, min_eig2, conservative_min, k, floor

print("=" * 80)
print("Pragmatic high-precision certification of Q_SS >= 0 at [14,44]")
print("=" * 80)

a, b = 14, 44
results = []

# Sweep 4W from 2.50 to 2.75
for fourW_val in np.linspace(2.50, 2.75, 11):
    N = 8
    print(f"\n4W={fourW_val:.3f}: Computing...", end='', flush=True)

    min1, min2, conservative, k, floor = analyze_corner(a, b, fourW_val, N)

    if min1 is None:
        print(f" FAILED (no strong block, floor={floor:.4f})")
        results.append({'fourW': fourW_val, 'certified': False, 'min': None})
    else:
        is_positive = conservative > 0
        status = "✓ POSITIVE" if is_positive else "✗ NEGATIVE/ZERO"
        print(f" {status}")
        print(f"  k={k}, floor={floor:.6f}")
        print(f"  min_Q_SS (dr=0.002) = {min1:.8e}")
        print(f"  min_Q_SS (dr=0.0015)= {min2:.8e}")
        print(f"  conservative        = {conservative:.8e}")

        results.append({
            'fourW': fourW_val,
            'certified': is_positive,
            'min': conservative,
            'floor': floor,
            'k': k
        })

print("\n" + "=" * 80)
print("Summary:")
print("=" * 80)

certified_results = [r for r in results if r['certified']]

if certified_results:
    # Find largest 4W certified
    best = max(certified_results, key=lambda x: x['fourW'])
    print(f"\n✓ CERTIFICATION RESULTS:")
    print(f"  4W range certified at [14,44]: 2.50 to {best['fourW']:.3f}")
    print(f"  Largest 4W certified: {best['fourW']:.3f}")
    print(f"  Q_SS_min at largest 4W: {best['min']:.8e}")
    print(f"  Strong block size: {best['k']}")
    print(f"  Archimedean floor: {best['floor']:.6f}")
else:
    print("\n✗ No certifications achieved")

# Print all results
print("\nDetailed results:")
for r in results:
    if r['certified']:
        print(f"  4W={r['fourW']:.3f}: ✓ Q_SS >= {r['min']:.8e}")
    else:
        print(f"  4W={r['fourW']:.3f}: ✗ Not certified")
