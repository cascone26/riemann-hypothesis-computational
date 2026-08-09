#!/usr/bin/env python3
"""Quick test: certify Q_SS at [14,44] 4W=2.5 (shallow corner)"""

import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from scipy.special import digamma
from scipy.linalg import eigh

# Float computation first
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

def build_float(a, b, W, N, dr=0.004):
    """Build arch+pole and prime matrices in float"""
    fourW = 4*W
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

    # Normalize
    gw, gv = eigh(G)
    B = gv @ np.diag(1/np.sqrt(gw)) @ gv.T

    return B @ AP @ B, B @ prime @ B

# Test [14,44] 4W=2.5
a, b = 14, 44
fourW = 2.5
W = fourW / 4
N = 8

print(f"Building [{a},{b}] 4W={fourW:.2f} N={N}")
AP_norm, prime_norm = build_float(a, b, W, N)

mu, Wv = eigh(prime_norm)
floor = eigh(AP_norm, eigvals_only=True)[0]

idx = np.argsort(mu)[::-1]
mu_sort = mu[idx]
Wv_sort = Wv[:, idx]

k_strong = int((mu_sort > floor).sum())

print(f"Floor={floor:.6f}, strong block k={k_strong}")
print(f"Prime eigenvalues (desc): {mu_sort[:5]}")

# Q in prime eigenbasis
Q = AP_norm - prime_norm
Qp = Wv.T @ Q @ Wv
Qp_sort = Qp[np.ix_(np.argsort(mu)[::-1], np.argsort(mu)[::-1])]

Q_SS = Qp_sort[:k_strong, :k_strong]
print(f"\nQ_SS shape: {Q_SS.shape}")
print(f"Q_SS min eigenvalue (float): {eigh(Q_SS, eigvals_only=True)[0]:.6e}")
print(f"Q_SS eigenvalues (float): {np.round(eigh(Q_SS, eigvals_only=True)[:5], 6)}")
