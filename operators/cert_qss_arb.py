#!/usr/bin/env python3
"""Interval arithmetic certification of Q_SS >= 0 at [14,44] using arb (python-flint).
Builds Q_SS in the prime eigenbasis and uses LDL^T with bisection to certify min eigenvalue."""

import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from scipy.special import digamma
from scipy.linalg import eigh
from flint import arb, ctx

import arb_fejer_cert2 as ac
import arb_zerofree_cert3 as zc

# Float basis building (same as before)
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

    gw, gv = eigh(G)
    B = gv @ np.diag(1/np.sqrt(gw)) @ gv.T

    return B @ AP @ B, B @ prime @ B, B

def cert_qss_arb(a, b, fourW, N, rc, prec=200):
    """Certify Q_SS >= 0 using arb interval arithmetic.

    Strategy:
      1. Build Q=arch+pole-prime and prime in float
      2. Identify strong block S (prime eigenvalues > floor)
      3. Build Q_SS in arb using arb entries
      4. Use LDL^T bisection to find certified min eigenvalue
    """
    ctx.prec = prec
    W = arb(fourW) / 4
    ac.W = W
    ac.TW = 2 * W
    ac.c_tri = arb.pi() / W

    # Build float forms to identify strong block
    W_float = fourW / 4
    AP_float, prime_float, B = build_float(a, b, W_float, N)

    # Identify strong block in float
    mu, Wv = eigh(prime_float)
    floor = eigh(AP_float, eigvals_only=True)[0]
    idx = np.argsort(mu)[::-1]
    k_strong = int((mu[idx] > floor).sum())

    if k_strong == 0:
        print(f"    [No strong block (floor={floor:.4f}), unable to certify]")
        return None, None, None

    # Extract strong-block eigenvectors
    Wv_sort = Wv[:, idx[:k_strong]]  # shape (N, k_strong)

    # Build in arb: set up basis
    CEN_float = list(np.linspace(a, b, N))
    cs_arb = [arb(c) for c in CEN_float]

    # Build full Q and G matrices in arb
    Q_arb = [[None]*N for _ in range(N)]
    G_arb = [[None]*N for _ in range(N)]
    max_rad = 0.0

    for jj in range(N):
        for kk in range(jj, N):
            # Q entry = arch + pole - prime
            q_entry = (ac.pole_entry(cs_arb[jj], cs_arb[kk]) +
                      zc.arch_entry_v3(cs_arb[jj], cs_arb[kk], arb(rc)) -
                      ac.prime_entry(cs_arb[jj], cs_arb[kk]))
            Q_arb[jj][kk] = Q_arb[kk][jj] = q_entry
            max_rad = max(max_rad, float(q_entry.rad()))

            # G entry
            g_entry = 2 * arb.pi() * ac.gpair_closed_arb(cs_arb[jj], cs_arb[kk], arb(0))
            G_arb[jj][kk] = G_arb[kk][jj] = g_entry

    # Project to strong subspace: Q_SS = W_S^T Q W_S
    Q_SS = [[arb(0)] * k_strong for _ in range(k_strong)]
    G_SS = [[arb(0)] * k_strong for _ in range(k_strong)]

    for i in range(k_strong):
        for j in range(k_strong):
            for ii in range(N):
                for jj in range(N):
                    Q_SS[i][j] += arb(Wv_sort[ii, i]) * Q_arb[ii][jj] * arb(Wv_sort[jj, j])
                    G_SS[i][j] += arb(Wv_sort[ii, i]) * G_arb[ii][jj] * arb(Wv_sort[jj, j])

    # Test PSD via generalized LDL^T
    def ok(lam):
        d = [None] * k_strong
        L = [[arb(0)] * k_strong for _ in range(k_strong)]
        for i in range(k_strong):
            sd = Q_SS[i][i] - lam * G_SS[i][i]
            for kk in range(i):
                sd -= L[i][kk] * L[i][kk] * d[kk]
            if not (sd > 0):
                return False
            d[i] = sd
            for j2 in range(i+1, k_strong):
                v = Q_SS[j2][i] - lam * G_SS[j2][i]
                for kk in range(i):
                    v -= L[j2][kk] * L[i][kk] * d[kk]
                L[j2][i] = v / d[i]
        return True

    # Bisection
    if not ok(arb(0)):
        return None, max_rad, floor

    lo, hi = 0.0, 1.0
    for _ in range(50):
        mid = (lo + hi) / 2
        if ok(arb(mid)):
            lo = mid
        else:
            hi = mid

    return lo, max_rad, floor

# Main
print("=" * 80)
print("Certifying Q_SS >= 0 at [14,44] using interval arithmetic")
print("=" * 80)

a, b = 14, 44

for fourW in [2.50, 2.60, 2.70, 2.75]:
    N = 8  # Fixed basis size
    rc = 180 + int(fourW * 30)  # Scale rc with 4W

    print(f"\n[{a},{b}] 4W={fourW:.2f} N={N} rc={rc}")
    print("  Computing (this may take 30-60s per iteration)...")

    cert_lam, max_rad, floor = cert_qss_arb(a, b, fourW, N, rc, prec=200)

    if cert_lam is not None:
        print(f"  ✓ CERTIFIED: Q_SS >= {cert_lam:.8f}")
        print(f"    Floor = {floor:.6f}, maxrad = {max_rad:.2e}")
    else:
        print(f"  ✗ FAILED to certify (maxrad = {max_rad:.2e}, floor = {floor:.6f})")

print("\n" + "=" * 80)
