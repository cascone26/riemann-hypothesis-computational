#!/usr/bin/env python3
"""Minimal arb certification: Q_SS at [14,44] 4W=2.50 with small basis (N=4) for speed"""

import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from scipy.special import digamma
from scipy.linalg import eigh
from flint import arb, ctx

import arb_fejer_cert2 as ac
import arb_zerofree_cert3 as zc

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

def build_float_minimal(a, b, W, N):
    """Minimal float build for quick strong block identification"""
    fourW = 4*W
    RGRID = np.arange(-(b+120), (b+120), 0.004)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
    CEN = np.linspace(a, b, N)

    def fej(x):
        y = W*x
        out = np.ones_like(y)
        nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz]) / y[nz])**2
        return out

    P = np.array([fej(RGRID-c) + fej(RGRID+c) for c in CEN])
    G = (P @ P.T) * 0.004
    arch = (P * OMEGA[None, :]) @ P.T * 0.004 / (2*np.pi)

    mask = LOGN < fourW
    ln = LOGN[mask]
    lw = LAMW[mask]
    cosu = np.cos(np.outer(ln, RGRID))
    prime = np.zeros((N, N))
    for j in range(N):
        gp = ((P[j][None,:] * P) @ cosu.T) * 0.004 / (2*np.pi)
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

def cert_qss_minimal(a, b, fourW, N, rc, prec):
    """Certify Q_SS with minimal setup"""
    print(f"Certifying Q_SS at [{a},{b}] 4W={fourW} N={N} (prec={prec})...", flush=True)

    # Float to identify strong block
    W_float = fourW / 4
    AP_float, prime_float = build_float_minimal(a, b, W_float, N)

    mu, Wv = eigh(prime_float)
    floor = eigh(AP_float, eigvals_only=True)[0]
    idx = np.argsort(mu)[::-1]
    k = int((mu[idx] > floor).sum())

    if k == 0:
        print(f"  ERROR: No strong block")
        return None, floor

    Wv_S = Wv[:, idx[:k]]
    print(f"  Float analysis: k={k}, floor={floor:.6f}, mu_S={mu[idx][:k]}")

    # Now do arb
    ctx.prec = prec
    W_arb = arb(fourW) / 4
    ac.W = W_arb
    ac.TW = 2 * W_arb
    ac.c_tri = arb.pi() / W_arb

    CEN_float = np.linspace(a, b, N).tolist()
    cs_arb = [arb(c) for c in CEN_float]

    # Build Q in arb
    print(f"  Building Q in arb...", flush=True)
    Q_arb = [[None]*N for _ in range(N)]
    G_arb = [[None]*N for _ in range(N)]

    for jj in range(N):
        for kk in range(jj, N):
            q = (ac.pole_entry(cs_arb[jj], cs_arb[kk]) +
                 zc.arch_entry_v3(cs_arb[jj], cs_arb[kk], arb(rc)) -
                 ac.prime_entry(cs_arb[jj], cs_arb[kk]))
            Q_arb[jj][kk] = Q_arb[kk][jj] = q

            g = 2 * arb.pi() * ac.gpair_closed_arb(cs_arb[jj], cs_arb[kk], arb(0))
            G_arb[jj][kk] = G_arb[kk][jj] = g

    # Project Q_SS
    print(f"  Projecting to strong subspace...", flush=True)
    Q_SS = [[arb(0)] * k for _ in range(k)]
    G_SS = [[arb(0)] * k for _ in range(k)]

    for i in range(k):
        for j in range(k):
            for ii in range(N):
                for jj in range(N):
                    Q_SS[i][j] += arb(Wv_S[ii, i]) * Q_arb[ii][jj] * arb(Wv_S[jj, j])
                    G_SS[i][j] += arb(Wv_S[ii, i]) * G_arb[ii][jj] * arb(Wv_S[jj, j])

    # LDL^T test
    print(f"  Testing PSD...", flush=True)

    def ok(lam):
        d = [None] * k
        L = [[arb(0)] * k for _ in range(k)]
        for i in range(k):
            sd = Q_SS[i][i] - lam * G_SS[i][i]
            for kk in range(i):
                sd -= L[i][kk] * L[i][kk] * d[kk]
            if not (sd > 0):
                return False
            d[i] = sd
            for j2 in range(i+1, k):
                v = Q_SS[j2][i] - lam * G_SS[j2][i]
                for kk in range(i):
                    v -= L[j2][kk] * L[i][kk] * d[kk]
                L[j2][i] = v / d[i]
        return True

    if not ok(arb(0)):
        print(f"  FAILED: Q_SS not PSD")
        return None, floor

    # Bisect
    print(f"  Bisecting...", flush=True)
    lo, hi = 0.0, 1.0
    for _ in range(50):
        mid = (lo + hi) / 2
        if ok(arb(mid)):
            lo = mid
        else:
            hi = mid

    print(f"  ✓ Certified Q_SS >= {lo:.8f}", flush=True)
    return lo, floor

# Main
print("=" * 80)
print("Minimal arb certification of Q_SS at [14,44] corner")
print("=" * 80)

a, b = 14, 44

# Start with 4W=2.50 (easiest - largest margin)
cert_lam, floor = cert_qss_minimal(a, b, 2.50, 8, 270, 180)
if cert_lam:
    print(f"\n✓ HANDLE: Q_SS >= {cert_lam:.8f} at [14,44] 4W=2.50")
else:
    print(f"\n✗ Failed to certify at 4W=2.50")

# Try progressively larger 4W
if cert_lam:
    print(f"\nTrying larger 4W values...")
    for fourW_val in [2.60, 2.70, 2.75]:
        cert_lam, floor = cert_qss_minimal(a, b, fourW_val, 4, 220, 200)
        if cert_lam:
            print(f"  ✓ 4W={fourW_val}: Q_SS >= {cert_lam:.8f}")
        else:
            print(f"  ✗ 4W={fourW_val}: Failed")
