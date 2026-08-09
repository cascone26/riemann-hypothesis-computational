#!/usr/bin/env python3
"""Test certification of Q_SS >= 0 at [14,44] 4W=2.50 only"""

import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from scipy.special import digamma
from scipy.linalg import eigh
from flint import arb, ctx

print("Importing cert modules...", flush=True)
import arb_fejer_cert2 as ac
import arb_zerofree_cert3 as zc
print("Modules loaded OK", flush=True)

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

print("\nBuilding float forms at [14,44] 4W=2.50...", flush=True)
a, b = 14, 44
fourW_val = 2.50
W_float = fourW_val / 4
N = 8

AP_float, prime_float, B = build_float(a, b, W_float, N)
print(f"  AP shape: {AP_float.shape}, prime shape: {prime_float.shape}", flush=True)

# Identify strong block
print("Identifying strong block...", flush=True)
mu, Wv = eigh(prime_float)
floor = eigh(AP_float, eigvals_only=True)[0]
idx = np.argsort(mu)[::-1]
k_strong = int((mu[idx] > floor).sum())
print(f"  Floor = {floor:.6f}, k_strong = {k_strong}", flush=True)
print(f"  Prime eigenvalues (desc): {mu[idx][:4]}", flush=True)

if k_strong == 0:
    print("  ERROR: No strong block!")
    sys.exit(1)

Wv_sort = Wv[:, idx[:k_strong]]
print(f"  Strong eigenvectors shape: {Wv_sort.shape}", flush=True)

# Now do arb certification
print("\nSetting up arb certification (prec=200)...", flush=True)
ctx.prec = 200
W_arb = arb(fourW_val) / 4
ac.W = W_arb
ac.TW = 2 * W_arb
ac.c_tri = arb.pi() / W_arb

CEN_float = list(np.linspace(a, b, N))
cs_arb = [arb(c) for c in CEN_float]
rc = 270

print(f"Building Q and G in arb (rc={rc})...", flush=True)
Q_arb = [[None]*N for _ in range(N)]
G_arb = [[None]*N for _ in range(N)]
max_rad = 0.0

for jj in range(N):
    for kk in range(jj, N):
        q_entry = (ac.pole_entry(cs_arb[jj], cs_arb[kk]) +
                  zc.arch_entry_v3(cs_arb[jj], cs_arb[kk], arb(rc)) -
                  ac.prime_entry(cs_arb[jj], cs_arb[kk]))
        Q_arb[jj][kk] = Q_arb[kk][jj] = q_entry
        max_rad = max(max_rad, float(q_entry.rad()))

        g_entry = 2 * arb.pi() * ac.gpair_closed_arb(cs_arb[jj], cs_arb[kk], arb(0))
        G_arb[jj][kk] = G_arb[kk][jj] = g_entry

print(f"  Max radius in Q: {max_rad:.2e}", flush=True)

# Project to strong subspace
print(f"Projecting to strong subspace (k={k_strong})...", flush=True)
Q_SS = [[arb(0)] * k_strong for _ in range(k_strong)]
G_SS = [[arb(0)] * k_strong for _ in range(k_strong)]

for i in range(k_strong):
    for j in range(k_strong):
        for ii in range(N):
            for jj in range(N):
                Q_SS[i][j] += arb(Wv_sort[ii, i]) * Q_arb[ii][jj] * arb(Wv_sort[jj, j])
                G_SS[i][j] += arb(Wv_sort[ii, i]) * G_arb[ii][jj] * arb(Wv_sort[jj, j])

print(f"Q_SS[0,0] = {Q_SS[0][0]}")
print(f"G_SS[0,0] = {G_SS[0][0]}")

# Test PSD via LDL^T
print("\nTesting PSD via generalized LDL^T...", flush=True)

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

print("  Testing lam=0...", flush=True)
if not ok(arb(0)):
    print("  FAILED: Q_SS - 0*G_SS is not PSD!", flush=True)
    sys.exit(1)

print("  Bisecting to find certified lower bound...")
lo, hi = 0.0, 1.0
for iter_num in range(50):
    mid = (lo + hi) / 2
    if ok(arb(mid)):
        lo = mid
        if iter_num % 10 == 0:
            print(f"    iter {iter_num}: lo={lo:.8f}, hi={hi:.8f}", flush=True)
    else:
        hi = mid

print(f"\n✓ CERTIFIED: Q_SS >= {lo:.8f} at [14,44] 4W=2.50", flush=True)
print(f"  strong block size: {k_strong}")
print(f"  floor(AP): {floor:.6f}")
print(f"  max radius: {max_rad:.2e}")
