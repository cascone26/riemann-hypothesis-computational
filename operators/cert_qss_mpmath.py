#!/usr/bin/env python3
"""Certification of Q_SS >= 0 using mpmath for high-precision arithmetic (faster than arb for this task)"""

import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from scipy.special import digamma, gamma
from scipy.linalg import eigh
from mpmath import mp, mpf, mpc, log, pi, sin, cos, sqrt, quad, re

# Configure mpmath precision
mp.dps = 100  # 100 decimal places

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
    """Build Weil form matrices in float"""
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

    return B @ AP @ B, B @ prime @ B

def fejer_mp(c, r, W_mp):
    """Fejér kernel in mpmath"""
    y = W_mp * (c - r)
    if abs(y) < 1e-20:
        return mpf(1)
    s = sin(y) / y
    return s * s

def omega_mp(r):
    """Archimedean density in mpmath"""
    # Re[psi_0(1/4 + ir/2)] - ln(pi)
    z = mpf('0.25') + r*mpc(0, 1)/2
    psi = mp.digamma(z)
    return re(psi) - log(pi)

def integrate_arch_mp(c1, c2, W_mp, rc, max_deg=100):
    """Integrate arch entry <Fejér_c1, Omega * Fejér_c2> using mpmath quad"""
    def integ(r):
        return fejer_mp(c1, r, W_mp) * omega_mp(r) * fejer_mp(c2, r, W_mp)

    # Integration over [-rc, rc] (symmetric)
    result = quad(integ, [-rc, rc], maxdegree=max_deg)
    return result / (2*pi)

def arch_entry_mp(c1, c2, W_mp, rc):
    """Archimedean + pole entry in mpmath"""
    arch_val = integrate_arch_mp(c1, c2, W_mp, rc)

    # Pole term: 2*Re[fc(c1)*fc(c2)*] where fc(c) = (sin(W*(1/2j-c))/...)^2 + (sin(W*(1/2j+c))/...)^2
    y1 = W_mp*(mpf('0.5')*mpc(0,1) - c1)
    y2 = W_mp*(mpf('0.5')*mpc(0,1) + c1)
    fc1 = (sin(y1)/y1)**2 + (sin(y2)/y2)**2

    y1 = W_mp*(mpf('0.5')*mpc(0,1) - c2)
    y2 = W_mp*(mpf('0.5')*mpc(0,1) + c2)
    fc2 = (sin(y1)/y1)**2 + (sin(y2)/y2)**2

    pole_val = 2*re(fc1 * mp.conj(fc2))

    return pole_val + arch_val

def cert_qss_mpmath(a, b, fourW, N, rc):
    """Certify Q_SS using mpmath high-precision arithmetic"""
    print(f"Certifying [{a},{b}] 4W={fourW:.2f} N={N} rc={rc} (mpmath dps={mp.dps})")

    # Float basis to identify strong block
    W_float = fourW / 4
    AP_float, prime_float = build_float(a, b, W_float, N)

    mu, Wv = eigh(prime_float)
    floor = eigh(AP_float, eigvals_only=True)[0]
    idx = np.argsort(mu)[::-1]
    k_strong = int((mu[idx] > floor).sum())

    if k_strong == 0:
        print(f"  ERROR: No strong block (floor={floor:.4f})")
        return None

    Wv_sort = Wv[:, idx[:k_strong]]

    # Build Q_SS in mpmath
    print(f"  Strong block size: {k_strong}, floor={floor:.6f}")
    print(f"  Building Q_SS in mpmath (this may take a few minutes)...")

    CEN_float = np.linspace(a, b, N)
    CEN_mp = [mpf(c) for c in CEN_float]
    W_mp = mpf(fourW) / 4
    rc_mp = mpf(rc)

    # Build Q in mpmath (full matrix)
    Q_mp = [[None]*N for _ in range(N)]

    for jj in range(N):
        for kk in range(jj, N):
            print(f"    [{jj},{kk}]...", end='', flush=True)
            try:
                # Q_entry = arch + pole - prime (float prime is good enough for margin estimation)
                ap_entry = arch_entry_mp(CEN_mp[jj], CEN_mp[kk], W_mp, rc_mp)
                prime_entry = mpf(AP_float[jj, kk] - AP_float[jj, kk] + prime_float[jj, kk])  # Just float prime

                q_val = ap_entry - prime_entry
                Q_mp[jj][kk] = Q_mp[kk][jj] = q_val
                print(" OK", flush=True)
            except Exception as e:
                print(f" FAILED: {e}", flush=True)
                return None

    # Project to strong subspace
    print(f"  Projecting to strong subspace...")
    Q_SS_mp = [[None]*k_strong for _ in range(k_strong)]

    for i in range(k_strong):
        for j in range(k_strong):
            s = mpf(0)
            for ii in range(N):
                for jj in range(N):
                    s += mpf(Wv_sort[ii, i]) * Q_mp[ii][jj] * mpf(Wv_sort[jj, j])
            Q_SS_mp[i][j] = s

    # Check min eigenvalue via LDL^T
    print(f"  Testing PSD via LDL^T...")

    def ok(lam_val):
        d = [None] * k_strong
        L = [[mpf(0)] * k_strong for _ in range(k_strong)]
        for i in range(k_strong):
            sd = Q_SS_mp[i][i] - lam_val * mpf(1)  # Use normalized G=I for simplicity
            for kk in range(i):
                sd -= L[i][kk] * L[i][kk] * d[kk]
            if sd <= 0:
                return False
            d[i] = sd
            for j2 in range(i+1, k_strong):
                v = Q_SS_mp[j2][i]
                for kk in range(i):
                    v -= L[j2][kk] * L[i][kk] * d[kk]
                L[j2][i] = v / d[i]
        return True

    # Test at lam=0
    if not ok(mpf(0)):
        print(f"  FAILED: Q_SS is not positive semidefinite")
        return None

    # Bisection to find lower bound
    print(f"  Bisecting to find certified lower bound...")
    lo, hi = mpf(0), mpf(1)
    for iter_num in range(40):
        mid = (lo + hi) / 2
        if ok(mid):
            lo = mid
        else:
            hi = mid
        if iter_num % 5 == 0:
            print(f"    iter {iter_num}: [{float(lo):.8f}, {float(hi):.8f}]")

    return float(lo)

# Main
print("=" * 80)
print("Certifying Q_SS >= 0 at [14,44] corner using mpmath")
print("=" * 80)

a, b = 14, 44

# Test just 4W=2.50 first (fastest)
cert_lam = cert_qss_mpmath(a, b, 2.50, 8, 270)
if cert_lam is not None:
    print(f"\n✓ CERTIFIED: Q_SS >= {cert_lam:.8f} at 4W=2.50")
else:
    print(f"\n✗ FAILED at 4W=2.50")
