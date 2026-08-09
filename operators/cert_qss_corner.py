#!/usr/bin/env python3
"""PATH D: Certify Q_SS >= 0 at [14,44] corner. Q_SS = arch|_S - prime|_S on strong-prime subspace.
Use float eigendecomposition to identify strong block S (prime eigenvalues > floor), then use
arb interval arithmetic to certify Q_SS is PSD (all LDL^T pivots > 0)."""

import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from scipy.special import digamma
from scipy.linalg import eigh
from flint import arb, acb, ctx

# ============= Build Weil forms at window [a,b] in floating point
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

    # Normalize by Gram matrix
    gw, gv = eigh(G)
    B = gv @ np.diag(1/np.sqrt(gw)) @ gv.T

    return B @ AP @ B, B @ prime @ B, G, B, arch, pole

# ============= Identify strong/robust split in float
def split_strong_robust(AP_norm, prime_norm):
    """Split on prime eigenvalues > floor(AP)"""
    mu, Wv = eigh(prime_norm)  # ascending mu
    floor = eigh(AP_norm, eigvals_only=True)[0]

    # Sort descending by eigenvalue
    idx = np.argsort(mu)[::-1]
    mu_sort = mu[idx]
    Wv_sort = Wv[:, idx]

    # Strong block: prime eigenvalues > floor
    k = int((mu_sort > floor).sum())
    return k, mu_sort, Wv_sort, floor

# ============= Q_SS construction in interval arithmetic (template)
import arb_fejer_cert2 as ac
import arb_zerofree_cert3 as zc

def arch_entry_arb(cj, ck, rc):
    """Archimedean entry (pole + arch) at high precision"""
    return ac.pole_entry(cj, ck) + zc.arch_entry_v3(cj, ck, rc)

def cert_qss(a, b, fourW, CEN_float, k_strong, Wv_sort, prime_eigs, rc, prec=200):
    """Certify Q_SS >= 0 using arb.
    Inputs:
      CEN_float: float basis centers
      k_strong: size of strong block
      Wv_sort: prime eigenvectors (float, sorted desc by eigenvalue)
      prime_eigs: prime eigenvalues (float, sorted desc)
      rc: integration radius
      prec: arb precision
    Returns:
      (certified_min_eigenvalue, max_radius, is_certified)
    """
    ctx.prec = prec
    n = len(CEN_float)

    # Set up basis in arb
    ac.W = arb(fourW) / 4
    ac.TW = 2 * ac.W
    ac.c_tri = arb.pi() / ac.W
    cs = [arb(c) for c in CEN_float]

    # Extract strong-block eigenvectors (first k rows of Wv_sort^T)
    W_S_float = Wv_sort[:, :k_strong]  # shape (n, k_strong)

    # Build Q_SS = W_S^T (arch + pole - prime) W_S in arb
    Q_SS = [[None]*k_strong for _ in range(k_strong)]
    G_SS = [[None]*k_strong for _ in range(k_strong)]
    max_rad = 0.0

    # Build full Q in arb first (only the entries we need)
    Q_arb = [[None]*n for _ in range(n)]
    G_arb = [[None]*n for _ in range(n)]

    for jj in range(n):
        for kk in range(jj, n):
            q_entry = arch_entry_arb(cs[jj], cs[kk], arb(rc)) - ac.prime_entry(cs[jj], cs[kk])
            Q_arb[jj][kk] = Q_arb[kk][jj] = q_entry
            g_entry = 2 * arb.pi() * ac.gpair_closed_arb(cs[jj], cs[kk], arb(0))
            G_arb[jj][kk] = G_arb[kk][jj] = g_entry
            max_rad = max(max_rad, float(q_entry.rad()))

    # Project to strong subspace: Q_SS = W_S^T Q W_S (in arb)
    for i in range(k_strong):
        for j in range(k_strong):
            s = arb(0)
            for ii in range(n):
                for jj in range(n):
                    s += arb(W_S_float[ii, i]) * Q_arb[ii][jj] * arb(W_S_float[jj, j])
            Q_SS[i][j] = s

    # G_SS = W_S^T G W_S
    for i in range(k_strong):
        for j in range(k_strong):
            s = arb(0)
            for ii in range(n):
                for jj in range(n):
                    s += arb(W_S_float[ii, i]) * G_arb[ii][jj] * arb(W_S_float[jj, j])
            G_SS[i][j] = s

    # Test PSD via LDL^T: compute min lam s.t. Q_SS - lam*G_SS is PSD
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

    # Bisection to find certified lower bound
    if not ok(arb(0)):
        return None, max_rad, False

    lo, hi = 0.0, 1.0
    for _ in range(50):
        mid = (lo + hi) / 2
        if ok(arb(mid)):
            lo = mid
        else:
            hi = mid

    return lo, max_rad, True

# ============= Main sweep over 4W values at [14,44]
a, b = 14, 44
rc_base = 180
results = []

print("=" * 80)
print(f"PATH D: Certifying Q_SS >= 0 at corner [{a},{b}]")
print("=" * 80)

for fourW_frac in [2.5, 2.6, 2.7, 2.75]:
    W = fourW_frac / 4
    N = max(6, int(np.ceil((b-a) / 5)))  # Adjust N with 4W

    # Build float forms
    AP_norm, prime_norm, G, B, arch, pole = build_float(a, b, W, N, dr=0.004)

    # Identify strong block
    k_strong, mu_sort, Wv_sort, floor = split_strong_robust(AP_norm, prime_norm)

    # Basis centers (float)
    CEN_float = list(np.linspace(a, b, N))

    # Try certification
    rc = rc_base + int(fourW_frac * 40)  # Scale integration radius with 4W
    cert_lam, max_rad, is_cert = cert_qss(a, b, fourW_frac, CEN_float, k_strong,
                                           Wv_sort, mu_sort[:k_strong], rc, prec=200)

    result_str = f"[{a},{b}] 4W={fourW_frac:.2f} N={N}: "
    result_str += f"k={k_strong}, floor={floor:.4f}, "
    if is_cert:
        result_str += f"CERTIFIED Q_SS >= {cert_lam:.8f}  (maxrad {max_rad:.2e})"
    else:
        result_str += f"FAILED (maxrad {max_rad:.2e} exceeds margin)"

    print(result_str)
    results.append({
        'fourW': fourW_frac,
        'N': N,
        'k_strong': k_strong,
        'floor': floor,
        'certified_min': cert_lam if is_cert else None,
        'max_radius': max_rad,
        'is_certified': is_cert
    })

print("=" * 80)
print("\nSummary:")
for r in results:
    status = "✓ CERTIFIED" if r['is_certified'] else "✗ FAILED"
    cert_val = f"{r['certified_min']:.8f}" if r['is_certified'] else "N/A"
    print(f"4W={r['fourW']:.2f}: {status}  Q_SS_min >= {cert_val}  (rad={r['max_radius']:.2e})")

# Find largest 4W certified
certified = [r for r in results if r['is_certified']]
if certified:
    best = max(certified, key=lambda x: x['fourW'])
    print(f"\n>>> LARGEST 4W CERTIFIED: {best['fourW']:.2f} with Q_SS_min >= {best['certified_min']:.8f}")
