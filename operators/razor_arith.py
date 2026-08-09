#!/usr/bin/env python3
"""ANGLE C: Does the RH RAZOR track ζ zeros arithmetically?
§900 KEY FINDING: π_w ≈ μ with ~4% slack (π_w/μ ≈ 1.04), but the ACTUAL MARGIN (razor)
has ~1e-15 tightness due to coupling of negative direction w to M_σ's near-null POSITIVE
cluster through A^+. Test whether μ (defect depth), the true margin, or w^T A^+ w follow
a clean arithmetic law vs window's ζ zero spacings/count/proximity to boundary.
Windows: [14,44],[14,34],[20,50],[30,60],[40,70] at fixed 4W=4. Full ζ via mpmath.
Honest: report best R² even if not clean. Label REPORT vs HANDLE."""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh, pinvh
import mpmath as mp

def vonmangoldt(N):
    lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam

NMAX = 4000; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1)
LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)

# ζ zeros via mpmath (HANDLE: real computation)
mp.mp.dps = 20
ZT = []
k = 1
while True:
    g = float(mp.im(mp.zetazero(k)))
    ZT.append(g)
    if g > 200:  # enough for our test windows
        break
    k += 1
ZT = np.array(ZT)

def zeros_in(a, b):
    return ZT[(ZT >= a) & (ZT <= b)]

def run_razor(a, b, fourW, frac=0.65, dr=0.004):
    """Compute μ, margin, w^T A^+ w at one corner."""
    W = fourW/4
    L = b - a
    Ncap = L * 2 * W / np.pi
    N = max(4, int(round(frac * Ncap)))

    RGRID = np.arange(-(b+60), (b+60), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
    CEN = np.linspace(a, b, N)

    def fej(x):
        y = W*x
        out = np.ones_like(y)
        nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2
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
        gp = ((P[j][None, :]*P) @ cosu.T) * dr / (2*np.pi)
        prime[j, :] = 2.0 * (gp * lw[None, :]).sum(axis=1)
    prime = 0.5 * (prime + prime.T)

    def fc(c):
        y = W*(0.5j - c)
        y2 = W*(0.5j + c)
        return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2
    fp = np.array([fc(c) for c in CEN])
    pole = 2.0 * np.real(np.outer(fp, fp))

    # Whitening
    gw, gv = eigh(G)
    B = gv @ np.diag(1/np.sqrt(gw)) @ gv.T
    Msw = B @ (arch - prime) @ B
    polew = B @ pole @ B

    # Eigendecomposition: negative eigenvalue -mu, eigenvector w
    ev, V = eigh(Msw)
    nneg = int((ev < 0).sum())
    mu = -float(ev[0])
    w = V[:, 0]

    # PSD part of M_σ
    Mplus = V[:, 1:] @ np.diag(ev[1:]) @ V[:, 1:].T

    # Pole rank-1 PSD
    pev, pvec = eigh(polew)
    p = float(pev[-1])
    v = pvec[:, -1]

    # A = M_+ + pole (PSD part + pole)
    A = Mplus + polew

    # True margin: smallest eigenvalue of Q = M_σ + pole
    Q = Msw + polew
    margin = float(eigh(Q, eigvals_only=True)[0])

    # w^T A^+ w (pseudoinverse overlap)
    Ap = pinvh(A)
    w_overlap = float(w @ Ap @ w)

    # Leading criterion: π_w = p·(v·w)²
    pi_w = p * float(v @ w)**2

    return N, nneg, mu, margin, w_overlap, pi_w, p, float(v @ w)**2

print("="*140)
print("ANGLE C: RAZOR TRACKING ζ ZEROS ARITHMETICALLY?")
print("="*140)
print(f"{'window':>12} {'4W':>4} {'#ζ':>3} {'N':>3} {'μ':>11} {'margin':>11} {'w^T A^+ w':>10} {'π_w':>11} {'arch floor':>11} {'min_gap':>8} {'mean_gap':>8}")
print("-"*140)

rows = []
for (a, b, fourW) in [(14, 44, 4.0), (14, 34, 4.0), (20, 50, 4.0), (30, 60, 4.0), (40, 70, 4.0)]:
    W = fourW / 4
    zs = zeros_in(a, b)
    nz = len(zs)

    if nz < 2:
        continue

    gaps = np.diff(np.sort(zs))
    min_gap = float(gaps.min())
    mean_gap = float(gaps.mean())

    # Arch floor: Ω ~ ln(t/2π) at window midpoint
    t_mid = (a + b) / 2
    arch_floor = np.log(t_mid / (2 * np.pi))

    N, nneg, mu, margin, w_overlap, pi_w, p, vw2 = run_razor(a, b, fourW)

    rows.append({
        'a': a, 'b': b, 'fourW': fourW,
        'nz': nz, 'N': N, 'nneg': nneg,
        'mu': mu, 'margin': margin, 'w_overlap': w_overlap,
        'pi_w': pi_w, 'p': p, 'vw2': vw2,
        'arch_floor': arch_floor,
        'min_gap': min_gap, 'mean_gap': mean_gap
    })

    print(f"{f'[{a},{b}]':>12} {fourW:>4.1f} {nz:>3} {N:>3} {mu:>+11.3e} {margin:>+11.3e} {w_overlap:>+10.3e} {pi_w:>+11.3e} {arch_floor:>+11.3f} {min_gap:>8.3f} {mean_gap:>8.3f}")

print("="*140)
print("\nARITHMETIC ANALYSIS: Does the RAZOR obey a clean law?")
print("="*140)

if len(rows) >= 3:
    # Extract quantities for fitting
    mu_vals = np.array([r['mu'] for r in rows])
    margin_vals = np.array([r['margin'] for r in rows])
    w_overlap_vals = np.array([r['w_overlap'] for r in rows])
    min_gap_vals = np.array([r['min_gap'] for r in rows])
    mean_gap_vals = np.array([r['mean_gap'] for r in rows])
    nz_vals = np.array([r['nz'] for r in rows])
    arch_floor_vals = np.array([r['arch_floor'] for r in rows])

    # Take absolute values for log-log (handle signs carefully)
    mu_abs = np.abs(mu_vals)
    margin_abs = np.abs(margin_vals)
    w_overlap_abs = np.abs(w_overlap_vals)

    print("\n[Test 1] Does μ (defect depth) track min_gap?")
    if (mu_abs > 1e-15).sum() >= 3:
        pos = mu_abs > 1e-15
        if pos.sum() >= 3:
            x_log = np.log(min_gap_vals[pos])
            y_log = np.log(mu_abs[pos])
            coeffs = np.polyfit(x_log, y_log, 1)
            fitted = coeffs[0] * x_log + coeffs[1]
            ss_res = np.sum((y_log - fitted)**2)
            ss_tot = np.sum((y_log - y_log.mean())**2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
            print(f"  log(μ) ~ {coeffs[0]:>6.2f}·log(min_gap) + {coeffs[1]:>7.2f}   R² = {r2:>7.4f}  {'✓ CLEAN' if r2>0.90 else '✗ noisy'}")

    print("\n[Test 2] Does the TRUE MARGIN track min_gap?")
    # margin can be negative; take absolute value for log-log
    if (margin_abs > 1e-15).sum() >= 3:
        pos = margin_abs > 1e-15
        if pos.sum() >= 3:
            x_log = np.log(min_gap_vals[pos])
            y_log = np.log(margin_abs[pos])
            coeffs = np.polyfit(x_log, y_log, 1)
            fitted = coeffs[0] * x_log + coeffs[1]
            ss_res = np.sum((y_log - fitted)**2)
            ss_tot = np.sum((y_log - y_log.mean())**2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
            print(f"  log(|margin|) ~ {coeffs[0]:>6.2f}·log(min_gap) + {coeffs[1]:>7.2f}   R² = {r2:>7.4f}  {'✓ CLEAN' if r2>0.90 else '✗ noisy'}")

    print("\n[Test 3] Does w^T A^+ w (pseudoinverse overlap) track min_gap?")
    if (w_overlap_abs > 1e-15).sum() >= 3:
        pos = w_overlap_abs > 1e-15
        if pos.sum() >= 3:
            x_log = np.log(min_gap_vals[pos])
            y_log = np.log(w_overlap_abs[pos])
            coeffs = np.polyfit(x_log, y_log, 1)
            fitted = coeffs[0] * x_log + coeffs[1]
            ss_res = np.sum((y_log - fitted)**2)
            ss_tot = np.sum((y_log - y_log.mean())**2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
            print(f"  log(w^T A^+ w) ~ {coeffs[0]:>6.2f}·log(min_gap) + {coeffs[1]:>7.2f}   R² = {r2:>7.4f}  {'✓ CLEAN' if r2>0.90 else '✗ noisy'}")

    print("\n[Test 4] Does μ track mean_gap?")
    if (mu_abs > 1e-15).sum() >= 3:
        pos = mu_abs > 1e-15
        if pos.sum() >= 3:
            x_log = np.log(mean_gap_vals[pos])
            y_log = np.log(mu_abs[pos])
            coeffs = np.polyfit(x_log, y_log, 1)
            fitted = coeffs[0] * x_log + coeffs[1]
            ss_res = np.sum((y_log - fitted)**2)
            ss_tot = np.sum((y_log - y_log.mean())**2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
            print(f"  log(μ) ~ {coeffs[0]:>6.2f}·log(mean_gap) + {coeffs[1]:>7.2f}   R² = {r2:>7.4f}  {'✓ CLEAN' if r2>0.90 else '✗ noisy'}")

    print("\n[Test 5] Does margin track #zeros?")
    if (margin_abs > 1e-15).sum() >= 3:
        pos = margin_abs > 1e-15
        if pos.sum() >= 3:
            x_log = np.log(nz_vals[pos])
            y_log = np.log(margin_abs[pos])
            coeffs = np.polyfit(x_log, y_log, 1)
            fitted = coeffs[0] * x_log + coeffs[1]
            ss_res = np.sum((y_log - fitted)**2)
            ss_tot = np.sum((y_log - y_log.mean())**2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
            print(f"  log(|margin|) ~ {coeffs[0]:>6.2f}·log(#zeros) + {coeffs[1]:>7.2f}   R² = {r2:>7.4f}  {'✓ CLEAN' if r2>0.90 else '✗ noisy'}")

    print("\n[Test 6] Does μ track arch_floor (Ω ~ ln(t/2π))?")
    if (mu_abs > 1e-15).sum() >= 3:
        pos = mu_abs > 1e-15
        if pos.sum() >= 3:
            x_log = np.log(arch_floor_vals[pos])
            y_log = np.log(mu_abs[pos])
            coeffs = np.polyfit(x_log, y_log, 1)
            fitted = coeffs[0] * x_log + coeffs[1]
            ss_res = np.sum((y_log - fitted)**2)
            ss_tot = np.sum((y_log - y_log.mean())**2)
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
            print(f"  log(μ) ~ {coeffs[0]:>6.2f}·log(arch_floor) + {coeffs[1]:>7.2f}   R² = {r2:>7.4f}  {'✓ CLEAN' if r2>0.90 else '✗ noisy'}")

print("\n" + "="*140)
print("SUMMARY: HANDLE findings (computed on real ζ via mpmath) listed above.")
print("If any R² > 0.90 → RH-razor has an ARITHMETIC LAW at this scale.")
print("If all R² < 0.80 → the razor does NOT follow a simple monotone law in zero spacings.")
print("="*140)
