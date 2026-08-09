#!/usr/bin/env python3
"""§853: is the corner razor-thinness INTRINSIC or a Fejér-BASIS artifact? Compare the Weil margin in a
GAUSSIAN test-function basis (g_c(r)=exp(-(r-c)^2/2s^2), the CCM Hermite/Gaussian family) vs the Fejér
basis, at the same window + matched effective bandwidth. If Gaussian margin is dramatically MORE robust
where Fejér is razor-thin → razor is basis-dependent (→ CCM Hermite is the 'right' basis). If equally
razor-thin → intrinsic. APPROXIMATE: compares the dominant arch−prime part (pole dropped, small ~1e-3);
Gaussian prime sum is Gaussian-truncated in log n (effective cutoff ~2/s), summed generously to n<400.
Matched bandwidth: Gaussian s = 1/(2W) so its u-decay ~ matches Fejér's 4W u-support."""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh

def vonmangoldt(N):
    lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam
NMAX = 400; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1)
LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)

def margins(a, b, W, N, dr=0.004):
    fourW = 4*W
    RGRID = np.arange(-(b+140), (b+140), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
    CEN = np.linspace(a, b, N)
    # ---- Fejér basis ----
    def fej(x):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out
    F = np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN])
    # ---- Gaussian basis (matched bandwidth s=1/(2W)) ----
    s = 1.0/(2*W)
    Gs = np.array([np.exp(-(RGRID-c)**2/(2*s**2)) + np.exp(-(RGRID+c)**2/(2*s**2)) for c in CEN])
    def wmargin(P, fourW_cut, gaussian):
        Gm = (P @ P.T)*dr
        arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)
        # prime cutoff: Fejér hard log n<4W; Gaussian generous (Gaussian-truncated), use n<400
        if gaussian: mask = np.ones_like(LOGN, dtype=bool)
        else: mask = LOGN < fourW_cut
        ln = LOGN[mask]; lw = LAMW[mask]
        cosu = np.cos(np.outer(ln, RGRID)); n = P.shape[0]; prime = np.zeros((n, n))
        for j in range(n):
            gp = ((P[j][None, :]*P) @ cosu.T)*dr/(2*np.pi)
            prime[j, :] = 2.0*(gp*lw[None, :]).sum(axis=1)
        prime = 0.5*(prime+prime.T)
        gw = np.linalg.eigvalsh(Gm); cond = float(gw[-1]/max(gw[0], 1e-300))
        mg = float(eigh(arch-prime, Gm, eigvals_only=True)[0])
        return mg, cond
    mF, cF = wmargin(F, fourW, False)
    mG, cG = wmargin(Gs, fourW, True)
    return mF, cF, mG, cG, s

print("="*100)
print("§853 Gaussian vs Fejér basis — is the corner razor-thinness intrinsic or basis-artifact? (arch−prime, approx)")
print("="*100)
print(f"{'window':>12} {'4W':>4} {'N':>3} {'s(gauss)':>8} {'margin_Fejer':>13} {'condF':>7} {'margin_Gauss':>13} {'condG':>7}")
for (a, b, fourW) in [(14,44,3.0),(14,44,4.0),(60,90,4.0),(120,150,4.0)]:
    W = fourW/4; L = b-a; N = max(4, round(0.75*L*2*W/np.pi))
    mF, cF, mG, cG, s = margins(a, b, W, N)
    print(f"{f'[{a},{b}]':>12} {fourW:>4.1f} {N:>3} {s:>8.3f} {mF:>+13.4e} {cF:>7.1f} {mG:>+13.4e} {cG:>7.1f}")
print("="*100)
print("Gaussian margin ≫ Fejér margin (robust where Fejér razor-thin) → razor is BASIS-DEPENDENT (CCM Hermite lead).")
print("Gaussian margin ~ Fejér (also razor-thin/→0) → the razor is INTRINSIC to Weil positivity.")
print("(Approximate: pole dropped; Gaussian prime Gaussian-truncated at n<400. Qualitative comparison.)")
print("="*100)
