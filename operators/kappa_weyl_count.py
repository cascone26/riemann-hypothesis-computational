#!/usr/bin/env python3
"""§965 IS THE POSITIVE-RANK A WEYL COUNT? — a counting-identity handle on κ≤1 (follows §964).

§964 found: as bandwidth grows, #pos(arch−prime) SATURATES at 15 while dim grows and κ stays 1; the growth
goes entirely into the near-null plunge. If that 15 is a bandwidth-INDEPENDENT geometric invariant of the
window [a,b] — specifically the number of critical-line zeros the window subtends (Riemann–von Mangoldt),
doubled by the ±c symmetrization — then κ is pinned by COUNTING, not alignment:
     dim = #neg + #pos + #zero,   #pos = Weyl(window) (fixed),   #zero = plunge (grows),   ⇒ #neg = deficit.
That would be a NON-§891 route to κ≤1. TEST quantitatively:
  (a) Riemann–von Mangoldt main term N(T)=(T/2π)ln(T/2π)−T/2π; zeros in [a,b] = N(b)−N(a); ×2 for ±symmetry.
  (b) inertia of ARCH ALONE (no primes) across bandwidth — does arch's own #pos also saturate near the same
      number? (arch is the 'signal' operator; the prime subtraction is what creates the single negative.)
  (c) compare observed saturated #pos to the Weyl prediction.
HANDLE (own code) — report the numbers honestly; the identity either holds numerically or it doesn't."""
import numpy as np
from scipy.special import digamma

def vonmangoldt(N):
    lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam
NMAX = 8000; LAM = vonmangoldt(NMAX)
NS = np.arange(2, NMAX+1); LOGN = np.log(NS); LAMv = LAM[2:]

def RvM(T):  # Riemann–von Mangoldt main term: # zeros with 0<Im<T
    if T <= 0: return 0.0
    return (T/(2*np.pi))*np.log(T/(2*np.pi)) - T/(2*np.pi) + 7/8

def build(cen_lo, cen_hi, fourW, frac, dr=0.004):
    W = fourW/4; L = cen_hi-cen_lo; N = max(4, round(frac*L*2*W/np.pi))
    RG = np.arange(-(cen_hi+120), (cen_hi+120), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RG)) - np.log(np.pi)
    CEN = np.linspace(cen_lo, cen_hi, N)
    def fej(x):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out
    P = np.array([fej(RG-c)+fej(RG+c) for c in CEN]); G = (P @ P.T)*dr
    arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)
    mask = LOGN < fourW; ln = LOGN[mask]; wts = LAMv[mask]/np.sqrt(NS[mask])
    cosu = np.cos(np.outer(ln, RG)); prime = np.zeros((N, N))
    for j in range(N):
        gp = ((P[j][None, :]*P) @ cosu.T)*dr/(2*np.pi)
        prime[j, :] = 2.0*(gp*wts[None, :]).sum(axis=1)
    return arch, 0.5*(prime+prime.T), G, N

def inertia(M, G):
    gw, gv = np.linalg.eigh(G); keep = gw > gw.max()*1e-12
    B = gv[:, keep]/np.sqrt(gw[keep]); ev = np.sort(np.linalg.eigvalsh(B.T @ M @ B))
    s = max(1.0, np.abs(ev).max()); tol = 1e-9*s
    return int((ev < -tol).sum()), int((ev > tol).sum()), int((np.abs(ev) <= tol).sum())

a, b = 14, 44
zeros_win = RvM(b) - RvM(a)
print("="*100)
print(f"§965 Weyl-count handle on κ. window [{a},{b}].")
print(f"Riemann–von Mangoldt zeros in ({a},{b}): N({b})-N({a}) = {RvM(b):.3f} - {RvM(a):.3f} = {zeros_win:.3f}")
print(f"  ×2 for ±c symmetrization = {2*zeros_win:.3f}   (compare to §964 saturated #pos = 15)")
print("="*100)
print(f"{'4W':>5} {'N':>4} | {'ARCH alone':>22} | {'ARCH−PRIME':>22}")
print(f"{'':>5} {'':>4} | {'#neg':>5}{'#pos':>6}{'#zero':>7} | {'#neg':>5}{'#pos':>6}{'#zero':>7}")
print("-"*100)
for fourW in [3.0, 4.0, 5.0, 6.0, 7.0, 8.0]:
    arch, prime, G, N = build(a, b, fourW, 0.8)
    an, ap, az = inertia(arch, G)
    dn, dp, dz = inertia(arch-prime, G)
    print(f"{fourW:>5.1f} {N:>4} | {an:>5}{ap:>6}{az:>7} | {dn:>5}{dp:>6}{dz:>7}")
print("-"*100)
print(f"Weyl prediction (2·RvM zeros in window) = {2*zeros_win:.1f}.")
print("READ: if ARCH-alone #pos and ARCH−PRIME #pos both saturate near 2·RvM, the positive rank is a WINDOW")
print("GEOMETRIC INVARIANT (a Weyl zero-count), independent of bandwidth. Then κ=#neg is pinned by the counting")
print("identity dim=#neg+#pos+#zero with #pos fixed by Weyl — a route to κ≤1 that does NOT pass through §891")
print("alignment. If arch-alone already has #neg>0, or #pos doesn't track 2·RvM, the identity fails — honest either way.")
print("="*100)
