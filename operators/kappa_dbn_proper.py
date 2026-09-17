#!/usr/bin/env python3
"""§971 PROPERLY CO-DEFORMED dBN flow: Q_t(g,g)=Q(g_t,g_t), ĝ_t(u)=e^{(t/2)u²}ĝ(u).  (HANDLE)

§970 caveat (3): deforming only the prime side blows up (unbalanced amplification). The RIGHT object: the dBN
flow evaluates the SAME Weil form on the DEFORMED test function g_t whose Fourier transform is e^{(t/2)u²}ĝ(u)
— so arch, prime, AND Gram all co-deform consistently. Our Fejér test functions have COMPACT u-support (FT of
sinc² is a triangle on [−2W,2W]), so the Gaussian multiplier is bounded on the support → NO blow-up. Implement
via FFT along the r-grid: φ_j → φ_j^t = Re IFFT[e^{(t/2)u²} FFT[φ_j]], rebuild everything from {φ_j^t}, measure
κ_t. Sanity: t=0 must reproduce κ=1. Then the sign & the t* where κ:1→0 are BOTH physically meaningful (dBN
convention pinned by which direction is 'toward real zeros'). HANDLE (own code, honest)."""
import numpy as np
from scipy.special import digamma

def vonmangoldt(N):
    lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam
NMAX = 4000; LAM = vonmangoldt(NMAX)
NS = np.arange(2, NMAX+1); LOGN = np.log(NS); LAMv = LAM[2:]

# fixed grid & base test functions
cen_lo, cen_hi, fourW, frac, dr = 14, 44, 6.0, 0.8, 0.004
W = fourW/4; L = cen_hi-cen_lo; N = max(4, round(frac*L*2*W/np.pi))
RG = np.arange(-(cen_hi+40), (cen_hi+40), dr); nR = len(RG)
OMEGA = np.real(digamma(0.25 + 0.5j*RG)) - np.log(np.pi)
CEN = np.linspace(cen_lo, cen_hi, N)
def fej(x):
    y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
    out[nz] = (np.sin(y[nz])/y[nz])**2; return out
P0 = np.array([fej(RG-c)+fej(RG+c) for c in CEN])       # base test functions φ_j(r)
UU = 2*np.pi*np.fft.fftfreq(nR, d=dr)                    # conjugate variable u
mask = LOGN < fourW; ln = LOGN[mask]; wts = LAMv[mask]/np.sqrt(NS[mask])
COSU = np.cos(np.outer(ln, RG))
Ucut = 4*W + 1.0                                          # zero FT beyond compact support (+margin) to kill noise

def deform(P, t):
    if t == 0.0: return P
    F = np.fft.fft(P, axis=1)
    mult = np.exp((t/2.0)*UU**2); mult[np.abs(UU) > Ucut] = 0.0   # bounded on support; noise outside killed
    return np.real(np.fft.ifft(F*mult[None, :], axis=1))

def kappa_t(t):
    Pt = deform(P0, t)
    G = (Pt @ Pt.T)*dr
    arch = (Pt*OMEGA[None, :]) @ Pt.T * dr/(2*np.pi)
    prime = np.zeros((N, N))
    for j in range(N):
        gp = ((Pt[j][None, :]*Pt) @ COSU.T)*dr/(2*np.pi)
        prime[j, :] = 2.0*(gp*wts[None, :]).sum(axis=1)
    M = arch - 0.5*(prime+prime.T)
    gw, gv = np.linalg.eigh(G); keep = gw > gw.max()*1e-12
    B = gv[:, keep]/np.sqrt(gw[keep]); ev = np.linalg.eigvalsh(B.T @ M @ B)
    s = max(1.0, np.abs(ev).max()); tol = 1e-9*s
    return int((ev < -tol).sum()), float(ev.min()), int(keep.sum())

print("="*90)
print(f"§971 PROPERLY co-deformed dBN flow. window[{cen_lo},{cen_hi}] 4W={fourW} N={N}. ĝ_t=e^(t/2 u²)ĝ.")
print("t=0 sanity must give κ=1. Sign pinned by which direction drives zeros→real (κ→0).")
print("="*90)
print(f"{'t':>8} | {'κ_t':>4} {'λ_min':>13} {'dim':>4} | note")
print("-"*90)
for t in [-0.5, -0.3, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2, 0.3, 0.5]:
    k, lmin, dim = kappa_t(t)
    note = "SANITY" if t == 0.0 else ("<-- PSD (κ=0)" if k == 0 else ("κ jumped >1" if k > 1 else ""))
    print(f"{t:>+8.2f} | {k:>4} {lmin:>+13.3e} {dim:>4} | {note}")
print("-"*90)
print("READ: with EVERYTHING co-deformed there is no artificial blow-up. The side where κ:1→0 is the true")
print("dBN 'toward real zeros' direction; the t* localizes where the single negative square dissolves under the")
print("balanced flow. If t=0 gives κ=1 and κ→0 only for t≥t* with t*≤0, the form is critical at/just above the")
print("true zeros (Λ≈0 picture). A rigorous PSD certificate at a specific t would bound Λ≤t (Polymath15 target).")
print("="*90)
