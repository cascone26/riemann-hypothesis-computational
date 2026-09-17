#!/usr/bin/env python3
"""§970 DE BRUIJN–NEWMAN FLOW on the Weil form — does the heat-flow dissolve the one negative square?  (HANDLE)

RH ⟺ Λ≤0 (de Bruijn–Newman); Rodgers–Tao Λ≥0, Polymath15 Λ≤0.22. H_t(z)=∫₀^∞ e^{tu²}Φ(u)cos(zu)du; zeros of
H_t all real ⟺ t≥Λ. In the explicit-formula/Weil picture the flow deforms the prime coefficients by the heat
weight b_n^t = exp((t/4)·log²n) (t>0 AMPLIFIES high n, t<0 DAMPS). Our arch−prime is Π₁ (κ=1) exactly at t=0.
QUESTION: track κ_t of the deformed form arch − Σ b_n^t·(Λ(n)/√n)·cos(r log n) across t (both signs). Which
direction drives κ:1→0 (the positivity margin opening), and is there a finite t* where the single negative
square dissolves? A clean κ_t curve connects our Π₁ picture to the dBN Λ program quantitatively. Prime-side
leading-order deformation (arch fixed) — honest about that scope. HANDLE (own code)."""
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

def build_pieces(cen_lo, cen_hi, fourW, frac=0.8, dr=0.005):
    """return arch, and per-term prime pieces {C_n} with base weight w_n and log n, so we can reweight by b_n^t."""
    W = fourW/4; L = cen_hi-cen_lo; N = max(4, round(frac*L*2*W/np.pi))
    RG = np.arange(-(cen_hi+40), (cen_hi+40), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RG)) - np.log(np.pi)
    CEN = np.linspace(cen_lo, cen_hi, N)
    def fej(x):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out
    P = np.array([fej(RG-c)+fej(RG+c) for c in CEN]); G = (P @ P.T)*dr
    arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)
    mask = LOGN < fourW; ln = LOGN[mask]; w = LAMv[mask]/np.sqrt(NS[mask])
    # precompute each prime term's matrix M_n = compression of cos(ln·r); prime(t)=Σ b_n^t w_n M_n
    cosu = np.cos(np.outer(ln, RG))
    Mns = np.empty((len(ln), N, N))
    for i in range(len(ln)):
        gp = ((P[None, :, :]*P[:, None, :]) @ cosu[i]) * dr/(2*np.pi)  # (N,N)
        Mns[i] = 2.0*gp
    return arch, Mns, w, ln, G, N

def kappa(M, G):
    gw, gv = np.linalg.eigh(G); keep = gw > gw.max()*1e-12
    B = gv[:, keep]/np.sqrt(gw[keep]); ev = np.linalg.eigvalsh(B.T @ M @ B)
    s = max(1.0, np.abs(ev).max()); tol = 1e-9*s
    return int((ev < -tol).sum()), float(ev.min())

arch, Mns, w, ln, G, N = build_pieces(14, 44, 6.0)
print("="*88)
print(f"§970 dBN heat-flow on the Weil form. window[14,44] 4W=6 N={N}, {len(ln)} prime terms.")
print("b_n^t = exp((t/4)log²n): t>0 amplifies high primes, t<0 damps. arch fixed (leading order).")
print("="*88)
print(f"{'t':>8} | {'κ_t=#neg':>9} {'λ_min':>13} | note")
print("-"*88)
for t in [-1.0, -0.5, -0.3, -0.2, -0.1, -0.05, 0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0]:
    bt = np.exp((t/4.0)*ln**2)
    prime_t = np.einsum('n,nij->ij', w*bt, Mns)
    M = arch - 0.5*(prime_t+prime_t.T)
    k, lmin = kappa(M, G)
    note = ""
    if k == 0: note = "<-- PSD! one negative square dissolved"
    print(f"{t:>+8.2f} | {k:>9} {lmin:>+13.3e} | {note}")
print("-"*88)
print("READ: the t-direction that drives κ:1→0 is where the flow opens the positivity margin. If κ→0 for t≥some")
print("t*≤0, that is consistent with RH (Λ≤0) — the form is critical at t=0 and positive just below. If κ→0 only")
print("for t<0 with the t>0 side AMPLIFYING (κ grows), that re-confirms the dBN obstruction with a number and")
print("localizes the critical t. If κ stays 1 both sides, the prime-only deformation is too crude — arch must")
print("co-deform. Honest either way; the κ_t curve is the deliverable.")
print("="*88)
