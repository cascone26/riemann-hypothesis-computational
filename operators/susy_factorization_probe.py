#!/usr/bin/env python3
"""ATTACK A12 probe: does tau admit a STRUCTURED (SUSY-like) factorization tau=A^dag A
where A is a near-local / first-order operator? If the matrix square root tau^{1/2} (a valid
symmetric A) is strongly BANDED and/or TOEPLITZ (convolution => a Fourier multiplier =
manifest square-root of the symbol), that's a real SUSY structure. If dense/structureless,
SUSY-in-this-basis is out. Also: how far is tau from TRIDIAGONAL (a true first-order A^dag A
would be tridiagonal)? And does adding the PRIME term make the square root more nonlocal?"""
import numpy as np, sys
sys.path.insert(0,"operators")
from critical_eigvec_id import build_tau_np, prime_powers_up_to

def sqrtm_psd(M):
    w,V=np.linalg.eigh(M); w=np.maximum(w,0)
    return V@np.diag(np.sqrt(w))@V.T

def band_profile(M):
    D=len(M); prof=[]
    for d in range(0,min(D,25)):
        vals=[abs(M[n,n+d]) for n in range(D-d)]
        prof.append(np.mean(vals))
    return np.array(prof)

def toeplitz_score(M):
    """for each offset d, std/mean of M[n,n+d] across n. ~0 => Toeplitz (convolution)."""
    D=len(M); sc=[]
    for d in range(1,min(D,15)):
        vals=np.array([M[n,n+d] for n in range(D-d)])
        m=np.mean(np.abs(vals))
        sc.append(np.std(vals)/m if m>1e-12 else np.nan)
    return np.array(sc)

def tridiag_fraction(M):
    D=len(M); tot=np.linalg.norm(M)
    tri=M.copy()
    for i in range(D):
        for j in range(D):
            if abs(i-j)>1: tri[i,j]=0
    return np.linalg.norm(tri)/tot

for c in [14]:
  for N in [48,64]:
    tau,ns=build_tau_np(N,c)
    # arch+boundary only (no prime) for comparison: rebuild with empty prime? easier: P = full - (A+B)
    # reconstruct A+B by zeroing the prime contribution is hard here; instead compare tau vs its own sqrt.
    St=sqrtm_psd(tau)
    prof=band_profile(np.abs(St)); prof=prof/prof[0]
    toe=toeplitz_score(St)
    trf=tridiag_fraction(tau)
    trf_sqrt=tridiag_fraction(St)
    # how much of tau^{1/2} energy is within band width b?
    D=len(St); frob=np.linalg.norm(St)
    def band_energy(b):
        m=St.copy()
        for i in range(D):
            for j in range(D):
                if abs(i-j)>b: m[i,j]=0
        return np.linalg.norm(m)/frob
    print(f"\n=== c={c} N={N} (D={D}) ===")
    print(f"  tau tridiagonal fraction: {trf:.3f}  (1.0=first-order A^dagA possible; <1 => not)")
    print(f"  tau^1/2 tridiagonal fraction: {trf_sqrt:.3f}")
    print(f"  tau^1/2 band-energy: within |d|<=1: {band_energy(1):.3f}  <=3: {band_energy(3):.3f}  <=7: {band_energy(7):.3f}")
    print(f"  tau^1/2 |offdiag| decay (d=0..8, norm to d0): {[f'{x:.3f}' for x in prof[:9]]}")
    print(f"  tau^1/2 Toeplitz score (d=1..6; ~0=convolution): {[f'{x:.2f}' for x in toe[:6]]}")
    # superpotential read: if A were bidiagonal L (Cholesky of tau+eps), is L banded?
    L=np.linalg.cholesky(tau+1e-12*np.eye(D))
    Lband=band_profile(np.abs(L)); Lband=Lband/Lband[0]
    print(f"  Cholesky L |offdiag| decay (d=0..8): {[f'{x:.3f}' for x in Lband[:9]]}")
