#!/usr/bin/env python3
"""§447: does the FE-closed Rosati (§446) detection SHARPEN as M (matrix size) and K (#zeros) grow? Find β*(M,K)
= smallest off-line displacement making the form indefinite. If β*→0 as M,K grow, the form sharply detects RH in
the limit (correct behavior for a truncation of the true infinite Rosati)."""
import numpy as np
G=np.load('zeros150.npy')
def rosati_fe(beta,M,K):
    g=G[:K]; rho=0.5+beta+1j*g; fe=1-rho; r=np.arange(M); A=np.zeros((M,M),dtype=complex)
    for k in range(K):
        er=np.exp(rho[k]*r); ef=np.exp(fe[k]*r); A+=np.outer(er,ef)+np.outer(ef,er)
    return np.real((A+A.conj().T)/2)
def bstar(M,K):
    if np.linalg.eigvalsh(rosati_fe(1.0,M,K))[0]>=0: return np.inf
    lo,hi=0.0,1.0
    for _ in range(26):
        mid=(lo+hi)/2
        if np.linalg.eigvalsh(rosati_fe(mid,M,K))[0]>=-1e-7*max(1,abs(np.linalg.eigvalsh(rosati_fe(mid,M,K))).max()): lo=mid
        else: hi=mid
    return (lo+hi)/2
print("§447 β*(M,K): off-line displacement threshold for FE-closed Rosati to go indefinite. β*→0 = sharp detection.")
print(f"{'M':>4}" + "".join(f"  K={K:<6}" for K in (20,40,80,150)))
for M in (8,10,12,16,20,24):
    row=[bstar(M,K) for K in (20,40,80,150)]
    print(f"{M:>4}" + "".join(f"  {b:8.4f}" if np.isfinite(b) else f"  {'inf':>8}" for b in row),flush=True)
print("""
reading: β* DECREASING as M and/or K grow => the FE-closed Rosati sharpens toward detecting ANY off-line zero (β*→0),
i.e. RH ⟺ the (limiting) form is PSD — the correct geometric-positivity behavior for ℚ, just weakly-conditioned at
finite truncation. β* flat/growing => the finite Rosati can't sharpen (the archimedean Γ-factor is needed to get
sharp detection). Either way quantifies the number-field detection strength.
""")
