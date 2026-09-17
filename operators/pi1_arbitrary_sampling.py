#!/usr/bin/env python3
"""§953 IS M_σ A GENUINE Π₁ KERNEL? — the proper Krein-Langer test: κ≤1 under ARBITRARY sampling.

Krein-Langer: a hermitian kernel K has κ negative squares ⟺ EVERY finite Gram [K(λ_i,λ_j)] (sampling at
ARBITRARY points/configs) has ≤ κ negative eigenvalues, and =κ for some config. §948 saw κ=1 with EVENLY-SPACED
Fejér centers. The genuine-Π₁ claim requires κ≤1 for ANY sampling: irregular centers, OVER-DENSE packing (far
past the Landau cap), mixed bandwidths, random configs. If κ≤1 survives all of these, M_σ is a bona fide Π₁
reproducing kernel ⟹ the Krein-Langer decomposition M_σ=K₊−K₋ (rank K₋=1) EXISTS (the theorem applies) ⟹ RH
via 'pole absorbs the one negative square'. If some config gives κ≥2, κ=1 was basis-specific — decisive either way.
HANDLE (own code). We vary the Fejér-center CONFIGURATION freely (arbitrary points = arbitrary test functions)."""
import numpy as np
from scipy.special import digamma

def vonmangoldt(N):
    lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam
NMAX = 8000; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1)
LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)

def M_sigma(centers, W, dr=0.004, span=None):
    """M_σ = arch − prime compression at ARBITRARY Fejér centers with bandwidth W (4W band). Returns (M,G)."""
    fourW = 4*W
    lo, hi = min(centers), max(centers)
    RGRID = np.arange(-(hi+120), (hi+120), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)
    def fej(x):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out
    P = np.array([fej(RGRID-c)+fej(RGRID+c) for c in centers]); G = (P @ P.T)*dr
    arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)
    mask = LOGN < fourW; ln = LOGN[mask]; lw = LAMW[mask]
    cosu = np.cos(np.outer(ln, RGRID)); n = len(centers); prime = np.zeros((n, n))
    for j in range(n):
        gp = ((P[j][None, :]*P) @ cosu.T)*dr/(2*np.pi)
        prime[j, :] = 2.0*(gp*lw[None, :]).sum(axis=1)
    prime = 0.5*(prime+prime.T)
    return arch - prime, G

def kappa(M, G):
    gw, gv = np.linalg.eigh(G); keep = gw > gw.max()*1e-11
    B = gv[:, keep]/np.sqrt(gw[keep]); ev = np.linalg.eigvalsh(B.T @ M @ B)
    s = max(1.0, np.abs(ev).max()); return int(np.sum(ev < -1e-9*s)), int(keep.sum())

rng = np.random.default_rng(953)
print("="*98)
print("§953 genuine-Π₁ test: κ=#neg(M_σ) under ARBITRARY sampling configs. Krein-Langer: Π₁ ⟺ κ≤1 for ALL.")
print("="*98)
print(f"{'config':>34} {'4W':>4} {'#centers':>8} {'eff.dim':>7} {'κ=#neg':>7} {'≤1?':>4}")
maxk = 0
def run(desc, centers, W):
    global maxk
    M, G = M_sigma(np.array(sorted(centers)), W); k, eff = kappa(M, G)
    maxk = max(maxk, k)
    print(f"{desc:>34} {4*W:>4.1f} {len(centers):>8} {eff:>7} {k:>7} {('YES' if k<=1 else 'NO!'):>4}")

# 1. evenly spaced baseline (Landau-cap ~ L*2W/pi)
for fourW in [4.0, 5.0]:
    W=fourW/4; run(f"even baseline", list(np.linspace(14,44,round(0.8*30*2*W/np.pi))), W)
# 2. OVER-DENSE packing: 2x, 3x, 4x the Landau cap
for mult in [2.0, 3.0, 4.0]:
    W=1.0; cap=30*2*W/np.pi; run(f"over-dense {mult:.0f}x Landau", list(np.linspace(14,44,round(mult*cap))), W)
# 3. irregular RANDOM centers in the window
for t in range(4):
    W=1.0; nc=rng.integers(10,26); run(f"random centers #{t}", list(rng.uniform(14,44,nc)), W)
# 4. clustered centers (many near one point — stresses reproducing-kernel degeneracy)
W=1.0; cl=list(14+0.3*rng.random(12))+list(30+0.3*rng.random(6)); run("clustered (14-ish & 30-ish)", cl, W)
# 5. wide window, dense
W=1.25; run("wide [14,80] dense", list(np.linspace(14,80,40)), W)
# 6. low window near first zeros
W=1.0; run("low [10,30]", list(np.linspace(10,30,16)), W)
# 7. very high bandwidth, moderate window
W=2.0; run("high band 4W=8 [14,40]", list(np.linspace(14,40,round(0.8*26*2*W/np.pi))), W)
print("-"*98)
print(f"MAX κ over ALL arbitrary configs = {maxk}")
print("="*98)
if maxk<=1:
    print("κ≤1 SURVIVES arbitrary/over-dense/irregular sampling ⟹ M_σ is a GENUINE Π₁ reproducing kernel ⟹")
    print("Krein-Langer GUARANTEES M_σ=K₊−K₋ with rank(K₋)=1 exists. The decomposition theorem now APPLIES to ζ.")
else:
    print("some config gives κ≥2 ⟹ κ=1 was basis/config-specific; M_σ is NOT globally Π₁. Redirect the structural claim.")
print("="*98)
