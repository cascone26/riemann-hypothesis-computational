#!/usr/bin/env python3
"""RH-CRUX probe: how does the Weil-form min-eigenvalue scale as the test-function basis grows?

Weil positivity criterion: W(φ)=arch+pole−prime ⪰ 0 ∀φ  ⟺  RH. On basis φ_a(x)=e^{−x²/2}cos(ax),
the form is the matrix Q_{jk}. Question that decides the PROOF STRUCTURE:
  - does min-eig(Q) stay robustly > 0 as we add basis functions (→ clean uniform lower bound → RH),
  - or does it decay toward 0 (→ marginal positivity, needs delicate limiting analysis)?
And crucially: does it ever go NEGATIVE at any size (→ RH FALSE on this basis — would be historic)?

Fast float build (rigorous Arb certifies fixed sizes separately, §731/§757/§767). Ω(r) is
BASIS-INDEPENDENT so precompute once on a grid; each arch entry = grid quadrature. Prime sum closed
form g_cross to p≤2e5 (float tail negligible). We sweep TWO families:
  (A) fixed spacing Δa=10, growing count N     — 'extensive' growth
  (B) fixed window [amin,amax], growing density — 'intensive' growth (resolves the inf over the window)
Report min-eig, its trend, and the pivot-margin, honestly."""
import numpy as np
import mpmath as mp
mp.mp.dps = 30

PI = np.pi
LOGPI = np.log(PI)

# --- Ω(r) = Re[ψ(1/4+ir/2)] − log π, basis-independent, precompute on grid ---
# R must exceed max basis freq + ~4 (φ̂ peak width ~1) or the arch integral loses the peak.
R = 240.0            # covers a up to ~215 safely
NG = 48000           # keep spacing 0.01
rg = np.linspace(-R, R, NG)
def Omega_grid():
    out = np.empty_like(rg)
    for i, r in enumerate(rg):
        z = mp.mpf(1)/4 + 1j*mp.mpf(r)/2
        out[i] = float(mp.re(mp.digamma(z))) - LOGPI
    return out
print("precomputing Ω(r) on grid (basis-independent, R=240)...", flush=True)
OM = Omega_grid()

def phihat(a, r):
    s2p = np.sqrt(2*PI)
    return s2p/2*(np.exp(-(r-a)**2/2) + np.exp(-(r+a)**2/2))

def arch_entry(aj, ak):
    integ = phihat(aj, rg)*phihat(ak, rg)*OM
    return np.trapz(integ, rg)/(2*PI) if hasattr(np,'trapz') else np.trapezoid(integ, rg)/(2*PI)

def phihat_ihalf(a):  # φ̂(a, i/2), real
    s2p = np.sqrt(2*PI)
    # e^{-(r-a)^2/2} at r=i/2: (i/2-a)^2 = -1/4 - i a; real part combine both terms
    z1 = np.exp(-((1j/2 - a)**2)/2); z2 = np.exp(-((1j/2 + a)**2)/2)
    return (s2p/2*(z1+z2))
def pole_entry(aj, ak):
    return 2*np.real(phihat_ihalf(aj)*phihat_ihalf(ak))

def g_cross(aj, ak, y):
    s=aj+ak; d=aj-ak
    return np.exp(-y**2/4)*(np.sqrt(PI)/2)*(np.exp(-s**2/4)*np.cos(d*y/2)+np.exp(-d**2/4)*np.cos(s*y/2))

from sympy import primerange
# von Mangoldt support: n = p^k ≤ NMAX, Λ(n)=log p. Build (n, Λ(n)) arrays ONCE.
NMAX = 200000
_p = list(primerange(2, NMAX+1))
_N = []; _LAM = []
for p in _p:
    pk = p
    lp = np.log(p)
    while pk <= NMAX:
        _N.append(pk); _LAM.append(lp)
        pk *= p
_N = np.array(_N, dtype=float); _LAM = np.array(_LAM)
_LOGN = np.log(_N); _INVSQRTN = 1.0/np.sqrt(_N)
def prime_entry(aj, ak):
    # 2 Σ_n Λ(n)/√n · g(log n)   over ALL prime powers n=p^k (von Mangoldt)
    return 2*np.sum(_LAM*_INVSQRTN*g_cross(aj, ak, _LOGN))

def weil_matrix(A):
    n=len(A); Q=np.zeros((n,n))
    for i in range(n):
        for j in range(i,n):
            v = arch_entry(A[i],A[j]) + pole_entry(A[i],A[j]) - prime_entry(A[i],A[j])
            Q[i,j]=Q[j,i]=v
    return Q

print(f"\n{'family':>10} {'N':>3} {'basis':>22} {'min-eig':>12} {'max-eig':>10} {'cond':>9}")
# (A) fixed spacing 10, growing count
for N in [2,3,4,5,6,8,10,12]:
    A=[15+10*k for k in range(N)]
    Q=weil_matrix(A); ev=np.linalg.eigvalsh(Q)
    print(f"{'A-space10':>10} {N:>3} {str(A[0])+'..'+str(A[-1]):>22} {ev.min():>12.6f} {ev.max():>10.3f} {ev.max()/max(ev.min(),1e-12):>9.1f}")
# (B) fixed window [15,55], growing density
for N in [3,5,7,9,12,16,20]:
    A=list(np.linspace(15,55,N))
    Q=weil_matrix(A); ev=np.linalg.eigvalsh(Q)
    print(f"{'B-dense':>10} {N:>3} {'[15,55]/'+str(N):>22} {ev.min():>12.6f} {ev.max():>10.3f} {ev.max()/max(ev.min(),1e-12):>9.1f}")
