#!/usr/bin/env python3
"""§955 IS ζ's WEYL m-FUNCTION A GENERALIZED NEVANLINNA N₁? (direct test of the §950/§953 route's crux)

The route (deep-extraction agent): RH via Krein-Langer needs ζ's Weyl/Nevanlinna m-function to be a GENERALIZED
NEVANLINNA function N_κ with κ=1 — its Nevanlinna kernel K_m(z,w)=(m(z)−conj(m(w)))/(z−conj(w)) has exactly ONE
negative square, the single negative-type pole ↔ ζ's pole at s=1 (arXiv:1011.2081: single simple pole ⟹ κ=1).
Here we build m(z)=vp(z)/vm(z) from the exact Suzuki canonical-system solutions (W_family: cplus,cminus from the
boundary data f_±=E_{±1} at s=1,s=0), sample K_m on a grid of UHP points, and COUNT its negative squares.
CONTROLS first (validate the Nevanlinna-Gram code): m(z)=z is Herglotz N₀ (κ=0); m(z)=z² is N₁ (κ=1).
Honest caveats printed. HANDLE (own code). If κ(K_m)=1 for ζ's Weyl function, that is numerical evidence the
ζ m-function is N₁ — the object the Krein-Langer route requires. If κ≠1, the naive ratio isn't the right m — redirect."""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from suzuki_exact_Q import build_Q, E_p, shat
from scipy.linalg import solve

def nevanlinna_kappa(mfun, zpts, tol=1e-9):
    """negative squares of the Nevanlinna kernel K(z,w)=(m(z)-conj(m(w)))/(z-conj(w)) sampled at zpts (UHP)."""
    m = np.array([mfun(z) for z in zpts]); Z = np.array(zpts)
    W = np.conj(Z)
    num = m[:,None] - np.conj(m)[None,:]
    den = Z[:,None] - W[None,:]
    # diagonal limit z->z: K(z,z) = Im m(z)/Im z
    K = np.empty((len(Z),len(Z)), dtype=complex)
    for i in range(len(Z)):
        for j in range(len(Z)):
            if abs(den[i,j]) < 1e-12:
                K[i,j] = np.imag(m[i])/np.imag(Z[i])
            else:
                K[i,j] = num[i,j]/den[i,j]
    K = 0.5*(K + K.conj().T)
    ev = np.linalg.eigvalsh(K); s = max(1.0, np.abs(ev).max())
    return int(np.sum(ev < -tol*s)), ev

# UHP sample grid (spread in x, moderate y>0)
rng = np.random.default_rng(955)
def uhp_grid(n=24, xlo=2, xhi=40, ylo=0.3, yhi=3.0):
    return [complex(x, y) for x, y in zip(rng.uniform(xlo,xhi,n), rng.uniform(ylo,yhi,n))]

zp = uhp_grid()
print("="*92)
print("§955 Nevanlinna negative-square count κ of ζ's Weyl m-function. Controls validate the kernel code.")
print("="*92)
# CONTROLS
kz,_ = nevanlinna_kappa(lambda z: z, zp);            print(f"  CONTROL m(z)=z   (Herglotz N₀, expect κ=0): κ = {kz}")
kz2,_ = nevanlinna_kappa(lambda z: z*z, zp);         print(f"  CONTROL m(z)=z²  (N₁,        expect κ=1): κ = {kz2}")
kz3,_ = nevanlinna_kappa(lambda z: -1.0/z, zp);      print(f"  CONTROL m(z)=-1/z(Herglotz N₀, expect κ=0): κ = {kz3}")
kz4,_ = nevanlinna_kappa(lambda z: z**3, zp);        print(f"  CONTROL m(z)=z³  (N₁,        expect κ=1): κ = {kz4}")
print("-"*92)

def weyl_m(a, K, lam_shift=0.2, NG=400):
    kk = np.arange(1, K+1)
    Q = build_Q(a, K, NG); lam_a = np.linalg.eigvalsh(Q)[0]/a
    lam = lam_a - lam_shift; Ta = Q - lam*a*np.eye(K)
    fplus = E_p(1.0, a, kk); fminus = E_p(-1.0, a, kk)
    cplus = solve(Ta, fplus, assume_a='sym'); cminus = solve(Ta, fminus, assume_a='sym')
    npl = np.sqrt(abs(cplus@fplus)); nmi = np.sqrt(abs(cminus@fminus)); cminus *= npl/nmi
    def m(z):
        vp = np.sum(cplus*shat(z, a, kk)); vm = np.sum(cminus*shat(z, a, kk))
        return vp/vm
    return m, lam_a

print("  ζ Weyl m-function m(z)=vp/vm from exact Suzuki canonical system:")
print(f"  {'a':>5} {'K':>4} | {'κ(Nevanlinna kernel of m)':>26}")
for a in [1.0, 1.5, 2.0]:
    for K in [70, 100]:
        m, _ = weyl_m(a, K)
        k, ev = nevanlinna_kappa(m, zp)
        print(f"  {a:>5.1f} {K:>4} | {k:>26}")
print("="*92)
print("READ: controls must give 0,1,0,1 (else kernel code is wrong — ignore ζ rows). If ζ's m gives κ=1, the ζ")
print("Weyl m-function is (numerically) a generalized Nevanlinna N₁ — one negative-type pole = ζ's s=1 pole =")
print("the object the Krein-Langer route requires (single pole ⟹ κ=1 ⟹ M_σ is Π₁ ⟹ pole absorbs it ⟹ RH).")
print("Honest caveat: m=vp/vm is a natural but not-uniquely-canonical Weyl normalization; κ is the robust readout.")
print("="*92)
