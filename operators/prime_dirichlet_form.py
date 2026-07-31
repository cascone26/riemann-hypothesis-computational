#!/usr/bin/env python3
"""RT loop-2 Leg 2 (positivity-preserving / leverage card), resolved rigorously.
The prime form V (V[psi](u)=sum_n Lambda(n)/sqrt(n)[psi(u+ln n)+psi(u-ln n)]) has a DIRICHLET decomposition:
    <psi, V psi> = 2 C_N ||psi||^2  -  E_jump(psi),
    C_N = sum_{n<=N} Lambda(n)/sqrt(n)  (divergent ~2 sqrt(N), = PNT),
    E_jump(psi) = sum_n Lambda(n)/sqrt(n) * integral [psi(u)-psi(u+ln n)]^2 du  >= 0  (a genuine Dirichlet jump form).
Hence  <psi, H_F psi> = <psi,L0 psi> + kappa E_jump(psi)  -  2 kappa C_N ||psi||^2
     = [SUM OF TWO POSITIVE FORMS]  -  [divergent constant].
So H_F + 2 kappa C_N = L0 + kappa E_jump >= 0: self-adjoint, bounded below, REAL spectrum at finite N (Friedrichs).
The RT's >=0-kernel / positivity-preservation intuition is CORRECT (E_jump is a real positive Dirichlet form).
BUT leverage is NOT obtained: the content (a<=1 = RH) is whether the renormalized N->inf operator (constant
2 kappa C_N ~ 2 kappa sqrt(N) subtracted) stays bounded below -- the balance of the positive forms vs the
divergent constant. That IS the razor edge, not a shortcut past it.
This script VERIFIES the identity numerically (both sides match, E_jump>=0) to make the structure concrete.
"""
import numpy as np

def vonmangoldt(N):
    lam = np.zeros(N+1)
    for p in range(2, N+1):
        if lam[p] == 0:  # p prime (not yet marked)
            lp = np.log(p)
            pk = p
            while pk <= N:
                lam[pk] = lp
                pk *= p
    return lam

# u-grid (log x). test functions psi(u): a few smooth bumps, compactly supported well inside the grid.
du = 0.002
u = np.arange(-18.0, 18.0, du)   # wide: shifts up to ln(400)~6 keep psi (supported near 0) fully in-domain
LAM = vonmangoldt(400)
ns = np.arange(2, 401)
w = np.array([LAM[n]/np.sqrt(n) for n in ns])           # Lambda(n)/sqrt(n)
logn = np.log(ns)
C_N = w.sum()

def shift(psi, s):
    """psi(u+s) via linear interpolation on the grid (0 outside)."""
    return np.interp(u + s, u, psi, left=0.0, right=0.0)

def Vpsi(psi):
    out = np.zeros_like(psi)
    for wi, ln in zip(w, logn):
        out += wi*(shift(psi, ln) + shift(psi, -ln))
    return out

def test(psi, name):
    nrm2 = np.sum(psi*psi)*du
    lhs = np.sum(psi*Vpsi(psi))*du                       # <psi, V psi>
    Ejump = 0.0
    for wi, ln in zip(w, logn):
        d = psi - shift(psi, ln)
        Ejump += wi*np.sum(d*d)*du
    rhs = 2*C_N*nrm2 - Ejump
    print(f"{name:18} ||psi||^2={nrm2:.5f}  <psi,Vpsi>={lhs:+.5f}  2C_N||psi||^2-E_jump={rhs:+.5f}  "
          f"rel.err={abs(lhs-rhs)/(abs(lhs)+1e-12):.2e}  E_jump={Ejump:.4f}(>=0:{Ejump>=0})")
    return lhs, rhs, Ejump

print(f"C_N = sum_2^400 Lambda(n)/sqrt(n) = {C_N:.4f}   (~2 sqrt(400)=40; PNT growth)\n")
print("Verify identity  <psi,V psi> = 2 C_N ||psi||^2 - E_jump(psi),  E_jump>=0 (Dirichlet jump form):")
for c,wd,nm in [(0.0,0.6,"gauss@0 w=.6"),(1.0,0.5,"gauss@1 w=.5"),(-1.5,0.8,"gauss@-1.5 w=.8")]:
    psi = np.exp(-0.5*((u-c)/wd)**2)
    test(psi, nm)
# an oscillatory one (closer to a real eigenfunction)
psi = np.exp(-0.5*(u/2.0)**2)*np.cos(4*u)
test(psi, "gauss*cos(4u)")
print("\n=> identity holds; E_jump>=0 always. The prime form = (divergent const 2C_N) - (positive Dirichlet form).")
print("   H_F + 2*kappa*C_N = L0 + kappa*E_jump >= 0 (two positive forms): self-adjoint/real spectrum at finite N.")
print("   Leverage REFRAMED not obtained: a<=1 (=RH) is the renormalized N->inf balance of E_jump+L0 vs 2 kappa C_N~sqrt(N).")
