#!/usr/bin/env python3
"""
lens_g_execute.py — LEDGER §129 (Lens G execution, cancellation-free): the DAMAGE an off-line zero does to λ_a, via
the ZERO-SIDE screw function (no difference-of-rank-1 cancellation, unlike §127).
Zero-side screw (§47, matches prime-side to 4e-5):  g(t) = -2 Σ_{γ>0}(1-cos γt)/γ² = -0.0462 + 2 Σ_γ cos(γt)/γ².
Move one zero off-line: the functional-equation QUADRUPLE {½±δ±iγ₀} replaces the on-line pair ±γ₀.
  on-line pair ±γ₀ contributes:  2 cos(γ₀ t)/γ₀².
  off-line quadruple contributes: 4 Re[ cos((γ₀-iδ)t)/(γ₀-iδ)² ]  (derived: Σ_{4-tuple} e^{izt}/z² = 4Re[cos(z₁t)/z₁²]).
Compute λ_a (smallest gen-eigenvalue) via the compact G_a (§110) for g_online vs g_offline; damage = λ_online−λ_offline.
"""
import numpy as np

G = np.loadtxt("/Users/scones/projects/riemann/data/zeros_100k.txt", max_rows=3000)   # γ>0
C0 = -2*np.sum(1.0/G**2)   # constant -0.0462...

def g_online(t):
    t=np.atleast_1d(t)
    return C0 + 2*np.sum(np.cos(np.outer(t,G))/G**2, axis=1)

def g_offline(t, delta, gamma0):
    # start from online, replace the nearest-zero pair (γ₀) with the off-line quadruple
    t=np.atleast_1d(t)
    j=int(np.argmin(np.abs(G-gamma0))); g0=G[j]
    base = g_online(t)
    base -= 2*np.cos(g0*t)/g0**2                              # remove on-line pair ±g0
    z=g0-1j*delta
    base += 4*np.real(np.cos(z*t)/z**2)                       # add off-line quadruple
    return base

def lam_compact(gfun, a, M=1200):
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    dif=np.arange(-M,M+1)*h; gv=gfun(dif); idx=np.arange(M+1)
    Gm_=gv[(idx[:,None]-idx[None,:])+M]
    X=x[:,None]; Y=x[None,:]; N=(X**2+Y**2)/(4*a)-np.abs(X-Y)/2+a/6
    sw=np.sqrt(w); Gm=sw[:,None]*Gm_*sw[None,:]; Gm=(Gm+Gm.T)/2; Nm=sw[:,None]*N*sw[None,:]; Nm=(Nm+Nm.T)/2
    c=sw/np.linalg.norm(sw); P=np.eye(M+1)-np.outer(c,c); Gp=P@Gm@P; Np=P@Nm@P
    ev_N,V_N=np.linalg.eigh((Np+Np.T)/2); keep=ev_N>1e-12*ev_N.max()
    S=V_N[:,keep]/np.sqrt(ev_N[keep]); Ar=S.T@((Gp+Gp.T)/2)@S
    return np.linalg.eigvalsh((Ar+Ar.T)/2)[0]

print("Lens G EXECUTE: λ_a(online) vs λ_a(offline, zero at ½+δ+iγ₀); damage = how far off-line zero pushes λ_a down.")
print(f"{'a':>4} {'λ_online':>12} " + " ".join(f"{'dmg δ='+str(d):>13}" for d in [0.05,0.1,0.2]) + f"  (γ₀=14.13)")
for a in [1.0,1.5,2.0,2.5,3.0]:
    lon=lam_compact(g_online, a)
    dmgs=[lon-lam_compact(lambda t,dd=d: g_offline(t,dd,14.134), a) for d in [0.05,0.1,0.2]]
    print(f"{a:>4.1f} {lon:>12.4e} " + " ".join(f"{d:>13.4e}" for d in dmgs), flush=True)
print("\n  scaling check (damage vs a, δ=0.1, γ₀=14.13): does it grow ~e^{2δa}? ratio per Δa=0.5 predicted e^{2·0.1·0.5}=e^{0.1}=1.105")
prev=None
for a in [1.0,1.5,2.0,2.5,3.0]:
    d=lam_compact(g_online,a)-lam_compact(lambda t: g_offline(t,0.1,14.134),a)
    r='' if prev is None else f'  ratio={d/prev:.3f}'; prev=d
    print(f"   a={a}: damage={d:.4e}{r}",flush=True)
print("\nREAD: off-line zero pushes λ_a NEGATIVE (damage>0). Growth in a = the amplification law; scaling in δ = strength.")
print("This is the per-zero damage for the δ(a) bound + the detector (measured λ_a<−damage ⟹ off-line zero present).")
