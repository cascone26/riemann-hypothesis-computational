#!/usr/bin/env python3
"""
subgoalH_homotopy_refine.py — LEDGER §155 (Lens 2, the crux): §154 found λ_min(t) is CONCAVE and peaks near the true
primes t=1 at a tangency to 0. PIN IT: is the argmax EXACTLY t=1, robustly across a? Is it a smooth max or a CUSP
(eigenvalue crossing at t=1)? Test the prime-neutrality equation dλ_min/dt|_{t=1} = v_1'·(prime form)·v_1 = 0.
For each a: fine t-grid near 1 → argmax t*, peak λ*; lowest-2 eigenvalues near t=1 (crossing?); ground-state prime-
Rayleigh at t=1 (the derivative / prime-neutrality residual).
"""
import numpy as np, sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, screw_g_no_primes, vonmangoldt_upto

def setup(a, M):
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    dif=np.arange(-M,M+1)*h
    g_arch=np.array([screw_g_no_primes(d,Lam,Nmax) for d in dif])
    g_prime=np.array([screw_g(d,Lam,Nmax) for d in dif])-g_arch
    idx=np.arange(M+1); X=x[:,None]; Y=x[None,:]
    N=(X**2+Y**2)/(4*a)-np.abs(X-Y)/2+a/6
    sw=np.sqrt(w); Nm=sw[:,None]*N*sw[None,:]; Nm=(Nm+Nm.T)/2
    c=sw/np.linalg.norm(sw); Pp=np.eye(M+1)-np.outer(c,c)
    Np=Pp@Nm@Pp; Np=(Np+Np.T)/2
    ev_N,V_N=np.linalg.eigh(Np); keep=ev_N>1e-12*ev_N.max()
    Sw=V_N[:,keep]/np.sqrt(ev_N[keep])
    def whiten(gv):
        G=gv[(idx[:,None]-idx[None,:])+M]; Gm=sw[:,None]*G*sw[None,:]; Gm=(Gm+Gm.T)/2
        Gp=Pp@Gm@Pp; return Sw.T@((Gp+Gp.T)/2)@Sw
    Aarch=whiten(g_arch); Aprime=whiten(g_prime)     # whitened arch & prime operators; Ar(t)=Aarch+t·Aprime
    return Aarch, Aprime

print("Sub-goal H §155: pin the homotopy peak — is argmax(λ_min(t)) EXACTLY t=1? smooth max or cusp? prime-neutral?\n")
for a in [1.5, 2.0, 2.5]:
    Aarch,Aprime = setup(a, 1600)
    ts=np.arange(0.90,1.101,0.01)
    lam=[]; lam2=[]
    for t in ts:
        ev=np.linalg.eigvalsh(Aarch+t*Aprime); lam.append(ev[0]); lam2.append(ev[1])
    lam=np.array(lam); lam2=np.array(lam2)
    it=np.argmax(lam); tstar=ts[it]; peak=lam[it]
    # ground state at t=1 and prime-Rayleigh (dλ/dt)
    ev1,U1=np.linalg.eigh(Aarch+1.0*Aprime); u1=U1[:,0]
    dldt = u1@Aprime@u1                    # envelope-thm derivative dλ_min/dt|_{t=1}
    # smooth vs cusp: gap between λ_min and λ_2 at t=1 (small gap ⟹ near-crossing/cusp)
    gap1 = ev1[1]-ev1[0]
    print(f"=== a={a} ===")
    print(f"  argmax t* = {tstar:.3f}   peak λ* = {peak:+.4e}   λ_min(t=1) = {lam[np.argmin(np.abs(ts-1.0))]:+.4e}")
    print(f"  dλ_min/dt|_(t=1) = v1'·(prime)·v1 = {dldt:+.4e}   (≈0 ⟹ prime-NEUTRAL ground state ⟹ peak at t=1)")
    print(f"  gap λ2−λ1 at t=1 = {gap1:.3e}   (tiny ⟹ CUSP/near-crossing at true primes; O(1) ⟹ smooth max)")
    # print the peak shape
    near=(ts>=0.96)&(ts<=1.04)
    print("  λ_min near t=1:  " + "  ".join(f"{t:.2f}:{l:+.3e}" for t,l in zip(ts[near],lam[near])),flush=True)
    print()
print("READ: t*≈1.00 across all a + dλ/dt|_1≈0 ⟹ the true primes sit EXACTLY at the concave-max tangency (robust,")
print("structural). Small λ2−λ1 gap ⟹ the peak is a CUSP = an eigenvalue crossing at t=1 (even sharper). This is the")
print("prime-neutrality characterization of RH-marginality on OUR operator — a new, concrete critical-point equation.")
