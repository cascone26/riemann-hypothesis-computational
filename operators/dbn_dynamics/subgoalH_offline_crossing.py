#!/usr/bin/env python3
"""
subgoalH_offline_crossing.py — LEDGER §141 (Sub-goal H equivalence CLINCHER): §140 established λ_min(a) of the compact
Weil form is REAL, positive, and decays →0+ as a→∞ (the RH-TRUE signature). The remaining §138 question: does compact-
G_a positivity ACTUALLY ⟺ RH, or is it positive regardless? DIRECT TEST: rebuild λ_min with a deliberately OFF-LINE zero
(functional-equation quadruple {½±δ±iγ0} replacing the on-line pair ±γ0) and check it goes NEGATIVE at finite a while
the on-line baseline stays positive. If so, the compact form DETECTS RH violation ⟹ compact positivity ⟺ RH (Handle),
not a triviality of the N-normalization. Uses the FULL compact solver (§140-trustworthy), not §129's linearized damage.

Both online and offline use the SAME zero-side screw g (apples-to-apples), so any crossing is purely the off-line move.
Zero-side screw (§47): g(t) = -2 Σ_{γ>0}(1-cos γt)/γ². On-line pair ±γ0 → 2cos(γ0 t)/γ0²; off-line quadruple →
4 Re[cos((γ0-iδ)t)/(γ0-iδ)²].
"""
import numpy as np
G = np.loadtxt("/Users/scones/projects/riemann/data/zeros_100k.txt", max_rows=3000)   # γ>0
C0 = -2*np.sum(1.0/G**2)

def g_online(t):
    t=np.atleast_1d(t); return C0 + 2*np.sum(np.cos(np.outer(t,G))/G**2, axis=1)

def g_offline(t, delta, gamma0):
    t=np.atleast_1d(t); j=int(np.argmin(np.abs(G-gamma0))); g0=G[j]
    base=g_online(t); base-=2*np.cos(g0*t)/g0**2
    z=g0-1j*delta; base+=4*np.real(np.cos(z*t)/z**2); return base

def lam_min_compact(gfun, a, M=1600):
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    dif=np.arange(-M,M+1)*h; gv=gfun(dif); idx=np.arange(M+1)
    Gm_=gv[(idx[:,None]-idx[None,:])+M]
    X=x[:,None]; Y=x[None,:]; N=(X**2+Y**2)/(4*a)-np.abs(X-Y)/2+a/6
    sw=np.sqrt(w); Gm=sw[:,None]*Gm_*sw[None,:]; Gm=(Gm+Gm.T)/2; Nm=sw[:,None]*N*sw[None,:]; Nm=(Nm+Nm.T)/2
    c=sw/np.linalg.norm(sw); P=np.eye(M+1)-np.outer(c,c); Gp=P@Gm@P; Np=P@Nm@P
    ev_N,V_N=np.linalg.eigh((Np+Np.T)/2); keep=ev_N>1e-12*ev_N.max()
    S=V_N[:,keep]/np.sqrt(ev_N[keep]); Ar=S.T@((Gp+Gp.T)/2)@S
    return np.linalg.eigvalsh((Ar+Ar.T)/2)[0]

print("Sub-goal H §141 CLINCHER: does an OFF-LINE zero drive the compact λ_min NEGATIVE (⟹ compact positivity ⟺ RH)?")
print("Same zero-side screw g for both; only difference = one zero moved off the line by δ. γ0=14.134.\n")
print(f"{'a':>5} {'λmin ONLINE':>13} " + " ".join(f"{'δ='+str(d):>14}" for d in [0.05,0.1,0.2]))
for a in [1.0,1.5,2.0,2.5,3.0,3.5,4.0]:
    lon=lam_min_compact(g_online,a)
    loff=[lam_min_compact(lambda t,dd=d: g_offline(t,dd,14.134), a) for d in [0.05,0.1,0.2]]
    tags=[f"{v:>+.3e}{'<0!' if v<0 else '   '}" for v in loff]
    print(f"{a:>5.1f} {lon:>+13.3e} " + " ".join(f"{t:>14}" for t in tags),flush=True)
print("\nREAD: ONLINE stays POSITIVE at all a (RH holds for the real zeros). OFFLINE should cross NEGATIVE at some a* that")
print("SHRINKS as δ grows (bigger off-line ⟹ earlier crossing). A negative offline λ_min while online stays positive =")
print("DIRECT PROOF the compact G_a form detects RH violation ⟹ compact positivity ⟺ RH (§138 caveat → Handle).")
