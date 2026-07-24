#!/usr/bin/env python3
"""
subgoalH_dbn_flow.py — LEDGER §157 (collective-native pass): how does OUR compact-Weil λ_min(a) evolve under the
de Bruijn–Newman flow? The dBN flow moves ALL zeros collectively (RT zero dynamics, §131 convention: forward s makes
zeros more real). RT proved Λ≥0 with a position-only energy; RH ⟺ Λ≤0 needs an ARITHMETIC-COUPLED monotone functional —
and λ_min(a) from the compact operator IS arithmetic-coupled. Probe: flow the real zeros γ_j(s)=γ_j+s·v_j (first order,
v_j=2Σ_{k≠j}1/(γ_j−γ_k)), rebuild the zero-side screw g_s, compute λ_min(a;s). QUESTIONS:
  (1) is s=0 (true ζ) an EXTREMUM of λ_min(a;s)? (like the prime homotopy peaked at t=1)
  (2) is λ_min(a;s) MONOTONE in s (a candidate flow-monotone functional)?
  (3) slope dλ_min/ds at 0 — sign & magnitude (the arithmetic coupling to the flow).
Zero-side screw: g_s(τ) = −2Σ_j (1−cos(γ_j(s)τ))/γ_j(s)² = C0(s) + 2Σ_j cos(γ_j(s)τ)/γ_j(s)².
"""
import numpy as np

G = np.loadtxt("/Users/scones/projects/riemann/data/zeros_100k.txt", max_rows=3000)   # first 3000 γ>0 (resolved g)
# dBN flow velocities v_j = 2 Σ_{k≠j} 1/(γ_j − γ_k)   (RT / §131 forward convention)
D = G[:,None] - G[None,:]
np.fill_diagonal(D, np.inf)                    # skip k=j
V = 2.0*np.sum(1.0/D, axis=1)
print(f"loaded {len(G)} zeros; velocity stats: |v| median={np.median(np.abs(V)):.3f} max={np.max(np.abs(V)):.3f}")

def screw_s(tau, s):
    g = G + s*V                                # flowed zeros (real, first order)
    tau=np.atleast_1d(tau)
    C0 = -2.0*np.sum(1.0/g**2)
    return C0 + 2.0*np.sum(np.cos(np.outer(tau,g))/g**2, axis=1)

def lam_min(a, s, M=1400):
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    dif=np.arange(-M,M+1)*h; gv=screw_s(dif,s); idx=np.arange(M+1)
    Gm_=gv[(idx[:,None]-idx[None,:])+M]
    X=x[:,None]; Y=x[None,:]; N=(X**2+Y**2)/(4*a)-np.abs(X-Y)/2+a/6
    sw=np.sqrt(w); Gm=sw[:,None]*Gm_*sw[None,:]; Gm=(Gm+Gm.T)/2; Nm=sw[:,None]*N*sw[None,:]; Nm=(Nm+Nm.T)/2
    c=sw/np.linalg.norm(sw); P=np.eye(M+1)-np.outer(c,c); Gp=P@Gm@P; Np=P@Nm@P
    ev_N,V_N=np.linalg.eigh((Np+Np.T)/2); keep=ev_N>1e-12*ev_N.max()
    S=V_N[:,keep]/np.sqrt(ev_N[keep]); Ar=S.T@((Gp+Gp.T)/2)@S
    return np.linalg.eigvalsh((Ar+Ar.T)/2)[0]

print("\ndBN flow of the compact λ_min(a=1.5;s)  (s<0 backward toward Λ, s>0 forward more-real):\n")
print(f"{'s':>8} {'λ_min':>13} {'Δ vs s=0':>12}")
l0=lam_min(1.5,0.0)
ss=[-0.05,-0.03,-0.02,-0.01,-0.005,0.0,0.005,0.01,0.02,0.03,0.05]
lams=[]
for s in ss:
    lm=lam_min(1.5,s); lams.append(lm)
    print(f"{s:>8.3f} {lm:>13.5e} {lm-l0:>+12.3e}",flush=True)
lams=np.array(lams); ss=np.array(ss)
# is s=0 an extremum? slope each side
i0=list(ss).index(0.0)
sl_left=(lams[i0]-lams[i0-1])/(ss[i0]-ss[i0-1]); sl_right=(lams[i0+1]-lams[i0])/(ss[i0+1]-ss[i0])
print(f"\n slope just LEFT of 0: {sl_left:+.4e}   slope just RIGHT: {sl_right:+.4e}")
print(f" s=0 is {'a MAX (peak)' if sl_left>0>sl_right else 'a MIN' if sl_left<0<sl_right else 'MONOTONE-through' if sl_left*sl_right>0 else 'flat/cusp'}")
print(f" λ_min(s=0) = {l0:.5e}   min over range = {lams.min():.3e} @s={ss[np.argmin(lams)]:.3f}")
print("\nREAD: if s=0 is a MAX (peak) ⟹ the flow (either direction) reduces λ_min = same tangency as the prime homotopy")
print("(RH marginal, zeros perturbed = less positive). If MONOTONE increasing ⟹ a candidate flow-monotone arithmetic")
print("functional (forward flow ↑ positivity) → backward-extrapolate the s where it crosses 0 = a Λ estimate from OUR operator.")
