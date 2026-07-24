#!/usr/bin/env python3
"""
subgoalH_ladder_id.py — LEDGER §143 (Lens 1 of the §142 prism): IDENTIFY the near-null ladder λ_k≈9.26e-7·k².
§137 found the compact Weil form's low spectrum is a discrete Weyl ladder with a nearly a-independent floor ~9.26e-7.
Lens 1: that k² ladder = a vibrating-string spectrum ⟹ a FIXED limiting Sturm-Liouville operator at the bottom of G_a.
If so, "λ_min(a)>0 ∀a" collapses to "one fixed operator is positive." This tool:
  (1) fits λ_k = C·k^p (log-log, k=1..12) to pin the exponent p and coefficient C, vs a — is p≈2 and C≈const?
  (2) DECISIVE geometric-vs-arithmetic test: recompute the SAME ladder for ARCH-ONLY g (screw_g_no_primes, primes
      removed). If ARCH gives the SAME C·k² ladder as FULL ⟹ the bottom ladder is GEOMETRIC (arch/N structure), and
      the primes/RH content lives only in the BULK/interface (deflates Lens 1, sharpens Lens 2). If primes CHANGE the
      ladder ⟹ the bottom IS arithmetic (Lens 1 alive — the universal operator carries the primes).
  (3) mode-shape check: correlate mode k with sin(kπ(x+a)/2a) (particle-in-box) — sines ⟹ Laplacian-like bottom.
  (4) constant hunt: print candidate closed forms for C (1/γ₁ⁿ, (Σ2/γ²)², arch-derivative scales) to eyeball a match.
"""
import numpy as np, sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, screw_g_no_primes, vonmangoldt_upto

GAMMA1 = 14.134725142
SUM2_OVER_G2 = None  # filled below

def weil_modes(a, M, gfun, ncore=16):
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    dif=np.arange(-M,M+1)*h; gv=np.array([gfun(d,Lam,Nmax) for d in dif]); idx=np.arange(M+1)
    G=gv[(idx[:,None]-idx[None,:])+M]
    X=x[:,None]; Y=x[None,:]; N=(X**2+Y**2)/(4*a)-np.abs(X-Y)/2+a/6
    sw=np.sqrt(w); Gm=sw[:,None]*G*sw[None,:]; Gm=(Gm+Gm.T)/2; Nm=sw[:,None]*N*sw[None,:]; Nm=(Nm+Nm.T)/2
    c=sw/np.linalg.norm(sw); P=np.eye(M+1)-np.outer(c,c); Gp=P@Gm@P; Np=P@Nm@P
    ev_N,V_N=np.linalg.eigh((Np+Np.T)/2); keep=ev_N>1e-12*ev_N.max()
    S=V_N[:,keep]/np.sqrt(ev_N[keep]); Ar=S.T@((Gp+Gp.T)/2)@S
    wv,U=np.linalg.eigh((Ar+Ar.T)/2)
    funcs=(S@U); funcs=funcs/np.where(sw[:,None]>0,sw[:,None],1.0)
    return wv[:ncore], funcs[:,:ncore], x

def fit_ladder(ev):
    ev=ev[ev>0]; k=np.arange(1,len(ev)+1)
    kk=k[:12]; ee=ev[:12]
    p,logC=np.polyfit(np.log(kk), np.log(ee), 1)
    return p, np.exp(logC)

def sine_overlap(funcs, x, a, kmodes=6):
    ovs=[]
    for j in range(kmodes):
        best=0
        for kk in range(1,kmodes+3):
            s=np.sin(kk*np.pi*(x+a)/(2*a)); s=s/np.linalg.norm(s)
            v=funcs[:,j]/np.linalg.norm(funcs[:,j])
            best=max(best, abs(np.dot(s,v)))
        ovs.append(best)
    return ovs

M=2000
print("Lens 1 §143: identify the near-null ladder λ_k≈C·k^p — geometric (arch) or arithmetic (primes)?\n")
print(f"{'a':>4} | {'FULL: p':>7} {'C':>11} | {'ARCH: p':>7} {'C':>11} | {'C_full/C_arch':>13} | {'sine-overlap FULL (modes1-6)':>30}")
for a in [1.5,2.0,2.5,3.0]:
    evF,fF,x = weil_modes(a,M,screw_g)
    evA,fA,_ = weil_modes(a,M,screw_g_no_primes)
    pF,CF = fit_ladder(evF); pA,CA = fit_ladder(evA)
    ov = sine_overlap(fF,x,a)
    print(f"{a:>4.1f} | {pF:>7.3f} {CF:>11.3e} | {pA:>7.3f} {CA:>11.3e} | {CF/CA:>13.3f} | {np.round(ov,3)}",flush=True)

# constant hunt using a=2.0 full ladder coefficient
evF,_,_ = weil_modes(2.0,M,screw_g); _,C2 = fit_ladder(evF)
G3000 = np.loadtxt("/Users/scones/projects/riemann/data/zeros_100k.txt", max_rows=3000)
S2 = np.sum(2.0/G3000**2)   # Σ 2/γ² = -C0 = 0.0462...
print(f"\nConstant hunt (a=2.0 full ladder C = {C2:.4e}):")
cands = {
 "1/γ₁²":1/GAMMA1**2, "1/γ₁⁴":1/GAMMA1**4, "1/γ₁⁶":1/GAMMA1**6,
 "(2/γ₁²)²":(2/GAMMA1**2)**2, "(Σ2/γ²)²":S2**2, "(Σ2/γ²)²/γ₁²":S2**2/GAMMA1**2,
 "1/γ₁⁴/π²":1/GAMMA1**4/np.pi**2, "4/γ₁⁶":4/GAMMA1**6, "(2/γ₁³)²":(2/GAMMA1**3)**2,
}
for name,val in sorted(cands.items(), key=lambda kv: abs(np.log(kv[1]/C2))):
    print(f"   {name:>14} = {val:.4e}   ratio to C = {val/C2:>8.3f}")
print("\nREAD: p≈2 confirms the k² string ladder. If C_full/C_arch ≈ 1 ⟹ ladder is GEOMETRIC (primes don't set the floor;")
print("RH content is in the bulk/interface — Lens 2). If C_full/C_arch clearly ≠ 1 ⟹ the floor IS arithmetic (Lens 1 alive:")
print("a fixed prime-carrying limiting operator). High sine-overlap ⟹ the bottom operator is Laplacian-like (particle-in-box).")
