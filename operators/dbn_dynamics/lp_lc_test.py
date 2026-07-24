#!/usr/bin/env python3
"""
lp_lc_test.py — LEDGER §115: Weyl LIMIT-POINT vs LIMIT-CIRCLE at ∞ for the Suzuki screw canonical system (§114's
decisive question). The Weyl disk radius r_a(z) ∝ 1/‖v_z‖²_a (v_z = the ±i deficiency solution). So:
  ‖v_z‖²_a → ∞  ⟹ disk → POINT ⟹ LIMIT-POINT ⟹ m_a→m_∞ unique ⟹ (Herglotz closure §113) RH from convergence alone.
  ‖v_z‖²_a → finite ⟹ disk → CIRCLE ⟹ LIMIT-CIRCLE ⟹ m_∞ is θ-dependent (arithmetic θ(a) essential).
Compute the UNNORMALIZED deficiency solution v_-(x) = T_a^{-1}[e^{-x}] (compact G_a: G u_- = q_-, v_-=∫u_-, λ=0) and
track its L² norm ∫_{-a}^a|v_-|²dx AND its T_a-norm ⟨e^{-x},v_-⟩=‖v_-‖²_{T_a} as a grows. Divergence rate ⟹ LP/LC.
"""
import numpy as np, sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

def deficiency_norms(a, M=1600, lam=0.0):
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    dif=np.arange(-M,M+1)*h; gv=np.array([screw_g(d,Lam,Nmax) for d in dif]); idx=np.arange(M+1)
    G=gv[(idx[:,None]-idx[None,:])+M]
    X=x[:,None]; Y=x[None,:]; N=(X**2+Y**2)/(4*a)-np.abs(X-Y)/2+a/6
    sw=np.sqrt(w); Gm=sw[:,None]*G*sw[None,:]; Gm=(Gm+Gm.T)/2; Nm=sw[:,None]*N*sw[None,:]; Nm=(Nm+Nm.T)/2
    c=sw/np.linalg.norm(sw); P=np.eye(M+1)-np.outer(c,c); A=P@(Gm-lam*Nm)@P; A=(A+A.T)/2
    qm=-(np.exp(-a)-np.exp(-x))                       # q_-(t)=-(e^{-a}-e^{-t})
    um=(np.linalg.pinv(A,rcond=1e-12)@(P@(sw*qm)))/sw
    v=np.zeros(M+1); v[1:]=np.cumsum((um[1:]+um[:-1])/2*h)   # v_-(x)=∫_{-a}^x u_-
    L2=np.sum(v*v*w)                                   # ∫|v_-|²
    Tnorm=np.sum(np.exp(-x)*v*w)                       # ⟨e^{-x},v_-⟩ = ‖v_-‖²_{T_a}
    return L2, Tnorm

print("Weyl LIMIT-POINT vs LIMIT-CIRCLE test: ‖v_-‖² vs a (does it DIVERGE→LP or SATURATE→LC?)")
print(f"{'a':>5} {'∫|v_-|² (L²)':>16} {'‖v_-‖²_T =⟨e^-x,v⟩':>20} {'L² ratio to prev':>16}")
prevL2=None
res=[]
for a in [1.0,1.4,1.8,2.2,2.6,3.0,3.4]:
    L2,Tn=deficiency_norms(a)
    r=(L2/prevL2 if prevL2 else np.nan); prevL2=L2; res.append((a,L2,Tn))
    print(f"{a:>5.1f} {L2:>16.4e} {Tn:>20.4e} {r:>16.3f}",flush=True)
# fit growth: is L2 ~ e^{c a} (diverges, LP) or → const (LC)? log-linear slope of last points
aa=np.array([r[0] for r in res]); LL=np.log(np.abs([r[1] for r in res])+1e-300)
slope=np.polyfit(aa[-4:],LL[-4:],1)[0]
print(f"\n log‖v_-‖² slope over last 4 a: {slope:+.3f}  ({'DIVERGING (e^{'+f'{slope:.2f}'+'a}) ⇒ LIMIT-POINT ⇒ RH from convergence alone' if slope>0.3 else 'flat/saturating ⇒ LIMIT-CIRCLE ⇒ θ(a) arithmetic essential'})")
print(" (LP: disk→point, m_∞ unique, Herglotz closure gives RH. LC: disk→circle, θ-dependent limit, harder.)")
