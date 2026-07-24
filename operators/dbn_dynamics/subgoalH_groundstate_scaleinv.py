#!/usr/bin/env python3
"""
subgoalH_groundstate_scaleinv.py — LEDGER §159 (RG Step 1): the §156 constant-lift sum rule, via first-order
perturbation, PREDICTS the ground state ψ_1 of the Weil form has ψ_1'(x) ~ e^{|x|/2} (derivative growing toward the
interval boundaries at rate EXACTLY 1/2 = the critical line). Derivation: adding prime p shifts λ_min by lift_p ≈
(log p/√p)·F(log p)/‖ψ_1‖², F(ℓ)=∫∫_{|x−y|≥ℓ}(|x−y|−ℓ)ψ_1'(x)ψ_1'(y); F''(ℓ)=2·A(ℓ) with A=autocorr of ψ_1'. Constant
lift ⟹ F(ℓ)∝e^{ℓ/2}/ℓ ⟹ A(τ)∝e^{τ/2} ⟹ ψ_1'(x)~e^{|x|/2}. TEST all links directly on the true-ζ ground state.
Sine-derivative form (correct positive λ_min at the full operator, clean ψ_1').
"""
import numpy as np, sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

a=2.5; M=2400; K=48
Nmax=int(np.floor(np.exp(2*a)))+1; Lam=vonmangoldt_upto(min(Nmax,30_000_000))
x=np.linspace(-a,a,M+1); h=x[1]-x[0]; w=np.full(M+1,h); w[0]=w[-1]=0.5*h
dif=np.arange(-M,M+1)*h; gv=np.array([screw_g(d,Lam,Nmax) for d in dif]); idx=np.arange(M+1)
G=gv[(idx[:,None]-idx[None,:])+M]
kk=np.arange(1,K+1); arg=np.pi*np.outer((x+a),kk)/(2*a)
Vsin=np.sin(arg); Vp=(kk*np.pi/(2*a))*np.cos(arg)          # basis values & derivatives
Wp=Vp*w[:,None]; Q=(Wp.T@G@Wp); Q=(Q+Q.T)/2
ev,C=np.linalg.eigh(Q/a)
print(f"a={a}: ground λ_min = {ev[0]:.4e}   (full ζ operator)")
c1=C[:,0]                                                   # ground-state coefficients
psi1  = Vsin@c1                                             # ψ_1(x)
psi1p = Vp@c1                                               # ψ_1'(x)

# (1) does |ψ_1'(x)| ~ e^{|x|/2|? fit log|ψ_1'| vs |x| on the outer region (avoid nodes near center)
xr=x[x>0]; dpr=np.abs(psi1p[x>0])
sel=(xr>0.4*a)&(dpr>1e-12)                                  # outer half, away from central nodes
slope_R=np.polyfit(xr[sel], np.log(dpr[sel]),1)[0]
xl=x[x<0]; dpl=np.abs(psi1p[x<0]); sell=(xl<-0.4*a)&(dpl>1e-12)
slope_L=np.polyfit(np.abs(xl[sell]), np.log(dpl[sell]),1)[0]
print(f"(1) log|ψ_1'| growth rate toward boundary:  RIGHT={slope_R:+.3f}  LEFT={slope_L:+.3f}   (predict +0.500)")

# (2) autocorrelation A(τ)=∫ψ_1'(x)ψ_1'(x+τ)dx ; does A(τ)~e^{τ/2}?
def autocorr(tau):
    n=int(round(tau/h))
    if n>=M: return 0.0
    return np.sum(psi1p[:M+1-n]*psi1p[n:]*w[:M+1-n])
taus=np.linspace(0.3,2.0,18); A=np.array([autocorr(t) for t in taus])
selA=np.abs(A)>1e-12
rateA=np.polyfit(taus[selA], np.log(np.abs(A[selA])),1)[0]
print(f"(2) ψ_1'-autocorrelation A(τ)~e^(rate·τ): rate = {rateA:+.3f}   (predict +0.500)")

# (3) lift prediction per prime: lift_p ≈ (log p/√p)·F(log p)/‖ψ_1‖²,  compare CONSTANCY vs p
nrm=np.sum(psi1**2*w)
def Ffun(ell):
    Kk=np.where(adif:=np.abs(dif)>=ell, np.abs(dif)-ell, 0.0)  # (|τ|-ℓ)·1[|τ|≥ℓ]
    Gk=Kk[(idx[:,None]-idx[None,:])+M]
    return (Wp.T@Gk@Wp)  # matrix; quadratic form on c1:
def lift_pred(p):
    ell=np.log(p); Kk=np.where(np.abs(dif)>=ell, np.abs(dif)-ell, 0.0)
    Gk=Kk[(idx[:,None]-idx[None,:])+M]
    Fval=(Wp.T@Gk@Wp); Fq=c1@Fval@c1
    return (np.log(p)/np.sqrt(p))*Fq/nrm
print("(3) first-order lift prediction per prime (should be ~CONSTANT, matching §156):")
for p in [2,3,5,7,11,13,17,23,31,41,53]:
    print(f"     p={p:>3}: lift_pred = {lift_pred(p):+.4e}",flush=True)
print("\nREAD: RIGHT/LEFT growth ≈ +0.5 AND autocorr rate ≈ +0.5 ⟹ ψ_1'~e^{|x|/2} CONFIRMED = the ground state is")
print("scale-invariant in the ½-line metric (the fixed-point signature). lift_pred ~ constant in p ⟹ the first-order")
print("theory REPRODUCES §156's constant lift ⟹ the sum rule is PROVEN-from-ψ_1-structure (RG Step 1 deposit).")
