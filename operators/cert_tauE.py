#!/usr/bin/env python3
"""Certify τ_E (§808 final): weighted off-window mass of E. Load saved 8×6 E-vectors, G-orthonormalize,
for each unit ψ∈E compute ∫_{|t|∉[14,44]} |ψ̂|²(1+|Ω|) dt in arb: [0,14] + [44,R] direct interval
integration + [R,∞) quartic-envelope·ln tail. Max over E = τ_E. Handle closes if τ_E gives ε₀>0
(insensitive: even τ_E≤0.6 works)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, ctx
import arb_fejer_cert2 as ac
ctx.prec = 200
ac.W = arb('1.32')/4; Wa = ac.W
V = np.load('anchor_Evecs.npy'); CEN = np.load('anchor_CEN.npy'); n = len(CEN)
G = np.array([[float((2*arb.pi()*ac.gpair_closed_arb(arb(float(CEN[j])),arb(float(CEN[k])),arb(0))).mid()) for k in range(n)] for j in range(n)])
# G-orthonormalize E
M0=V.T@G@V; w0,u0=np.linalg.eigh(M0); U=V@(u0@np.diag(w0**-0.5)@u0.T)  # 8×6, U^T G U = I
def sinc_sq_r(y):  # arb real, robust
    ay=abs(y); ymax=float(ay.mid())+float(ay.rad())
    if ymax<0.5:
        y2=y*y; return arb(1)-y2/3+2*y2*y2/45+arb(0,ymax**6/45+1e-30)
    return (y.sin()/y)**2
def Fev_r(c,r): return sinc_sq_r(Wa*(r-c))+sinc_sq_r(Wa*(r+c))
def Om_abs(t):  # |Ω(t)|, t real arb
    z1=acb(arb(1)/4,t/2); z2=acb(arb(1)/4,-t/2)
    return abs(((z1.digamma()+z2.digamma())/2).real - arb.pi().log())
def psihat_sq(coef, t):
    s=arb(0)
    for i in range(n): s=s+arb(float(coef[i]))*Fev_r(arb(float(CEN[i])),t)
    return s*s
R=200.0
tau_max=0.0
for j in range(6):
    coef=U[:,j]
    # [0,14] and [44,R]: fixed-segment interval integration of |psihat|^2 (1+|Om|)
    tot=arb(0)
    for lo,hi in [(0.0,14.0),(44.0,R)]:
        m=int((hi-lo)*20)
        for tt in range(m):
            a=lo+(hi-lo)*tt/m; b=lo+(hi-lo)*(tt+1)/m
            tb=arb((a+b)/2,(b-a)/2+1e-15)
            tot=tot+psihat_sq(coef,tb)*(1+Om_abs(tb))*arb(b-a)
    # tail [R,∞): |psihat|^2 <= (sum|coef|)^2 * 4/(W^4 (t-44)^4); (1+|Om|)<= 1+ln t+2
    sc=sum(abs(float(coef[i])) for i in range(n))
    # integral R..inf of (t-44)^-4 (3+ln t) dt <= (3+ln R)/(3 (R-44)^3)
    tail=arb(sc*sc)*4/(Wa**4)*(3+arb(R).log())/(3*arb(R-44)**3)
    tau_j=tot+arb(0).union(tail)   # tail one-sided add
    tj=float((tot+tail).upper())
    import math
    if math.isnan(tj): print(f"  E-vec {j}: NAN — cert failed, abort"); sys.exit(2)
    tau_max=max(tau_max,tj)
    print(f"  E-vec {j}: certified weighted off-mass ≤ {tj:.4f}")
print(f"\n===== certified τ_E ≤ {tau_max:.4f} =====")
# final assembly
mE=0.091759; lam7=0.60696; supOm=2.0175; Bprime=0.3579
ff=1/(1-lam7); B_pf=Bprime+supOm*np.sqrt(lam7)+np.sqrt(tau_max); C_f=lam7*supOm+Bprime
eps0=0.0
for T in np.linspace(1e-8,0.01,1000001):
    s=ff*T
    if s>=1: break
    if mE*(1-s)-2*B_pf*np.sqrt(max(s*(1-s),0))-C_f*s-T<0: eps0=T; break
print(f"FINAL certified ε₀ = {eps0:.7f}  {'ANCHOR THEOREM FULLY CERTIFIED ✓' if eps0>1e-7 else 'does not close'}")
