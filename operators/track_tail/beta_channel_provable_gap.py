#!/usr/bin/env python3
"""
Rie: the beta-channel bound is the exponentially-sensitive gate on M*, and it is NOT proven.

VERIFIED SOLID (Isa's result): the Cauchy corner coupling ||B_C|| <= (1/2)sqrt(log((2M+2)/(M+25)))
<= (1/2)sqrt(log2) = 0.4163, near-tight (op 0.354 vs bound 0.356 at M=117). Real, clean, done.

THE GAP: beta-channel = sum-channel + diff-channel of the Loewner form.
  sum-channel  = beta*sin(t(k+l)/2)cos(t(k-l)/2)/(k+l)   -- modulated Cauchy,  provable <= (pi/2)beta = 0.980
  diff-channel = beta*cos(t(k+l)/2)sin(t(k-l)/2)/(k-l)   -- modulated Toeplitz T[j]=sin(t j/2)/j
    T's symbol theta/2 + 2 Sigma sin(tj/2)cos(jx)/j has SUP = pi (peak at x=0). So ||T|| = pi,
    and diff = (beta/2)(D_{t/2} T D_{t/2} + D_{-t/2} T D_{-t/2}), triangle -> ||diff|| <= beta*pi = 1.96.
  => PROVABLE-TODAY beta-channel bound (triangle) = 0.980 + 1.96 = 2.94.
  MEASURED whole beta-channel = 0.464 (Ein, converged): the two twisted Toeplitz copies CANCEL
  6x. That cancellation is the deferred "interference" -- unproven.

WHY IT MATTERS EXPONENTIALLY: the only slack is d0(k) = log(kb/2pi) - ... which grows LOGARITHMICALLY.
The three-term test needs d_C >= (pi/2 + beta_bnd + 0.03) + (0.4163 + beta_bnd + 0.03)^2 / lam_A.
So M* ~ exp(need), and a 6x gap in beta_bnd becomes a ~2500x gap in M*:
   beta_bnd = 0.464 (measured, UNPROVEN)      -> M* ~ 40      (7 tail modes, trivial)
   beta_bnd = 0.980 (Isa assumed, NO derivation)-> M* ~ 118    (46 tail modes, laptop)
   beta_bnd = 2.94  (CLEAN PROVABLE triangle)  -> M* ~ 98,000  (~49k tail modes, INFEASIBLE)

CONCLUSION: "4' closes with crude classical bounds, laptop tonight" is premature. Under
currently-provable constants M* is infeasible. Feasibility rests entirely on proving the
difference-channel cancellation ||(1/2)(D_{t/2} T D_{t/2} + D_{-t/2} T D_{-t/2})|| ~ 0.46
(vs pi triangle). Target any provable <= 0.6 -> M* ~ 50 (laptop). Routes: (a) symbol of the
SUM of the two twisted Toeplitz copies (does the twist cancel the x=0 peak?); (b) Cotlar-Stein
almost-orthogonality between the +t/2 and -t/2 copies; (c) direct beta-channel bound via the
hhat Fourier structure. Well-posed estimate with a numerical answer key -- not a wall.
"""
import numpy as np, math
b=math.pi/(2*0.549); log2=math.log(2); theta=b*log2; beta=4*(log2/math.sqrt(2))/math.pi

# Toeplitz symbol sup
x=np.linspace(0,2*math.pi,4001); j=np.arange(1,4000)
sym=np.array([theta/2 + 2*np.sum(np.sin(theta*j/2)/j*np.cos(j*xi)) for xi in x])
Tsup=np.max(np.abs(sym))
prov=math.pi/2*beta + beta*Tsup
# measured whole beta-channel
tm=np.arange(25,4002,2).astype(float); K,L=np.meshgrid(tm,tm,indexing='ij'); nz=(K!=L)
whole=np.where(nz, beta*(K*np.sin(theta*L)-L*np.sin(theta*K))/((K-L)*(K+L)), 0)
meas=np.linalg.norm(whole,2)

def Mstar(bb):
    need=(math.pi/2+bb+0.032)+(0.4163+bb+0.032)**2/2.0
    M=4
    while math.log(M*b/(2*math.pi))-0.3614-0.9/M < need: M+=2
    return M
if __name__=="__main__":
    print(f"Toeplitz symbol sup ||T|| = {Tsup:.3f} (= pi)")
    print(f"PROVABLE beta-channel (triangle) = {prov:.3f}   MEASURED = {meas:.3f}   gap {prov/meas:.1f}x")
    for bb,tag in [(meas,'measured UNPROVEN'),(math.pi/2*beta,'Isa 0.98 underived'),(prov,'clean provable')]:
        print(f"  beta<= {bb:.3f} ({tag:20s}) -> M* ~ {Mstar(bb)}")
