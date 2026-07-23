#!/usr/bin/env python3
"""§191: are the 2 arch negative modes RECOGNIZABLE closed-form functions? If u₁,u₂ are (say) low-degree
polynomials / prolates, the p=2 alignment inequality becomes explicit & provable. Extract u₁,u₂ at a=0.5 as
functions v(x) on (−a,a) (v'=Σc_k sin_k basis), report: parity, Legendre-poly fit residuals, and the p=2
involution overlap structure. mpmath. Reuses weil_matrices."""
import sys, math
import mpmath as mp
sys.path.insert(0,'/Users/scones/projects/riemann/operators/dbn_dynamics')
from mpmath_weil_hiprec import weil_matrices, gen_eig_sym
mp.mp.dps=40; a=mp.mpf('0.5'); K=24; M=400
Qf,Qa,Mass,Vv,xs,w = weil_matrices(a,K,M)
Ea,Va = gen_eig_sym(Qa,Mass)
print(f"§191 arch negative-mode structure @ a={a} K={K}")
print(f"arch spectrum lowest 4: {[mp.nstr(e,5) for e in Ea[:4]]}")
# reconstruct v_i'(x) = Σ_k Va[k,i]*sin(kπ(x+a)/2a); v_i(x)=integral. Just examine v'(x) shape + parity of v'.
def vprime(i):
    return [mp.fsum(Va[k,i]*mp.sin((k+1)*mp.pi*(xs[j]+a)/(2*a)) for k in range(K)) for j in range(M+1)]
for i in (0,1):
    vp=vprime(i)
    # parity of v' about x=0 (index M/2): even v'⟺ v odd; report
    mid=M//2
    sym=mp.fsum(vp[j]*vp[M-j]*w[j] for j in range(M+1))/mp.fsum(vp[j]*vp[j]*w[j] for j in range(M+1))
    par='EVEN(v'"'"' )' if sym>0.9 else ('ODD(v'"'"')' if sym<-0.9 else f'mix({mp.nstr(sym,3)})')
    # how many sine coeffs dominate? (low-mode = smooth/polynomial-like)
    cs=[abs(Va[k,i]) for k in range(K)]; tot=mp.fsum(cs); cum=mp.mpf(0); ndom=0
    for k in range(K):
        cum+=cs[k]
        if cum>0.99*tot: ndom=k+1; break
    top=sorted(range(K),key=lambda k:-abs(Va[k,i]))[:5]
    print(f"\narch mode {i} (α={mp.nstr(Ea[i],6)}): v' parity {par}; {ndom} sine-modes hold 99% energy")
    print(f"   dominant sine indices(k+1): {[t+1 for t in top]}  coeffs {[mp.nstr(Va[t,i],4) for t in top]}")
print("\nREAD: few dominant low sine-modes + clean parity ⟹ negative modes are SMOOTH/low-degree ⟹ closed-form-attackable.")
