import sys,os
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
import mpmath as mp
mp.mp.dps=18
from dh_arch_numeric import make_psiR, zeta_closed_WR
c=14; N=6; L=float(mp.log(c))
psiR=make_psiR(L, mp.mpf(1)/4, 1)
h=mp.mpf(10)**(-5)
print("DIAGONAL (tail-corrected, dps=18): numeric psi_R'(n) vs Groskin closed diag", flush=True)
for n in [1,2,3]:
    d=float((psiR(n+h)-psiR(n-h))/(2*h))
    clo=zeta_closed_WR(n,n,N,c)
    print(f"  n={n}: numeric={d:+.6e} closed={clo:+.6e} rel.err={abs(d-clo)/abs(clo):.2e}", flush=True)
print("OFF-DIAG:", flush=True)
vals={k:psiR(k) for k in range(4)}
for (m,n) in [(1,0),(2,1),(3,1)]:
    num=float((vals[m]-vals[n])/(m-n)); clo=zeta_closed_WR(n,m,N,c)
    print(f"  (m={m},n={n}): numeric={num:+.6e} closed={clo:+.6e} rel.err={abs(num-clo)/(abs(clo)+1e-30):.2e}", flush=True)
