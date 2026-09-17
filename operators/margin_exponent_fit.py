import numpy as np, sys, os, math
sys.path.insert(0,'/Users/scones/projects/riemann/operators'); os.chdir('/Users/scones/projects/riemann/operators')
from bs_schur_threshold import analyze
# evaluate M(a) at PRIME THRESHOLDS a=0.5 log p (cleaner), two K to check convergence, then fit exponent
prs=[3,5,7,11,13,17,19,23,29,37,53,73,101]
ath=[0.5*math.log(p) for p in prs]
rows=[]
print(f"{'p':>4} {'a':>6} {'M(K=100)':>11} {'M(K=140)':>11} {'rel.diff':>9}")
for p,a in zip(prs,ath):
    Mgrid=min(3000,max(1600,int(500/a)))
    _,_,_,_,_,L1,R1=analyze(a,M=Mgrid,K=100); m1=L1-R1
    _,_,_,_,_,L2,R2=analyze(a,M=Mgrid,K=140); m2=L2-R2
    rd=abs(m1-m2)/max(abs(m2),1e-30)
    rows.append((a,m2))
    print(f"{p:>4} {a:6.3f} {m1:11.4e} {m2:11.4e} {rd:9.2%}",flush=True)
# fit log M = log C - p log a  (use K=140 values, converged region a>0.7)
A=np.array([a for a,m in rows]); Mv=np.array([m for a,m in rows])
good=(A>0.7)&(Mv>0)
if good.sum()>=3:
    slope,inter=np.polyfit(np.log(A[good]), np.log(Mv[good]), 1)
    print(f"\nFIT (a>0.7, K=140): M(a) ~ {math.exp(inter):.3e} * a^({slope:.3f})   [exponent p={-slope:.3f}]")
    # also fit M*a^2 trend
    print(f"M*a^2 values: {np.round(Mv[good]*A[good]**2*1e6,3).tolist()} (x1e-6) -- flat if p=2")
