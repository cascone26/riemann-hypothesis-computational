import numpy as np, sys, os
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
os.chdir('/Users/scones/projects/riemann/operators')
from bs_schur_threshold import analyze
import math
# prime thresholds a=0.5 log p
prs=[2,3,5,7,11,13]; thr={0.5*math.log(p):p for p in prs}
print("a       thr    lambda_a     LHS-1      beta      margin M=LHS-RHS")
agrid=np.arange(0.36,1.30,0.03)
for a in agrid:
    M=min(2600,max(1500,int(400/a)))
    la,mu1,p1sq,beta,s,L,R=analyze(a,M=M,K=80)
    # nearest threshold marker
    near=min(thr.keys(), key=lambda t:abs(t-a)); mark=f"~p{thr[near]}" if abs(near-a)<0.02 else ""
    print(f"{a:5.3f}  {mark:>5}  {la:10.3e}  {L-1:9.4e}  {beta:9.4e}  {L-R:11.4e}",flush=True)
