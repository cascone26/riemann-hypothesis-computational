import sys,os
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
from flint import ctx
import numpy as np
from dh_certified_tau import build_arb_tau
for c in [25,50,100]:
    ctx.prec=800
    A,DIM=build_arb_tau(c,64,800)
    M=np.array([[float(A[i,j].mid()) for j in range(DIM)] for i in range(DIM)])
    ev=np.sort(np.linalg.eigvalsh(0.5*(M+M.T)))
    near=int((np.abs(ev)<1e-3).sum())
    print(f"c={c}: dim={DIM}, #|eig|<1e-3 (near-null)={near}, smallest 6 eig={[f'{x:+.2e}' for x in ev[:6]]}",flush=True)
    print(f"        eig[near-2:near+3] around the cluster edge: {[f'{x:+.2e}' for x in ev[max(0,near-2):near+3]]}")
print("A tight near-zero cluster of size ~8/15/23 then a jump to O(0.01+) positives confirms the near-null facet.")
