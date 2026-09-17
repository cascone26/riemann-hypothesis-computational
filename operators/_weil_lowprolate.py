import numpy as np
from scipy.linalg import eigh
import sys; sys.path.insert(0,"operators")
from coincidence_prolate import build

print("Weil-form MIN eigenvector: cumulative overlap with top-k MOST-CONCENTRATED prolate modes")
print("(CCM/Connes IR claim: the Weil-min-vector k_λ is a LOW-prolate object -> overlap saturates fast)")
for (a,b,fourW,N) in [(14,44,4.0,13),(14,44,5.0,17),(14,54,5.0,21),(14,64,6.0,27)]:
    archw,primew,Cw = build(a,b,fourW,N)
    tau = archw - primew
    ev,V = eigh(tau); w = V[:,0]           # min eigenvector (the razor / zero-interpolating dir)
    cev,cV = eigh(Cw)                        # concentration; high eigenvalue = most concentrated = low prolate
    order = np.argsort(cev)[::-1]           # most-concentrated first
    cum=[]
    for k in [1,2,3,5,8]:
        P = cV[:,order[:k]]
        cum.append(float(np.sum((P.T@w)**2)))
    print(f"[{a},{b}] N={N}: overlap² top-k = "
          + " ".join(f"k{k}:{c:.3f}" for k,c in zip([1,2,3,5,8],cum))
          + f"  (μ_min={ev[0]:+.2e})", flush=True)
print("\nread: if k3 overlap →~1, the Weil-min-vector IS a low-prolate object (CCM k_λ IR claim supported);")
print("      if it stays ~0.5 like §906's plunge measure, the low-prolate approximation is only partial.")
