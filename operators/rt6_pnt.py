#!/usr/bin/env python3
"""§416 (ATTACK): decompose prime_push (§415c) into PNT-SMOOTH vs ARITHMETIC-FLUCTUATION.
prime term = Σ_{prime powers n} Λ(n)/√n · g_pair(c_j,c_k; log n). PNT says the von Mangoldt measure
Λ(n)dn → dx (density 1). Build prime_SMOOTH = ∫ (1/√x) g_pair(c_j,c_k; log x) dx (x≥2) — the 'expected'
prime term for ANY density-1 measure (no arithmetic). Then:
  push_true   = λmin(archpole) − λmin(archpole − prime_true)
  push_smooth = λmin(archpole) − λmin(archpole − prime_smooth)
If push_smooth << archmargin (smooth is safely covered) and push_true → archmargin exactly at the edge, the
RH-content is the ARITHMETIC EXCESS push_true − push_smooth = the fluctuation Σ(Λ(n)−1)... — quantified."""
import numpy as np
from scipy import integrate
import l2_edge_sweep as L
BAND=(14.0,44.0); w=3.0
def g_pair_scalar(cj,ck,w,logn):  # matches L.g_pair signature (array logn) — wrap for scalar
    return L.g_pair(cj,ck,w,np.array([logn]))[0]
def prime_smooth_mat(cs,w):
    n=len(cs); M=np.zeros((n,n))
    for j in range(n):
        for k in range(j,n):
            f=lambda x: (1.0/np.sqrt(x))*L.g_pair(cs[j],cs[k],w,np.log(x))  # measure dx, weight 1/√x
            val,_=integrate.quad(f,2.0,np.exp(np.log(max(L.NS))+1),limit=200)
            M[j,k]=M[k,j]=2.0*val
    return M
print(f"§416 PNT-smooth vs arithmetic prime push. band{BAND} w={w}.")
print(f"{'s':>5} {'n':>3} {'archmargin':>11} {'push_true':>10} {'push_smooth':>11} {'arith excess':>12} {'true/arch':>9}")
for s in np.arange(6.0,2.9,-0.4):
    cs=L.probes(BAND,s)
    if len(cs)<4: continue
    f=L.build_forms(cs,w); AP=f['archpole']; PT=f['prime']; PS=prime_smooth_mat(cs,w)
    la=float(np.linalg.eigvalsh(AP)[0])
    pt=la-float(np.linalg.eigvalsh(AP-PT)[0])
    ps=la-float(np.linalg.eigvalsh(AP-PS)[0])
    print(f"{s:5.1f} {len(cs):3d} {la:+11.4e} {pt:+10.4e} {ps:+11.4e} {pt-ps:+12.4e} {pt/la if la>1e-15 else np.nan:9.3f}",flush=True)
print("""
reading: push_smooth < archmargin comfortably => generic density is always safe; the RH-risk is ONLY the
  arithmetic excess (pt−ps). If (pt−ps) is what closes the gap to archmargin at the edge, RH-on-this-band =
  'the arithmetic fluctuation excess never exceeds archmargin − push_smooth' = a sharpened, smaller target
  than raw Weil positivity (the PNT main term is already paid for). If push_smooth ALSO → archmargin, the
  degeneracy is density-generic and arithmetic is a rider (matches §411 random-jitter-safe).
""")
