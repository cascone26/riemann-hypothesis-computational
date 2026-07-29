#!/usr/bin/env python3
"""§423 (ATTACK, second-space construction, step 1): build the PLACE-INDEXED cross-prime pairing directly.
The explicit formula's archimedean term pairs two prime-frequencies via M_arch[p,q] = ∫ Ω(r) cos(r log p)
cos(r log q) dr, Ω(r)=Re ψ(1/4+ir/2)−log π. This is the arithmetic-surface intersection number of the fibers
at p and q (archimedean/Arakelov contribution) — a REAL, derivable object (Handle-earning). Question: does the
place-indexed form (diagonal = prime self-energies, off-diagonal = M_arch cross terms) carry a HODGE-INDEX
signature (1, n−1) — the structural positivity the analytic instrument (§414) LACKED? If the raw archimedean
place-form is (1,n−1) or (k,1,·), the second-space Hodge structure lives HERE (unlike the probe-space §414).
Then: what off-diagonal correction (the §332 ⟨δ_p,δ_q⟩ / a twin-prime 𝔖 term) is needed to pin signature exactly?"""
import numpy as np
from scipy.special import digamma
from scipy import integrate
PRIMES=np.array([2,3,5,7,11,13,17,19,23,29])
def Om(r): return np.real(digamma(0.25+0.5j*r))-np.log(np.pi)
def Mcos(a,b,wgt):  # ∫ wgt(r) cos(a r) cos(b r) dr over a damped range (Gaussian reg to converge)
    f=lambda r: wgt(r)*np.cos(a*r)*np.cos(b*r)*np.exp(-(r/40.0)**2)
    v,_=integrate.quad(f,0,120,limit=400); return 2*v
n=len(PRIMES); L=np.log(PRIMES)
M=np.zeros((n,n))
for i in range(n):
    for j in range(i,n):
        M[i,j]=M[j,i]=Mcos(L[i],L[j],Om)
ev=np.linalg.eigvalsh(M)
def sig(e,tol=1e-9):
    m=max(1,abs(e).max()); return int((e>tol*m).sum()),int((abs(e)<=tol*m).sum()),int((e<-tol*m).sum())
print(f"§423 archimedean place-form M_arch[p,q]=∫Ω cos(log p·r)cos(log q·r)dr, primes={list(PRIMES)}")
print(f"eigenvalues: {np.array2string(ev,precision=3,max_line_width=200)}")
print(f"signature (+,0,−) = {sig(ev)}")
print(f"diag(M) = {np.array2string(np.diag(M),precision=3)}")
# compare: is it dominated by 1 positive (ample-like) direction? ratio top/second
print(f"top eig / |2nd|: {ev[-1]/abs(ev[-2]) if abs(ev[-2])>1e-12 else np.inf:.2f}  | most-negative: {ev[0]:+.4f}")
print("""
reading: signature (1,·,n−1) [one +, rest −] => the archimedean place-form IS a Hodge-index intersection form
  (ample anchor = the + direction) — the structure §414 couldn't find in probe-space lives in PLACE-space. Then
  Weil positivity of the FULL form (with primes) = Hodge index on this lattice. signature mixed/other => the
  raw archimedean pairing isn't yet Hodge; the §332 cross-prime ⟨δ_p,δ_q⟩ correction is what supplies it — solve
  for it next.""")
