#!/usr/bin/env python3
"""§201b: the two exponents for the rigorous tail lemma — (1) diagonal growth Q[k,k] vs log(k); (2) low-mode
coupling decay |Q[1,l]|,|Q[3,l]| vs l (determines Schur correction & K0). Reuse tail_structure engine."""
import sys; sys.path.insert(0,'/Users/scones/projects/riemann/operators/dbn_dynamics')
from tail_structure import Q
import mpmath as mp
print("(1) diagonal vs log k:  Q[k,k] / log(k*b) -> const?  (arch symbol Omega(r)~log r prediction)")
b=mp.pi/(2*mp.mpf('0.549'))
for k in [11,21,31,41,51]:
    print(f"   k={k:>3} Q[k,k]={mp.nstr(Q(k,k),6)}  Q/log(k*b)={mp.nstr(Q(k,k)/mp.log(k*b),5)}")
print("\n(2) low-mode coupling decay |Q[1,l]| and |Q[3,l]| vs l (fit l^-p):")
prev1=prev3=None
for l in [9,15,21,31,41,51]:
    q1=abs(Q(1,l)); q3=abs(Q(3,l))
    p1=p3=''
    if prev1: p1=mp.log(prev1[1]/q1)/mp.log(mp.mpf(l)/prev1[0]); p3=mp.log(prev3[1]/q3)/mp.log(mp.mpf(l)/prev3[0])
    print(f"   l={l:>3} |Q[1,l]|={mp.nstr(q1,4)} (p~{mp.nstr(p1,3) if p1 else '-'})  |Q[3,l]|={mp.nstr(q3,4)} (p~{mp.nstr(p3,3) if p3 else '-'})",flush=True)
    prev1=(mp.mpf(l),q1); prev3=(mp.mpf(l),q3)
print("\nREAD: (1) const ⟹ Q[k,k]~C log k (rigorous via arch symbol). (2) p>=2 ⟹ Σ_{l>K0}|Q[1,l]|^2 ~ K0^{-(2p-1)}")
print("tiny ⟹ Schur correction to the near-null << margin at modest K0. Both feed the explicit-constant tail lemma.")
