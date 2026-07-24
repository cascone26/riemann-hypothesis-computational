#!/usr/bin/env python3
"""§201c: RIGOROUS coupling bound engine. Analytic claim (from the D_kl closed form, k fixed, l->inf, odd):
   D_kl(t) -> -2kβ·sin(lβt)   [subleading O(1/l)]
hence  Q[k,l]=∫_0^{2a} g(t)D_kl(t)dt ~ -2kβ·∫_0^{2a} g sin(lβt)dt
      = -2kβ·{[-g cos(lβt)/(lβ)]_0^{2a} + O(1/l^2)}   (I.B.P.; lβ·2a=lπ ⟹ cos=(-1)^l)
      = -2k·[g(0) - (-1)^l g(2a)]/l  + O(1/l^2).
For ODD l (even sector): Q[k,l] ~ -2k[g(0)+g(2a)]/l. VERIFY both: (A) D_kl vs -2kβ sin(lβt) pointwise;
(B) l·Q[k,l] -> const = -2k[g(0)+g(2a)]. Get g(0),g(2a) so the coupling constant is EXPLICIT."""
import sys; sys.path.insert(0,'/Users/scones/projects/riemann/operators/dbn_dynamics')
import mpmath as mp
from dkl_closedform import Dcf as D_kl
from tail_structure import Q, gval   # gval(t) = screw kernel value if present
a=mp.mpf('0.549'); b=mp.pi/(2*a)
mp.mp.dps=40
# (A) pointwise D_kl vs asymptote
print("(A) D_kl(t) vs -2kβ sin(lβt)  [k=1, t=0.4]:")
t=mp.mpf('0.4')
for l in [11,21,41,81,161]:
    d=D_kl(1,l,t,a); asy=-2*1*b*mp.sin(l*b*t)
    print(f"   l={l:>4} D={mp.nstr(d,6)}  asy={mp.nstr(asy,6)}  ratio={mp.nstr(d/asy,6)}")
# g(0), g(2a)
try:
    g0=gval(mp.mpf('1e-9')); g2a=gval(2*a-mp.mpf('1e-9'))
    print(f"\n g(0+)={mp.nstr(g0,6)}  g(2a-)={mp.nstr(g2a,6)}  g0+g2a={mp.nstr(g0+g2a,6)}")
    print(f" predicted l·Q[1,l] -> -2[g0+g2a] = {mp.nstr(-2*(g0+g2a),6)}")
except Exception as e:
    print(f"\n gval unavailable ({e}); infer const from l·Q directly")
# (B) l·Q[k,l] -> const
print("\n(B) l·Q[1,l] and l·Q[3,l] (odd l) -> const?  (=coupling constant per low mode):")
for l in [21,41,81,121]:
    print(f"   l={l:>4}  l·Q[1,l]={mp.nstr(l*Q(1,l),6)}   l·Q[3,l]={mp.nstr(l*Q(3,l),6)}",flush=True)
print("\nREAD: (A) ratio->1 confirms D_kl~-2kβ sin(lβt). (B) l·Q[k,l]->const ⟹ |Q[k,l]|<=C_k/l, C_k~2k[g0+g2a]")
print("EXPLICIT. Then Σ_{l>K0}|Q[k,l]|^2 <= C_k^2/(2K0) ⟹ Schur correction bound(K0) ⟹ pick K0. Tail lemma constants LOCKED.")
