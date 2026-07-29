#!/usr/bin/env python3
"""§402 (angle 10): Krein/Schur one-probe-at-a-time. B_full PSD <=> every successive Schur complement
s_k = B[k,k] - B[k,<k] inv(B[<k,<k]) B[<k,k] > 0 (Verblunsky/orthogonal-polynomial recursion wrt the
zero measure). Track {s_k} across the grid: which step is BINDING (s_k->0 first), and what is its
arithmetic content? If the binding step is a LOW-order scalar inequality, the corner decomposes into
weaker pieces. Compare the s_k sequence for TRUE vs bent primes (does one step carry the sign-sensitivity?)."""
import numpy as np
import l2_edge_sweep as L
w = 3.0
def schur_seq(B):
    n = B.shape[0]; s = []
    for k in range(n):
        if k == 0: s.append(B[0,0]); continue
        A = B[:k,:k]; b = B[:k,k]
        s.append(float(B[k,k] - b @ np.linalg.solve(A, b)))
    return np.array(s)
# bent primes (Lambda(2)x1.05)
scale = np.ones_like(L.LAMW)
for i,nn in enumerate(L.NS):
    m=nn
    while m%2==0: m//=2
    if m==1: scale[i]=1.05
print("Schur-complement sequence s_k (each must be >0 for PSD). band [14,44], w=3, order = ascending c.")
print(f"{'s':>5} {'n':>3} {'min s_k':>10} {'argmin':>7} {'last s_k':>10} {'min s_k(bent)':>13} {'bent argmin':>11}")
for s in np.arange(6.0, 2.9, -0.3):
    cs = L.probes(L.BAND, s)
    f = L.build_forms(cs, w); fb = L.build_forms(cs, w, lam_scale=scale)
    sk = schur_seq(f['full']); skb = schur_seq(fb['full'])
    print(f"{s:5.1f} {len(cs):3d} {sk.min():+10.3e} {int(sk.argmin()):7d} {sk[-1]:+10.3e} "
          f"{skb.min():+13.3e} {int(skb.argmin()):11d}", flush=True)
print("\ndetail at near-critical s=3.6 (the binding-step arithmetic):")
cs = L.probes(L.BAND, 3.6); f = L.build_forms(cs, w)
sk = schur_seq(f['full'])
print("  s_k =", " ".join(f"{x:+.2e}" for x in sk))
print("\nreading: if min s_k is at a LOW index k (early step binds) => corner reduces to a small-block")
print("scalar inequality (weaker target). If bent primes flip a SPECIFIC s_k negative while others stay")
print("positive => that Schur step CARRIES the arithmetic sign-sensitivity = the one inequality to attack.")
