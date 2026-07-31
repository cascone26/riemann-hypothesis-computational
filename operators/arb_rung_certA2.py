#!/usr/bin/env python3
"""CERT A v2 — rigorous route via the explicit-formula identity (no arch integrals).
CLAIM: at 4W=1.32 (support past log3), the p2-only Weil form (n=3 term removed) is INDEFINITE.
IDENTITY: full(=pole+arch-prime_all) = zero-side = 2 sum_gamma psi(gamma)^2  (Guinand-Weil, theorem).
Hence  Q_no3(v) = 2 sum_gamma psi(gamma)^2 + 2*(Lambda(3)/sqrt(3))*ghat_vv(log 3)
(only n=3 matters: log 9 = 2.197 > 4W = 1.32 excludes higher powers).
Both terms rigorously computable WITHOUT integrals:
  - ghat_vv(log3) = sum_jk v_j v_k gpair(c_j,c_k,log3): CLOSED FORM (arb_fejer_cert2.gpair_closed_arb).
  - zero side: psi(r) = sum v_j F_j(r), Fejer closed form; sum over zeros gamma <= G with mpmath zeros
    (dps 30, ball radius 1e-18 conservative -- provenance: mpmath zetazero, cross-checked literature
    tables; all gamma < 1000 are certified on the line by Platt-Trudgian 3e12), plus a rigorous tail:
    for gamma > G >= 2*44: F_j(gamma) <= 2/(W(gamma-c_j))^2 <= 2/(W(gamma-44))^2, so
    psi(gamma)^2 <= (sum|v_j|)^2 * 4/(W(gamma-44))^4, and with the zero-count bound
    N(t+1)-N(t) <= 2 ln t (crude, valid t>=2 -- e.g. Backlund/Trudgian give ln t; use 2 ln t safe),
    tail <= sum_{m>=G} 2 ln(m+1) * (sum|v_j|)^2 * 4/(W(m-44))^4  (rigorous arb sum + integral bound).
RESULT: rigorous enclosure of Q_no3(v). Negative upper bound => CERTIFIED indefinite.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, ctx
import arb_fejer_cert2 as ac
import l2_fejer as fj
import mpmath as mp

fourW = 1.32; S = 2.0
ctx.prec = 200

# --- float witness on the same probe set ---
fj.W = fourW/4.0
cs_f = fj.probes((14.0, 44.0), S)
scale = np.ones_like(fj.LAMW); scale[np.isin(fj.NS, [3, 9, 27, 81])] = 0.0
fmat = fj.build(cs_f, lam_scale=scale)
FU = 0.5*(fmat['full']+fmat['full'].T)
ev, V = np.linalg.eigh(FU)
print(f"float: n={len(cs_f)} probes (s={S}), min eig(no-3) = {ev[0]:.4e}")
v = V[:, 0]

# --- arb setup ---
ac.W = arb(fourW)/4; ac.TW = 2*ac.W; ac.c_tri = arb.pi()/ac.W
Wa = ac.W
cs = [arb(float(c)) for c in cs_f]
va = [arb(float(x)) for x in v]
n = len(cs)

# term 1: 2*(Lambda(3)/sqrt3)*sum_jk v_j v_k gpair(cj,ck,log3)  [closed form]
u3 = arb(3).log()
g3 = arb(0)
for j in range(n):
    for k in range(n):
        g3 = g3 + va[j]*va[k]*ac.gpair_closed_arb(cs[j], cs[k], u3)
term_prime3 = 2*arb(3).log()/arb(3).sqrt()*g3
print(f"prime-3 term (closed form): {term_prime3}")

# term 2: zero side 2 sum psi(gamma)^2, gamma <= G, + rigorous tail
G = 800
ZC = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'zeros500.npy'))
mp.mp.dps = 30
def fej_arb(c, r):
    out = arb(0)
    for cc in (c, -c):
        y = Wa*(r-cc)
        # sinc^2 with removable singularity: y is arb; use arb sinc via sin(y)/y guarded
        if y.contains(arb(0)):
            out = out + arb(1).union((y.sin()/y)**2 if not (y.mid()==0 and y.rad()==0) else arb(1))
        else:
            out = out + (y.sin()/y)**2
    return out

zero_sum = arb(0); nz = 0
for g in ZC:
    if g > G: break
    gb = arb(float(g), 1e-9)    # ball around cached mpmath zero (conservative radius; PT on-line cert)
    psi = arb(0)
    for j in range(n):
        psi = psi + va[j]*fej_arb(cs[j], gb)
    zero_sum = zero_sum + psi*psi
    nz += 1
term_zero = 2*zero_sum
print(f"zero side (gamma<={G}, {nz} zeros): {term_zero}")

# rigorous tail (sharpened): psi(gamma) <= sum|v_j| * 2/(W(gamma-44))^2  (fej(y)<=1/y^2, both mirror
# terms bounded by the near one), so psi^2 <= |v|_1^2 * 4/(W^4 (gamma-44)^4).
# Zero density: N(t+1)-N(t) <= 2 ln t for t >= 100 (safe cover of (1/2pi)ln(t/2pi)+O; Backlund/Trudgian
# explicit bounds give < ln t here, factor-2 safety). Then
#   sum_{gamma>G} (gamma-44)^{-4} <= int_G^inf 2 ln x (x-44)^{-4} dx <= 2 ln'(..)-free bound:
#   <= 2 [ ln G /(3 (G-44)^3) + int_G^inf (x-44)^{-4} dx/ (x ln-deriv) ] <= 2 (ln G + 1/2)/(3 (G-44)^3).
sv = arb(0)
for j in range(n): sv = sv + abs(va[j])
Gm = arb(G-44)
tail = sv*sv * 4/(Wa**4) * 2*(arb(G).log() + arb(1)/2)/(3*Gm**3)
term_tail = 2*arb(0).union(tail)
print(f"|v|_1 = {sv}")
print(f"rigorous zero-tail bound (gamma>{G}): +/- {2*tail}")

Q = term_zero + term_tail + term_prime3
print(f"\nRIGOROUS Q_no3(v) = {Q}")
if Q < 0:
    print("CERT A v2 PASSED: Q_no3(v) < 0 RIGOROUSLY => the p2-only Weil form is INDEFINITE past log3.")
    print("(Provenance: Guinand-Weil identity [theorem]; cached zeros gamma<=800 [mpmath, ball r=1e-9,")
    print(" on-line certified by Platt-Trudgian 3e12]; closed-form arb prime term; rigorous density tail.)")
    print(f"explicit bounds: Q in [{float(Q.mid())-float(Q.rad()):.3e}, {float(Q.mid())+float(Q.rad()):.3e}]")
elif Q > 0:
    print("Q > 0 rigorously -- float witness not confirmed; investigate (identity check or witness quality).")
else:
    print("INCONCLUSIVE: enclosure straddles 0 (raise prec / tighten zero balls / larger G).")
