#!/usr/bin/env python3
"""p⊕r assembly item (1): the JOINT probe-span certificate in PENCIL (L²) normalization (§802).
Zero-free Weil Gram Q on {14,20,26,32,38,44} at 4W=1.32 (the §788-addendum span, corner+complement
together per §801's p⊕r form) vs the L² overlap Gram G (closed form, §801-verified):
certified generalized λ_min(Q, G) — the class-metric margin Gauss's normalization rule requires.
Hypothesis set: ∅ (zero-free route; G closed form)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx
import arb_fejer_cert2 as ac
import arb_zerofree_cert3 as zc

fourW = 1.32; RC = 400
ALL = [14, 20, 26, 32, 38, 44]
ctx.prec = 220
ac.W = arb(fourW)/4; ac.TW = 2*ac.W; ac.c_tri = arb.pi()/ac.W
cs = [arb(c) for c in ALL]; n = len(cs)
print(f"joint pencil cert: Q and G on {ALL}, 4W={fourW}...", flush=True)
Q = [[None]*n for _ in range(n)]
for j in range(n):
    for k in range(j, n):
        Q[j][k] = Q[k][j] = ac.pole_entry(cs[j], cs[k]) + zc.arch_entry_v3(cs[j], cs[k], arb(RC)) - ac.prime_entry(cs[j], cs[k])
    print(f"  Q row {j} done", flush=True)
G = [[None]*n for _ in range(n)]
for j in range(n):
    for k in range(j, n):
        G[j][k] = G[k][j] = 2*arb.pi()*ac.gpair_closed_arb(cs[j], cs[k], arb(0))
def ok(lam):
    d = [None]*n; L = [[arb(0)]*n for _ in range(n)]
    for i in range(n):
        sd = Q[i][i] - lam*G[i][i]
        for k in range(i): sd -= L[i][k]*L[i][k]*d[k]
        if not (sd > 0): return False
        d[i] = sd
        for j2 in range(i+1, n):
            v = Q[j2][i] - lam*G[j2][i]
            for k in range(i): v -= L[j2][k]*L[i][k]*d[k]
            L[j2][i] = v/d[i]
    return True
if not ok(arb(0)):
    print("Q not PD (should not happen — §788 certified)"); sys.exit(1)
lo, hi = 0.0, 2.0
for _ in range(40):
    mid = (lo+hi)/2
    if ok(arb(mid)): lo = mid
    else: hi = mid
print(f"\n===== JOINT PENCIL CERT (§802) =====")
print(f"Q(ψ) ≥ {lo:.6e} · ‖ψ‖²_L2   for ALL ψ in span{{F_14..F_44}} at the open rung  [hypothesis: ∅]")
print(f"(this is the p-block margin of the p⊕r anchor assembly, class-metric, per §801)")
