#!/usr/bin/env python3
"""§184 Yoshida break-mining: extend arch_floor_vs_full to large c to characterize the archimedean floor
μ_arch (smallest POSITIVE pivot of A+B) asymptotically — is it bounded below by μ_0>0 as c→∞ (a robust
structural fact = candidate unconditional sub-lemma), or does it →0? A+B is INDEFINITE (negative boundary
modes n=±1,±2), so μ_arch is the smallest positive pivot, not a min-eigenvalue; P both definitizes A+B into
PD τ AND critically cancels the smallest pivot to e^{−4πc}. Reuses build_split + min_pivot from the original
tool (additive; original untouched). prec scaled to resolve the collapsing full-τ margin."""
import sys, os, json, math
sys.path.insert(0, os.path.dirname(__file__))
from arch_floor_vs_full import build_split, min_pivot

N = 64
cs = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else [300, 400, 500, 700]
print(f"§184 arch floor asymptotics  N={N}\n{'c':>5} {'prec':>7} {'mu_arch(A+B floor)':>20} {'n_neg(A+B)':>11} {'margin_full':>13} {'floor/margin':>13}")
out = []
for c in cs:
    prec = max(4000, int(5.46 * c * 3.33) + 2500)   # bits > log2(10)*|log10 margin| + headroom
    Aa, Af, DIM = build_split(c, N, prec)
    ma, na, ua = min_pivot(Aa, DIM)
    mf, nf, uf = min_pivot(Af, DIM)
    ratio = (ma / mf) if (ma and mf and mf > 0) else float('nan')
    print(f"{c:>5} {prec:>7} {('%.4e'%ma) if ma else 'n/a':>20} {na:>11} {('%.4e'%mf) if mf else 'n/a':>13} {ratio:>13.3e}", flush=True)
    out.append(dict(c=c, prec=prec, mu_arch=ma, n_neg_AB=na, margin_full=mf, floor_over_margin=ratio))
json.dump(out, open("results/arch_floor_ext.json", "w"), indent=2)
mus = [r["mu_arch"] for r in out if r["mu_arch"]]
print(f"\nREAD: mu_arch over c={cs}: min={min(mus):.4e} max={max(mus):.4e}. Bounded below (≥~0.05) ⇒ robust floor,")
print("candidate unconditional sub-lemma; the RH content concentrates in P definitizing+critically-cancelling A+B.")
