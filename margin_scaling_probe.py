#!/usr/bin/env python3
"""§972c (orchestrator's own non-colliding probe): how does the POLE-ABSORPTION MARGIN scale with resolution?

The exact Suzuki form (§954): Q=arch+pole−prime has κ=0; Q−pole has κ=1 (the single negative, min-eig μ<0);
the pole's positive part lifts it, margin = (pole quadratic form on the negative eigenvector) + μ ≥ 0 (absorbed).
At a=1.6,K=100 the margin was +0.028 — a HAIR above zero. The RH crux: as resolution K grows (and across the
half-width a), does this margin stay comfortably positive, or does it collapse toward 0⁺ (RH true but critical)?
A margin bounded away from 0 = comfortable positivity; margin→0⁺ = the razor is genuinely marginal. HANDLE."""
import numpy as np
import sys
sys.path.insert(0, 'operators')
from suzuki_exact_Q import build_Q

def analyze(a, K):
    Q, P = build_Q(a, K, return_pieces=True)
    Q = 0.5*(Q+Q.T); pole = 0.5*(P["pole"]+P["pole"].T)
    ev_full = np.linalg.eigvalsh(Q)
    Qmp = Q - pole
    w, V = np.linalg.eigh(0.5*(Qmp+Qmp.T))
    mu = w[0]; wvec = V[:, 0]                       # the single negative direction of Q−pole
    lift = wvec @ pole @ wvec                        # pole's quadratic form on it
    margin = lift + mu                               # ≥0 ⟺ pole absorbs the negative square
    # normalized margin (vs the size of the negative it must cover)
    rel = margin/abs(mu) if mu != 0 else float('nan')
    kneg_full = int((ev_full < -1e-9*max(1,abs(ev_full).max())).sum())
    return kneg_full, mu, lift, margin, rel

print("="*96)
print("§972c pole-absorption MARGIN scaling (exact Suzuki form). margin=lift+μ≥0 ⟺ pole absorbs the negative.")
print("="*96)
print(f"{'a':>4} {'K':>4} | {'κ(Q)':>5} {'μ=min(Q−pole)':>14} {'pole-lift':>11} {'margin':>11} {'margin/|μ|':>10}")
print("-"*96)
for a in [0.8, 1.2, 1.6, 2.0]:
    for K in [60, 100, 160, 220]:
        try:
            kn, mu, lift, margin, rel = analyze(a, K)
            flag = "" if margin >= 0 else "  <-- NEGATIVE (RH-violating!)"
            print(f"{a:>4.1f} {K:>4} | {kn:>5} {mu:>+14.4e} {lift:>+11.4e} {margin:>+11.4e} {rel:>+10.4f}{flag}")
        except Exception as e:
            print(f"{a:>4.1f} {K:>4} | ERROR {e}")
    print("-"*96)
print("READ: if margin stays > 0 and margin/|μ| is bounded away from 0 as K grows, the pole COMFORTABLY absorbs")
print("the one negative square (robust positivity). If margin/|μ|→0⁺ with K, the absorption is CRITICAL — the")
print("pole barely covers the negative, which is the exact quantitative face of 'RH is true but marginally'.")
print("Either way it quantifies the §891 razor in the compression-free exact form. Honest: report the trend.")
print("="*96)
