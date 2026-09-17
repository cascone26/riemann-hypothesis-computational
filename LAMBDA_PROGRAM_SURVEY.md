# de Bruijn–Newman Λ program — frontier survey (2026-08-10)

Deep-research harness (97 agents, 15 sources fetched, 25 claims 3-vote adversarially verified, 21 confirmed).
Commissioned to test whether the Λ program is a good compute-heavy lane toward RH. **Honest conclusion: no.**

## The facts (all HIGH confidence, cited)
- **RH ⟺ Λ ≤ 0.** de Bruijn (1950), Newman (1976). Λ = infimum of t such that all zeros of the heat flow
  H_t(z)=∫₀^∞ e^{tu²}Φ(u)cos(zu)du are real.
- **Λ ≥ 0** — Rodgers–Tao, arXiv:1801.05914, Forum of Math Pi 2020 (settled Newman's conjecture).
- **Λ ≤ 0.22** — Polymath15, arXiv:1904.12438, Res. Math. Sci. 2019 (improved prior Λ<1/2 of Ki–Kim–Lee).
- **Λ ≤ 0.2** — Platt–Trudgian, arXiv:2004.09765 (2020) via extended RH verification to higher height.
- So the rigorous interval is **0 ≤ Λ ≤ 0.2**. No improvement published 2020→2025 (5 years).

## Method (Polymath15)
Three components: (1) numerical zero-free regions for H_t(x+iy) in explicit rectangles; (2) rigorous
large-argument asymptotics of H_t; (3) the monotonicity **½σ_max(t)² + t is non-increasing** (∂_t(½σ_max²)≤−1
while σ_max>0), giving **Λ ≤ t + ½σ_max(t)²**. σ_max(t)=sup of imaginary parts of non-real H_t zeros, controlled
by verifying RH-zeros are on the line up to height T. Verify RH to higher T ⟹ take smaller t ⟹ push the bound down.

## Why this lane is NOT a good bet for us (the honest kill)
1. **It structurally cannot prove RH.** The method only pushes the UPPER bound down from above (0.22→0.2→…). Reaching
   Λ ≤ 0 (= RH) would require RH verified to *infinite* height. It approaches a wall it cannot cross by construction.
2. **The bottleneck is specialized RH-verification** (zeros on the line to height T) — a domain owned by dedicated
   teams (Platt; Gourdon–Demichel; Odlyzko) with decades-optimized Riemann–Siegel + Odlyzko–Schönhage + rigorous
   interval arithmetic. Our general-purpose compute has **no edge** here; we'd be reimplementing specialists' work.
3. **Diminishing returns**: no improvement below 0.2 in 5 years. The sources genuinely disagree (split vote) on
   whether further progress is "just more compute" or hits obstacles; claims of clean exponential scaling were
   *refuted* in verification — i.e. even the cost-to-improve is uncertain, not a clear "throw compute at it" win.

## Bottom line
The Λ program is real, cited, and gradable — but it (a) cannot reach a proof, and (b) its compute bottleneck is a
specialist lane where we have no advantage. My earlier recommendation of it was too optimistic; the survey corrects
it. It is the single most compute-amenable frontal-ish RH lane, and even it does not give a real shot at a proof.
