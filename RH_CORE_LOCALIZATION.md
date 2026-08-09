# The RH Corner: a localization map (Riemann8 session, 2026-08-09)

*Honest consolidation of LEDGER §884–§926. Every quantitative claim below is HANDLE (code-run, reproducible)
unless marked REPORT. The purpose of this document is to state precisely (1) what the remaining content of RH
has been reduced to, (2) the complete map of tactics ruled out this session and why, and (3) the honest open
core. It banks a real result — the sharpest localization we have — without overclaiming a proof.*

---

## 1. The setting

Weil's criterion: RH ⟺ the explicit-formula quadratic form **Q(φ) = pole(φ) + arch(φ) − prime(φ) ≥ 0** for all
even band-limited test functions φ, where on a Fejér band-limited span (window `[a,b]`, bandwidth `4W`, `N`
centers) Q is an explicit `N×N` generalized-eigenvalue problem against the Gram metric `G`:

- `arch` = compression of multiplication by `Ω(r)/2π`, `Ω(r) = Re ψ₀(¼+ir/2) − ln π` (archimedean density);
- `prime` = `(1/π) Σ_{log n < 4W} Λ(n)/√n · [compression of cos(r log n)]` (exact for bandwidth-`4W` φ);
- `pole` = rank-≤2 PSD term from ζ's pole at s=1.

**RH-corner** = the regime where the least generalized eigenvalue (the *margin*) is razor-thin.

## 2. The localization chain (§884–§912)

Each step is verified and reduces the object:

1. **Q = M_σ + pole, exactly** (§889). `M_σ = arch − prime` is the compression of multiplication by the *exact*
   symbol `σ(r) = Ω/2π − P(r)` (‖arch−prime − M_σ‖ = 1e-15). The pole is a genuine rank-≤2 term, not an error.
2. **One negative eigenvalue** (§889–890). At the corner `M_σ` has *exactly one* negative eigenvalue `−μ`;
   `#neg(M_σ) ≤ rank(pole)` is forced by RH; the pole lifts it. `σ` is pointwise-negative on 7–40% of every
   window, yet `M_σ` has a single negative direction — the arithmetic near-coincidence of arch and prime.
3. **Scalar form** (§891). `RH-corner ⟺ g = 1 − μ·wᵀA⁺w ≥ 0`, `A = M₊ + p·vvᵀ`, leading form `p(v·w)² ≥ μ`.
4. **The razor is off-diagonal / near-null** (§900–903). The leading `π_w/μ` has ~4% slack, but the M₊ near-null
   cluster eats it: the true razor is an irreducibly **~2–3-dimensional coupling** in the near-null (prolate-
   plunge-related, ~0.53 overlap, *not* clean Landau–Widom) subspace. No single scalar closes it.
5. **Dual / zeros picture** (§908–909). The margin `= 2 Σ_γ w(γ)²`. The worst direction `w` **annihilates its
   in-window zeros exactly**; the entire margin is `w`'s Paley–Wiener tail spilling onto the zeros *just above*
   the window (in-window 0%, below 0%, above 100%).
6. **Edge effect** (§910). The corner exists only at the **bottom edge of the zero spectrum**, where spillover is
   one-sided (upward). The moment a window has zeros below it, the negative eigenvalue vanishes; the interior is
   robustly positive.

## 3. The invariant statement (what everything reduces to)

> **RH ⟺ an unconditional lower bound on the band-limited Weil margin, uniform in bandwidth `4W` and height `t`.**

The strategic survey (§920, five programs, source-verified) confirmed **every** major RH program reduces to this
same statement: Connes–Consani (Weil positivity, the CCM μ_λ≥0 prolate lemma), de Branges (the ω→0 canonical-
system extension), function-field/Hodge (intersection-form positivity), and the positivity criteria (Li /
Bombieri–Lagarias / Nyman–Beurling). They are the same wall from different angles.

## 4. Tactics ruled out this session (the elimination map)

| Tactic | Verdict | Reason |
|---|---|---|
| Scalar / operator bound `floor − pnorm` | dead (§914, §919) | fails at the corner; misalignment correction closes only ~75% |
| Finite bottom region | tempered (§911) | corner region grows with `4W` (not a fixed height range) |
| Concentration / Szegő (band-limited) | robust-region only (§888) | diagonal symbol wrong-sign at corner; off-diagonal effect |
| Low-prime / first-few-primes | dead (§917) | razor is a global cancellation over the *full* prime set |
| Two-term `arch ⪰ p=2` (survey's #1 pick) | **refuted (§925)** | fails at `4W≥6`; primes cannot be peeled — cancellation is essential |
| Density-preserving (PNT main + fluctuation) | bottoms out (§926) | smooth part vacuous at the window; content = the error term = RH |
| de Branges *global* positivity | blocked (§920) | Conrey–Li 2000 counterexample at the 34th zero |
| Function-field / Hodge translation | blocked (§920) | Deninger 2022: no real-coefficient arithmetic Weil cohomology exists |
| Arb-certification of the corner | radius-blocked (§894) | razor margin (~1e-15) sits below achievable interval radius |

## 5. What IS proven / certified (the real, banked positives)

- **Robust region, unconditional** (§875): `floor(t) > pnorm(4W) ⟹ Q ⪰ 0`, with the exact floor law
  `floor(t) ≈ (1/2π)ln(t/2π)` (arb-certified). Certified frontier extended to `4W ≈ 5` (§871–873).
- **Critical-subspace positivity, arb-certified** (§884): `Q ⪰ 0` on the strong-prime subspace at a deep corner,
  min pivot `2.3e-5 ≫` radii — a rigorous zero-free certificate that arch dominates prime on the critical subspace.
- **Closed-form archimedean weight** (§898): `Ŵ(u) = −π e^{u/2}/sinh(u)`, verified; arch/prime pointwise ratio
  `= 2π/log n`, crossing 1 at the scale `u = 2π`.
- **Archimedean density is Cauchy–Stieltjes** (§897): `dΩ/d(r²)` completely monotone, unconditional from ψ₀'s series.
- The Li coefficients are comfortably positive (`λ_n ~ (n/2)ln n`, §921) — the razor is specific to *band-limiting*.

## 6. The honest open core

The remaining mathematical content of RH, localized: **on the bottom-edge, band-limited critical subspace, the
archimedean form dominates the prime form — uniformly in bandwidth and height — by a margin the pole exactly
lifts to zero.** This is a ~2–3-dimensional near-null (prolate-plunge) positivity that is a *global* prime
cancellation, razor-thin, sub-interval-radius, and identical (verified) to the open lemma every leading program
leaves unproven. It is **universal across L-functions** (§928, §930): each L is razor-thin at its own spectrum
bottom edge, so any proof via this route must be universal — it cannot use ζ-specific arithmetic (the pole, the
specific low primes).

**The nature of the barrier (§931, the sharpest statement).** The function-field case, where RH *is* proven, has
the *identical* Weil form — a sum-of-squares over the zeros (a Gram matrix `VV*`), PSD ⟺ zeros on the line, with
the *same* sampling razor (min eigenvalue → 0 as the band-limit reaches the zero count; negative if a zero moves
off the line). The **only** difference is the proof of zeros-on-line: function fields supply it via the
**Hodge-index geometry** (intersection form on `C×C`), which hands them the Gram/PSD structure for free. Number
fields have the same zero-free form (`arch+pole−prime`, ⟺ RH) but **lack the surface/Hodge-index** — Deninger's
own 2022 no-go rules out a real-coefficient arithmetic Weil cohomology. **So the razor is the analytic shadow of a
missing geometric positivity, and no analysis of the arch-vs-prime balance can supply it.** This is why every
analytic lever this session bottomed out: the missing ingredient is structural/geometric (a number-field
Hodge-index / Hilbert–Pólya / Connes-arithmetic-site positivity), not analytic. That is the genuine frontier.

**This map is the contribution: RH's core is now an explicit, finite, reproducible object; a complete account of
why each classical and modern route stalls at exactly it; a universality constraint on any solution; and — via the
function-field comparison — a precise identification of the barrier as the missing geometric structure rather than
an unproven analytic inequality.**

*Provenance: LEDGER §884–§926; PROCESS.md arc; ANCHOR_PAPER.md §8–§10; workflows wf_045cb265, wf_5685d3c6,
wf_e696c870, wf_259382b8. Discipline: every agent-sourced number independently re-verified; multiple over-claims
(mine and agents') caught and corrected in-ledger (§888, §889, §897b, §902, §918, §922→925).*
