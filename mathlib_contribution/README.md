# mathlib contribution: Bessel functions J₀, J₁ (draft)

`BesselJ01.lean` is a self-contained, general-purpose file extracted from the RiemannLean project's
`Bessel.lean` (all RH-specific lemmas removed). It **compiles clean against our pinned mathlib**
(Lean v4.29.0) with only the trusted kernel axioms `{propext, Classical.choice, Quot.sound}`.

## What it contains
- `besselJ0`, `besselJ1` — power-series definitions.
- `besselJ0_series_summable`, `besselJ0_zero`, `besselJ0_even`.
- `besselJ0_pos` (`0<J₀` on `|x|<2`); two-sided bounds `besselJ0_ge_one_sub`, `besselJ0_le_one`,
  the sharp `besselJ0_two_term` (`|J₀−(1−(x/2)²)|≤(x/2)⁴/4`) + corollary `besselJ0_le_two_term`.
- `besselJ1_pos` (`0<J₁` on `0<x`, `x²<8`), `besselJ1_le_half`.

## Status / honest process notes
mathlib currently has **no Bessel functions** — this is a genuine gap. But mathlib strongly prefers a
**Zulip design discussion before a PR for a new area**, to agree on the canonical definition (order-`0,1`
series vs general `Jᵥ` vs the ODE/integral rep). This file deliberately covers only orders 0 and 1 with
elementary bounds; a reviewer will likely want a generalization plan.

**Recommended next step (needs Jacob):**
1. Post the draft below to the mathlib Zulip `#mathlib4` (or `#new members`) stream to propose the addition
   and get design guidance.
2. `git clone` mathlib master, `lake exe cache get`, drop this file under
   `Mathlib/Analysis/SpecialFunctions/Bessel.lean`, fix any API drift from v4.29.0 → master, `lake build`.
3. Sign the mathlib CLA (once), open the PR (draft first), iterate with the reviewer.

I (Claude) did not open the PR unilaterally: it attaches Jacob's GitHub identity to a public contribution
and the design conversation should be his. The file is ready; the submission is a one-word go.

## Zulip proposal draft

> **Proposing: Bessel functions J₀, J₁ (elementary bounds).** mathlib has no Bessel functions. I have a
> self-contained file defining `J₀`/`J₁` by their power series and proving elementary results below the
> first zero: `J₀(0)=1`, evenness, `J₀>0` on `|x|<2`, the two-sided bounds `1−(x/2)² ≤ J₀ ≤ 1` and the
> sharp `|J₀−(1−(x/2)²)| ≤ (x/2)⁴/4`, and `0<J₁≤x/2` on `0<x, x²<8` — all via
> `alternating_series_error_bound`. Before a PR I'd like design input: is the preferred entry point
> order-`0,1` power series (as here), or should this be `Jᵥ` for general order from the start (with these
> as instances)? And is `Mathlib/Analysis/SpecialFunctions/Bessel.lean` the right home? Happy to
> generalize; the order-`0,1` proofs are self-contained and need no new infrastructure. File compiles
> clean (trusted axioms only) and I can push a branch immediately.
