# First-X Meta-Pass — re-mining the walled basins with §184's newly-discovered heuristic (2026-07-23)

## The heuristic, stated precisely
§184's breakthrough was NOT a tool — it was a *pattern* that didn't exist before it: Yoshida/Bombieri prove
Weil positivity **unconditionally** for support ≤ log 2; it **breaks at the first prime (p=2)**, and that first
obstruction *sets the criticality* (p=2 alone collapses the positivity margin ~340×). The whole first-prime
program came from asking, of an "exhausted" basin: **"where does the known-unconditional partial result stop, and
what is the FIRST obstruction beyond that boundary?"** That question was never asked of the *other* walled basins,
because it didn't exist until §184 — each of them was walled as a *whole basin*, not probed for an interior
boundary. This pass re-asks it, once, against every closed basin. Output is Report-tier (ranked NEW angles for
ATTACK_QUEUE), not results.

**Sharper form of the heuristic (the generalization that matters):** every RH-equivalent basin has an
*unconditionally-proven partial version* with a *boundary*. Progress lives at that boundary + the first thing that
crosses it. The lever is real (not circular) only if the first obstruction is on the **prime/analytic side**
(something we can bound), not the **zero side** (using zeros to constrain zeros = circular).

---

## Basin-by-basin

### 1. Crystalline rigidity / Lee-Yang (§183) — HIGHEST re-mining promise
- **Unconditional partial result:** Lee-Yang proves RH-for-a-partition-function when reflection symmetry **plus
  coefficient positivity** both hold; ξ has the reflection (functional equation) but the positivity is exactly the
  Weil A+B≥P wall. In the quasicrystal/structure-factor picture, the prime diffraction is *decisive* (physics
  deepdive 2026-07-18).
- **Boundary + first obstruction:** the first prime's **structure-factor / Bragg contribution**. §184 already used
  "first prime" on the *Weil* face of this basin — but NOT on the *diffraction* face. The first-X lever here is:
  does the first prime's diffraction peak have its own criticality threshold (an analogue of the 340× margin
  collapse) expressed in structure-factor language? That would be a **second, independent handle on the same
  criticality** — potentially sharper because diffraction intensity is a positivity-manifest quantity.
- **Verdict:** re-mine. Concrete next probe — compute the single-prime (p=2) structure-factor contribution to the
  unfolded-zero diffraction and look for a threshold in a mirroring the Weil margin collapse. Prime-side ⟹ not
  circular. **Rank #1.**

### 2. Zero-density (Guth-Maynard) — reframes into a LIVE frontier
- **Unconditional partial result:** the classical **zero-free region** Re(s) > 1 − c/log t (and Vinogradov-Korobov
  widening). This is the honest unconditional "partial RH."
- **Boundary + first obstruction:** the edge of the zero-free region. The first obstruction to widening it is the
  first prime-power/character-sum term that limits the classical method (the same Λ(n)/√n prime norms that appear
  in our §189 norm-budget relay!). Guth-Maynard was walled as "density axis, not RH-path" — but the first-X reframe
  points not at *density* (how many can be off) but at the **zero-free-region boundary** (whether ANY can be off
  near σ=1), where the first-prime budget is the actual limiting resource.
- **Verdict:** the reframe connects the zero-free-region frontier to OUR §189 prime-norm budget — the same
  Λ(p)/√p ledger that covers arch-negativity growth may be the resource that bounds the zero-free region. Worth a
  focused check: is the first-prime norm budget the shared currency between our Weil program and zero-free-region
  widening? If yes, a result in one transfers. Prime-side ⟹ not circular. **Rank #2.**

### 3. de Bruijn-Newman (§35-47e) — novel angle, but circularity risk
- **Unconditional partial result:** Newman Λ ≥ 0; current unconditional upper bound Λ ≤ ~0.2 (and RH ⟺ Λ ≤ 0, so
  Λ = 0 under RH). The known-good region is t ≥ Λ (all zeros real).
- **Boundary + first obstruction:** as t decreases toward Λ⁺, the FIRST (lowest) zero-pair to approach collision
  sets the near-boundary dynamics. First-X = the lowest zeros / first spectral gap of the heat flow.
- **Verdict:** the first-X here lives on the **zero side** (lowest zeros drive the flow), so using it to *force*
  Λ ≤ 0 risks circularity (§35-47e already flagged DBN as detector-only). BUT there's a non-circular sub-angle:
  the heat-flow generator near t=0 is a prime-side object via the explicit formula — if the first-prime budget
  controls the t→0⁺ flow's coercivity, that's a real lever. Lower confidence. **Rank #3, conditional on finding a
  prime-side (not zero-side) formulation of the near-t=0 flow.**

### 4. LeClair spectral-flow (§125 = Lagarias) — no clean non-circular first-X
- **Unconditional partial result:** the smooth counting-function approximation N(T) ≈ (T/2π)log(T/2π) − T/2π + 7/8;
  the remainder S(T) is unconditionally O(log T).
- **Boundary + first obstruction:** the first large excursion of S(T). But S(T) is unbounded (grows like
  √(log T log log T) even under RH — Montgomery/Selberg), so there is no "first obstruction" threshold that a
  finite budget crosses; the criticality isn't set by a *first* instance.
- **Verdict:** no clean first-X lever that isn't already the whole Lagarias sum. **Rank #4 — pass (do not re-mine).**

---

## Output (for ATTACK_QUEUE)
The meta-pass is productive on **2 of 4** basins, both prime-side (non-circular):
1. **Crystalline/Lee-Yang, diffraction face:** hunt the first-prime structure-factor criticality threshold — a
   second independent handle on the same margin-collapse §184 found on the Weil face.
2. **Zero-free-region ↔ §189 prime-norm budget:** test whether the first-prime Λ(p)/√p budget is the shared
   limiting resource between our Weil program and zero-free-region widening; if so, results transfer.
Both are genuinely NEW angles (the question that generated them didn't exist before §184) and both keep the lever
on the analytic/prime side, avoiding the zero-side circularity that walled DBN and LeClair. The meta-move itself —
"locate the boundary of the unconditional partial result and probe the first obstruction beyond it" — is now a
standing tool to apply to any future 'walled' verdict before accepting it.
