import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.Tactic

/-!
# §1774 — When the arithmetic Hodge form has a positive direction: the sharp scalar criterion

The Arakelov ℝ-model (forge §1713) carries a positive-definite pairing `⟪·,·⟫` and a degree functional
`deg = ⟪g, ·⟫` (Riesz representative `g`).  The Hodge form is `H(x,y) = deg(x)·deg(y) − ⟪x,y⟫`
(`= deg⊗deg − pairing`).  The arithmetic Hodge-index signature `(1, dim−1)` needs EXACTLY one positive
direction — i.e. a *hodge-positive vector* (the opaque hypothesis of forge §1719/§1720/§1721).

This brick replaces that opaque existence hypothesis with a **sharp concrete scalar criterion**:

    (∃ w, 0 < H(w,w))  ↔  1 < ⟪g, g⟫    (= ‖deg‖² in the pairing-dual metric).

So the positive direction exists iff the degree functional has pairing-dual norm `> 1` — a single
number depending only on the model's normalization.  Proof: Cauchy–Schwarz `deg(w)² = ⟪g,w⟫² ≤ ⟪g,g⟫⟪w,w⟫`
gives `H(w,w) ≤ (⟪g,g⟫−1)⟪w,w⟫`, and `w = g` realizes it (`H(g,g) = ⟪g,g⟫(⟪g,g⟫−1)`).  Over any real
inner-product space.  EXIT0 / std-axioms.
-/

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

open scoped InnerProductSpace

/-- **Sharp positive-direction criterion for the arithmetic Hodge form.**  With `deg = ⟪g,·⟫`, the Hodge
form `H(x,y) = ⟪g,x⟫⟪g,y⟫ − ⟪x,y⟫` has a positive vector iff `1 < ⟪g,g⟫`. -/
theorem hodge_positive_vector_iff (g : E) :
    (∃ w : E, 0 < ⟪g, w⟫_ℝ * ⟪g, w⟫_ℝ - ⟪w, w⟫_ℝ) ↔ 1 < ⟪g, g⟫_ℝ := by
  constructor
  · rintro ⟨w, hw⟩
    have hcs : ⟪g, w⟫_ℝ * ⟪g, w⟫_ℝ ≤ ⟪g, g⟫_ℝ * ⟪w, w⟫_ℝ := real_inner_mul_inner_self_le g w
    have hww_pos : 0 < ⟪w, w⟫_ℝ := by
      rcases (real_inner_self_nonneg (x := w)).lt_or_eq with h | h
      · exact h
      · exfalso
        have hw0 : w = 0 := inner_self_eq_zero.mp h.symm
        rw [hw0] at hw
        simp at hw
    nlinarith [hw, hcs, hww_pos]
  · intro hg
    refine ⟨g, ?_⟩
    nlinarith [hg]
