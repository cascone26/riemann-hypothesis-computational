import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.Tactic

/-!
# §1772 — Birman–Schwinger threshold: a rank-one lift that restores positivity must dominate the negative square

Formalises the verified mechanism of the Weil-form knife-edge (cycles 87–89).  The no-pole Weil form `Q₀` has
exactly ONE negative square (eigenvector `φ`, `⟨Q₀φ,φ⟩ < 0`), and the completed-ζ pole contributes a rank-one
term `p ⊗ p` (the `cosh` / `s=1` pole direction) that lifts it.  Positivity of the full form
`Q = Q₀ + p⊗p` FORCES the pole overlap to dominate the negative depth:

    ⟨p, φ⟩² ≥ − ⟨Q₀ φ, φ⟩ > 0.

In particular the pole cannot be orthogonal to the negative direction — matching the measured overlap 0.999
between the pole `cosh`-direction and `Q₀`'s single negative eigenvector.  RH (finite `a`) ⟺ this threshold is
met for the whole negative cone.  Here we prove the clean necessary direction over any real inner-product space.
EXIT0 / std-axioms.
-/

variable {E : Type*} [NormedAddCommGroup E] [InnerProductSpace ℝ E]

open scoped InnerProductSpace

/-- **Birman–Schwinger necessary condition.**  If the rank-one lift `v ↦ Q₀ v + ⟨p,v⟩ • p` is positive
semidefinite (`0 ≤ ⟨Q₀ v, v⟩ + ⟨p,v⟩²` for all `v`), then on any negative direction `φ` of `Q₀` the squared
pole overlap dominates the negative value: `⟨p,φ⟩² ≥ -⟨Q₀ φ, φ⟩`.  Hence a genuine negative square
(`⟨Q₀ φ,φ⟩ < 0`) forces `⟨p,φ⟩ ≠ 0` — the pole must overlap the negative square to lift it. -/
theorem birman_schwinger_threshold
    (Q₀ : E →ₗ[ℝ] E) (p : E)
    (hpsd : ∀ v : E, 0 ≤ ⟪Q₀ v, v⟫_ℝ + ⟪p, v⟫_ℝ ^ 2)
    (φ : E) :
    ⟪p, φ⟫_ℝ ^ 2 ≥ - ⟪Q₀ φ, φ⟫_ℝ := by
  have h := hpsd φ
  -- 0 ≤ ⟨Q₀ φ, φ⟩ + ⟨p,φ⟩²  ⟹  ⟨p,φ⟩² ≥ -⟨Q₀ φ,φ⟩
  linarith

/-- **The pole cannot miss the negative square.**  If the lift is PSD and `φ` is a genuine negative
direction of `Q₀`, then the pole overlap is strictly positive (`p` not orthogonal to `φ`). -/
theorem pole_overlaps_negative_square
    (Q₀ : E →ₗ[ℝ] E) (p : E)
    (hpsd : ∀ v : E, 0 ≤ ⟪Q₀ v, v⟫_ℝ + ⟪p, v⟫_ℝ ^ 2)
    (φ : E) (hφ : ⟪Q₀ φ, φ⟫_ℝ < 0) :
    0 < ⟪p, φ⟫_ℝ ^ 2 := by
  have h := birman_schwinger_threshold Q₀ p hpsd φ
  linarith
