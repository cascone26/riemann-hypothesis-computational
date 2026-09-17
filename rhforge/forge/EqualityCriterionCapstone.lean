import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.Tactic

/-!
# §1761 — Equality-criterion capstone: equality on `V°(0)` ⟹ strict Weil positivity (RH-direction of Cor 1.5)

Assembles the equality-criterion skeleton (§1754/§1756/§1757/§1760) into the single statement that IS the
RH-direction of Suzuki's Corollary 1.5 (arXiv:2301.00421).  The Weil form has the Gram representation
`W(x,y) = ⟨L x, L y⟩` (`L = P_d D`, §1756).  Suzuki's equality criterion supplies the hypothesis
`⟨L x, L x⟩ = c⟨x, x⟩` with `c = 2 > 0` on `V°(0)`.  From this ALONE we get, on that space:

1. `W` equals the scaled inner product as a full bilinear form (`⟨L x, L y⟩ = c⟨x, y⟩`)          — §1760 polarization
2. `W` is positive semidefinite (`0 ≤ ⟨L x, L x⟩`)                                                — Gram
3. `W` is positive **definite** (`⟨L x, L x⟩ = 0 ↔ x = 0`), given `c > 0`                          — the strict step
4. `L` is injective (`L x = 0 ⟹ x = 0`)

So "the Weil form is `2‖·‖²` on `V°(0)`" (RH, via Cor 1.5) forces the Weil form to be strictly positive there,
with `L = P_d D` a `√c`-similarity onto its image — "the nonnegativity of Weil's form is explained by the
nonnegativity of the L²-norm."  EXIT0 / std-axioms.
-/

variable {V W : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
  [NormedAddCommGroup W] [InnerProductSpace ℝ W]

/-- Polarization (from §1760): diagonal Gram equality ⟹ full Gram equality. -/
theorem gram_full_of_diag (L : V →ₗ[ℝ] W) (c : ℝ)
    (h : ∀ x, (inner ℝ (L x) (L x) : ℝ) = c * (inner ℝ x x : ℝ)) (x y : V) :
    (inner ℝ (L x) (L y) : ℝ) = c * (inner ℝ x y : ℝ) := by
  have hxy := h (x + y)
  rw [map_add, real_inner_add_add_self, real_inner_add_add_self] at hxy
  nlinarith [hxy, h x, h y]

/-- **PSD.**  The Gram/Weil form is nonnegative. -/
theorem weil_psd (L : V →ₗ[ℝ] W) (x : V) : 0 ≤ (inner ℝ (L x) (L x) : ℝ) :=
  real_inner_self_nonneg

/-- **Strict positivity (the capstone).**  Under the equality criterion with `c > 0`, the Weil form is
positive-definite: it vanishes only at `0`. -/
theorem equality_criterion_posDef (L : V →ₗ[ℝ] W) (c : ℝ) (hc : 0 < c)
    (h : ∀ x, (inner ℝ (L x) (L x) : ℝ) = c * (inner ℝ x x : ℝ)) (x : V) :
    (inner ℝ (L x) (L x) : ℝ) = 0 ↔ x = 0 := by
  rw [h x]
  constructor
  · intro hh
    have hx : (inner ℝ x x : ℝ) = 0 := by
      rcases mul_eq_zero.mp hh with h1 | h2
      · exact absurd h1 (ne_of_gt hc)
      · exact h2
    exact inner_self_eq_zero.mp hx
  · intro hh; rw [hh]; simp

/-- **`L` injective on `V°(0)`.**  Under the equality criterion with `c > 0`, `L = P_d D` has trivial kernel. -/
theorem equality_criterion_injective (L : V →ₗ[ℝ] W) (c : ℝ) (hc : 0 < c)
    (h : ∀ x, (inner ℝ (L x) (L x) : ℝ) = c * (inner ℝ x x : ℝ)) (x : V) (hx : L x = 0) :
    x = 0 := by
  refine (equality_criterion_posDef L c hc h x).mp ?_
  rw [hx, inner_zero_left]
