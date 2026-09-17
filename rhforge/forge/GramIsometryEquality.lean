import Mathlib.Analysis.InnerProductSpace.Basic
import Mathlib.Tactic

/-!
# §1760 — The equality criterion's geometric core: diagonal Gram equality ⟹ L is a similarity

Suzuki's RH equality criterion (Cor 1.5, arXiv:2301.00421) says RH ⟺ the Weil form equals `2‖ψ‖²_{L²}` on
`V°(0)`, where the Weil form has the Gram representation `⟨ψ,ψ⟩_W = ‖L ψ‖²` with `L = P_d D` (Hardy projection
of the momentum operator, §1756).  The geometric content of "`‖Lψ‖² = c‖ψ‖²` on the diagonal" is that `L`
is a *similarity* of factor `√c` on the subspace — it preserves inner products up to the scalar `c`:

`⟨L x, L x⟩ = c⟨x, x⟩  ∀x  ⟹  ⟨L x, L y⟩ = c⟨x, y⟩  ∀x,y.`

This is polarization applied to the Gram form `⟨L·,L·⟩ − c⟨·,·⟩` (symmetric, diagonal-zero ⟹ zero).  It closes
the gap from the Gram representation (§1756) to the Cor-1.5 equality: under RH, `P_d D` restricted to `V°(0)`
is exactly a `√2`-isometry into `L²`, which is why "Weil positivity is explained by the nonnegativity of the
L²-norm."  Combined with §1753/§1754/§1757 this is the formal skeleton of the equality criterion.
EXIT0 / std-axioms.
-/

variable {V W : Type*} [NormedAddCommGroup V] [InnerProductSpace ℝ V]
  [NormedAddCommGroup W] [InnerProductSpace ℝ W]

/-- **Diagonal Gram equality ⟹ full Gram equality (L is a `√c`-similarity).** -/
theorem gram_similarity_of_diag (L : V →ₗ[ℝ] W) (c : ℝ)
    (h : ∀ x, (inner ℝ (L x) (L x) : ℝ) = c * (inner ℝ x x : ℝ)) (x y : V) :
    (inner ℝ (L x) (L y) : ℝ) = c * (inner ℝ x y : ℝ) := by
  have hxy := h (x + y)
  rw [map_add, real_inner_add_add_self, real_inner_add_add_self] at hxy
  have hx := h x
  have hy := h y
  nlinarith [hxy, hx, hy]

/-- **Corollary (Cor 1.5 shape).**  If the Gram/Weil form equals `c·I` on the diagonal, then it equals the
scaled inner product as a full bilinear form: `⟨L x, L y⟩ = c⟨x,y⟩`.  With `c = 2` and `L = P_d D` this is
exactly Suzuki's `⟨ψ,ψ'⟩_W = 2⟨ψ,ψ'⟩_{L²}` on `V°(0)` — the equality criterion, polarized. -/
theorem weil_eq_scaled_inner_of_diag (L : V →ₗ[ℝ] W) (c : ℝ)
    (h : ∀ x, (inner ℝ (L x) (L x) : ℝ) = c * (inner ℝ x x : ℝ)) :
    ∀ x y, (inner ℝ (L x) (L y) : ℝ) = c * (inner ℝ x y : ℝ) :=
  gram_similarity_of_diag L c h
