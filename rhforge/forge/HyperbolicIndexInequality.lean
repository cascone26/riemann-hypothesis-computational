import Mathlib.LinearAlgebra.BilinearForm.Basic
import Mathlib.LinearAlgebra.BilinearForm.Properties
import Mathlib.Tactic

set_option maxHeartbeats 1000000

noncomputable section

open LinearMap

variable {F M : Type*} [Field F] [LinearOrder F] [IsStrictOrderedRing F]
  [AddCommGroup M] [Module F M]

/-!
# §1725 ★★★★★ — WEIL'S HYPERBOLIC-PLANE INEQUALITY (the connector's algebraic core)

The one piece of the function-field template our chain had NOT yet formalized (flagged honestly in
THE_BRIDGE.md as the "correspondence connector"): the inequality that turns Hodge-index negativity
into the Frobenius **trace bound** `t² ≤ 4q` feeding §1724.

Weil's move on `C×C` uses the two rulings `f₁ = [pt×C]`, `f₂ = [C×pt]`: a **hyperbolic pair**,
`f₁² = f₂² = 0`, `f₁·f₂ = 1`, carrying the single positive direction of the Hodge index; the form is
negative on their orthogonal complement (the primitive part, §1720/§1713). For ANY class `v`,
    `B(v,v) ≤ 2 · (v·f₁) · (v·f₂)`.
Applied to the graph of Frobenius (with `v·f₁ = 1`, `v·f₂ = q`, `B(v,v) = 2·(#fixed pts) − …`) this
is exactly the inequality that becomes `t² ≤ 4q`. Here it is proved abstractly, over any ordered
field, from the hyperbolic-pair relations plus negativity-on-the-complement — no curve needed. The
curve-specific numerics (`v·f₁=1, v·f₂=q, t = trace`) are the only remaining Report-tier input.

Proof: `u := v − (v·f₂)f₁ − (v·f₁)f₂` is orthogonal to both `f₁,f₂` (hence in the negative
complement), so `B(u,u) ≤ 0`; and `B(u,u) = B(v,v) − 2(v·f₁)(v·f₂)` by bilinearity.
-/

/-- **Weil's hyperbolic-plane index inequality.** Given a symmetric bilinear form `B`, a hyperbolic
pair `f₁, f₂` (`B f₁ f₁ = 0`, `B f₂ f₂ = 0`, `B f₁ f₂ = 1`), and negativity of `B` on the orthogonal
complement of the pair, every vector `v` satisfies `B v v ≤ 2 · (B v f₁) · (B v f₂)`. -/
theorem hyperbolic_index_ineq (B : LinearMap.BilinForm F M)
    (hsymm : ∀ x y, B x y = B y x)
    (f₁ f₂ : M) (h11 : B f₁ f₁ = 0) (h22 : B f₂ f₂ = 0) (h12 : B f₁ f₂ = 1)
    (hnp : ∀ z : M, B z f₁ = 0 → B z f₂ = 0 → B z z ≤ 0)
    (v : M) :
    B v v ≤ 2 * (B v f₁) * (B v f₂) := by
  set u : M := v - (B v f₂) • f₁ - (B v f₁) • f₂ with hu
  have h21 : B f₂ f₁ = 1 := by rw [hsymm]; exact h12
  -- u is orthogonal to f₁ and f₂
  have huf1 : B u f₁ = 0 := by
    simp only [hu, map_sub, map_smul, LinearMap.sub_apply, LinearMap.smul_apply, smul_eq_mul,
      h11, h21, mul_zero, mul_one]
    ring
  have huf2 : B u f₂ = 0 := by
    simp only [hu, map_sub, map_smul, LinearMap.sub_apply, LinearMap.smul_apply, smul_eq_mul,
      h22, h12, mul_zero, mul_one]
    ring
  -- B u u = B v v − 2 (B v f₁)(B v f₂) by full bilinear expansion
  have huu : B u u = B v v - 2 * (B v f₁) * (B v f₂) := by
    simp only [hu, map_sub, map_smul, LinearMap.sub_apply, LinearMap.smul_apply, smul_eq_mul,
      h11, h22, h12, h21, hsymm f₁ v, hsymm f₂ v]
    ring
  have hle : B u u ≤ 0 := hnp u huf1 huf2
  rw [huu] at hle
  linarith

#print axioms hyperbolic_index_ineq
