import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Tactic

/-!
# §1780 — Screw-increment identity: on mean-zero vectors the increment form equals the kernel form

The analytic input feeding §1779 (`PSDCodimOneNegOne`).  Suzuki proves the ζ-screw function `g_ξ` lies in
Kreĭn's class `𝒢_∞`: the increment kernel `G_g(t,u)=g(t-u)-g(t)-g(-u)+g(0)` is nonnegative definite.  To
turn that into "the plain Weil–CS form `∑ g(t_i-t_j) v_i v_j` is `≥ 0` on the mean-zero subspace" one needs
the elementary identity that the rank part `s(i)+sm(j)-g0` (with `s(i)=g(t_i)`, `sm(j)=g(-t_j)`, `g0=g(0)`)
contributes **zero** against any vector with `∑ v = 0`:

    ∑ᵢ∑ⱼ (K i j - s i - sm j + g0) vᵢ vⱼ = ∑ᵢ∑ⱼ K i j vᵢ vⱼ    when ∑ v = 0.

Combined with the `𝒢_∞` nonnegativity of the increment form and §1779, this yields `#neg ≤ 1` for the Weil–CS
form — the CS "simple isolated minimum" hypothesis.  Pure finite algebra (holds for ANY `K`, `s`, `sm`, `g0`).
EXIT0 / std-axioms.
-/

open Finset

variable {ι : Type*} [Fintype ι]

/-- **Mean-zero kills the rank part.**  The separable/rank piece `s i + sm j - g0` pairs to `0` against a
mean-zero vector `v` (each summand carries a factor `∑ v = 0`). -/
theorem meanzero_kills_rank_terms (s sm : ι → ℝ) (g0 : ℝ) (v : ι → ℝ)
    (hv : ∑ i, v i = 0) :
    ∑ i, ∑ j, (s i + sm j - g0) * v i * v j = 0 := by
  set T := ∑ j, sm j * v j with hT
  have inner : ∀ i, ∑ j, (s i + sm j - g0) * v i * v j = v i * T := by
    intro i
    have step : ∀ j, (s i + sm j - g0) * v i * v j
        = v i * ((s i - g0) * v j + sm j * v j) := by intro j; ring
    rw [Finset.sum_congr rfl (fun j _ => step j), ← Finset.mul_sum, Finset.sum_add_distrib,
      ← Finset.mul_sum, hv, mul_zero, zero_add, hT]
  rw [Finset.sum_congr rfl (fun i _ => inner i), ← Finset.sum_mul, hv, zero_mul]

/-- **Screw-increment identity.**  On a mean-zero vector, the increment-kernel quadratic form
`∑ (K i j - s i - sm j + g0) vᵢ vⱼ` equals the kernel form `∑ K i j vᵢ vⱼ`.  With `K i j = g(tᵢ-tⱼ)`,
`s i = g tᵢ`, `sm j = g(-tⱼ)`, `g0 = g 0`, the left side is the `G_g`-form (nonneg definite for `g ∈ 𝒢_∞`),
so the plain Weil–CS form is `≥ 0` on the mean-zero subspace. -/
theorem incrementForm_eq_kernelForm_on_meanZero (K : ι → ι → ℝ) (s sm : ι → ℝ) (g0 : ℝ)
    (v : ι → ℝ) (hv : ∑ i, v i = 0) :
    ∑ i, ∑ j, (K i j - s i - sm j + g0) * v i * v j = ∑ i, ∑ j, K i j * v i * v j := by
  have h := meanzero_kills_rank_terms s sm g0 v hv
  have split : ∀ i j, (K i j - s i - sm j + g0) * v i * v j
      = K i j * v i * v j - (s i + sm j - g0) * v i * v j := by intro i j; ring
  simp only [split, Finset.sum_sub_distrib]
  rw [h, sub_zero]
