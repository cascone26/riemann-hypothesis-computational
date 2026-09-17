import OneLeDivisorZeta
import ZetaDivisorFinite

/-!
# §1353 — the kept DLVP zero lies in the divisor support (§1336's `hρ₀`)

`§1336`'s master theorem requires the kept zero `β+iγ` to be an element of the finite divisor-support Finset `S`
(so its `1/(σ-β)` term is actually present in the partial-fractions sum).  From §1351 the divisor at a genuine zero
is `≥ 1`, hence `≠ 0`, so the zero lies in `Function.support (divisor ζ (ball c R))`, i.e. in
`(divisor_riemannZeta_ball_support_finite c R).toFinset` — precisely §1349's index set.  EXIT0 / std-axioms.
-/

open Complex Metric

theorem zeta_zero_mem_divisor_support {c : ℂ} {R : ℝ} (hball : (1 : ℂ) ∉ ball c R)
    {ρ z₀ : ℂ} (hρ : ρ ∈ ball c R) (hζρ : riemannZeta ρ = 0)
    (hz₀ : z₀ ∈ ball c R) (hz₀ne : riemannZeta z₀ ≠ 0) :
    ρ ∈ (divisor_riemannZeta_ball_support_finite c R).toFinset := by
  rw [Set.Finite.mem_toFinset, Function.mem_support]
  have h1 : 1 ≤ MeromorphicOn.divisor riemannZeta (ball c R) ρ :=
    one_le_divisor_riemannZeta_at_zero hball hρ hζρ hz₀ hz₀ne
  omega

#print axioms zeta_zero_mem_divisor_support
