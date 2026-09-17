import HankelConverseCore
import RealLinearSurj

/-!
# §1799 — Full Hermite converse for a single conjugate pair (end-to-end composition)

Composes §1798 (`real_linear_surj_at_nonreal`, build the isolating real vector) with §1797
(`hankel_isolated_nonreal_neg`, an isolating vector drives the form to `-2`) to prove the **complete** converse
for the two-root case: for a non-real `r₀`, the complex Hankel form of the conjugate pair `{r₀, r̄₀}` is NOT
positive semidefinite — there is a real coefficient vector making it `-2 < 0`.  This verifies the converse
machinery composes end to end; the general case adds only the product `Q` vanishing on the remaining roots.
EXIT0 / std-axioms.
-/

open Finset

/-- **Converse, single-pair case.**  A non-real conjugate pair's Hankel form is not PSD (a real vector makes
it `-2`). -/
theorem hankel_pair_not_psd (r₀ : ℂ) (hr : r₀.im ≠ 0) :
    ∃ c : Fin 2 → ℝ,
      ((r₀ ::ₘ (starRingEnd ℂ r₀) ::ₘ 0).map
        (fun r => (∑ i : Fin 2, (c i : ℂ) * r ^ (i : ℕ)) ^ 2)).sum = -2 := by
  obtain ⟨a, b, hab⟩ := real_linear_surj_at_nonreal r₀ hr Complex.I
  refine ⟨![a, b], ?_⟩
  apply hankel_isolated_nonreal_neg (0 : Multiset ℂ) ![a, b] r₀
  · rw [show (∑ i : Fin 2, ((![a, b] i : ℝ) : ℂ) * r₀ ^ (i : ℕ)) = (a : ℂ) + (b : ℂ) * r₀ from by
        simp [Fin.sum_univ_two]]
    exact hab
  · intro s hs
    simp at hs
