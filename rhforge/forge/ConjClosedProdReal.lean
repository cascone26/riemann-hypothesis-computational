import Mathlib.Algebra.Polynomial.Roots
import Mathlib.Tactic

/-!
# §1800 — A conjugation-closed root set gives a real-coefficiented product

General-case converse, step 1.  For the Hermite converse's isolating polynomial we need `Q = ∏_{s∈rest}(X−s)`
to have REAL coefficients so that (multiplied by §1798's real linear factor) it yields a real coefficient
vector for §1797.  If the "rest" multiset `S` is closed under conjugation (`S.map conj = S`, true for the
non-`{r₀,r̄₀}` roots of a real polynomial), then `Polynomial.map conj (∏_{s∈S}(X − C s)) = ∏_{s∈S}(X − C s)`,
i.e. every coefficient is fixed by conjugation — hence real.  EXIT0 / std-axioms.
-/

open Polynomial

/-- **Conj-closed ⟹ real product.**  `∏_{s∈S}(X − C s)` is fixed by coefficientwise conjugation when `S` is
conjugation-closed. -/
theorem conj_closed_prod_real (S : Multiset ℂ) (hS : S.map (starRingEnd ℂ) = S) :
    Polynomial.map (starRingEnd ℂ) (S.map (fun s => X - C s)).prod
      = (S.map (fun s => X - C s)).prod := by
  rw [Polynomial.map_multiset_prod, Multiset.map_map]
  conv_rhs => rw [← hS, Multiset.map_map]
  refine congrArg Multiset.prod (Multiset.map_congr rfl (fun s _ => ?_))
  simp only [Function.comp_apply, Polynomial.map_sub, Polynomial.map_X, Polynomial.map_C]
