import Mathlib.Analysis.SpecialFunctions.Complex.Circle
import Mathlib.Analysis.Meromorphic.Divisor

/-!
# §1381 — splitting the zero-product by a predicate (half-disk vs annulus)

The polynomial half of the two-radii factorization: over the (finite) zero set `T`, the product factors as
`∏ = (∏ over p) · (∏ over ¬p)`.  With `p = (· ∈ closedBall c (R/2))` this splits `P = P_half · P_rest`, where
`P_half` carries the inner-disk zeros (bounded below by §1376 on the sphere) and `P_rest` the annulus zeros (folded
into the nonvanishing factor `g_full = P_rest · g`).  Direct from `Finset.prod_filter_mul_prod_filter_not`.
EXIT0 / std-axioms.
-/

open Classical in
theorem prod_split_pred {T : Finset ℂ} {s : ℂ} (m : ℂ → ℕ) (p : ℂ → Prop) :
    ∏ ρ ∈ T, (s - ρ) ^ (m ρ)
      = (∏ ρ ∈ T.filter p, (s - ρ) ^ (m ρ))
        * (∏ ρ ∈ T.filter (fun ρ => ¬ p ρ), (s - ρ) ^ (m ρ)) :=
  (Finset.prod_filter_mul_prod_filter_not T p _).symm

#print axioms prod_split_pred
