import SmulFactorizedEval
import ZpowToNat
import ProdSplitPred

/-!
# §1386 — split-evaluation of the factorized-rational smul into `P_half · (P_rest · g)`

Combining §1338 (evaluate `((∏ᶠ (·-u)^{d u}) • g) s = ∏_{support} (s-u)^{d u} · g s`), §1380 (`zpow→ℕ` on the
nonneg divisor) and §1381 (split the Finset product by a predicate), the factorized-rational value splits into the
inner/outer product form `§1378` consumes:

  `((∏ᶠ (·-u)^{d u}) • g) s = P_half(s) · (P_rest(s) · g s)`,

with `P_half = ∏_{support, p} (s-u)^{(d u).toNat}` (inner/half-disk zeros) and `P_rest` the complement.  Setting
`p = (· ∈ closedBall c (R/2))` gives exactly `P_half` and `g_full = P_rest · g`.  EXIT0 / std-axioms.
-/

open Function

open Classical in
theorem smul_split_eval_nat {d : ℂ → ℤ} (hd : (support d).Finite) (hdnn : ∀ u, 0 ≤ d u)
    {g : ℂ → ℂ} (s : ℂ) (p : ℂ → Prop) :
    ((∏ᶠ u, (· - u) ^ (d u)) • g) s
      = (∏ u ∈ hd.toFinset.filter p, (s - u) ^ (d u).toNat)
        * ((∏ u ∈ hd.toFinset.filter (fun u => ¬ p u), (s - u) ^ (d u).toNat) * g s) := by
  rw [smul_factorized_rational_eval hd]
  have hconv : (∏ u ∈ hd.toFinset, (s - u) ^ (d u))
      = ∏ u ∈ hd.toFinset, (s - u) ^ (d u).toNat :=
    Finset.prod_congr rfl (fun u _ => zpow_eq_pow_toNat (hdnn u))
  rw [hconv, prod_split_pred (fun u => (d u).toNat) p]
  ring

#print axioms smul_split_eval_nat
