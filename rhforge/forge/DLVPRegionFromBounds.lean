import DLVPGapFromBounds

/-!
# §1311 ★ Parameter-level quantitative zero-free region in ONE theorem: `β ≤ 1 - 1/(14K)`

Composing §1310 (four bounds → gap) with the optimal choice `σ = 1 + 1/(2K)`, `K = 3A+4B+C`.  The gap RHS
collapses (`4·(1/(2K))/(3+K/(2K)) = 4/(7K)`), giving `4/(7K) ≤ 1+1/(2K)-β`, i.e. `β ≤ 1 - 1/(14K)`.  With the
Borel–Carathéodory growth `B ≍ log|t|` (so `K ≍ log|t|`), this is the quantitative de la Vallée Poussin region
`β ≤ 1 - c/log|t|`, `c = 1/14`, ENTIRELY assembled at the parameter level.  The ONLY remaining ζ-specific input
is `hzero` — the upper bound `Re(-ζ'/ζ(σ+iγ)) ≤ B - 1/(σ-β)` (the Hadamard / partial-fractions bridge).  This is
the analytic analogue of §1277's `heig`: the region is now ONE named lemma away.  EXIT0 / std-axioms.
-/

theorem dlvp_region_of_four_bounds {β T₀ T₁ T₂ A B C : ℝ}
    (hK : 0 < 3 * A + 4 * B + C) (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C)
    (hβ : β < 1 + 1 / (2 * (3 * A + 4 * B + C)))
    (hdlvp : 0 ≤ 3 * T₀ + 4 * T₁ + T₂)
    (hpole : T₀ ≤ 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1) + A)
    (hzero : T₁ ≤ B - 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - β))
    (hgrowth : T₂ ≤ C) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  set K := 3 * A + 4 * B + C with hKdef
  have hKne : K ≠ 0 := hK.ne'
  have hσ : (1 : ℝ) < 1 + 1 / (2 * K) := by
    have : 0 < 1 / (2 * K) := by positivity
    linarith
  have hgap := dlvp_gap_of_four_bounds hσ hβ hA hB hC hdlvp hpole hzero hgrowth
  have hσ1 : (1 + 1 / (2 * K)) - 1 = 1 / (2 * K) := by ring
  rw [hσ1] at hgap
  have hcollapse : 4 * (1 / (2 * K)) / (3 + K * (1 / (2 * K))) = 4 / (7 * K) := by
    field_simp; ring
  rw [hcollapse] at hgap
  have harith : 1 + 1 / (2 * K) - 4 / (7 * K) = 1 - 1 / (14 * K) := by
    field_simp; ring
  linarith [hgap, harith]

#print axioms dlvp_region_of_four_bounds
