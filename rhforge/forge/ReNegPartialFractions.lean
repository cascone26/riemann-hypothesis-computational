import DropZeros

/-!
# §1312 — hzero-shape from the partial-fractions decomposition (abstract)

The DLVP `hzero` input `Re(-ζ'/ζ(s)) ≤ B - 1/(σ-β)` follows from the Hadamard / partial-fractions form
`-ζ'/ζ(s) = -((∑_{ρ∈S} 1/(s-ρ)) + g)` (finite disk-zero set `S ∋ ρ₀`, bounded remainder `g` with `Re g ≥ -B`)
by: drop all zeros but `ρ₀` (§1147, each term has `Re ≥ 0`) and bound `-Re g ≤ B`.  Pure complex-analysis
lemma — it isolates the ONE genuinely-open analytic fact (the decomposition itself) from the algebra around it,
exactly as §1311 isolated `hzero` from the region.  EXIT0 / std-axioms.
-/

open Complex

theorem re_neg_partial_fractions_le {s ρ₀ g : ℂ} {S : Finset ℂ} {B : ℝ}
    (hρ₀ : ρ₀ ∈ S) (hSleft : ∀ ρ ∈ S, ρ.re < s.re) (hgB : -B ≤ g.re) :
    (-((∑ ρ ∈ S, 1 / (s - ρ)) + g)).re ≤ B - (1 / (s - ρ₀)).re := by
  have hdrop := re_sum_one_div_sub_ge_single S ρ₀ hρ₀ hSleft
  rw [Complex.neg_re, Complex.add_re]
  linarith [hdrop, hgB]

#print axioms re_neg_partial_fractions_le
