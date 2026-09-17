import DLVPAssemble
import DLVPGap

/-!
# §1310 — Four DLVP bounds → explicit zero gap (composing §1148 + §1144)

`dlvp_assemble` (§1148) turns the four analytic inputs into `4/(σ-β) ≤ 3/(σ-1) + (3A+4B+C)`, and
`dlvp_zero_gap` (§1144) turns that shape into the explicit gap `4(σ-1)/(3+K(σ-1)) ≤ σ-β`.  They were never
composed.  This lemma does it in one step: from the trig positivity, the pole bound, the zero bound, and the
growth bound (with the three constants `A,B,C ≥ 0`), the gap `σ-β` is bounded below explicitly with
`K = 3A+4B+C`.  Feeding `σ = 1+1/(2K)` (the optimal choice, §1145) then yields `β ≤ 1 - 1/(14K)`, `K ≍ log|t|`.
This closes the parameter-level region derivation to a single composition; only the ζ-analytic `hzero`
(the Hadamard / partial-fractions bridge) remains open.  EXIT0 / std-axioms.
-/

theorem dlvp_gap_of_four_bounds {σ β T₀ T₁ T₂ A B C : ℝ}
    (hσ : 1 < σ) (hβ : β < σ) (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C)
    (hdlvp : 0 ≤ 3 * T₀ + 4 * T₁ + T₂)
    (hpole : T₀ ≤ 1 / (σ - 1) + A)
    (hzero : T₁ ≤ B - 1 / (σ - β))
    (hgrowth : T₂ ≤ C) :
    4 * (σ - 1) / (3 + (3 * A + 4 * B + C) * (σ - 1)) ≤ σ - β :=
  dlvp_zero_gap hσ hβ (by positivity) (dlvp_assemble hdlvp hpole hzero hgrowth)

#print axioms dlvp_gap_of_four_bounds
