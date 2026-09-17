import DropZeros

/-!
# §1334 — Multiplicity drop-zeros: `m ρ₀·Re(1/(s-ρ₀)) ≤ Re(∑ m ρ/(s-ρ))`

The multiplicity-weighted version of §1147/§1312.  For a finite zero set `S ∋ ρ₀` all left of `Re s`, with
non-negative integer multiplicities `m ρ ≥ 0`, every term `m ρ·Re(1/(s-ρ)) ≥ 0` (§1142), so the sum dominates the
single kept term.  At `s = σ+iγ`, `ρ₀ = β+iγ` this yields `m ρ₀/(σ-β) ≥ 1/(σ-β)` (`m ρ₀ ≥ 1`) — the negative
contribution driving the region, for a zero of any multiplicity.  EXIT0 / std-axioms.
-/

open Complex

theorem re_sum_zpow_div_ge_single {s : ℂ} {S : Finset ℂ} {m : ℂ → ℤ} {ρ₀ : ℂ} (hρ₀ : ρ₀ ∈ S)
    (hm : ∀ ρ ∈ S, 0 ≤ m ρ) (hleft : ∀ ρ ∈ S, ρ.re < s.re) :
    (m ρ₀ : ℝ) * (1 / (s - ρ₀)).re ≤ (∑ ρ ∈ S, (m ρ : ℂ) / (s - ρ)).re := by
  have hterm : ∀ ρ ∈ S, ((m ρ : ℂ) / (s - ρ)).re = (m ρ : ℝ) * (1 / (s - ρ)).re := by
    intro ρ _
    rw [div_eq_mul_one_div, ← Complex.ofReal_intCast, Complex.re_ofReal_mul]
  rw [Complex.re_sum, ← hterm ρ₀ hρ₀]
  refine Finset.single_le_sum (f := fun ρ => ((m ρ : ℂ) / (s - ρ)).re) (fun ρ hρ => ?_) hρ₀
  rw [hterm ρ hρ]
  exact mul_nonneg (by exact_mod_cast hm ρ hρ) (re_one_div_sub_pos s ρ (hleft ρ hρ)).le

#print axioms re_sum_zpow_div_ge_single
