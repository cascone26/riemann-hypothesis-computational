import ReSumZpowGe
import ReOneDivSameIm

/-!
# §1335 — Multiplicity `hzero` from the decomposition: `Re(-F) ≤ B - 1/(σ-β)`

The multiplicity version of §1314.  Given `F = ∑_{ρ∈S} m ρ/((σ+iγ)-ρ) + g` (the decomposition §1333 produces from
`ζ = P·h`, with `m` the divisor multiplicities and `g = logDeriv h`, `Re g ≥ -B`), the kept zero `ρ₀ = β+iγ` with
`m ρ₀ ≥ 1` contributes `m ρ₀/(σ-β) ≥ 1/(σ-β)` (§1334 drop-zeros + §1313 value), so `Re(-F) ≤ B - 1/(σ-β)`.  This is
the exact `hzero`/`T₁` input the region assembly §1326 consumes — now for a zero of any multiplicity.
EXIT0 / std-axioms.
-/

open Complex

theorem re_neg_le_of_zpow_decomp {σ β γ B : ℝ} {g F : ℂ} {S : Finset ℂ} {m : ℂ → ℤ}
    (hρ₀ : ((β : ℂ) + I * (γ : ℝ)) ∈ S) (hleft : ∀ ρ ∈ S, ρ.re < σ)
    (hm : ∀ ρ ∈ S, 0 ≤ m ρ) (hm₀ : 1 ≤ m ((β : ℂ) + I * (γ : ℝ)))
    (hσβ : β < σ) (hgB : -B ≤ g.re)
    (hdecomp : F = (∑ ρ ∈ S, (m ρ : ℂ) / (((σ : ℂ) + I * (γ : ℝ)) - ρ)) + g) :
    (-F).re ≤ B - 1 / (σ - β) := by
  have hsre : (((σ : ℂ) + I * (γ : ℝ))).re = σ := by simp
  have hleft' : ∀ ρ ∈ S, ρ.re < (((σ : ℂ) + I * (γ : ℝ))).re :=
    fun ρ hρ => by rw [hsre]; exact hleft ρ hρ
  have hdrop := re_sum_zpow_div_ge_single (s := (σ : ℂ) + I * (γ : ℝ)) hρ₀ hm hleft'
  rw [re_one_div_sub_same_im σ β γ] at hdrop
  have hσβ0 : (0 : ℝ) < σ - β := by linarith
  have hpos : (0 : ℝ) ≤ 1 / (σ - β) := (div_pos one_pos hσβ0).le
  have hstep : (1 : ℝ) / (σ - β) ≤ (m ((β : ℂ) + I * (γ : ℝ)) : ℝ) * (1 / (σ - β)) :=
    le_mul_of_one_le_left hpos (by exact_mod_cast hm₀)
  rw [hdecomp, Complex.neg_re, Complex.add_re]
  linarith [hdrop, hstep, hgB]

#print axioms re_neg_le_of_zpow_decomp
