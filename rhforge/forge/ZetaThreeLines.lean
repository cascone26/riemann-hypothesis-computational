import ZetaEEdges

/-!
# S2984 — Hadamard three-lines for E(s) = (s−1)ζ(s) with a Gaussian localiser (lane 1: convexity route)

Complex-power normalisers (s+2)^{−(3−s)/2} carry an e^{−arg·t/2} ≈ e^{−π|t|/4} factor that is awkward on both half-strips.
Instead localise at height T with the entire Gaussian G_T(s) = exp((s − iT)²), |G_T(σ+it)| = e^{σ² − (t−T)²}, and use
Mathlib's Hadamard three-lines theorem (`Complex.HadamardThreeLines.norm_le_interp_of_mem_verticalClosedStrip₀₁'`):
* `norm_gauss` — ‖exp((s − iT)²)‖ = exp(Re(s)² − (Im(s) − T)²).
* `zetaE_three_lines` — if F = E·G_T is bounded on the closed strip 0 ≤ Re ≤ 1 and the localised edges satisfy
  ‖E(it)‖e^{−(t−T)²} ≤ a, ‖E(1+it)‖e^{1−(t−T)²} ≤ b (all t), then for 0 ≤ σ ≤ 1:
  ‖E(σ + iT)‖·e^{σ²} ≤ a^{1−σ} b^{σ}.
With a ≍ T^{3/2} log T (S2982 left edge) and b ≍ T log T (S2982 right edge) this gives |ζ(σ+iT)| ≪ T^{(1−σ)/2} log T.
Open input (stated as hypothesis `hB`): boundedness of E·G_T on the closed strip — needs a polynomial bound for ζ near
Re s = 0 (functional equation + a Gamma-ratio lower bound); next lane-1 brick.
-/

open Complex Set

namespace ZetaThreeLines

theorem norm_gauss (s : ℂ) (T : ℝ) :
    ‖Complex.exp ((s - (T : ℂ) * I) ^ 2)‖ = Real.exp (s.re ^ 2 - (s.im - T) ^ 2) := by
  rw [Complex.norm_exp]
  congr 1
  simp only [sq, Complex.mul_re, Complex.sub_re, Complex.sub_im, Complex.mul_im,
    Complex.ofReal_re, Complex.ofReal_im, Complex.I_re, Complex.I_im]
  ring

theorem zetaE_three_lines (T a b σ : ℝ) (hσ0 : 0 ≤ σ) (hσ1 : σ ≤ 1)
    (hB : BddAbove ((norm ∘ fun s => ZetaEntire.zetaE s * Complex.exp ((s - (T : ℂ) * I) ^ 2)) ''
      Complex.HadamardThreeLines.verticalClosedStrip 0 1))
    (ha : ∀ t : ℝ, ‖ZetaEntire.zetaE ((t : ℂ) * I)‖ * Real.exp (-(t - T) ^ 2) ≤ a)
    (hb : ∀ t : ℝ, ‖ZetaEntire.zetaE (1 + (t : ℂ) * I)‖ * Real.exp (1 - (t - T) ^ 2) ≤ b) :
    ‖ZetaEntire.zetaE ((σ : ℂ) + (T : ℂ) * I)‖ * Real.exp (σ ^ 2) ≤ a ^ (1 - σ) * b ^ σ := by
  set F : ℂ → ℂ := fun s => ZetaEntire.zetaE s * Complex.exp ((s - (T : ℂ) * I) ^ 2) with hF
  have hd : DiffContOnCl ℂ F (Complex.HadamardThreeLines.verticalStrip 0 1) := by
    have : Differentiable ℂ F := by
      rw [hF]
      exact ZetaEntire.zetaE_differentiable.mul
        ((((differentiable_id.sub_const _)).pow 2).cexp)
    exact this.diffContOnCl
  have hz : ((σ : ℂ) + (T : ℂ) * I) ∈ Complex.HadamardThreeLines.verticalClosedStrip 0 1 := by
    simp [Complex.HadamardThreeLines.verticalClosedStrip, hσ0, hσ1]
  have hedge0 : ∀ z ∈ re ⁻¹' {(0 : ℝ)}, ‖F z‖ ≤ a := by
    intro z hz0
    have hre : z.re = 0 := hz0
    have hzeq : z = ((z.im : ℝ) : ℂ) * I := by
      apply Complex.ext <;> simp [hre]
    rw [hF]; simp only
    rw [norm_mul, norm_gauss, hre, hzeq]
    have := ha z.im
    simpa [hre] using this
  have hedge1 : ∀ z ∈ re ⁻¹' {(1 : ℝ)}, ‖F z‖ ≤ b := by
    intro z hz1
    have hre : z.re = 1 := hz1
    have hzeq : z = 1 + ((z.im : ℝ) : ℂ) * I := by
      apply Complex.ext <;> simp [hre]
    rw [hF]; simp only
    rw [norm_mul, norm_gauss, hre, hzeq]
    have := hb z.im
    simpa [hre] using this
  have key := Complex.HadamardThreeLines.norm_le_interp_of_mem_verticalClosedStrip₀₁' F hz hd hB
    hedge0 hedge1
  have hre : ((σ : ℂ) + (T : ℂ) * I).re = σ := by simp
  have him : ((σ : ℂ) + (T : ℂ) * I).im = T := by simp
  rw [hF] at key; simp only at key
  rw [norm_mul, norm_gauss, hre, him, sub_self] at key
  simpa using key

end ZetaThreeLines

#print axioms ZetaThreeLines.norm_gauss
#print axioms ZetaThreeLines.zetaE_three_lines
