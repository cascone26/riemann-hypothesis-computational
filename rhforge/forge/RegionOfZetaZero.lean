import ZetaDecompAt
import ZetaDivisorNonneg
import OneLeDivisorZeta
import ZetaZeroMemSupport
import DivisorSupportReLeOne
import RegionZpow

/-!
# §1355 ★ Form-match: concrete quantitative region for an actual ζ zero

Threading the machine-checked discrete inputs into the multiplicity master theorem §1336
(`region_of_zpow_decomp_and_bounds`), instantiated at the real Riemann zeta:

  * `hdecomp`  ← §1349 (`logDeriv_riemannZeta_decomp_at`), evaluated at the DLVP test point `z = σ_K + iγ`,
                 `σ_K = 1 + 1/(2K)`;
  * `hm`       ← §1350 (`divisor_riemannZeta_nonneg`) on the pole-free ball;
  * `hm₀`      ← §1351 (`one_le_divisor_riemannZeta_at_zero`) at the kept zero;
  * `hρ₀`      ← §1353 (`zeta_zero_mem_divisor_support`);
  * `hleft`    ← §1354 (`divisor_support_re_le_one`), `Re ρ ≤ 1 < σ_K`.

Everything structural is discharged internally.  What remains is supplied as three clean hypotheses:
`hgB` (the Borel–Carathéodory bound `Re (logDeriv g) ≥ -B` for *any* analytic nonvanishing factor `g` — exactly
§1352's output once `‖ζ/P‖ ≤ M'` is furnished by the Jensen sub-arc) and the DLVP `T₀/T₂` positivity/pole/growth
bounds (`hdlvp`/`hpole`/`hgrowth`, from §1133/§1138).  Conclusion: any zero `β+iγ` of ζ in the ball obeys
`β ≤ 1 - 1/(14K)`.  This crystallizes the concrete zero-free region into a single theorem over the actual ζ.
EXIT0 / std-axioms.
-/

open Complex Metric

theorem region_of_riemannZeta_zero {c : ℂ} {R γ β A B C T₀ T₂ : ℝ}
    (hball : (1 : ℂ) ∉ ball c R)
    (hK : 0 < 3 * A + 4 * B + C) (hA : 0 ≤ A) (hB : 0 ≤ B) (hC : 0 ≤ C)
    (hβ : β < 1 + 1 / (2 * (3 * A + 4 * B + C)))
    (hρ_ball : ((β : ℂ) + I * (γ : ℝ)) ∈ ball c R)
    (hζρ : riemannZeta ((β : ℂ) + I * (γ : ℝ)) = 0)
    (hz_test : (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ)) ∈ ball c R)
    {z₀ : ℂ} (hz₀ : z₀ ∈ ball c R) (hz₀1 : z₀ ≠ 1) (hz₀ne : riemannZeta z₀ ≠ 0)
    (hgB : ∀ g : ℂ → ℂ, AnalyticOnNhd ℂ g (ball c R) → (∀ u : ball c R, g u ≠ 0) →
      -B ≤ (logDeriv g (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ))).re)
    (hdlvp : 0 ≤ 3 * T₀ + 4 * (-logDeriv riemannZeta
        (((1 + 1 / (2 * (3 * A + 4 * B + C)) : ℝ) : ℂ) + I * (γ : ℝ))).re + T₂)
    (hpole : T₀ ≤ 1 / ((1 + 1 / (2 * (3 * A + 4 * B + C))) - 1) + A)
    (hgrowth : T₂ ≤ C) :
    β ≤ 1 - 1 / (14 * (3 * A + 4 * B + C)) := by
  set K := 3 * A + 4 * B + C with hKdef
  set σK : ℝ := 1 + 1 / (2 * K) with hσK
  set z : ℂ := ((σK : ℝ) : ℂ) + I * (γ : ℝ) with hzdef
  have hzre : 1 < z.re := by
    have hre : z.re = σK := by
      simp [hzdef, Complex.add_re, Complex.ofReal_re, Complex.ofReal_im]
    rw [hre, hσK]
    have : 0 < 1 / (2 * K) := by positivity
    linarith
  obtain ⟨g, hg_ana, hg_ne, hdecomp_eq⟩ :=
    logDeriv_riemannZeta_decomp_at hz_test hzre hz₀ hz₀1 hz₀ne
  have hρ₀ : ((β : ℂ) + I * (γ : ℝ)) ∈ (divisor_riemannZeta_ball_support_finite c R).toFinset :=
    zeta_zero_mem_divisor_support hball hρ_ball hζρ hz₀ hz₀ne
  have hleft : ∀ ρ ∈ (divisor_riemannZeta_ball_support_finite c R).toFinset, ρ.re < σK := by
    intro ρ hρmem
    have hle := divisor_support_re_le_one hρmem
    rw [hσK]
    have : 0 < 1 / (2 * K) := by positivity
    linarith
  have hm : ∀ ρ ∈ (divisor_riemannZeta_ball_support_finite c R).toFinset,
      0 ≤ MeromorphicOn.divisor riemannZeta (ball c R) ρ := by
    intro ρ hρmem
    rw [Set.Finite.mem_toFinset, Function.mem_support] at hρmem
    exact divisor_riemannZeta_nonneg hball
      ((MeromorphicOn.divisor riemannZeta (ball c R)).supportWithinDomain hρmem)
  have hm₀ : 1 ≤ MeromorphicOn.divisor riemannZeta (ball c R) ((β : ℂ) + I * (γ : ℝ)) :=
    one_le_divisor_riemannZeta_at_zero hball hρ_ball hζρ hz₀ hz₀ne
  exact region_of_zpow_decomp_and_bounds hK hA hB hC hβ hρ₀ hleft hm hm₀
    (hgB g hg_ana hg_ne) hdecomp_eq hdlvp hpole hgrowth

#print axioms region_of_riemannZeta_zero
