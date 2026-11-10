import ExistsHoloBranchBall
import ReEqLogNormOfExp
import Mathlib.Analysis.Complex.Schwarz

/-!
# S3018 — SHARP Borel–Carathéodory derivative bound via the Schwarz lemma (lane 1, structural lever)

The pipeline's normalized BC (§1367–§1369) gives `-(4·A/R) ≤ Re f'(0)`, A = M − Re f(0).  The sharp constant is 2:
with g = f − f(0) (Re g ≤ A) the Cayley-type map φ = g/(2A − g) sends `ball 0 R` into the closed unit disk
(|g| ≤ |2A − g| ⟺ Re g ≤ A), φ(0) = 0, φ'(0) = f'(0)/(2A); Mathlib's Schwarz lemma `norm_deriv_le_div_of_mapsTo_ball`
gives |φ'(0)| ≤ 1/R, i.e. |f'(0)| ≤ 2A/R.
* `re_deriv_ge_of_re_le_sharp`        — −(2A/R) ≤ Re f'(0)   (vs 4A/R in §1367).
* `re_logDeriv_ge_of_bounded_sharp`   — −(2(log M' − log‖h 0‖)/R) ≤ Re logDeriv h 0  (vs §1368).
* `re_logDeriv_ge_of_bounded_sharp_at` — centre-c version (vs §1369).
Propagating the 2 through the §1378→§1444 chain halves every growth term: κ 4 → 2, c₀ ≈ 7.16e-3 → ≈ 1.43e-2 (target).
-/

open Complex Metric Set

namespace SharpBC

theorem re_deriv_ge_of_re_le_sharp {f : ℂ → ℂ} {R M : ℝ} (hR : 0 < R)
    (hM' : 0 < M - (f 0).re)
    (hf : DifferentiableOn ℂ f (ball 0 R))
    (hre : ∀ z ∈ ball (0 : ℂ) R, (f z).re ≤ M) :
    -(2 * (M - (f 0).re) / R) ≤ (deriv f 0).re := by
  set A : ℝ := M - (f 0).re with hA
  set g : ℂ → ℂ := fun z => f z - f 0 with hg
  have hgre : ∀ z ∈ ball (0 : ℂ) R, (g z).re ≤ A := by
    intro z hz; simp only [hg, Complex.sub_re, hA]; linarith [hre z hz]
  have hden : ∀ z ∈ ball (0 : ℂ) R, (2 * (A : ℂ) - g z) ≠ 0 := by
    intro z hz h
    have := congrArg Complex.re h
    simp at this
    linarith [hgre z hz]
  set φ : ℂ → ℂ := fun z => g z / (2 * (A : ℂ) - g z) with hφ
  have hgd : DifferentiableOn ℂ g (ball 0 R) := hf.sub_const _
  have hφd : DifferentiableOn ℂ φ (ball 0 R) :=
    hgd.div ((differentiableOn_const _).sub hgd) hden
  have h0 : (0 : ℂ) ∈ ball (0 : ℂ) R := mem_ball_self hR
  have hg0 : g 0 = 0 := by simp [hg]
  have hφ0 : φ 0 = 0 := by simp [hφ, hg0]
  have hmaps : MapsTo φ (ball 0 R) (closedBall (φ 0) 1) := by
    intro z hz
    rw [hφ0, mem_closedBall, dist_zero_right]
    simp only [hφ]
    rw [norm_div, div_le_one (norm_pos_iff.mpr (hden z hz))]
    have hre' := hgre z hz
    set w := g z
    have hA0 : 0 < A := hM'
    rw [← sq_le_sq₀ (norm_nonneg _) (norm_nonneg _), Complex.sq_norm, Complex.sq_norm,
      Complex.normSq_apply, Complex.normSq_apply]
    simp
    nlinarith
  have hS := Complex.norm_deriv_le_div_of_mapsTo_ball hφd hmaps hR
  -- derivative of φ at 0
  have hfd0 : HasDerivAt f (deriv f 0) 0 :=
    ((hf 0 h0).differentiableAt (isOpen_ball.mem_nhds h0)).hasDerivAt
  have hgd0 : HasDerivAt g (deriv f 0) 0 := hfd0.sub_const _
  have hdd0 : HasDerivAt (fun z => 2 * (A : ℂ) - g z) (-(deriv f 0)) 0 := hgd0.const_sub _
  have hden0 : (2 * (A : ℂ) - g 0) ≠ 0 := hden 0 h0
  have hφd0 := hgd0.div hdd0 hden0
  have hAc : (2 * (A : ℂ)) ≠ 0 := by
    have : (0 : ℝ) < A := hM'
    exact_mod_cast (by positivity : (2 * A : ℝ) ≠ 0)
  have hderiv : deriv φ 0 = deriv f 0 / (2 * (A : ℂ)) := by
    rw [show φ = (g / fun z => 2 * (A : ℂ) - g z) from rfl, hφd0.deriv, hg0]
    have hA0c : (A : ℂ) ≠ 0 := by exact_mod_cast (ne_of_gt hM')
    field_simp
    ring
  rw [hderiv, norm_div] at hS
  have hnA : ‖(2 * (A : ℂ))‖ = 2 * A := by
    have : (0 : ℝ) < A := hM'
    rw [show (2 * (A : ℂ)) = ((2 * A : ℝ) : ℂ) by push_cast; ring, Complex.norm_real, Real.norm_eq_abs,
      abs_of_pos (by positivity)]
  rw [hnA, div_le_div_iff₀ (by have : (0:ℝ) < A := hM'; positivity) hR] at hS
  have hre_le : -‖deriv f 0‖ ≤ (deriv f 0).re := by
    have := Complex.abs_re_le_norm (deriv f 0); rw [abs_le] at this; linarith
  have : ‖deriv f 0‖ ≤ 2 * A / R := by rw [le_div_iff₀ hR]; linarith
  linarith

theorem re_logDeriv_ge_of_bounded_sharp {R M' : ℝ} (hR : 0 < R) {h : ℂ → ℂ}
    (hh : DifferentiableOn ℂ h (ball 0 R)) (hh0 : ∀ z ∈ ball 0 R, h z ≠ 0)
    (hb : ∀ z ∈ ball 0 R, ‖h z‖ ≤ M') (hh0lt : ‖h 0‖ < M') :
    -(2 * (Real.log M' - Real.log ‖h 0‖) / R) ≤ (logDeriv h 0).re := by
  obtain ⟨f, hfd, hfe⟩ := exists_holo_branch_on_ball hR hh hh0
  have h0ball : (0 : ℂ) ∈ ball 0 R := mem_ball_self hR
  have hh00 : h 0 ≠ 0 := hh0 0 h0ball
  have hfdiff : DifferentiableOn ℂ f (ball 0 R) :=
    fun z hz => (hfd z hz).differentiableAt.differentiableWithinAt
  have hf0 : (f 0).re = Real.log ‖h 0‖ := re_eq_log_norm_of_exp (hfe 0 h0ball)
  have hre : ∀ z ∈ ball (0 : ℂ) R, (f z).re ≤ Real.log M' := by
    intro z hz
    rw [re_eq_log_norm_of_exp (hfe z hz)]
    exact Real.log_le_log (norm_pos_iff.mpr (hh0 z hz)) (hb z hz)
  have hM0' : 0 < Real.log M' - (f 0).re := by
    rw [hf0]
    have := Real.log_lt_log (norm_pos_iff.mpr hh00) hh0lt
    linarith
  have hbc := re_deriv_ge_of_re_le_sharp hR hM0' hfdiff hre
  rw [hf0] at hbc
  rwa [(hfd 0 h0ball).deriv] at hbc

theorem re_logDeriv_ge_of_bounded_sharp_at {c : ℂ} {R M' : ℝ} (hR : 0 < R) {h : ℂ → ℂ}
    (hh : DifferentiableOn ℂ h (ball c R)) (hh0 : ∀ z ∈ ball c R, h z ≠ 0)
    (hb : ∀ z ∈ ball c R, ‖h z‖ ≤ M') (hh0lt : ‖h c‖ < M') :
    -(2 * (Real.log M' - Real.log ‖h c‖) / R) ≤ (logDeriv h c).re := by
  set sh : ℂ → ℂ := fun w => c + w with hsh
  have hmaps : Set.MapsTo sh (ball 0 R) (ball c R) := by
    intro w hw
    rw [mem_ball, dist_zero_right] at hw
    rw [mem_ball, Complex.dist_eq]
    simpa [hsh] using hw
  have hshdiff : ∀ x, DifferentiableAt ℂ sh x := by
    intro x; rw [hsh]; fun_prop
  have hshon : DifferentiableOn ℂ sh (ball (0 : ℂ) R) :=
    fun x _ => (hshdiff x).differentiableWithinAt
  have hh' : DifferentiableOn ℂ (h ∘ sh) (ball 0 R) := hh.comp hshon hmaps
  have hh0' : ∀ z ∈ ball (0 : ℂ) R, (h ∘ sh) z ≠ 0 := fun z hz => hh0 _ (hmaps hz)
  have hb' : ∀ z ∈ ball (0 : ℂ) R, ‖(h ∘ sh) z‖ ≤ M' := fun z hz => hb _ (hmaps hz)
  have h00 : (h ∘ sh) 0 = h c := by simp [hsh]
  have hb0lt : ‖(h ∘ sh) 0‖ < M' := by rw [h00]; exact hh0lt
  have hdc : DifferentiableAt ℂ h c :=
    (hh c (mem_ball_self hR)).differentiableAt (isOpen_ball.mem_nhds (mem_ball_self hR))
  have hld : logDeriv (h ∘ sh) 0 = logDeriv h c := by
    have hcomp := logDeriv_comp (by simpa [hsh] using hdc) (hshdiff 0)
    have hsh0 : sh 0 = c := by simp [hsh]
    have hderiv : deriv sh 0 = 1 := by
      simp [hsh]
    rw [hcomp, hsh0, hderiv, mul_one]
  have hbc := re_logDeriv_ge_of_bounded_sharp hR hh' hh0' hb' hb0lt
  rw [h00, hld] at hbc
  exact hbc

end SharpBC

#print axioms SharpBC.re_deriv_ge_of_re_le_sharp
#print axioms SharpBC.re_logDeriv_ge_of_bounded_sharp
#print axioms SharpBC.re_logDeriv_ge_of_bounded_sharp_at
