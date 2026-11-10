import ZetaImagAxisBound
import ConvexSphere

/-!
# S3006 — ζ is O(log t) on 1 ≤ Re z ≤ 2 (lane 1: the Re z > 1 part of the convexity sphere bound)

* `divdiff_le_log` — x > 1, N ≥ 1 ⇒ (N^{1−x} − 1)/(1 − x) ≤ log N   (N^{−u} ≥ 1 − u·log N).
* `norm_zeta_le_log_right` — 1 < Re z ≤ 2, |Im z| ≥ 2 ⇒ ‖ζ(z)‖ ≤ log|Im z| + 6
  (forge two-sided N-term bound S1303 with N = ⌊|Im z|⌋).
Together with S2977 (Re z = 1) and S3004 (0 ≤ Re z ≤ 1): on the pipeline sphere, log|ζ| ≤ ((1−σ₀+R)/2)·log|γ| + O(log log|γ|)
everywhere — the input for growth (8/R)·log Ms ≈ 4·log|γ|.
-/

open Complex

namespace ZetaReGeOne

theorem divdiff_le_log {N : ℕ} (hN : 1 ≤ N) {x : ℝ} (hx : 1 < x) :
    ((N : ℝ) ^ (1 - x) - 1) / (1 - x) ≤ Real.log N := by
  have hNR : (1 : ℝ) ≤ N := by exact_mod_cast hN
  have hNpos : (0 : ℝ) < N := by linarith
  set u := x - 1 with hu
  have hu0 : 0 < u := by linarith
  have hpow : (N : ℝ) ^ (1 - x) = Real.exp (-(u * Real.log N)) := by
    rw [Real.rpow_def_of_pos hNpos]; congr 1; rw [hu]; ring
  have hexp : 1 - u * Real.log N ≤ Real.exp (-(u * Real.log N)) := by
    have := Real.add_one_le_exp (-(u * Real.log N)); linarith
  rw [hpow, show (1 - x) = -u by rw [hu]; ring, div_neg, ← neg_div, div_le_iff₀ hu0]
  linarith

theorem norm_zeta_le_log_right {z : ℂ} (h1 : 1 < z.re) (h2 : z.re ≤ 2) (ht : 2 ≤ |z.im|) :
    ‖riemannZeta z‖ ≤ Real.log |z.im| + 6 := by
  set t := |z.im| with htdef
  have hN : 1 ≤ ⌊t⌋₊ := Nat.le_floor (by norm_num; linarith)
  have hb := norm_riemannZeta_le_offline hN (by linarith : 0 < z.re) (by linarith : z.re ≠ 1)
  have hNx : (⌊t⌋₊ : ℝ) ≤ t := Nat.floor_le (by linarith)
  have hNlt : t < (⌊t⌋₊ : ℝ) + 1 := Nat.lt_floor_add_one t
  have hNpos : (0 : ℝ) < (⌊t⌋₊ : ℝ) := by linarith
  have hN1 : (1 : ℝ) ≤ (⌊t⌋₊ : ℝ) := by exact_mod_cast hN
  have hA := divdiff_le_log hN h1
  have hlog : Real.log (⌊t⌋₊ : ℝ) ≤ Real.log t := Real.log_le_log hNpos hNx
  have hp1 : (⌊t⌋₊ : ℝ) ^ (1 - z.re) ≤ 1 := Real.rpow_le_one_of_one_le_of_nonpos hN1 (by linarith)
  have hp0 : 0 ≤ (⌊t⌋₊ : ℝ) ^ (1 - z.re) := Real.rpow_nonneg hNpos.le _
  have hz1 : t ≤ ‖z - 1‖ := by
    have := Complex.abs_im_le_norm (z - 1); simpa [htdef] using this
  have hB : (⌊t⌋₊ : ℝ) ^ (1 - z.re) / ‖z - 1‖ ≤ 1 / 2 := by
    rw [div_le_iff₀ (by linarith)]; nlinarith
  have hzn : ‖z‖ ≤ 2 + t := by
    have := Complex.norm_le_abs_re_add_abs_im z
    rw [abs_of_pos (by linarith : (0 : ℝ) < z.re)] at this; linarith
  have hq : (⌊t⌋₊ : ℝ) ^ (-z.re) ≤ (⌊t⌋₊ : ℝ)⁻¹ := by
    rw [← Real.rpow_neg_one]
    exact Real.rpow_le_rpow_of_exponent_le hN1 (by linarith)
  have hq0 : 0 ≤ (⌊t⌋₊ : ℝ) ^ (-z.re) := Real.rpow_nonneg hNpos.le _
  have hC : ‖z‖ * ((⌊t⌋₊ : ℝ) ^ (-z.re) / z.re) ≤ 4 := by
    have hdiv : (⌊t⌋₊ : ℝ) ^ (-z.re) / z.re ≤ (⌊t⌋₊ : ℝ)⁻¹ := by
      refine le_trans (div_le_self hq0 (by linarith)) hq
    have h4 : (2 + t) * (⌊t⌋₊ : ℝ)⁻¹ ≤ 4 := by
      rw [← div_eq_mul_inv, div_le_iff₀ hNpos]; linarith
    calc ‖z‖ * ((⌊t⌋₊ : ℝ) ^ (-z.re) / z.re) ≤ (2 + t) * (⌊t⌋₊ : ℝ)⁻¹ :=
          mul_le_mul hzn hdiv (by positivity) (by linarith)
      _ ≤ 4 := h4
  linarith

end ZetaReGeOne

#print axioms ZetaReGeOne.divdiff_le_log
#print axioms ZetaReGeOne.norm_zeta_le_log_right
