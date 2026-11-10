import Mathlib

/-!
# S3025 — detector sensitivity: certified residual spread ⇒ cap on any extra (non-Pick) term (lane 2)

S3023 certifies, for every grid step j, D_j := s_j − P_j(T) ∈ [d, d + η] (η = spread, ≤ 1.4e-8) where s_j are the
Arb-certified arithmetic-side steps and P_j(T) the certified critical-line zero sum below T.  If the truth is
D_j = t_j + w·e_j with t_j the (Pick) zero tail above T, varying across the grid by at most θ, and w·e_j the step profile of
one extra term (e.g. an off-line zero pair) with known profile e, then any two grid steps i, k pin the weight:
* `weight_mul_le` — |w|·|e_i − e_k| ≤ η + θ;
* `weight_le`     — e_i ≠ e_k ⇒ |w| ≤ (η + θ)/|e_i − e_k|.
* `weight_le_of_all` — the best pair over the grid (any i, k) gives the bound.
-/

namespace CCMDetectorSensitivity

theorem weight_mul_le {D t e : ℕ → ℝ} {w d η τ θ : ℝ} {i k : ℕ}
    (hD : ∀ j, D j = t j + w * e j)
    (hDi : d ≤ D i ∧ D i ≤ d + η) (hDk : d ≤ D k ∧ D k ≤ d + η)
    (hti : τ ≤ t i ∧ t i ≤ τ + θ) (htk : τ ≤ t k ∧ t k ≤ τ + θ) :
    |w| * |e i - e k| ≤ η + θ := by
  rw [← abs_mul, mul_sub]
  have h1 := hD i
  have h2 := hD k
  rw [abs_le]
  constructor <;> linarith [hDi.1, hDi.2, hDk.1, hDk.2, hti.1, hti.2, htk.1, htk.2]

theorem weight_le {D t e : ℕ → ℝ} {w d η τ θ : ℝ} {i k : ℕ}
    (hD : ∀ j, D j = t j + w * e j)
    (hDi : d ≤ D i ∧ D i ≤ d + η) (hDk : d ≤ D k ∧ D k ≤ d + η)
    (hti : τ ≤ t i ∧ t i ≤ τ + θ) (htk : τ ≤ t k ∧ t k ≤ τ + θ) (hne : e i ≠ e k) :
    |w| ≤ (η + θ) / |e i - e k| := by
  have hpos : 0 < |e i - e k| := abs_pos.mpr (sub_ne_zero.mpr hne)
  rw [le_div_iff₀ hpos]
  exact weight_mul_le hD hDi hDk hti htk

theorem weight_le_of_all {D t e : ℕ → ℝ} {w d η τ θ : ℝ} {n : ℕ}
    (hD : ∀ j, D j = t j + w * e j)
    (hDall : ∀ j, j ≤ n → d ≤ D j ∧ D j ≤ d + η)
    (htall : ∀ j, j ≤ n → τ ≤ t j ∧ t j ≤ τ + θ) :
    ∀ i k, i ≤ n → k ≤ n → |w| * |e i - e k| ≤ η + θ := by
  intro i k hi hk
  exact weight_mul_le hD (hDall i hi) (hDall k hk) (htall i hi) (htall k hk)

end CCMDetectorSensitivity

#print axioms CCMDetectorSensitivity.weight_mul_le
#print axioms CCMDetectorSensitivity.weight_le
#print axioms CCMDetectorSensitivity.weight_le_of_all
