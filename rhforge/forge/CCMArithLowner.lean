import CCMW02Rank2
import CCMLownerDisplacement

/-!
# S2970 — the whole CCM arithmetic matrix A = W02 − W_R − W_P is ONE Löwner matrix (lane 2)

S2966: `W02 = c(uuᵀ − vvᵀ)`.  S2968: W_R, W_P are Löwner (`(g_j − g_i)/(π(ν_i − ν_j))` off-diagonal).
New: W02 is Löwner too, with generator `g₀ = −(c/4)·v`, because
  `(c/4)(v_i − v_j)/(π(ν_i − ν_j)) = c(L² − 16π²ν_iν_j)/(D_iD_j)`   (D_i = L² + 16π²ν_i²),
via `ν_iD_j − ν_jD_i = (ν_i − ν_j)(L² − 16π²ν_iν_j)`.
* `W02_eq_lowner` — for injective ν: `W02 = lowner (diag W02) (−(c/4)•v) ν`.
* `lowner_neg`, `lowner_sub` — Löwner matrices form a linear family in (d, g).
* `ccm_arith_lowner` — `W02 − lowner d_R g_R ν − lowner d_P g_P ν = lowner (d₀ − d_R − d_P) (g₀ − g_R − g_P) ν`.
* `ccm_arith_displacement` — `[diag ν, A] = (1/π)(𝟙 g_Aᵀ − g_A 𝟙ᵀ)` with `g_A = g₀ − g_R − g_P`.
Numeric cross-check (python parts(20,8)): max |W02 − lowner(g₀)| = 8.3e-17 (scale 7.2).
Reading (Report): the CCM Weil matrix is the divided-difference matrix of ONE function g_A of the frequency, plus a
diagonal — its off-diagonal part is entirely determined by g_A. Positivity questions reduce to (d_A, g_A).
-/

open Matrix Finset Real

namespace CCMArith

variable {ι : Type*} [DecidableEq ι]

omit [DecidableEq ι] in
theorem W02_offdiag {L : ℝ} (hL : 0 < L) (ν : ι → ℝ) {i j : ι} (hij : ν i ≠ ν j) :
    CCMW02.W02 L ν i j
      = ((-(CCMW02.cW L / 4) • CCMW02.vW L ν) j - (-(CCMW02.cW L / 4) • CCMW02.vW L ν) i)
        / (π * (ν i - ν j)) := by
  have hi := (CCMW02.den_pos hL (ν i)).ne'
  have hj := (CCMW02.den_pos hL (ν j)).ne'
  have hne : ν i - ν j ≠ 0 := sub_ne_zero.mpr hij
  have hpi : (π : ℝ) ≠ 0 := Real.pi_ne_zero
  simp only [CCMW02.W02, CCMW02.vW, CCMW02.cW, Pi.smul_apply, smul_eq_mul]
  rw [div_eq_div_iff (mul_ne_zero hi hj) (mul_ne_zero hpi hne)]
  field_simp
  ring

theorem W02_eq_lowner {L : ℝ} (hL : 0 < L) {ν : ι → ℝ} (hν : Function.Injective ν) :
    CCMW02.W02 L ν
      = CCMLowner.lowner (fun i => CCMW02.W02 L ν i i) (-(CCMW02.cW L / 4) • CCMW02.vW L ν) ν := by
  ext i j
  simp only [CCMLowner.lowner]
  by_cases h : i = j
  · subst h; simp
  · rw [if_neg h]
    exact W02_offdiag hL ν (fun e => h (hν e))

theorem lowner_neg (d g ν : ι → ℝ) : -CCMLowner.lowner d g ν = CCMLowner.lowner (-d) (-g) ν := by
  ext i j
  simp only [Matrix.neg_apply, CCMLowner.lowner, Pi.neg_apply]
  by_cases h : i = j
  · simp [h]
  · rw [if_neg h, if_neg h]; ring

theorem lowner_sub (d g d' g' ν : ι → ℝ) :
    CCMLowner.lowner d g ν - CCMLowner.lowner d' g' ν = CCMLowner.lowner (d - d') (g - g') ν := by
  rw [sub_eq_add_neg, lowner_neg, CCMLowner.lowner_add, sub_eq_add_neg, sub_eq_add_neg]

theorem ccm_arith_lowner {L : ℝ} (hL : 0 < L) {ν : ι → ℝ} (hν : Function.Injective ν)
    (dR gR dP gP : ι → ℝ) :
    CCMW02.W02 L ν - CCMLowner.lowner dR gR ν - CCMLowner.lowner dP gP ν
      = CCMLowner.lowner ((fun i => CCMW02.W02 L ν i i) - dR - dP)
          (-(CCMW02.cW L / 4) • CCMW02.vW L ν - gR - gP) ν := by
  conv_lhs => rw [W02_eq_lowner hL hν]
  rw [lowner_sub, lowner_sub]

theorem ccm_arith_displacement [Fintype ι] {L : ℝ} (hL : 0 < L) {ν : ι → ℝ}
    (hν : Function.Injective ν) (dR gR dP gP : ι → ℝ) :
    let A := CCMW02.W02 L ν - CCMLowner.lowner dR gR ν - CCMLowner.lowner dP gP ν
    let gA := -(CCMW02.cW L / 4) • CCMW02.vW L ν - gR - gP
    diagonal ν * A - A * diagonal ν
      = (1 / π) • (vecMulVec (fun _ => (1 : ℝ)) gA - vecMulVec gA (fun _ => (1 : ℝ))) := by
  intro A gA
  simp only [A, gA]
  rw [ccm_arith_lowner hL hν]
  exact CCMLowner.lowner_displacement _ _ hν

end CCMArith

#print axioms CCMArith.W02_offdiag
#print axioms CCMArith.W02_eq_lowner
#print axioms CCMArith.lowner_sub
#print axioms CCMArith.ccm_arith_lowner
#print axioms CCMArith.ccm_arith_displacement
