import ZetaLogSharp

/-!
# §1299 — `|t|`-explicit ζ log bound: `‖ζ(σ+it)‖ ≤ e·log(|t|+1) + 2 + e + e/σ`

Instantiating the sharp bound §1298 at `N = ⌈|t|⌉` (`t = Im s`).  For `0 < σ < 1`, `|t| ≥ 1`, and the region
hypothesis `(1-σ)·log(|t|+1) ≤ 1` (i.e. `σ ≥ 1 − 1/log(|t|+1)`):

  `‖ζ(σ+it)‖ ≤ e·log(|t|+1) + 2 + e + e/σ`.

Uses `|t| ≤ ⌈|t|⌉ < |t|+1` (so `log⌈|t|⌉ ≤ log(|t|+1)`), `‖s-1‖ ≥ |Im s| = |t| ≥ 1` (so `e/‖s-1‖ ≤ e`).  This is
the genuine, unconditional `‖ζ(σ+it)‖ ≲ log|t|` growth bound in the critical strip — the standard input to the
de la Vallée Poussin / Borel–Carathéodory sharp zero-free region argument.  EXIT0 / std-axioms.
-/

theorem norm_riemannZeta_le_log_im {s : ℂ} (hs0 : 0 < s.re) (hs1 : s.re < 1)
    (hst : 1 ≤ |s.im|) (hreg : (1 - s.re) * Real.log (|s.im| + 1) ≤ 1) :
    ‖riemannZeta s‖ ≤ Real.exp 1 * Real.log (|s.im| + 1) + 2 + Real.exp 1 + Real.exp 1 / s.re := by
  have ht0 : (0 : ℝ) < |s.im| := lt_of_lt_of_le one_pos hst
  have hN1 : 1 ≤ ⌈|s.im|⌉₊ := Nat.one_le_ceil_iff.mpr ht0
  have hNR : (1 : ℝ) ≤ (⌈|s.im|⌉₊ : ℝ) := by exact_mod_cast hN1
  have hNpos : (0 : ℝ) < (⌈|s.im|⌉₊ : ℝ) := lt_of_lt_of_le one_pos hNR
  have hleN : |s.im| ≤ (⌈|s.im|⌉₊ : ℝ) := Nat.le_ceil _
  have hNlt : (⌈|s.im|⌉₊ : ℝ) < |s.im| + 1 := Nat.ceil_lt_add_one ht0.le
  have hlogN : Real.log (⌈|s.im|⌉₊ : ℝ) ≤ Real.log (|s.im| + 1) := Real.log_le_log hNpos hNlt.le
  have h1ms : (0 : ℝ) ≤ 1 - s.re := by linarith
  have hregN : (1 - s.re) * Real.log (⌈|s.im|⌉₊ : ℝ) ≤ 1 :=
    le_trans (mul_le_mul_of_nonneg_left hlogN h1ms) hreg
  have hsharp := norm_riemannZeta_le_log_sharp hN1 hs0 hs1 hregN hleN
  have hd1 : |s.im| ≤ ‖s - 1‖ := by
    have h := Complex.abs_im_le_norm (s - 1)
    simpa [Complex.sub_im, Complex.one_im] using h
  have hdpos : (0 : ℝ) < ‖s - 1‖ := lt_of_lt_of_le ht0 hd1
  have h1d : (1 : ℝ) ≤ ‖s - 1‖ := le_trans hst hd1
  have hterm3 : Real.exp 1 / ‖s - 1‖ ≤ Real.exp 1 := by
    rw [div_le_iff₀ hdpos]; nlinarith [Real.exp_pos 1, h1d]
  have hterm1 : Real.exp 1 * Real.log (⌈|s.im|⌉₊ : ℝ) ≤ Real.exp 1 * Real.log (|s.im| + 1) :=
    mul_le_mul_of_nonneg_left hlogN (Real.exp_pos 1).le
  linarith [hsharp, hterm3, hterm1]

#print axioms norm_riemannZeta_le_log_im
