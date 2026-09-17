import Mathlib.Analysis.SpecialFunctions.Complex.Circle
import Mathlib.Tactic

/-!
# §1781 — Joukowski on-circle criterion: z + z⁻¹ real in [-2,2] ⟺ |z| = 1

The finite real-algebra substitution under the Chebyshev / Schur-Cohn on-circle test.  For the
self-inversive (palindromic) eigenvector polynomial of §1776, the map `z ↦ z + z⁻¹` sends the unit
circle onto the real segment `[-2,2]`; this brick proves the exact equivalence

    (z + z⁻¹ ∈ ℝ  ∧  |z + z⁻¹| ≤ 2)   ⟺   |z| = 1        (z ≠ 0).

It reduces "a zero lies on the unit circle" to "z + 1/z is real with magnitude ≤ 2" — the change of
variable that turns the on-circle question (§1777's owed input) into a real-rootedness question,
avoiding the (mathlib-absent) Herglotz/trigonometric-moment theory.  Stated via `normSq` (so the
magnitude bound is `normSq (z + z⁻¹) ≤ 4` and `|z|=1` is `normSq z = 1`), giving a pure real-algebra
proof.  EXIT0 / std-axioms.
-/

open Complex

/-- **Joukowski on-circle criterion.**  For `z ≠ 0`, the quantity `z + z⁻¹` is real with
`normSq (z + z⁻¹) ≤ 4` (i.e. `|z + z⁻¹| ≤ 2`) if and only if `normSq z = 1` (i.e. `|z| = 1`). -/
theorem joukowski_on_circle (z : ℂ) (hz : z ≠ 0) :
    ((z + z⁻¹).im = 0 ∧ normSq (z + z⁻¹) ≤ 4) ↔ normSq z = 1 := by
  have hspos : 0 < normSq z := normSq_pos.mpr hz
  have hs0 : normSq z ≠ 0 := ne_of_gt hspos
  have hns : normSq z = z.re ^ 2 + z.im ^ 2 := by rw [normSq_apply]; ring
  -- clean identities for the imaginary part and the (denominator-cleared) modulus
  have key_im : (z + z⁻¹).im = z.im * (normSq z - 1) / normSq z := by
    rw [add_im, inv_im]; field_simp; ring
  have key_ns : normSq (z + z⁻¹) * normSq z ^ 2
      = z.re ^ 2 * (normSq z + 1) ^ 2 + z.im ^ 2 * (normSq z - 1) ^ 2 := by
    rw [normSq_apply, add_re, add_im, inv_re, inv_im]; field_simp; ring
  constructor
  · rintro ⟨h1, h2⟩
    -- imaginary part zero ⟹ z.im (normSq z - 1) = 0
    rw [key_im, div_eq_zero_iff] at h1
    have h1' : z.im * (normSq z - 1) = 0 := by rcases h1 with h | h; exacts [h, absurd h hs0]
    rcases mul_eq_zero.mp h1' with hb | hs1
    · -- z.im = 0 ⟹ normSq z = z.re², and the modulus bound forces (normSq z - 1)² ≤ 0
      have hsre : normSq z = z.re ^ 2 := by rw [hns, hb]; ring
      have hmul : normSq (z + z⁻¹) * normSq z ^ 2 ≤ 4 * normSq z ^ 2 := by
        have := mul_le_mul_of_nonneg_right h2 (sq_nonneg (normSq z))
        linarith [this]
      rw [key_ns, hb] at hmul
      -- hmul : z.re²(normSq z+1)² + 0 ≤ 4 normSq z² ; with normSq z = z.re²  ⟹  (normSq z-1)² ≤ 0
      have hsq : (normSq z - 1) ^ 2 ≤ 0 := by nlinarith [hmul, hsre, hspos]
      have hz1 : (normSq z - 1) ^ 2 = 0 := le_antisymm hsq (sq_nonneg _)
      have hzz : normSq z - 1 = 0 := by
        by_contra hne; exact (pow_ne_zero 2 hne) hz1
      linarith [hzz]
    · linarith [hs1]
  · intro h
    refine ⟨?_, ?_⟩
    · rw [key_im, h]; ring
    · -- normSq z = 1 ⟹ normSq (z + z⁻¹) = 4 z.re² ≤ 4
      have hb : z.re ^ 2 ≤ 1 := by nlinarith [hns, h, sq_nonneg z.im]
      have hmul : normSq (z + z⁻¹) * normSq z ^ 2
          = z.re ^ 2 * (normSq z + 1) ^ 2 + z.im ^ 2 * (normSq z - 1) ^ 2 := key_ns
      rw [h] at hmul
      -- hmul : normSq(z+z⁻¹) * 1 = z.re² * 4 + 0
      nlinarith [hmul, hb]
