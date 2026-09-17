import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Tactic

/-!
# §1776 — Palindromic (even) eigenvector ⟹ self-inversive polynomial (reciprocal-paired zeros)

The next ingredient of the Connes–van Suijlekom / Carathéodory-Fejér engine.  By §1775+§1771, the minimal
eigenvector `v` of the finite Weil-Toeplitz form is EVEN, i.e. PALINDROMIC: `v (Fin.rev i) = v i`.  Its
associated polynomial `P(z) = ∑ v_i z^i` is then SELF-INVERSIVE:

    ∑ v_i ρ^{(n-1)-i}  =  ∑ v_i ρ^i      (i.e. `P(ρ) = ρ^{n-1} P(1/ρ)`).

Hence the zeros come in reciprocal pairs `{ρ, 1/ρ}` — symmetric about the unit circle `|z|=1`.  (The PSD /
singular Toeplitz structure — Carathéodory — upgrades "symmetric about" to "on" the circle; that is the
remaining half.)  This brick establishes the reciprocal-pairing.  EXIT0 / std-axioms.
-/

open Finset

variable {𝕜 : Type*} [Field 𝕜] {n : ℕ}

/-- **Self-inversive identity.**  For a palindromic coefficient vector, the reversed polynomial equals the
polynomial: `∑ v_i ρ^{(n-1)-i} = ∑ v_i ρ^i`. -/
theorem self_inversive_identity (v : Fin n → 𝕜) (hpal : ∀ i, v (Fin.rev i) = v i) (ρ : 𝕜) :
    (∑ i : Fin n, v i * ρ ^ ((n - 1) - (i : ℕ))) = ∑ i : Fin n, v i * ρ ^ (i : ℕ) := by
  refine Fintype.sum_equiv (Equiv.mk Fin.rev Fin.rev Fin.rev_rev Fin.rev_rev) _ _ (fun i => ?_)
  simp only [Equiv.coe_fn_mk]
  have h1 : v (Fin.rev i) = v i := hpal i
  have h2 : ((Fin.rev i : Fin n) : ℕ) = (n - 1) - (i : ℕ) := by
    simp only [Fin.val_rev]; have := i.isLt; omega
  rw [h1, h2]

/-- **Reciprocal-paired zeros.**  If `ρ ≠ 0` is a zero of `P(z) = ∑ v_i z^i` (palindromic `v`), then so is
`1/ρ`: the zeros are symmetric about the unit circle. -/
theorem self_inversive_reciprocal_zero (v : Fin n → 𝕜) (hpal : ∀ i, v (Fin.rev i) = v i)
    (ρ : 𝕜) (hρ : ρ ≠ 0) (hz : (∑ i : Fin n, v i * ρ ^ (i : ℕ)) = 0) :
    (∑ i : Fin n, v i * (ρ⁻¹) ^ (i : ℕ)) = 0 := by
  -- `∑ v_i ρ^{n-1-i} = 0` from the identity + hypothesis
  have hrev : (∑ i : Fin n, v i * ρ ^ ((n - 1) - (i : ℕ))) = 0 :=
    (self_inversive_identity v hpal ρ).trans hz
  -- multiply the target by `ρ^{n-1}` and match termwise
  have key : ρ ^ (n - 1) * (∑ i : Fin n, v i * (ρ⁻¹) ^ (i : ℕ))
      = ∑ i : Fin n, v i * ρ ^ ((n - 1) - (i : ℕ)) := by
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro i _
    have hsplit : ρ ^ (n - 1) = ρ ^ ((n - 1) - (i : ℕ)) * ρ ^ (i : ℕ) := by
      rw [← pow_add]; congr 1; have := i.isLt; omega
    have hterm : ρ ^ (n - 1) * (ρ⁻¹) ^ (i : ℕ) = ρ ^ ((n - 1) - (i : ℕ)) := by
      rw [inv_pow, hsplit, mul_assoc, mul_inv_cancel₀ (pow_ne_zero _ hρ), mul_one]
    calc ρ ^ (n - 1) * (v i * (ρ⁻¹) ^ (i : ℕ))
        = v i * (ρ ^ (n - 1) * (ρ⁻¹) ^ (i : ℕ)) := by ring
      _ = v i * ρ ^ ((n - 1) - (i : ℕ)) := by rw [hterm]
  rw [hrev] at key
  exact (mul_eq_zero.mp key).resolve_left (pow_ne_zero _ hρ)
