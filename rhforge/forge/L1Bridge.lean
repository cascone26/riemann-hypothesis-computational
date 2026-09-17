import Mathlib.NumberTheory.NumberField.Completion.FinitePlace
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Data.Int.WithZero
open IsDedekindDomain NumberField WithZeroMulInt

noncomputable section
variable {K : Type*} [Field K] [NumberField K]

/-! L1Bridge: Connection between FinitePlace and HeightOneSpectrum valuations

This module bridges finite places to their underlying height-one-spectrum valuations
and ideal norms, enabling the product formula's decomposition across places.

Central mathematical connection (Mathlib's FinitePlace.norm_embedding'):
  ‖embedding v x‖ = toNNReal(absNorm_ne_zero v)(v.valuation K x)

where:
  - v : HeightOneSpectrum (𝓞 K) is a prime ideal
  - v.valuation K x : ℤᵐ⁰ is the adic valuation (WithZero (Multiplicative ℤ))
  - toNNReal (e : ℝ≥0) : ℤᵐ⁰ → ℝ≥0 computes e^(-(valuation))
  - absNorm v.asIdeal : ℕ is the ideal norm

This proves the bridge without sorries by leveraging Mathlib's existing lemmas.
-/

/-- Core bridge: FinitePlace value equals embedding norm into the completion. -/
theorem finitePlace_eq_embedding_norm (w : FinitePlace K) (x : K) :
    w x = ‖FinitePlace.embedding (FinitePlace.maximalIdeal w) x‖ :=
  (FinitePlace.norm_embedding_eq w x).symm

/-- FinitePlace.mk realizes the embedding norm. -/
theorem finitePlace_mk_is_embedding_norm (v : HeightOneSpectrum (𝓞 K)) (x : K) :
    FinitePlace.mk v x = ‖FinitePlace.embedding v x‖ :=
  FinitePlace.mk_apply v x

/-- Absolute value positivity characterizes nonzero elements. -/
theorem finitePlace_norm_pos_iff (w : FinitePlace K) (x : K) :
    0 < w x ↔ x ≠ 0 :=
  FinitePlace.pos_iff

/-- HeightOneSpectrum ideal norm is strictly > 1. -/
theorem ideal_absNorm_gt_one (v : HeightOneSpectrum (𝓞 K)) :
    1 < (Ideal.absNorm v.asIdeal : ℝ) :=
  mod_cast HeightOneSpectrum.one_lt_absNorm v

/-- Ideal norm is nonzero. -/
theorem ideal_absNorm_ne_zero (v : HeightOneSpectrum (𝓞 K)) :
    (Ideal.absNorm v.asIdeal : ℝ) ≠ 0 := by
  linarith [ideal_absNorm_gt_one v]

/-- FinitePlace to HeightOneSpectrum connection. -/
theorem finitePlace_maximalIdeal_is_embedding_norm (w : FinitePlace K) (x : K) :
    w x = ‖FinitePlace.embedding w.maximalIdeal x‖ :=
  finitePlace_eq_embedding_norm w x

/-- The mk/maximalIdeal equivalence. -/
theorem finitePlace_mk_maximalIdeal_symm (w : FinitePlace K) :
    FinitePlace.mk w.maximalIdeal = w :=
  FinitePlace.mk_maximalIdeal w

/-- The adic absolute value equals embedding norm. -/
theorem adicAbv_eq_embedding_norm (v : HeightOneSpectrum (𝓞 K)) (x : K) :
    HeightOneSpectrum.adicAbv K v x = ‖FinitePlace.embedding v x‖ :=
  (FinitePlace.norm_embedding v x).symm

/-- The L1 bridge lemma (from Mathlib FinitePlace.norm_embedding'):
    embedding norm equals toNNReal applied to the valuation.
    toNNReal(e)(n) computes e^n for n : ℤᵐ⁰. -/
theorem norm_embedding_is_toNNReal (v : HeightOneSpectrum (𝓞 K)) (x : K) :
    ‖FinitePlace.embedding v x‖ = toNNReal
      (HeightOneSpectrum.absNorm_ne_zero v) (v.valuation K x) :=
  FinitePlace.norm_embedding' v x

/-- The bridge: FinitePlace value via toNNReal. -/
theorem finitePlace_via_toNNReal (w : FinitePlace K) (x : K) :
    w x = toNNReal
      (HeightOneSpectrum.absNorm_ne_zero w.maximalIdeal)
      (w.maximalIdeal.valuation K x) := by
  calc w x
      = ‖FinitePlace.embedding w.maximalIdeal x‖ := finitePlace_eq_embedding_norm w x
    _ = toNNReal (HeightOneSpectrum.absNorm_ne_zero w.maximalIdeal)
        (w.maximalIdeal.valuation K x) := norm_embedding_is_toNNReal w.maximalIdeal x

/-- The adicAbv is the canonical bridge formula. -/
theorem adicAbv_equals_toNNReal (v : HeightOneSpectrum (𝓞 K)) (x : K) :
    HeightOneSpectrum.adicAbv K v x =
      toNNReal (HeightOneSpectrum.absNorm_ne_zero v) (v.valuation K x) := rfl

/-- For nonzero elements, support is finite. -/
theorem finitePlace_valuation_support (x : 𝓞 K) (hx : x ≠ 0) :
    (fun w : FinitePlace K => w ↑x).HasFiniteMulSupport :=
  FinitePlace.hasFiniteMulSupport_int hx
