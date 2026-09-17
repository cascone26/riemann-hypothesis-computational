import ArakelovVecRHodge
import Mathlib.NumberTheory.NumberField.Completion.FinitePlace

set_option maxHeartbeats 1000000

open IsDedekindDomain NumberField

noncomputable section

variable (K : Type*) [Field K] [NumberField K]

namespace ArakelovVecR

/-!
# §1714 — THE HYPERBOLIC PAIR in the ℝ-model (real-place case)

With a REAL place `w₀` (`mult = 1`) and any finite prime `v` (`L := log N(v) > 1·0`), the Hodge
form has an explicit hyperbolic pair: `F := (0, (1/L)·δ_{w₀})` is isotropic for EVERY value at a
real place (`t² − t²·1 = 0`), `S := (δ_v, b·δ_{w₀})` with `b := (1−L²)/(2L)` is isotropic
(`(L+b)² − (1+b²) = L² + 2Lb − 1 = 0`), `B(F,S) = (1/L)(L+b) − b/L = 1` exactly, and they are
independent (only `S` has finite support).  Together with §1713's negdef-on-`ker deg`, these are
ALL the pointwise hypotheses of §1659 `hodge_signature_of_adjunction` — what remains for the
focused session is the finite-dimensional-subspace restriction and the orthogonal-complement
bookkeeping (and the no-real-place case via the `2a²+2b²+8ab = 0` conic).
-/

theorem exists_hyperbolic_pair (w₀ : InfinitePlace K) (hw : w₀.IsReal)
    (v : HeightOneSpectrum (RingOfIntegers K)) :
    ∃ F S : ArakelovVecR K, hodgeBilin K F F = 0 ∧ hodgeBilin K S S = 0 ∧
      hodgeBilin K F S = 1 ∧ LinearIndependent ℝ ![F, S] := by
  classical
  have hm1 : w₀.mult = 1 := by simp [InfinitePlace.mult, hw]
  have hm1R : ((w₀.mult : ℕ) : ℝ) = 1 := by rw [hm1]; norm_num
  set L : ℝ := Real.log (Ideal.absNorm v.asIdeal : ℝ) with hLdef
  have hL0 : 0 < L := by
    rw [hLdef]
    apply Real.log_pos
    exact_mod_cast HeightOneSpectrum.one_lt_absNorm v
  set b : ℝ := (1 - L ^ 2) / (2 * L) with hbdef
  set F : ArakelovVecR K := (0, Pi.single w₀ (1 / L)) with hFdef
  set S : ArakelovVecR K := (Finsupp.single v 1, Pi.single w₀ b) with hSdef
  -- archimedean sum evaluators
  have hdeg_single : ∀ t : ℝ,
      (∑ w : InfinitePlace K, (w.mult : ℝ) * (Pi.single w₀ t : InfinitePlace K → ℝ) w)
        = t := by
    intro t
    rw [Finset.sum_eq_single_of_mem w₀ (Finset.mem_univ _)]
    · rw [Pi.single_eq_same, hm1R, one_mul]
    · intro w _ hne
      rw [Pi.single_eq_of_ne hne, mul_zero]
  have hpair_single : ∀ s t : ℝ,
      (∑ w : InfinitePlace K, (Pi.single w₀ s : InfinitePlace K → ℝ) w
          * (Pi.single w₀ t : InfinitePlace K → ℝ) w * (w.mult : ℝ)) = s * t := by
    intro s t
    rw [Finset.sum_eq_single_of_mem w₀ (Finset.mem_univ _)]
    · rw [Pi.single_eq_same, Pi.single_eq_same, hm1R, mul_one]
    · intro w _ hne
      rw [Pi.single_eq_of_ne hne, zero_mul, zero_mul]
  -- the six evaluations
  have hdegF : degFun K F = 1 / L := by
    unfold degFun
    rw [hFdef]
    rw [Finsupp.sum_zero_index]
    rw [hdeg_single]
    ring
  have hdegS : degFun K S = L + b := by
    unfold degFun
    rw [hSdef]
    rw [Finsupp.sum_single_index (zero_mul _)]
    rw [hdeg_single]
    rw [← hLdef]
    ring
  have hpFF : pairingFun K F F = (1 / L) * (1 / L) := by
    unfold pairingFun
    rw [hFdef]
    rw [Finsupp.sum_zero_index, hpair_single]
    ring
  have hpSS : pairingFun K S S = 1 + b * b := by
    unfold pairingFun
    rw [hSdef]
    rw [Finsupp.sum_single_index (zero_mul _), Finsupp.single_eq_same, hpair_single]
    ring
  have hpFS : pairingFun K F S = (1 / L) * b := by
    unfold pairingFun
    rw [hFdef, hSdef]
    rw [Finsupp.sum_zero_index, hpair_single]
    ring
  refine ⟨F, S, ?_, ?_, ?_, ?_⟩
  · rw [hodgeBilin_apply, hdegF, hpFF]
    ring
  · rw [hodgeBilin_apply, hdegS, hpSS, hbdef]
    field_simp
    ring
  · rw [hodgeBilin_apply, hdegF, hdegS, hpFS]
    field_simp
    ring
  · rw [LinearIndependent.pair_iff]
    intro s t h
    have h1 := congrArg Prod.fst h
    have h2 := congrArg Prod.snd h
    have hfst : (s • F + t • S).1 = s • F.1 + t • S.1 := rfl
    have hb0 : t = 0 := by
      have hv := DFunLike.congr_fun (by rw [← hfst]; exact h1 :
        s • F.1 + t • S.1 = (0 : HeightOneSpectrum (RingOfIntegers K) →₀ ℝ)) v
      rw [hFdef, hSdef] at hv
      simp only [Finsupp.coe_add, Finsupp.coe_smul, Pi.add_apply, Pi.smul_apply,
        Finsupp.coe_zero, Pi.zero_apply, smul_eq_mul] at hv
      rw [Finsupp.single_eq_same] at hv
      simpa using hv
    have ha0 : s = 0 := by
      have hsnd : (s • F + t • S).2 = s • F.2 + t • S.2 := rfl
      have hw₀eval := congrFun (by rw [← hsnd]; exact h2 :
        s • F.2 + t • S.2 = (0 : InfinitePlace K → ℝ)) w₀
      rw [hFdef, hSdef] at hw₀eval
      simp only [Pi.add_apply, Pi.smul_apply, Pi.zero_apply, smul_eq_mul,
        Pi.single_eq_same] at hw₀eval
      rw [hb0] at hw₀eval
      have hLne : (1 : ℝ) / L ≠ 0 := by positivity
      have : s * (1 / L) = 0 := by linarith [hw₀eval]
      rcases mul_eq_zero.mp this with h' | h'
      · exact h'
      · exact absurd h' hLne
    exact ⟨ha0, hb0⟩

#print axioms exists_hyperbolic_pair
