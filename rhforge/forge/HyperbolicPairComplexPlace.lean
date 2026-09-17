import HyperbolicPairRealPlace

set_option maxHeartbeats 1000000

open IsDedekindDomain NumberField

noncomputable section

variable (K : Type*) [Field K] [NumberField K]

namespace ArakelovVecR

/-!
# §1715 — THE HYPERBOLIC PAIR, complex-place case (+ the every-field corollary)

At a complex place (`mult = 2`) pure-archimedean isotropics don't exist (`(2t)² − 2t² = 2t² ≠ 0`),
so mix with a prime: on `D = (δ_v, t·δ_{w₀})` the Hodge form is `Q(t) = 2t² + 4Lt + (L²−1)`, an
indefinite conic with roots `t₁,₂ = (−2L ± √(2L²+2))/2`.  Both roots give isotropic vectors; Vieta
gives `B(F₀,S) = −(2L²+2) ≠ 0`, and rescaling `F := (1/c)·F₀` normalizes the cross term to 1.
With §1714 this covers EVERY infinite place, hence every number field.
-/

theorem exists_hyperbolic_pair_of_complex (w₀ : InfinitePlace K) (hw : ¬w₀.IsReal)
    (v : HeightOneSpectrum (RingOfIntegers K)) :
    ∃ F S : ArakelovVecR K, hodgeBilin K F F = 0 ∧ hodgeBilin K S S = 0 ∧
      hodgeBilin K F S = 1 ∧ LinearIndependent ℝ ![F, S] := by
  classical
  have hm2 : w₀.mult = 2 := by simp [InfinitePlace.mult, hw]
  have hm2R : ((w₀.mult : ℕ) : ℝ) = 2 := by rw [hm2]; norm_num
  set L : ℝ := Real.log (Ideal.absNorm v.asIdeal : ℝ) with hLdef
  have hL0 : 0 < L := by
    rw [hLdef]
    apply Real.log_pos
    exact_mod_cast HeightOneSpectrum.one_lt_absNorm v
  set s : ℝ := Real.sqrt (2 * L ^ 2 + 2) with hsdef
  have hs0 : 0 < s := by
    rw [hsdef]
    apply Real.sqrt_pos.mpr
    positivity
  have hs2 : s ^ 2 = 2 * L ^ 2 + 2 := by
    rw [hsdef]
    exact Real.sq_sqrt (by positivity)
  set t₁ : ℝ := (-2 * L + s) / 2 with ht₁def
  set t₂ : ℝ := (-2 * L - s) / 2 with ht₂def
  have hroot₁ : 2 * t₁ ^ 2 + 4 * L * t₁ + (L ^ 2 - 1) = 0 := by
    rw [ht₁def]
    have : ((-2 * L + s) / 2) ^ 2 = (4 * L ^ 2 - 4 * L * s + s ^ 2) / 4 := by ring
    rw [this, hs2]
    ring
  have hroot₂ : 2 * t₂ ^ 2 + 4 * L * t₂ + (L ^ 2 - 1) = 0 := by
    rw [ht₂def]
    have : ((-2 * L - s) / 2) ^ 2 = (4 * L ^ 2 + 4 * L * s + s ^ 2) / 4 := by ring
    rw [this, hs2]
    ring
  have htne : t₁ ≠ t₂ := by
    rw [ht₁def, ht₂def]
    intro h
    have : s = 0 := by linarith
    exact hs0.ne' this
  set F₀ : ArakelovVecR K := (Finsupp.single v 1, Pi.single w₀ t₁) with hF₀def
  set S : ArakelovVecR K := (Finsupp.single v 1, Pi.single w₀ t₂) with hSdef
  -- archimedean sum evaluators (mult = 2)
  have hdeg_single : ∀ t : ℝ,
      (∑ w : InfinitePlace K, (w.mult : ℝ) * (Pi.single w₀ t : InfinitePlace K → ℝ) w)
        = 2 * t := by
    intro t
    rw [Finset.sum_eq_single_of_mem w₀ (Finset.mem_univ _)]
    · rw [Pi.single_eq_same, hm2R]
    · intro w _ hne
      rw [Pi.single_eq_of_ne hne, mul_zero]
  have hpair_single : ∀ a b : ℝ,
      (∑ w : InfinitePlace K, (Pi.single w₀ a : InfinitePlace K → ℝ) w
          * (Pi.single w₀ b : InfinitePlace K → ℝ) w * (w.mult : ℝ)) = 2 * (a * b) := by
    intro a b
    rw [Finset.sum_eq_single_of_mem w₀ (Finset.mem_univ _)]
    · rw [Pi.single_eq_same, Pi.single_eq_same, hm2R]
      ring
    · intro w _ hne
      rw [Pi.single_eq_of_ne hne, zero_mul, zero_mul]
  -- evaluations
  have hdegF₀ : degFun K F₀ = L + 2 * t₁ := by
    unfold degFun
    rw [hF₀def]
    rw [Finsupp.sum_single_index (zero_mul _), hdeg_single, ← hLdef]
    ring
  have hdegS : degFun K S = L + 2 * t₂ := by
    unfold degFun
    rw [hSdef]
    rw [Finsupp.sum_single_index (zero_mul _), hdeg_single, ← hLdef]
    ring
  have hpF₀F₀ : pairingFun K F₀ F₀ = 1 + 2 * (t₁ * t₁) := by
    unfold pairingFun
    rw [hF₀def]
    rw [Finsupp.sum_single_index (zero_mul _), Finsupp.single_eq_same, hpair_single]
    ring
  have hpSS : pairingFun K S S = 1 + 2 * (t₂ * t₂) := by
    unfold pairingFun
    rw [hSdef]
    rw [Finsupp.sum_single_index (zero_mul _), Finsupp.single_eq_same, hpair_single]
    ring
  have hpF₀S : pairingFun K F₀ S = 1 + 2 * (t₁ * t₂) := by
    unfold pairingFun
    rw [hF₀def, hSdef]
    rw [Finsupp.sum_single_index (zero_mul _), Finsupp.single_eq_same, hpair_single]
    ring
  -- isotropy and the cross term
  have hBF₀F₀ : hodgeBilin K F₀ F₀ = 0 := by
    rw [hodgeBilin_apply, hdegF₀, hpF₀F₀]
    linear_combination hroot₁
  have hBSS : hodgeBilin K S S = 0 := by
    rw [hodgeBilin_apply, hdegS, hpSS]
    linear_combination hroot₂
  have hsum : t₁ + t₂ = -2 * L := by
    rw [ht₁def, ht₂def]
    ring
  have hprod : t₁ * t₂ = (L ^ 2 - 1) / 2 := by
    rw [ht₁def, ht₂def]
    have : (-2 * L + s) / 2 * ((-2 * L - s) / 2) = (4 * L ^ 2 - s ^ 2) / 4 := by ring
    rw [this, hs2]
    ring
  set c : ℝ := -(2 * L ^ 2 + 2) with hcdef
  have hc0 : c ≠ 0 := by
    rw [hcdef]
    have : (0:ℝ) < 2 * L ^ 2 + 2 := by positivity
    linarith
  have hBF₀S : hodgeBilin K F₀ S = c := by
    rw [hodgeBilin_apply, hdegF₀, hdegS, hpF₀S, hcdef]
    have hexp : (L + 2 * t₁) * (L + 2 * t₂)
        = L ^ 2 + 2 * L * (t₁ + t₂) + 4 * (t₁ * t₂) := by ring
    rw [hexp, hsum, hprod]
    ring
  -- smul lemmas for the bilinear form
  have hsmul_left : ∀ (a : ℝ) (D E : ArakelovVecR K),
      hodgeBilin K (a • D) E = a * hodgeBilin K D E := by
    intro a D E
    rw [hodgeBilin_apply, hodgeBilin_apply, degFun_smul, pairingFun_smul_left]
    ring
  have hsmul_right : ∀ (a : ℝ) (D E : ArakelovVecR K),
      hodgeBilin K D (a • E) = a * hodgeBilin K D E := by
    intro a D E
    rw [hodgeBilin_apply, hodgeBilin_apply, degFun_smul, pairingFun_smul_right]
    ring
  refine ⟨(1 / c) • F₀, S, ?_, hBSS, ?_, ?_⟩
  · rw [hsmul_left, hsmul_right, hBF₀F₀]
    ring
  · rw [hsmul_left, hBF₀S]
    field_simp
  · rw [LinearIndependent.pair_iff]
    intro a u h
    have h1 := congrArg Prod.fst h
    have hfst : (a • ((1 / c) • F₀) + u • S).1
        = (a * (1 / c)) • F₀.1 + u • S.1 := by
      rw [smul_smul]
      rfl
    have hveval : (a * (1 / c)) * 1 + u * 1 = 0 := by
      have hv := DFunLike.congr_fun (by rw [← hfst]; exact h1 :
        (a * (1 / c)) • F₀.1 + u • S.1
          = (0 : HeightOneSpectrum (RingOfIntegers K) →₀ ℝ)) v
      rw [hF₀def, hSdef] at hv
      simpa [Finsupp.single_eq_same] using hv
    have h2 := congrArg Prod.snd h
    have hsnd : (a • ((1 / c) • F₀) + u • S).2
        = (a * (1 / c)) • F₀.2 + u • S.2 := by
      rw [smul_smul]
      rfl
    have hweval : (a * (1 / c)) * t₁ + u * t₂ = 0 := by
      have hw' := congrFun (by rw [← hsnd]; exact h2 :
        (a * (1 / c)) • F₀.2 + u • S.2 = (0 : InfinitePlace K → ℝ)) w₀
      rw [hF₀def, hSdef] at hw'
      simpa [Pi.single_eq_same] using hw'
    have hu0 : u = 0 := by
      have hac : a * (1 / c) = -u := by linarith [hveval]
      rw [hac] at hweval
      have : u * (t₂ - t₁) = 0 := by linarith [hweval]
      rcases mul_eq_zero.mp this with h' | h'
      · exact h'
      · exact absurd (by linarith [h'] : t₁ = t₂) htne
    have ha0 : a = 0 := by
      rw [hu0] at hveval
      have : a * (1 / c) = 0 := by linarith [hveval]
      rcases mul_eq_zero.mp this with h' | h'
      · exact h'
      · exact absurd h' (by simpa using hc0)
    exact ⟨ha0, hu0⟩

/-- **Every number field admits a hyperbolic pair** for the Hodge form on the ℝ-model. -/
theorem exists_hyperbolic_pair_any (v : HeightOneSpectrum (RingOfIntegers K)) :
    ∃ F S : ArakelovVecR K, hodgeBilin K F F = 0 ∧ hodgeBilin K S S = 0 ∧
      hodgeBilin K F S = 1 ∧ LinearIndependent ℝ ![F, S] := by
  obtain ⟨w₀⟩ : Nonempty (InfinitePlace K) := inferInstance
  by_cases hw : w₀.IsReal
  · exact exists_hyperbolic_pair K w₀ hw v
  · exact exists_hyperbolic_pair_of_complex K w₀ hw v

#print axioms exists_hyperbolic_pair_of_complex
#print axioms exists_hyperbolic_pair_any
