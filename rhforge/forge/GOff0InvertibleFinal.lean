import CauchyVectorsLinearIndependent
import QOffKMonicPoly
import STOffRealSymm
import Eq7p17DistinctnessBasis

open scoped MatrixOrder Matrix Polynomial
open Polynomial

-- STEP A: S_T_off_sum's entries are purely real (im = 0), reusing STOffRealValued.lean's
-- S_T_off_sum_eq_ofReal_re.
lemma S_T_off_sum_entries_real (N : ℕ) (L : ℝ) (K : ℕ) (z : Fin K → ℂ) (i j : Fin (2 * N + 1)) :
    (S_T_off_sum N L K z i j).im = 0 := by
  rw [S_T_off_sum_eq_ofReal_re N L K z i j]
  simp

-- STEP B: real kernel condition complexifies.
lemma complexify_kernel (N : ℕ) (L : ℝ) (K : ℕ) (z : Fin K → ℂ)
    (u : Fin (2 * N + 1) → ℝ)
    (hu : (S_T_off_real N L K z).mulVec u = 0) :
    (S_T_off_sum N L K z).mulVec (fun n => (u n : ℂ)) = 0 := by
  funext i
  have hui : (S_T_off_real N L K z).mulVec u i = 0 := congrFun hu i
  show ∑ j, S_T_off_sum N L K z i j * (u j : ℂ) = 0
  have heq : ∀ j, S_T_off_sum N L K z i j * (u j : ℂ)
      = ((S_T_off_real N L K z i j * u j : ℝ) : ℂ) := by
    intro j
    rw [Complex.ofReal_mul]
    congr 1
    exact S_T_off_sum_eq_ofReal_re N L K z i j
  rw [Finset.sum_congr rfl (fun j _ => heq j), ← Complex.ofReal_sum]
  have : ∑ j, S_T_off_real N L K z i j * u j = (S_T_off_real N L K z).mulVec u i := rfl
  rw [this, hui]
  simp

-- STEP C: expand the rank-1 structure of S_T_off_sum.mulVec.
lemma vecMulVec_mulVec_comm {n : Type*} [Fintype n] (a w : n → ℂ) :
    (Matrix.vecMulVec a a).mulVec w = (a ⬝ᵥ w) • a := by
  rw [Matrix.vecMulVec_mulVec]
  simp

lemma H_L_G_C_mulVec (N : ℕ) (L : ℝ) (w : ℂ) (x : Fin (2 * N + 1) → ℂ) :
    (H_L_G_C N L w).mulVec x
      = (vPlusLC N L w ⬝ᵥ x) • vPlusLC N L w + (vMinusLC N L w ⬝ᵥ x) • vMinusLC N L w := by
  unfold H_L_G_C
  rw [Matrix.add_mulVec, vecMulVec_mulVec_comm (vPlusLC N L w) x,
    vecMulVec_mulVec_comm (vMinusLC N L w) x]

lemma S_T_off_sum_expansion (N : ℕ) (L : ℝ) (K : ℕ) (z : Fin K → ℂ)
    (u_c : Fin (2 * N + 1) → ℂ) :
    (S_T_off_sum N L K z).mulVec u_c =
    ∑ k : Fin K,
      ((vPlusLC N L (z k) ⬝ᵥ u_c) • (vPlusLC N L (z k)) +
       (vMinusLC N L (z k) ⬝ᵥ u_c) • (vMinusLC N L (z k)) +
       (vPlusLC N L (starRingEnd ℂ (z k)) ⬝ᵥ u_c) • (vPlusLC N L (starRingEnd ℂ (z k))) +
       (vMinusLC N L (starRingEnd ℂ (z k)) ⬝ᵥ u_c) • (vMinusLC N L (starRingEnd ℂ (z k)))) := by
  unfold S_T_off_sum
  rw [Matrix.sum_mulVec]
  apply Finset.sum_congr rfl
  intro k _
  rw [Matrix.add_mulVec, H_L_G_C_mulVec, H_L_G_C_mulVec]
  ring

#print axioms S_T_off_sum_entries_real
#print axioms complexify_kernel
#print axioms S_T_off_sum_expansion

-- STEP D: quartet indexing and injectivity.
def quartet_a (K : ℕ) (z : Fin K → ℂ) (p : Fin 4 × Fin K) : ℂ :=
  match p.1 with
  | 0 => z p.2
  | 1 => -(z p.2)
  | 2 => star (z p.2)
  | 3 => -star (z p.2)

lemma quartet_a_mem_f (K : ℕ) (z : Fin K → ℂ) (p : Fin 4 × Fin K) :
    quartet_a K z p ∈ [z p.2, -(z p.2), star (z p.2), -star (z p.2)] := by
  obtain ⟨p1, p2⟩ := p
  fin_cases p1 <;> simp [quartet_a]

lemma pairwise_disjoint_of_mem {α β : Type*} {f : α → List β} :
    ∀ {l₁ : List α}, List.Pairwise (Function.onFun List.Disjoint f) l₁ →
      ∀ a ∈ l₁, ∀ b ∈ l₁, a ≠ b → List.Disjoint (f a) (f b)
  | [], _, a, ha, _, _, _ => by simp at ha
  | x :: xs, h, a, ha, b, hb, hab => by
      cases h with
      | cons hr ht =>
        simp only [List.mem_cons] at ha hb
        rcases ha with ha | ha <;> rcases hb with hb | hb
        · exact absurd (ha.trans hb.symm) hab
        · subst ha; exact hr b hb
        · subst hb; exact (hr a ha).symm
        · exact pairwise_disjoint_of_mem ht a ha b hb hab

lemma quartet_a_injective (K : ℕ) (z : Fin K → ℂ)
    (hdist : quartets_distinct_nodup K z) :
    Function.Injective (quartet_a K z) := by
  intro p q hpq
  have hnd := (List.nodup_flatMap).mp hdist
  obtain ⟨p1, p2⟩ := p
  obtain ⟨q1, q2⟩ := q
  by_cases hk : p2 = q2
  · subst hk
    have hnodup4 : List.Nodup [z p2, -(z p2), star (z p2), -star (z p2)] :=
      hnd.1 p2 (by simp)
    simp only [List.nodup_cons, List.mem_cons, List.mem_singleton, List.not_mem_nil,
      or_false, List.nodup_singleton, and_true, not_or] at hnodup4
    obtain ⟨⟨h01, h02, h03⟩, ⟨h12, h13⟩, h23, -, -⟩ := hnodup4
    fin_cases p1 <;> fin_cases q1 <;> simp only [quartet_a] at hpq <;>
      first
        | rfl
        | exact absurd hpq h01 | exact absurd hpq.symm h01
        | exact absurd hpq h02 | exact absurd hpq.symm h02
        | exact absurd hpq h03 | exact absurd hpq.symm h03
        | exact absurd hpq h12 | exact absurd hpq.symm h12
        | exact absurd hpq h13 | exact absurd hpq.symm h13
        | exact absurd hpq h23 | exact absurd hpq.symm h23
  · exfalso
    have hmemuniv1 : p2 ∈ (Finset.univ : Finset (Fin K)).toList := by simp
    have hmemuniv2 : q2 ∈ (Finset.univ : Finset (Fin K)).toList := by simp
    have hdisj : List.Disjoint ([z p2, -(z p2), star (z p2), -star (z p2)])
        ([z q2, -(z q2), star (z q2), -star (z q2)]) :=
      pairwise_disjoint_of_mem hnd.2 p2 hmemuniv1 q2 hmemuniv2 hk
    have hmemp := quartet_a_mem_f K z (p1, p2)
    have hmemq := quartet_a_mem_f K z (q1, q2)
    exact hdisj hmemp (hpq ▸ hmemq)

#print axioms quartet_a_injective

lemma quartet_a_im_ne_zero (K : ℕ) (z : Fin K → ℂ) (hz : ∀ k, (z k).im ≠ 0)
    (p : Fin 4 × Fin K) : (quartet_a K z p).im ≠ 0 := by
  obtain ⟨p1, p2⟩ := p
  have him := hz p2
  fin_cases p1 <;> simp only [quartet_a] <;> simp <;> [skip; skip; skip; skip] <;>
    first
      | exact him
      | (intro h; apply him; linarith)
      | (simpa using him)

#print axioms quartet_a_im_ne_zero

-- STEP E: package the quartet indexing for Fin (4*K), and apply cauchy_vectors_linear_independent.
noncomputable def quartet_a_fin (K : ℕ) (z : Fin K → ℂ) : Fin (4 * K) → ℂ :=
  quartet_a K z ∘ (finProdFinEquiv (m := 4) (n := K)).symm

lemma quartet_a_fin_injective (K : ℕ) (z : Fin K → ℂ) (hdist : quartets_distinct_nodup K z) :
    Function.Injective (quartet_a_fin K z) :=
  (quartet_a_injective K z hdist).comp (Equiv.injective _)

lemma quartet_a_fin_im_ne_zero (K : ℕ) (z : Fin K → ℂ) (hz : ∀ k, (z k).im ≠ 0) (j : Fin (4 * K)) :
    (quartet_a_fin K z j).im ≠ 0 :=
  quartet_a_im_ne_zero K z hz _


#print axioms quartet_a_fin_injective
#print axioms quartet_a_fin_im_ne_zero

-- STEP F: the nonzero quartet constants.
lemma sin_ne_zero_of_im_ne_zero (w : ℂ) (h : w.im ≠ 0) : Complex.sin w ≠ 0 := by
  intro hs
  rw [Complex.sin_eq_zero_iff] at hs
  obtain ⟨n, hn⟩ := hs
  apply h
  rw [hn]
  simp

noncomputable def quartet_const (K : ℕ) (L : ℝ) (z : Fin K → ℂ) (p : Fin 4 × Fin K) : ℂ :=
  match p.1 with
  | 0 => (Real.sqrt 2 : ℂ) * Complex.sin (L * z p.2 / 2) / L
  | 1 => -((Real.sqrt 2 : ℂ) * Complex.sin (L * z p.2 / 2) / L)
  | 2 => (Real.sqrt 2 : ℂ) * Complex.sin (L * (starRingEnd ℂ (z p.2)) / 2) / L
  | 3 => -((Real.sqrt 2 : ℂ) * Complex.sin (L * (starRingEnd ℂ (z p.2)) / 2) / L)

lemma half_mul_ne_zero_im (L : ℝ) (hL : L ≠ 0) (w : ℂ) (hw : w.im ≠ 0) :
    ((L : ℂ) * w / 2).im ≠ 0 := by
  have heq : ((L : ℂ) * w / 2).im = L * w.im / 2 := by
    simp [Complex.div_im, Complex.mul_im]
  rw [heq]
  exact div_ne_zero (mul_ne_zero hL hw) (by norm_num)

lemma quartet_const_ne_zero (K : ℕ) (L : ℝ) (hL : L ≠ 0) (z : Fin K → ℂ)
    (hz : ∀ k, (z k).im ≠ 0) (p : Fin 4 × Fin K) : quartet_const K L z p ≠ 0 := by
  obtain ⟨p1, p2⟩ := p
  have himz : (z p2).im ≠ 0 := hz p2
  have himstar : (starRingEnd ℂ (z p2)).im ≠ 0 := by simpa using himz
  have hsin1 : Complex.sin ((L:ℂ) * z p2 / 2) ≠ 0 :=
    sin_ne_zero_of_im_ne_zero _ (half_mul_ne_zero_im L hL (z p2) himz)
  have hsin2 : Complex.sin ((L:ℂ) * (starRingEnd ℂ (z p2)) / 2) ≠ 0 :=
    sin_ne_zero_of_im_ne_zero _ (half_mul_ne_zero_im L hL _ himstar)
  have hL' : (L : ℂ) ≠ 0 := Complex.ofReal_ne_zero.mpr hL
  have hsqrt2 : (Real.sqrt 2 : ℂ) ≠ 0 := by
    simp [Real.sqrt_eq_zero']
  fin_cases p1 <;> simp only [quartet_const] <;>
    simp [div_eq_zero_iff, mul_eq_zero, hsin1, hsin2, hsqrt2, hL']


#print axioms sin_ne_zero_of_im_ne_zero
#print axioms half_mul_ne_zero_im
#print axioms quartet_const_ne_zero

-- STEP G: entrywise Cauchy-vector form.
noncomputable def freq_node_N (N : ℕ) (L : ℝ) (n : Fin (2 * N + 1)) : ℝ :=
  2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L

lemma vPlusLC_eq_cauchy (N : ℕ) (L : ℝ) (w : ℂ) (n : Fin (2 * N + 1)) :
    vPlusLC N L w n = ((Real.sqrt 2 : ℂ) * Complex.sin (L * w / 2) / L) / (w + (freq_node_N N L n : ℂ)) := by
  unfold vPlusLC freq_node_N
  rw [div_div]
  push_cast
  ring_nf

lemma vMinusLC_eq_cauchy (N : ℕ) (L : ℝ) (w : ℂ) (n : Fin (2 * N + 1)) :
    vMinusLC N L w n = (-((Real.sqrt 2 : ℂ) * Complex.sin (L * w / 2) / L))
      / (-w + (freq_node_N N L n : ℂ)) := by
  unfold vMinusLC freq_node_N
  have hdenom : (-w + ((2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℝ) : ℂ))
      = -(w - ((2 * Real.pi * (((n : ℕ) : ℝ) - (N : ℝ)) / L : ℝ) : ℂ)) := by ring
  rw [hdenom, neg_div_neg_eq, div_div]
  push_cast
  ring_nf

lemma quartet_vec_eq_cauchy (K : ℕ) (L : ℝ) (z : Fin K → ℂ) (p : Fin 4 × Fin K)
    (n : Fin (2 * (2 * K) + 1)) :
    (quartet_const K L z p) / (quartet_a K z p + (freq_node_N (2 * K) L n : ℂ)) =
      (match p.1 with
        | 0 => vPlusLC (2 * K) L (z p.2) n
        | 1 => vMinusLC (2 * K) L (z p.2) n
        | 2 => vPlusLC (2 * K) L (starRingEnd ℂ (z p.2)) n
        | 3 => vMinusLC (2 * K) L (starRingEnd ℂ (z p.2)) n) := by
  obtain ⟨p1, p2⟩ := p
  fin_cases p1 <;> simp only [quartet_const, quartet_a]
  · exact (vPlusLC_eq_cauchy (2 * K) L (z p2) n).symm
  · rw [show (-(z p2) : ℂ) + (freq_node_N (2 * K) L n : ℂ) = -(z p2) + ↑(freq_node_N (2 * K) L n)
      from rfl]
    exact (vMinusLC_eq_cauchy (2 * K) L (z p2) n).symm
  · exact (vPlusLC_eq_cauchy (2 * K) L (starRingEnd ℂ (z p2)) n).symm
  · exact (vMinusLC_eq_cauchy (2 * K) L (starRingEnd ℂ (z p2)) n).symm


#print axioms vPlusLC_eq_cauchy
#print axioms vMinusLC_eq_cauchy
#print axioms quartet_vec_eq_cauchy

-- STEP H: the main kernel-vanishing theorem.
lemma index_eq (K : ℕ) : 2 * (2 * K) + 1 = 4 * K + 1 := by ring

noncomputable def freq_node_fin (K : ℕ) (L : ℝ) : Fin (4 * K + 1) → ℝ :=
  freq_node_N (2 * K) L ∘ (finCongr (index_eq K)).symm

lemma freq_node_fin_injective (K : ℕ) (L : ℝ) (hL : L ≠ 0) :
    Function.Injective (freq_node_fin K L) := by
  have hinj : Function.Injective (freq_node_N (2 * K) L) := by
    intro a b hab
    unfold freq_node_N at hab
    have hpi : (2 : ℝ) * Real.pi ≠ 0 := by positivity
    have h1 : 2 * Real.pi * (((a:ℕ):ℝ) - ((2*K:ℕ):ℝ)) = 2 * Real.pi * (((b:ℕ):ℝ) - ((2*K:ℕ):ℝ)) := by
      have hmul := congrArg (· * L) hab
      simp only [div_mul_cancel₀ _ hL] at hmul
      exact hmul
    have h2 : ((a:ℕ):ℝ) = ((b:ℕ):ℝ) := by
      have := mul_left_cancel₀ hpi h1
      linarith
    exact Fin.ext (by exact_mod_cast h2)
  exact hinj.comp (Equiv.injective _)

noncomputable def quartet_D (K : ℕ) (L : ℝ) (z : Fin K → ℂ) (u_c : Fin (2 * (2 * K) + 1) → ℂ)
    (p : Fin 4 × Fin K) : ℂ :=
  match p.1 with
  | 0 => vPlusLC (2 * K) L (z p.2) ⬝ᵥ u_c
  | 1 => vMinusLC (2 * K) L (z p.2) ⬝ᵥ u_c
  | 2 => vPlusLC (2 * K) L (starRingEnd ℂ (z p.2)) ⬝ᵥ u_c
  | 3 => vMinusLC (2 * K) L (starRingEnd ℂ (z p.2)) ⬝ᵥ u_c

lemma quartet_term_eq (K : ℕ) (L : ℝ) (z : Fin K → ℂ) (u_c : Fin (2 * (2 * K) + 1) → ℂ)
    (n : Fin (2 * (2 * K) + 1)) (k : Fin K) :
    (vPlusLC (2 * K) L (z k) ⬝ᵥ u_c) * vPlusLC (2 * K) L (z k) n
      + (vMinusLC (2 * K) L (z k) ⬝ᵥ u_c) * vMinusLC (2 * K) L (z k) n
      + (vPlusLC (2 * K) L (starRingEnd ℂ (z k)) ⬝ᵥ u_c) * vPlusLC (2 * K) L (starRingEnd ℂ (z k)) n
      + (vMinusLC (2 * K) L (starRingEnd ℂ (z k)) ⬝ᵥ u_c) * vMinusLC (2 * K) L (starRingEnd ℂ (z k)) n
    = ∑ r : Fin 4, quartet_D K L z u_c (r, k) * quartet_const K L z (r, k)
        / (quartet_a K z (r, k) + (freq_node_N (2 * K) L n : ℂ)) := by
  rw [Fin.sum_univ_four]
  simp only [quartet_D, mul_div_assoc, quartet_vec_eq_cauchy]

theorem quartet_coeffs_vanish (K : ℕ) (L : ℝ) (hL : 0 < L) (z : Fin K → ℂ)
    (hz : ∀ k, (z k).im ≠ 0) (hdist : quartets_distinct_nodup K z) (hK : 1 ≤ K)
    (u_c : Fin (2 * (2 * K) + 1) → ℂ)
    (hzero : (S_T_off_sum (2 * K) L K z).mulVec u_c = 0) :
    ∀ p : Fin 4 × Fin K, quartet_D K L z u_c p = 0 := by
  have hLne : L ≠ 0 := ne_of_gt hL
  have hexp := S_T_off_sum_expansion (2 * K) L K z u_c
  rw [hexp] at hzero
  have hpoint : ∀ n : Fin (2 * (2 * K) + 1),
      ∑ p : Fin 4 × Fin K, quartet_D K L z u_c p * quartet_const K L z p
        / (quartet_a K z p + (freq_node_N (2 * K) L n : ℂ)) = 0 := by
    intro n
    have hn := congrFun hzero n
    simp only [Finset.sum_apply, Pi.add_apply, Pi.smul_apply, smul_eq_mul,
      Pi.zero_apply] at hn
    have heq2 : ∑ p : Fin 4 × Fin K, quartet_D K L z u_c p * quartet_const K L z p
          / (quartet_a K z p + (freq_node_N (2 * K) L n : ℂ))
        = ∑ k : Fin K, ((vPlusLC (2 * K) L (z k) ⬝ᵥ u_c) * vPlusLC (2 * K) L (z k) n
            + (vMinusLC (2 * K) L (z k) ⬝ᵥ u_c) * vMinusLC (2 * K) L (z k) n
            + (vPlusLC (2 * K) L (starRingEnd ℂ (z k)) ⬝ᵥ u_c) * vPlusLC (2 * K) L (starRingEnd ℂ (z k)) n
            + (vMinusLC (2 * K) L (starRingEnd ℂ (z k)) ⬝ᵥ u_c) * vMinusLC (2 * K) L (starRingEnd ℂ (z k)) n) := by
      rw [Fintype.sum_prod_type_right (f := fun p => quartet_D K L z u_c p * quartet_const K L z p
        / (quartet_a K z p + (freq_node_N (2 * K) L n : ℂ)))]
      exact Finset.sum_congr rfl (fun k _ => (quartet_term_eq K L z u_c n k).symm)
    rw [heq2]
    convert hn using 2
  have hpoint_fin : ∀ n' : Fin (4 * K + 1),
      ∑ j : Fin (4 * K),
        (quartet_D K L z u_c ((finProdFinEquiv (m := 4) (n := K)).symm j)
          * quartet_const K L z ((finProdFinEquiv (m := 4) (n := K)).symm j))
        / (quartet_a_fin K z j + (freq_node_fin K L n' : ℂ)) = 0 := by
    intro n'
    have hthis := hpoint ((finCongr (index_eq K)).symm n')
    have heq := Fintype.sum_equiv (finProdFinEquiv (m := 4) (n := K))
      (fun p => quartet_D K L z u_c p * quartet_const K L z p
        / (quartet_a K z p + (freq_node_N (2 * K) L ((finCongr (index_eq K)).symm n') : ℂ)))
      (fun j => quartet_D K L z u_c ((finProdFinEquiv (m := 4) (n := K)).symm j)
        * quartet_const K L z ((finProdFinEquiv (m := 4) (n := K)).symm j)
        / (quartet_a_fin K z j + (freq_node_fin K L n' : ℂ)))
      (fun p => by simp [quartet_a_fin, freq_node_fin])
    exact heq.symm.trans hthis
  have hvanish := cauchy_vectors_linear_independent (4 * K) (by omega) (quartet_a_fin K z)
    (freq_node_fin K L) (quartet_a_fin_injective K z hdist) (quartet_a_fin_im_ne_zero K z hz)
    (freq_node_fin_injective K L hLne)
    (fun j => quartet_D K L z u_c ((finProdFinEquiv (m := 4) (n := K)).symm j)
      * quartet_const K L z ((finProdFinEquiv (m := 4) (n := K)).symm j))
    hpoint_fin
  intro p
  have hj := hvanish ((finProdFinEquiv (m := 4) (n := K)) p)
  simp only [Equiv.symm_apply_apply] at hj
  have hc := quartet_const_ne_zero K L hLne z hz p
  exact (mul_eq_zero.mp hj).resolve_right hc

#print axioms quartet_coeffs_vanish


-- STEP I: convert quartet_D vanishing into Cauchy-divided form.
lemma quartet_cauchy_dot_vanish (K : ℕ) (L : ℝ) (hL : 0 < L) (z : Fin K → ℂ)
    (hz : ∀ k, (z k).im ≠ 0) (hdist : quartets_distinct_nodup K z) (hK : 1 ≤ K)
    (u_c : Fin (2 * (2 * K) + 1) → ℂ)
    (hzero : (S_T_off_sum (2 * K) L K z).mulVec u_c = 0) (p : Fin 4 × Fin K) :
    ∑ n : Fin (2 * (2 * K) + 1), u_c n / (quartet_a K z p + (freq_node_N (2 * K) L n : ℂ)) = 0 := by
  have hD := quartet_coeffs_vanish K L hL z hz hdist hK u_c hzero p
  have hne := quartet_const_ne_zero K L (ne_of_gt hL) z hz p
  have hkey : quartet_D K L z u_c p
      = quartet_const K L z p
        * ∑ n : Fin (2 * (2 * K) + 1), u_c n / (quartet_a K z p + (freq_node_N (2 * K) L n : ℂ)) := by
    rw [Finset.mul_sum]
    have heq : ∀ n : Fin (2 * (2 * K) + 1),
        quartet_const K L z p * (u_c n / (quartet_a K z p + (freq_node_N (2 * K) L n : ℂ)))
          = (match p.1 with
              | 0 => vPlusLC (2 * K) L (z p.2) n
              | 1 => vMinusLC (2 * K) L (z p.2) n
              | 2 => vPlusLC (2 * K) L (starRingEnd ℂ (z p.2)) n
              | 3 => vMinusLC (2 * K) L (starRingEnd ℂ (z p.2)) n) * u_c n := by
      intro n
      have hqc := quartet_vec_eq_cauchy K L z p n
      rw [← hqc]
      ring
    rw [Finset.sum_congr rfl (fun n _ => heq n)]
    obtain ⟨p1, p2⟩ := p
    unfold quartet_D dotProduct
    fin_cases p1 <;> rfl
  rw [hkey] at hD
  exact (mul_eq_zero.mp hD).resolve_left hne

#print axioms quartet_cauchy_dot_vanish

-- STEP J: R_vec and the eval-at-negated-pole fact (mirrors the internal argument of
-- cauchy_vectors_linear_independent).
noncomputable def R_vec (K : ℕ) (L : ℝ) (v : Fin (2 * (2 * K) + 1) → ℂ) : Polynomial ℂ :=
  ∑ n : Fin (2 * (2 * K) + 1), C (v n) * ∏ m ∈ Finset.univ.erase n, (X - C (freq_node_N (2 * K) L m : ℂ))

lemma R_vec_linear (K : ℕ) (L : ℝ) (v w : Fin (2 * (2 * K) + 1) → ℂ) (c : ℂ) :
    R_vec K L (fun n => v n - c * w n) = R_vec K L v - c • R_vec K L w := by
  unfold R_vec
  rw [Finset.smul_sum, ← Finset.sum_sub_distrib]
  apply Finset.sum_congr rfl
  intro n _
  show C (v n - c * w n) * ∏ m ∈ Finset.univ.erase n, (X - C (freq_node_N (2 * K) L m : ℂ))
      = C (v n) * ∏ m ∈ Finset.univ.erase n, (X - C (freq_node_N (2 * K) L m : ℂ))
        - c • (C (w n) * ∏ m ∈ Finset.univ.erase n, (X - C (freq_node_N (2 * K) L m : ℂ)))
  rw [C_sub, C_mul, smul_eq_C_mul]
  ring

lemma eval_neg_eq_zero_of_dot_zero (K : ℕ) (L : ℝ) (a : ℂ) (ha : a.im ≠ 0)
    (v : Fin (2 * (2 * K) + 1) → ℂ)
    (hv : ∑ n : Fin (2 * (2 * K) + 1), v n / (a + (freq_node_N (2 * K) L n : ℂ)) = 0) :
    (R_vec K L v).eval (-a) = 0 := by
  have hne : ∀ n : Fin (2 * (2 * K) + 1), a + (freq_node_N (2 * K) L n : ℂ) ≠ 0 := by
    intro n heq
    apply ha
    have : (a + (freq_node_N (2 * K) L n : ℂ)).im = 0 := by rw [heq]; simp
    simpa using this
  have hmul : (∑ n, v n / (a + (freq_node_N (2 * K) L n : ℂ)))
      * ∏ m, (a + (freq_node_N (2 * K) L m : ℂ)) = 0 := by rw [hv]; ring
  have hexpand : (∑ n, v n / (a + (freq_node_N (2 * K) L n : ℂ)))
      * ∏ m, (a + (freq_node_N (2 * K) L m : ℂ))
      = ∑ n, v n * ∏ m ∈ Finset.univ.erase n, (a + (freq_node_N (2 * K) L m : ℂ)) := by
    rw [Finset.sum_mul]
    apply Finset.sum_congr rfl
    intro n _
    rw [← Finset.mul_prod_erase Finset.univ (fun m => a + (freq_node_N (2 * K) L m : ℂ))
      (Finset.mem_univ n), ← mul_assoc, div_mul_cancel₀ _ (hne n)]
  rw [hexpand] at hmul
  unfold R_vec
  rw [eval_finsetSum]
  rw [← hmul]
  apply Finset.sum_congr rfl
  intro n _
  simp only [eval_mul, eval_C, eval_prod, eval_sub, eval_X]
  congr 1
  have hcard : (Finset.univ.erase n).card = 4 * K := by
    rw [Finset.card_erase_of_mem (Finset.mem_univ n), Finset.card_univ, Fintype.card_fin]
    omega
  have heven : Even (4 * K) := ⟨2 * K, by ring⟩
  have hstep : ∏ m ∈ Finset.univ.erase n, (-a - (freq_node_N (2 * K) L m : ℂ))
      = ∏ m ∈ Finset.univ.erase n, (-(a + (freq_node_N (2 * K) L m : ℂ))) := by
    apply Finset.prod_congr rfl; intro m _; ring
  rw [hstep, Finset.prod_neg, hcard, heven.neg_one_pow, one_mul]

#print axioms R_vec_linear
#print axioms eval_neg_eq_zero_of_dot_zero

-- STEP K: R_vec(u_c) vanishes at all 4K quartet points.
lemma quartet_a_neg_eq (K : ℕ) (z : Fin K → ℂ) (p : Fin 4 × Fin K) :
    ∃ p' : Fin 4 × Fin K, -(quartet_a K z p) = quartet_a K z p' := by
  obtain ⟨p1, p2⟩ := p
  fin_cases p1
  · exact ⟨(1, p2), by simp [quartet_a]⟩
  · exact ⟨(0, p2), by simp [quartet_a]⟩
  · exact ⟨(3, p2), by simp [quartet_a]⟩
  · exact ⟨(2, p2), by simp [quartet_a]⟩

theorem R_vec_vanishes_at_quartet (K : ℕ) (L : ℝ) (hL : 0 < L) (z : Fin K → ℂ)
    (hz : ∀ k, (z k).im ≠ 0) (hdist : quartets_distinct_nodup K z) (hK : 1 ≤ K)
    (u_c : Fin (2 * (2 * K) + 1) → ℂ)
    (hzero : (S_T_off_sum (2 * K) L K z).mulVec u_c = 0) (p : Fin 4 × Fin K) :
    (R_vec K L u_c).eval (quartet_a K z p) = 0 := by
  obtain ⟨p', hp'⟩ := quartet_a_neg_eq K z p
  have hdot := quartet_cauchy_dot_vanish K L hL z hz hdist hK u_c hzero p'
  have := eval_neg_eq_zero_of_dot_zero K L (quartet_a K z p') (quartet_a_im_ne_zero K z hz p') u_c hdot
  rwa [← hp', neg_neg] at this

#print axioms quartet_a_neg_eq
#print axioms R_vec_vanishes_at_quartet

-- STEP L: degree bound, matching R_vec(u_c) to c*Q_off_K_poly via the quartet finset.
lemma R_vec_natDegree_le (K : ℕ) (L : ℝ) (v : Fin (2 * (2 * K) + 1) → ℂ) :
    (R_vec K L v).natDegree ≤ 4 * K := by
  unfold R_vec
  apply Polynomial.natDegree_sum_le_of_forall_le
  intro n _
  calc (C (v n) * ∏ m ∈ Finset.univ.erase n, (X - C (freq_node_N (2 * K) L m : ℂ))).natDegree
      ≤ (C (v n)).natDegree
        + (∏ m ∈ Finset.univ.erase n, (X - C (freq_node_N (2 * K) L m : ℂ))).natDegree :=
        Polynomial.natDegree_mul_le
    _ = (∏ m ∈ Finset.univ.erase n, (X - C (freq_node_N (2 * K) L m : ℂ))).natDegree := by
        simp [Polynomial.natDegree_C]
    _ = (Finset.univ.erase n).card := by
        rw [Polynomial.natDegree_prod_of_monic _ _ (fun m _ => Polynomial.monic_X_sub_C _)]
        simp
    _ = 4 * K := by
        rw [Finset.card_erase_of_mem (Finset.mem_univ n), Finset.card_univ, Fintype.card_fin]
        omega

noncomputable def quartet_finset (K : ℕ) (z : Fin K → ℂ) : Finset ℂ :=
  Finset.image (quartet_a K z) Finset.univ

lemma quartet_finset_card (K : ℕ) (z : Fin K → ℂ) (hdist : quartets_distinct_nodup K z) :
    (quartet_finset K z).card = 4 * K := by
  unfold quartet_finset
  rw [Finset.card_image_of_injective _ (quartet_a_injective K z hdist), Finset.card_univ]
  simp [Fintype.card_prod]

#print axioms R_vec_natDegree_le
#print axioms quartet_finset_card

-- Connect Q_off_K_poly (mapped to ℂ) to Q_off_K_complex, and show vanishing at quartet points.
lemma Q_off_K_poly_map_eval (K : ℕ) (z : Fin K → ℂ) (w : ℂ) :
    ((Q_off_K_poly K z).map (algebraMap ℝ ℂ)).eval w = Q_off_K_complex K z w := by
  unfold Q_off_K_poly Q_off_K_complex
  rw [Polynomial.map_prod, Polynomial.eval_prod]
  apply Finset.prod_congr rfl
  intro k _
  simp only [Polynomial.map_add, Polynomial.map_sub, Polynomial.map_mul, Polynomial.map_pow,
    Polynomial.map_C, Polynomial.map_X, Polynomial.eval_add, Polynomial.eval_sub,
    Polynomial.eval_mul, Polynomial.eval_pow, Polynomial.eval_C, Polynomial.eval_X,
    Complex.coe_algebraMap]
  push_cast
  have h1 : (z k) ^ 2 + (starRingEnd ℂ (z k)) ^ 2
      = 2 * (((z k).re : ℂ) ^ 2 - ((z k).im : ℂ) ^ 2) := by
    have e1 := Complex.add_conj (z k)
    have e2 := Complex.mul_conj (z k)
    have hsplit : (z k) ^ 2 + (starRingEnd ℂ (z k)) ^ 2
        = (z k + starRingEnd ℂ (z k)) ^ 2 - 2 * (z k * starRingEnd ℂ (z k)) := by ring
    rw [hsplit, e1, e2, Complex.normSq_apply]
    push_cast
    ring
  have h2 : (((z k).re : ℂ) ^ 2 + ((z k).im : ℂ) ^ 2) ^ 2
      = ((z k) * (starRingEnd ℂ (z k))) ^ 2 := by
    rw [Complex.mul_conj, Complex.normSq_apply]
    push_cast
    ring
  rw [← h1, h2]
  ring

lemma Q_off_K_poly_map_vanishes_at_quartet (K : ℕ) (z : Fin K → ℂ) (p : Fin 4 × Fin K) :
    ((Q_off_K_poly K z).map (algebraMap ℝ ℂ)).eval (quartet_a K z p) = 0 := by
  rw [Q_off_K_poly_map_eval]
  unfold Q_off_K_complex quartet_a
  obtain ⟨p1, p2⟩ := p
  apply Finset.prod_eq_zero (Finset.mem_univ p2)
  fin_cases p1 <;> simp <;> ring_nf <;> simp

#print axioms Q_off_K_poly_map_eval
#print axioms Q_off_K_poly_map_vanishes_at_quartet

theorem R_vec_eq_smul_Q_off_K (K : ℕ) (L : ℝ) (hL : 0 < L) (z : Fin K → ℂ)
    (hz : ∀ k, (z k).im ≠ 0) (hdist : quartets_distinct_nodup K z) (hK : 1 ≤ K)
    (u_c : Fin (2 * (2 * K) + 1) → ℂ)
    (hzero : (S_T_off_sum (2 * K) L K z).mulVec u_c = 0) :
    R_vec K L u_c = (R_vec K L u_c).leadingCoeff • (Q_off_K_poly K z).map (algebraMap ℝ ℂ) := by
  set R := R_vec K L u_c with hRdef
  set Q := (Q_off_K_poly K z).map (algebraMap ℝ ℂ) with hQdef
  have hQmonic : Q.Monic := (Q_off_K_poly_monic K z).map _
  have hQdeg : Q.natDegree = 4 * K := by
    rw [hQdef, Polynomial.natDegree_map_eq_of_injective
      (algebraMap ℝ ℂ).injective, Q_off_K_poly_natDegree]
  have hRdeg : R.natDegree ≤ 4 * K := R_vec_natDegree_le K L u_c
  have hcard : (quartet_finset K z).card = 4 * K := quartet_finset_card K z hdist
  have hRvan : ∀ x ∈ quartet_finset K z, R.eval x = 0 := by
    intro x hx
    unfold quartet_finset at hx
    simp only [Finset.mem_image, Finset.mem_univ, true_and] at hx
    obtain ⟨p, hp⟩ := hx
    rw [← hp, hRdef]
    exact R_vec_vanishes_at_quartet K L hL z hz hdist hK u_c hzero p
  have hQvan : ∀ x ∈ quartet_finset K z, Q.eval x = 0 := by
    intro x hx
    unfold quartet_finset at hx
    simp only [Finset.mem_image, Finset.mem_univ, true_and] at hx
    obtain ⟨p, hp⟩ := hx
    rw [← hp, hQdef]
    exact Q_off_K_poly_map_vanishes_at_quartet K z p
  by_cases hRd : R.natDegree = 4 * K
  · set c := R.leadingCoeff with hcdef
    have hQc4K : Q.coeff (4 * K) = 1 := by
      rw [← hQdeg]; exact hQmonic
    have hRc4K : R.coeff (4 * K) = c := by rw [hcdef, Polynomial.leadingCoeff, hRd]
    set D := R - c • Q with hDdef
    have hD4K : D.coeff (4 * K) = 0 := by
      rw [hDdef]
      simp only [Polynomial.coeff_sub, Polynomial.coeff_smul, smul_eq_mul]
      rw [hRc4K, hQc4K, mul_one, sub_self]
    have hDdeg_le : D.degree ≤ (4 * K : ℕ) := by
      rw [hDdef]
      calc (R - c • Q).degree ≤ max R.degree (c • Q).degree := Polynomial.degree_sub_le _ _
        _ ≤ (4 * K : ℕ) := by
            apply max_le
            · exact le_trans Polynomial.degree_le_natDegree (by exact_mod_cast hRdeg)
            · calc (c • Q).degree ≤ Q.degree := Polynomial.degree_smul_le c Q
                _ ≤ (4 * K : ℕ) := by
                    rw [← hQdeg]; exact Polynomial.degree_le_natDegree
    have hDvan : ∀ x ∈ quartet_finset K z, D.eval x = 0 := by
      intro x hx
      rw [hDdef]
      simp only [Polynomial.eval_sub, Polynomial.eval_smul, smul_eq_mul]
      rw [hRvan x hx, hQvan x hx, mul_zero, sub_zero]
    have hDeq0 : D = 0 := by
      by_cases hD0 : D = 0
      · exact hD0
      · exfalso
        have hDdeg_lt : D.degree < (4 * K : ℕ) := by
          rcases lt_or_eq_of_le hDdeg_le with h | h
          · exact h
          · exfalso
            have : D.natDegree = 4 * K := by
              have := Polynomial.degree_eq_natDegree hD0
              rw [this] at h
              exact_mod_cast h
            have : D.coeff (4 * K) = D.leadingCoeff := by rw [Polynomial.leadingCoeff, this]
            rw [hD4K] at this
            exact hD0 (Polynomial.leadingCoeff_eq_zero.mp this.symm)
        have := Polynomial.eq_zero_of_degree_lt_of_eval_finset_eq_zero (quartet_finset K z)
          (by rw [hcard]; exact hDdeg_lt) hDvan
        exact hD0 this
    exact sub_eq_zero.mp hDeq0
  · have hRdeg_lt : R.degree < (4 * K : ℕ) := by
      rcases lt_or_eq_of_le hRdeg with h | h
      · exact lt_of_le_of_lt Polynomial.degree_le_natDegree (by exact_mod_cast h)
      · exact absurd h hRd
    have hR0 : R = 0 :=
      Polynomial.eq_zero_of_degree_lt_of_eval_finset_eq_zero (quartet_finset K z)
        (by rw [hcard]; exact hRdeg_lt) hRvan
    rw [hR0]
    simp

#print axioms R_vec_eq_smul_Q_off_K

-- FINAL STEP: connect to c_off and conclude.
lemma R_vec_eval_node (K : ℕ) (L : ℝ) (v : Fin (2 * (2 * K) + 1) → ℂ) (j : Fin (2 * (2 * K) + 1)) :
    (R_vec K L v).eval (freq_node_N (2 * K) L j : ℂ)
      = v j * ∏ m ∈ Finset.univ.erase j, ((freq_node_N (2 * K) L j : ℂ) - (freq_node_N (2 * K) L m : ℂ)) := by
  unfold R_vec
  rw [eval_finsetSum]
  rw [Finset.sum_eq_single j]
  · simp only [eval_mul, eval_C, eval_prod, eval_sub, eval_X]
  · intro n _ hn
    have : j ∈ Finset.univ.erase n := Finset.mem_erase.mpr ⟨Ne.symm hn, Finset.mem_univ j⟩
    simp only [eval_mul, eval_C, eval_prod, eval_sub, eval_X]
    rw [Finset.prod_eq_zero this]
    · ring
    · ring
  · intro h; exact absurd (Finset.mem_univ j) h

lemma freq_node_N_injective (K : ℕ) (L : ℝ) (hL : L ≠ 0) :
    Function.Injective (freq_node_N (2 * K) L) := by
  intro a b hab
  unfold freq_node_N at hab
  have hpi : (2 : ℝ) * Real.pi ≠ 0 := by positivity
  have h1 : 2 * Real.pi * (((a : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ))
      = 2 * Real.pi * (((b : ℕ) : ℝ) - ((2 * K : ℕ) : ℝ)) := by
    have hmul := congrArg (· * L) hab
    simp only [div_mul_cancel₀ _ hL] at hmul
    exact hmul
  have h2 : ((a : ℕ) : ℝ) = ((b : ℕ) : ℝ) := by
    have := mul_left_cancel₀ hpi h1
    linarith
  exact Fin.ext (by exact_mod_cast h2)

lemma R_vec_eq_zero_imp_eq_zero (K : ℕ) (L : ℝ) (hL : L ≠ 0)
    (v : Fin (2 * (2 * K) + 1) → ℂ) (h0 : R_vec K L v = 0) : v = 0 := by
  funext j
  have heval := R_vec_eval_node K L v j
  rw [h0] at heval
  simp only [eval_zero] at heval
  have hprod_ne : (∏ m ∈ Finset.univ.erase j,
      ((freq_node_N (2 * K) L j : ℂ) - (freq_node_N (2 * K) L m : ℂ))) ≠ 0 := by
    apply Finset.prod_ne_zero_iff.mpr
    intro m hm
    simp only [Finset.mem_erase] at hm
    intro heq
    apply hm.1
    have hre : (freq_node_N (2 * K) L j : ℝ) = (freq_node_N (2 * K) L m : ℝ) := by
      have := sub_eq_zero.mp heq
      exact_mod_cast this
    exact (freq_node_N_injective K L hL hre).symm
  exact (mul_eq_zero.mp heval.symm).resolve_right hprod_ne

lemma R_vec_c_off_eq_Q_off_K_poly_map (K : ℕ) (L : ℝ) (hL : L ≠ 0) (z : Fin K → ℂ) :
    R_vec K L (fun n => (c_off K L z n : ℂ)) = (Q_off_K_poly K z).map (algebraMap ℝ ℂ) := by
  set lhs := R_vec K L (fun n => (c_off K L z n : ℂ)) with hlhs
  set rhs := (Q_off_K_poly K z).map (algebraMap ℝ ℂ) with hrhs
  have hlhs_deg : lhs.natDegree ≤ 4 * K := R_vec_natDegree_le K L _
  have hrhs_deg : rhs.natDegree = 4 * K := by
    rw [hrhs, Polynomial.natDegree_map_eq_of_injective (algebraMap ℝ ℂ).injective,
      Q_off_K_poly_natDegree]
  have hrhs_ne_zero : rhs ≠ 0 :=
    Polynomial.Monic.ne_zero ((Q_off_K_poly_monic K z).map _)
  have heq_nodes : ∀ j : Fin (2 * (2 * K) + 1), lhs.eval (freq_node_N (2 * K) L j : ℂ)
      = rhs.eval (freq_node_N (2 * K) L j : ℂ) := by
    intro j
    rw [hlhs, R_vec_eval_node]
    have hprod_ne : (∏ m ∈ Finset.univ.erase j,
        (freq_node_N (2 * K) L j - freq_node_N (2 * K) L m)) ≠ 0 := by
      apply Finset.prod_ne_zero_iff.mpr
      intro m hm
      simp only [Finset.mem_erase] at hm
      intro heq
      exact hm.1 (freq_node_N_injective K L hL (sub_eq_zero.mp heq)).symm
    have hc : c_off K L z j = Q_off_K K z (freq_node_N (2 * K) L j)
        / ∏ m ∈ Finset.univ.erase j, (freq_node_N (2 * K) L j - freq_node_N (2 * K) L m) :=
      c_off_eq_div_prod K L z j
    rw [hc]
    have hcast : rhs.eval (freq_node_N (2 * K) L j : ℂ)
        = (((Q_off_K_poly K z).eval (freq_node_N (2 * K) L j) : ℝ) : ℂ) := by
      rw [hrhs, Q_off_K_poly_map_eval, Q_off_K_complex_eq_ofReal, Q_off_K_poly_eval]
    rw [hcast, Q_off_K_poly_eval]
    push_cast
    rw [div_mul_cancel₀ _ (by exact_mod_cast hprod_ne : (∏ m ∈ Finset.univ.erase j,
      ((freq_node_N (2 * K) L j : ℂ) - (freq_node_N (2 * K) L m : ℂ))) ≠ 0)]
  have hD0 : lhs - rhs = 0 := by
    by_contra hne
    have hdeg_lt : (lhs - rhs).degree < (4 * K + 1 : ℕ) := by
      calc (lhs - rhs).degree ≤ max lhs.degree rhs.degree := Polynomial.degree_sub_le _ _
        _ < (4 * K + 1 : ℕ) := by
            apply max_lt
            · exact lt_of_le_of_lt Polynomial.degree_le_natDegree (by exact_mod_cast
                (Nat.lt_succ_of_le hlhs_deg))
            · rw [Polynomial.degree_eq_natDegree hrhs_ne_zero, hrhs_deg]
              exact_mod_cast Nat.lt_succ_self _
    have hvanish : ∀ x ∈ (Finset.image (fun j => (freq_node_N (2 * K) L j : ℂ)) Finset.univ),
        (lhs - rhs).eval x = 0 := by
      intro x hx
      simp only [Finset.mem_image, Finset.mem_univ, true_and] at hx
      obtain ⟨j, hj⟩ := hx
      rw [← hj, Polynomial.eval_sub, heq_nodes j, sub_self]
    have hcard : (Finset.image (fun j => (freq_node_N (2 * K) L j : ℂ))
        (Finset.univ : Finset (Fin (2 * (2 * K) + 1)))).card = 4 * K + 1 := by
      rw [Finset.card_image_of_injective _ (fun a b hab =>
        freq_node_N_injective K L hL (by exact_mod_cast hab)),
        Finset.card_univ, Fintype.card_fin]
      omega
    exact hne (Polynomial.eq_zero_of_degree_lt_of_eval_finset_eq_zero _
      (by rw [hcard]; exact hdeg_lt) hvanish)
  exact sub_eq_zero.mp hD0

theorem ker_S_T_off_real_eq_span_e_off (K : ℕ) (L : ℝ) (hL : 0 < L) (z : Fin K → ℂ)
    (hz : ∀ k, (z k).im ≠ 0) (hdist : quartets_distinct_nodup K z) (hK : 1 ≤ K)
    (u : Fin (2 * (2 * K) + 1) → ℝ)
    (hu : (S_T_off_real (2 * K) L K z).mulVec u = 0) :
    ∃ c : ℝ, u = c • (e_off K L z) := by
  have hLne : L ≠ 0 := ne_of_gt hL
  set u_c : Fin (2 * (2 * K) + 1) → ℂ := fun n => (u n : ℂ) with hu_c
  have hzero : (S_T_off_sum (2 * K) L K z).mulVec u_c = 0 := complexify_kernel (2 * K) L K z u hu
  have hmain := R_vec_eq_smul_Q_off_K K L hL z hz hdist hK u_c hzero
  have hcoff := R_vec_c_off_eq_Q_off_K_poly_map K L hLne z
  set c : ℂ := (R_vec K L u_c).leadingCoeff with hcdef
  have hcombine : R_vec K L (fun n => u_c n - c * (c_off K L z n : ℂ)) = 0 := by
    rw [R_vec_linear, hmain, hcoff]
    simp
  have hveq := R_vec_eq_zero_imp_eq_zero K L hLne _ hcombine
  have hpointwise : ∀ n, u_c n = c * (c_off K L z n : ℂ) := by
    intro n
    have := congrFun hveq n
    simpa [sub_eq_zero] using this
  have hex : ∃ j, c_off K L z j ≠ 0 := by
    by_contra hcon
    push_neg at hcon
    have hsum := c_off_sum_eq_one K L hLne z
    rw [Finset.sum_eq_zero (fun j _ => hcon j)] at hsum
    exact absurd hsum (by norm_num)
  obtain ⟨j, hj⟩ := hex
  have hc_real : ∃ cr : ℝ, c = (cr : ℂ) := by
    refine ⟨(u j) / (c_off K L z j), ?_⟩
    have hp := hpointwise j
    simp only [hu_c] at hp
    have hjc : (c_off K L z j : ℂ) ≠ 0 := by exact_mod_cast hj
    push_cast
    exact (eq_div_iff hjc).mpr hp.symm
  obtain ⟨cr, hcr⟩ := hc_real
  refine ⟨cr / Real.sqrt L, ?_⟩
  funext n
  have hp := hpointwise n
  rw [hu_c, hcr] at hp
  simp only [Pi.smul_apply, smul_eq_mul]
  unfold e_off
  have : (u n : ℂ) = ((cr * c_off K L z n : ℝ) : ℂ) := by push_cast; linear_combination hp
  have hreal : u n = cr * c_off K L z n := by exact_mod_cast this
  rw [hreal]
  field_simp

#print axioms R_vec_eval_node
#print axioms R_vec_eq_zero_imp_eq_zero
#print axioms R_vec_c_off_eq_Q_off_K_poly_map
#print axioms ker_S_T_off_real_eq_span_e_off
