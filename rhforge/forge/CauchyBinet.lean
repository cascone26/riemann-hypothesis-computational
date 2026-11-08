/-
CauchyBinet.lean (HP RH center, LEDGER S2917, 2026-11-08)
The Cauchy–Binet formula (absent from Mathlib), built on CauchyBinetExpand (S2916):
  det (A * B) = Σ_{S ⊆ Fin n, |S| = k} det (A[:, S]) * det (B[S, :]),
with S enumerated increasingly by `orderEmbOfFin`.
-/
import Mathlib
import CauchyBinetExpand

open Matrix Finset Equiv

namespace CauchyBinet

variable {R : Type*} [CommRing R] {k n : ℕ}

/-- The increasing enumeration of a `k`-subset. -/
noncomputable abbrev emb (S : {S : Finset (Fin n) // S.card = k}) : Fin k ↪o Fin n :=
  S.1.orderEmbOfFin S.2

lemma term_expand (A : Matrix (Fin k) (Fin n) R) (B : Matrix (Fin n) (Fin k) R)
    (S : {S : Finset (Fin n) // S.card = k}) :
    det (A.submatrix id (emb S)) * det (B.submatrix (emb S) id)
      = ∑ τ : Perm (Fin k), det (A.submatrix id (emb S ∘ τ)) * ∏ i, B (emb S (τ i)) i := by
  rw [det_apply' (B.submatrix _ id), Finset.mul_sum]
  refine sum_congr rfl fun τ _ => ?_
  have h : A.submatrix id (emb S ∘ τ) = (A.submatrix id (emb S)).submatrix id τ := by
    ext; rfl
  rw [h, det_permute']
  simp only [submatrix_apply, id]
  ring

lemma range_emb_comp (S : {S : Finset (Fin n) // S.card = k}) (τ : Perm (Fin k)) :
    Set.range (emb S ∘ τ) = (S.1 : Set (Fin n)) := by
  rw [Set.range_comp, τ.surjective.range_eq, Set.image_univ, Finset.range_orderEmbOfFin]

/-- **Cauchy–Binet.** -/
theorem cauchy_binet (A : Matrix (Fin k) (Fin n) R) (B : Matrix (Fin n) (Fin k) R) :
    det (A * B) = ∑ S : {S : Finset (Fin n) // S.card = k},
      det (A.submatrix id (emb S)) * det (B.submatrix (emb S) id) := by
  rw [det_mul_expand_inj]
  simp_rw [term_expand A B]
  rw [← Fintype.sum_sigma (fun x : (Σ _ : {S : Finset (Fin n) // S.card = k}, Perm (Fin k)) =>
      det (A.submatrix id (emb x.1 ∘ x.2)) * ∏ i, B (emb x.1 (x.2 i)) i)]
  symm
  refine Finset.sum_bij (fun x _ => emb x.1 ∘ x.2) ?_ ?_ ?_ ?_
  · intro x _
    simp only [mem_filter, mem_univ, true_and]
    exact (emb x.1).injective.comp x.2.injective
  · rintro ⟨S, τ⟩ _ ⟨S', τ'⟩ _ h
    have hS : S = S' := by
      apply Subtype.ext
      have := congrArg Set.range h
      simp only at this
      rw [range_emb_comp, range_emb_comp] at this
      exact_mod_cast this
    subst hS
    have : τ = τ' := by
      ext i
      have := congrFun h i
      simp only [Function.comp_apply] at this
      exact congrArg Fin.val ((emb S).injective this)
    subst this; rfl
  · intro p hp
    simp only [mem_filter, mem_univ, true_and] at hp
    have hcard : (univ.image p).card = k := by
      rw [card_image_of_injective _ hp, card_univ, Fintype.card_fin]
    let S : {S : Finset (Fin n) // S.card = k} := ⟨univ.image p, hcard⟩
    have hmem : ∀ i, p i ∈ Set.range (emb S) := by
      intro i; rw [Finset.range_orderEmbOfFin]; simp [S]
    let f : Fin k → Fin k := fun i => (hmem i).choose
    have hf : ∀ i, emb S (f i) = p i := fun i => (hmem i).choose_spec
    have hfi : Function.Injective f := by
      intro a b hab; apply hp; rw [← hf a, ← hf b, hab]
    refine ⟨⟨S, Equiv.ofBijective f (Finite.injective_iff_bijective.mp hfi)⟩, mem_univ _, ?_⟩
    funext i; exact hf i
  · intro x _; rfl

end CauchyBinet

#print axioms CauchyBinet.term_expand
#print axioms CauchyBinet.cauchy_binet
