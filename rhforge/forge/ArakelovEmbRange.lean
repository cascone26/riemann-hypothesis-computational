import ArakelovEmbFiniteHom

set_option maxHeartbeats 1000000

open IsDedekindDomain NumberField

namespace ArakelovDivisor

variable (K : Type*) [Field K] [NumberField K]

/-!
# §1675 — the finite-divisor subgroup: kernel and range of `embFiniteHom`

Completes the §1650 bundling: the embedding has trivial kernel (it is injective), and its range —
the subgroup of finite Arakelov divisors — is characterized exactly as the divisors with vanishing
archimedean part. The finite divisors sit inside the Arakelov group as `{D | D.infinite = 0}`,
the kernel of the archimedean projection.
-/

theorem ker_embFiniteHom : (embFiniteHom K).ker = ⊥ :=
  (AddMonoidHom.ker_eq_bot_iff _).mpr (embFiniteHom_injective K)

theorem mem_range_embFiniteHom_iff {x : ArakelovDivisor K} :
    x ∈ (embFiniteHom K).range ↔ x.infinite = 0 := by
  constructor
  · rintro ⟨D, rfl⟩
    rfl
  · intro hx
    refine ⟨x.finite, ?_⟩
    apply ext
    · rfl
    · funext w
      rw [hx]
      rfl

end ArakelovDivisor

#print axioms ArakelovDivisor.ker_embFiniteHom
#print axioms ArakelovDivisor.mem_range_embFiniteHom_iff
