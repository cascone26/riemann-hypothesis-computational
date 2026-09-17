import Mathlib.LinearAlgebra.Matrix.Hermitian
import Mathlib.Tactic

/-!
# §1771 — Simple eigenvalue + commuting parity involution ⟹ definite parity of the eigenvector

Suzuki 2606.09096 §5.3 (parity of the lowest eigenfunction) in abstract form, completing the class-C lane
(§1763–1768).  Once the ground state is SIMPLE (§5.2 positivity-improving Dirichlet form; our §1764 Kramers /
§1768 class-C give the structural simplicity), a commuting parity involution `P` (`P² = 1`, `A ∘ P = P ∘ A`)
forces the eigenvector to have DEFINITE parity `P v = ± v`.  This is exactly the even/odd dichotomy Yoshida [17]
and Suzuki use to split RH into an even sector (RH up to real zeros) and an odd sector (full RH).

Statement: over any field `𝕜` and `𝕜`-vector space `E`, if `A` commutes with an involution `P`, and `μ` is a
SIMPLE eigenvalue of `A` (its eigenspace is spanned by `v`), then `P v = c • v` with `c = ±1`.  EXIT0 / std-axioms.
-/

variable {𝕜 : Type*} [Field 𝕜]
variable {E : Type*} [AddCommGroup E] [Module 𝕜 E] [NoZeroSMulDivisors 𝕜 E]

/-- **Simple eigenvalue ⟹ definite parity.**  `A` commutes with the involution `P` (`P∘P = id`); `μ` is a
simple eigenvalue with eigenvector `v` (every `μ`-eigenvector is a scalar multiple of `v`).  Then `v` is an
eigenvector of `P` with eigenvalue `±1` — i.e. it is even (`P v = v`) or odd (`P v = -v`). -/
theorem simple_eigen_definite_parity
    (A P : E →ₗ[𝕜] E)
    (hcomm : A ∘ₗ P = P ∘ₗ A)
    (hinv : P ∘ₗ P = LinearMap.id)
    (v : E) (hv : v ≠ 0) (μ : 𝕜) (heig : A v = μ • v)
    (hsimple : ∀ u, A u = μ • u → ∃ c : 𝕜, u = c • v) :
    ∃ c : 𝕜, P v = c • v ∧ (c = 1 ∨ c = -1) := by
  -- (1) `P v` is again a `μ`-eigenvector, since `A (P v) = P (A v) = P (μ • v) = μ • (P v)`.
  have hPeig : A (P v) = μ • (P v) := by
    have h := LinearMap.congr_fun hcomm v
    simp only [LinearMap.comp_apply] at h
    rw [h, heig, map_smul]
  -- (2) simplicity gives `P v = c • v`.
  obtain ⟨c, hc⟩ := hsimple (P v) hPeig
  -- (3) `P (P v) = v` from the involution law.
  have hPP : P (P v) = v := by
    have h := LinearMap.congr_fun hinv v
    simpa using h
  -- (4) hence `c² • v = v`.
  have hc2 : (c ^ 2) • v = v := by
    have h1 : P (P v) = (c ^ 2) • v := by
      rw [hc, map_smul, hc, smul_smul, ← sq]
    rw [← h1, hPP]
  -- (5) `(c² - 1) • v = 0`, and `v ≠ 0`, so `c² = 1`.
  have hzero : (c ^ 2 - 1) • v = 0 := by
    rw [sub_smul, hc2, one_smul, sub_self]
  have hc21 : c ^ 2 - 1 = 0 := by
    rcases smul_eq_zero.mp hzero with h | h
    · exact h
    · exact absurd h hv
  -- (6) `c² = 1 ⟹ (c-1)(c+1) = 0 ⟹ c = ±1`.
  have hfac : (c - 1) * (c + 1) = 0 := by
    have hexp : (c - 1) * (c + 1) = c ^ 2 - 1 := by ring
    rw [hexp, hc21]
  refine ⟨c, hc, ?_⟩
  rcases mul_eq_zero.mp hfac with h | h
  · exact Or.inl (sub_eq_zero.mp h)
  · exact Or.inr (by rwa [add_eq_zero_iff_eq_neg] at h)
