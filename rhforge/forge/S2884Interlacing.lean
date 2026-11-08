/-
S2884Interlacing.lean (HP RH center, LEDGER S2905b, 2026-11-07)
CAPSTONE: S2884 (one-step interlacing of Shi's unshifted contrast pencil) in Lean.
Given the pencil data (S, ν odd, displacement, nested even/odd families E, O on m+1 levels as in PencilTri),
the pencil block B_ij := E_iᵀ (S D) O_j has zero top-right column, so the squared pencil eigenvalues at level m
interlace with those at level m+1:   θ²_{k+1}(m+1) ≤ θ²_k(m) ≤ θ²_k(m+1).
-/
import Mathlib
import PencilLowerTriangular
import PencilThetaInterlace

open Matrix

namespace S2884

theorem pencil_interlacing {N m : ℕ} (S : Matrix (Fin (2 * N + 1)) (Fin (2 * N + 1)) ℝ)
    (ν : Fin (2 * N + 1) → ℝ) (hodd : ∀ i, ν i.rev = -ν i) (c : Fin (2 * N + 1)) (hc : ν c = 0)
    (b : Fin (2 * N + 1) → ℝ)
    (hdisp : diagonal ν * S - S * diagonal ν = vecMulVec b (fun _ => 1) - vecMulVec (fun _ => 1) b)
    (K : Fin (m + 1) → Set (Fin (2 * N + 1)))
    (E O : Fin (m + 1) → Fin (2 * N + 1) → ℝ)
    (hE : ∀ i, E i ∈ Submodule.span ℝ (ParityDisp.evc c '' K i))
    (hOmem : ∀ j, ∃ Kj : Set (Fin (2 * N + 1)), O j ∈ Submodule.span ℝ (ParityDisp.od '' Kj))
    (hO : ∀ i j, i < j → ∀ w ∈ Submodule.span ℝ (ParityDisp.od '' K i), w ⬝ᵥ (S *ᵥ O j) = 0)
    (hBig : (toEuclideanLin ((Matrix.of fun i j => E i ⬝ᵥ ((S * diagonal ν) *ᵥ O j)) *
      (Matrix.of fun i j => E i ⬝ᵥ ((S * diagonal ν) *ᵥ O j))ᵀ)).IsSymmetric)
    (hSmall : (toEuclideanLin (((Matrix.of fun i j => E i ⬝ᵥ ((S * diagonal ν) *ᵥ O j)).submatrix
      Fin.castSucc Fin.castSucc) * ((Matrix.of fun i j => E i ⬝ᵥ ((S * diagonal ν) *ᵥ O j)).submatrix
      Fin.castSucc Fin.castSucc)ᵀ)).IsSymmetric) (k : Fin m) :
    hBig.eigenvalues (MatrixInterlace.finrank_e (m + 1)) k.succ
        ≤ hSmall.eigenvalues (MatrixInterlace.finrank_e m) k ∧
      hSmall.eigenvalues (MatrixInterlace.finrank_e m) k
        ≤ hBig.eigenvalues (MatrixInterlace.finrank_e (m + 1)) k.castSucc := by
  apply ThetaInterlace.theta_sq_interlace
  intro i
  simp only [Matrix.of_apply]
  exact PencilTri.lower_triangular S ν hodd c hc b hdisp K E O hE hOmem hO (Fin.castSucc_lt_last i)

end S2884

#print axioms S2884.pencil_interlacing
