import Mathlib

open Finset

/-- **Eq (5.21) (Shi, arXiv:2609.04908)**: the difference-basis diagonal identity underlying
Theorem 5.2's pencil identification. In coordinates: for a finite index set `s` with a
distinguished node `z ∈ s` where `ν z = 0`, barycentric weights `c` on `s` summing to `1`
(§2162's `barycentric_weights_sum_eq_one`, taken here as an explicit hypothesis rather than
re-derived), and `A m n := ν n * (if m = n then 1 else 0) - ν n * c m` (Eq 5.20's
`(δ_mn - c_{m,N,L})·ν_{n,L}`), the diagonal action on the difference basis vector `q_n := u_n - u_z`
minus the pencil action `∑_{m≠z} A m n • q_m` equals `ν n • c` (the paper's `(ν_n/√L)·e_{N,L}`,
with the `√L` rescaling of `e := √L·c` absorbed away since it cancels on both sides — proving the
mathematically essential content directly in terms of the barycentric weights `c`). Fully general
over ℝ, no zeta/RH dependence beyond the stated hypotheses. -/
theorem difference_basis_diagonal_identity {ι : Type*} [DecidableEq ι]
    (s : Finset ι) (z : ι) (hz : z ∈ s) (ν c : ι → ℝ) (hνz : ν z = 0)
    (hc : ∑ i ∈ s, c i = 1) (n : ι) (hn : n ∈ s) (hnz : n ≠ z) (i : ι) (hi : i ∈ s) :
    ν i * ((if i = n then (1 : ℝ) else 0) - (if i = z then 1 else 0))
      - ∑ m ∈ s.erase z, (ν n * (if m = n then (1 : ℝ) else 0) - ν n * c m)
          * ((if i = m then (1 : ℝ) else 0) - (if i = z then 1 else 0))
      = ν n * c i := by
  have hnerase : n ∈ s.erase z := Finset.mem_erase.mpr ⟨hnz, hn⟩
  have hcerase : ∑ m ∈ s.erase z, c m = 1 - c z := by
    have hadd := Finset.sum_erase_add s c hz
    linarith [hadd, hc]
  have hsum_ite_n : ∑ m ∈ s.erase z, (if m = n then (1 : ℝ) else 0) = 1 := by
    rw [Finset.sum_ite_eq' (s.erase z) n (fun _ => (1:ℝ)), if_pos hnerase]
  have hsplit : ∑ m ∈ s.erase z, (ν n * (if m = n then (1 : ℝ) else 0) - ν n * c m)
        * ((if i = m then (1 : ℝ) else 0) - (if i = z then 1 else 0))
      = (∑ m ∈ s.erase z, (if i = m then (1:ℝ) else 0)
            * (ν n * (if m = n then (1:ℝ) else 0) - ν n * c m))
        - (if i = z then (1:ℝ) else 0)
            * ∑ m ∈ s.erase z, (ν n * (if m = n then (1:ℝ) else 0) - ν n * c m) := by
    rw [show (fun m => (ν n * (if m = n then (1 : ℝ) else 0) - ν n * c m)
          * ((if i = m then (1 : ℝ) else 0) - (if i = z then 1 else 0)))
        = (fun m => (if i = m then (1:ℝ) else 0) * (ν n * (if m = n then (1:ℝ) else 0) - ν n * c m)
            - (if i = z then (1:ℝ) else 0) * (ν n * (if m = n then (1:ℝ) else 0) - ν n * c m))
        from by funext m; ring]
    rw [Finset.sum_sub_distrib, Finset.mul_sum]
  rw [hsplit]
  have hpart1raw : ∑ m ∈ s.erase z, (if i = m then (1:ℝ) else 0)
      * (ν n * (if m = n then (1:ℝ) else 0) - ν n * c m)
      = ∑ m ∈ s.erase z, (if m = i then (ν n * (if m = n then (1:ℝ) else 0) - ν n * c m) else 0) := by
    apply Finset.sum_congr rfl
    intro m _
    by_cases h : m = i
    · rw [if_pos h, if_pos h.symm]; ring
    · rw [if_neg h, if_neg (Ne.symm h)]; ring
  have hpart1 : ∑ m ∈ s.erase z, (if i = m then (1:ℝ) else 0)
      * (ν n * (if m = n then (1:ℝ) else 0) - ν n * c m)
      = if i ∈ s.erase z then ν n * (if i = n then (1:ℝ) else 0) - ν n * c i else 0 := by
    rw [hpart1raw]
    exact Finset.sum_ite_eq' (s.erase z) i (fun m => ν n * (if m = n then (1:ℝ) else 0) - ν n * c m)
  have hpart2 : ∑ m ∈ s.erase z, (ν n * (if m = n then (1:ℝ) else 0) - ν n * c m) = ν n * c z := by
    have hexp : ∑ m ∈ s.erase z, (ν n * (if m = n then (1:ℝ) else 0) - ν n * c m)
        = ν n * (∑ m ∈ s.erase z, (if m = n then (1:ℝ) else 0)) - ν n * (∑ m ∈ s.erase z, c m) := by
      rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_sub_distrib]
    rw [hexp, hsum_ite_n, hcerase]
    ring
  rw [hpart1, hpart2]
  by_cases hiz : i = z
  · have hnotin : i ∉ s.erase z := by
      rw [hiz]; simp
    rw [if_neg (by rw [hiz]; exact hnz.symm : i ≠ n), if_pos hiz, if_neg hnotin]
    rw [hiz, hνz]
    ring
  · have hiine : i ∈ s.erase z := Finset.mem_erase.mpr ⟨hiz, hi⟩
    rw [if_neg hiz, if_pos hiine]
    by_cases hin : i = n
    · subst hin; ring
    · rw [if_neg hin]; ring

#print axioms difference_basis_diagonal_identity
