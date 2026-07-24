# SUZUKI GAP MAP — unconditional vs under-RH in the screw-function program (2026-07-20)

Companion to OBSTRUCTION_MAP.md, for the CURRENT frontier: Suzuki arXiv:2301.00421 ("[14]", Canad. J. Math.)
+ arXiv:2606.09096 (June 2026). Sources extracted locally: `membrane_pulls/suzuki_2301.00421.txt`,
`membrane_pulls/suzuki_2606.09096.txt`. Purpose: pin down EXACTLY what is proven unconditionally, where RH
enters, and the precise open gap — so effort goes at the gap, not at re-demonstrating the machine (LEDGER §46b).

## A. UNCONDITIONAL (proven, no RH anywhere)
1. **The screw function g(t)** (2606 eq 1.3) = our A+B−P as one continuous kernel (verified §42). Near-origin
   expansion g = |t|log|t| + A|t| + ½Σ_{n≤e^{|t|}}Λ(n)/√n(|t|−log n) + O(t²), A = ½(log2π−ψ(2)) (2606 eq 2.2).
2. **Operator chain** A_a = Friedrichs ext of D*G_aD (2606 Thm 1.1); spectrum discrete, bounded below; λ_a
   CONTINUOUS in a (Thm 1.3); λ_a > 0 for small a with λ_a = log(1/a)+μ₁−log2π+ψ(2)−1+O(a) (Thm 1.4).
   [Our addition: μ₁ = M ≈ 0.233, ground state of the prime-free confined H^{1/2} functional — §42c.]
3. **Yoshida/Suzuki criterion**: RH ⟺ λ_a ≥ 0 for ALL a > 0. (Failure of RH ⟺ some λ_a < 0, by continuity.)
4. **Thm 1.5 (2606)**: the extensions D_{a,θ} of i d/dx on ℋ(T_a) have characteristic function W(a,θ;z) with
   ALL REAL zeros — proven unconditionally (needs only finitely many primes at fixed a).
5. **The family S_t(z)** ([14] eq 1.5-1.6): explicit, zero-free (poles + Σ_{n≤e^t}Λ(n)-part + Γ-part + Lerch
   tail), S_t ∈ L²(ℝ) for every t (Prop 1.2). The spaces H₀, K₀ built from {S_t} — unconditional ([14] §3.3).
6. **THE EQUALITY CRITERIA** (the sharpest unconditional statements — [14]):
   - Thm 1.4: RH ⟺ ‖P̂_{Dψ}‖²_{L²} = π⟨ψ,ψ⟩_W for all ψ ∈ Cc^∞ (countably many suffice).
   - Cor 1.5: RH ⟺ 2‖ψ‖²_{L²} = ⟨ψ,ψ⟩_W on the explicit subspace V°(0).
   - **Cor 4.3 (single-variable form): RH ⟺ (1/2π)‖S_t‖²_{L²(ℝ)} = −g(t) for all t ≥ t₀ (some t₀ ≥ 0).**
   These convert Weil's INEQUALITY criterion into IDENTITIES ("nonnegativity of the Weil form is explained by
   the nonnegativity of the L² norm").

## B. UNDER RH (where the hypothesis enters, and how)
1. **Thm 1.1 [14]**: H_W ≅ H(E_ξ) (de Branges) ≅ K(Θ_ξ) — the isomorphism REQUIRES RH even to define H_W
   (the Weil form must be ≥ 0 to be an inner product). Zeros Γ = eigenvalues of the self-adjoint extensions of
   multiplication — the Hilbert–Pólya realization — lives on top of this, hence under RH.
2. **Eq (1.2)**: lim c_a·v̂_a(z) = ξ(1/2+iz) — the a→∞ eigenfunction limit; formulated/expected under RH
   (§46/§46b: our numerics confirm its content at finite a, which is explicit-formula-forced).
3. **Cor 1.6 (2606)**: the θ(a)-tuned limit e^φ W(a,θ;z) → z²ξ/ξ′ ⟹ RH — conjectural bridge, direction
   "limit exists ⟹ RH"; proving the limit IS the open problem.
4. **Proof mechanics of the ⟹ direction of the equality criteria**: uses the zero-representation of S_t
   ((3.2)) + Parseval over Γ — i.e., under RH the L² norm picks up exactly the on-line zero masses.

## C. THE GAP, stated precisely
**RH ⟺ D(t) := (1/2π)‖S_t‖²_{L²(ℝ)} + g(t) ≡ 0 on some half-line [t₀,∞).**
- D is ONE real function of ONE real variable, built unconditionally from primes+Γ-factor data (no zeros).
- Mechanism: by the zero-representation, D(t) = the excess L² mass contributed by any OFF-LINE zeros (each
  off-line quadruple contributes a strictly positive, explicitly-shaped function of t). On-line zeros contribute
  exactly −g(t) (Parseval). So D ≥ 0 always, and D ≡ 0 ⟺ no off-line zeros.
- HONEST: this is the same global positivity, RELOCATED (the May-2026 warning stands — no reformulation is
  progress per se). What makes THIS relocation valuable:
  (a) **Dimension collapse**: from a quadratic form over an infinite-dim test class (∀g) to ONE function of ONE
      variable on a half-line. The quantifier structure changed (cf. meta-rule M2, SUPERPRISM_2026-07-20b).
  (b) Each hypothetical off-line zero contributes a D-term with a DISTINCT analytic signature in t (growth
      e^{(1/2−σ)t}-type + oscillation e^{iγt}) — vanishing of D on a half-line vs. the rigidity of such sums is
      a Dirichlet-series/almost-periodicity question, a DIFFERENT toolkit (uniqueness of almost-periodic
      functions, Turán-type) than the positivity toolkit that OBSTRUCTION_MAP exhausted.
  (c) It is the natural home of the equality-not-inequality insight: proofs of IDENTITIES can go through
      analytic continuation / uniqueness arguments unavailable for inequalities.
- Immediate executables (next):
  1. Transcribe S_t exactly from the arXiv LaTeX source (OCR too mangled) → build the unconditional D(t)
     machine (compute D(t) with rigorous quadrature; validates our reading; gives the defect object for
     analysis). [bounded, ~1-2h]
  2. Derive the per-off-line-zero contribution D_ρ(t) in closed form (from (3.2)) and the almost-periodicity
     rigidity question it poses: can Σ_ρ D_ρ(t) ≡ 0 on a half-line without being empty? If a Turán/uniqueness
     argument applies to THIS specific class, that would be the proof. [the analytic frontier]
  3. Cross-connection: D(t)-vanishing vs our §38 margin-maximizer landscape (both say "ℤ sits exactly at the
     degenerate point"); and Cor 4.3's t₀ vs the verified-height T (does partial D-vanishing on [t₀,T'] give
     finite theorems à la RH_RESULTS? — likely yes and possibly new). [scoping]

## D. What NOT to do (per §46b + May doc)
- No more finite-a demonstrations of eq (1.2)/Cor 1.6 content — explicit-formula-forced, confirmation only.
- No θ-tuned matching runs — superseded (§46) and doubly superseded by the identity criteria.
- No re-derivation of the equivalences in yet another basis (May doc: provably circular).
