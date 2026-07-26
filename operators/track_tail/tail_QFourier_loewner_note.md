# Q-tail Fourier–Loewner structure — input note for tail-lemma v2 obligations 1–2
*(Hawk, 2026-07-24 roundtable. ADDITIVE to tail_lemma_v2_nearband_identity.md — does not supersede it.
Verification artifacts: loewner_cocycle_extract.py, results/loewner_cocycle_extract_out.txt.)*

## 1. The v2 identity collapses to a one-sequence form (exact algebra, no new analysis)

Product-to-sum on the v2 identity's two sines: sin((k+l)bt/2)·sin((k−l)bt/2) = ½[cos(lbt) − cos(kbt)].
Therefore, exactly:

    Q[k,l] = (2kl / (A(k²−l²))) · [ĥ(l) − ĥ(k)],     ĥ(m) := ∫₀^{2a} H′(t) cos(mbt) dt

i.e. the ENTIRE off-diagonal tail is driven by one scalar sequence — the half-range Fourier cosine
coefficients of H′ on [0,2a] (b·2a = π, so cos(mbt) is the standard basis). Equivalent Loewner form:
Q[k,l] = (2/A)·(k f(l) − l f(k))/(k²−l²) with f(m) = m·ĥ(m).

**Verified against the Q51 matrix (float, NL=4000):** R[k,l] := Q[k,l](k²−l²)/(kl) must satisfy the
cocycle R[k,l]+R[l,m]+R[m,k]=0 for every triple iff the form holds. Measured over all 2600 triples:
max violation 4.23e-14 absolute, 2.24e-13 relative. Exact to machine precision.

## 2. Measured ĥ behavior (extracted via φ(m) := R[1,m], φ = (2/A)(ĥ(m)−ĥ(1)))

- φ(m) → φ_∞ ≈ +0.0719 (window m=41..51 mean); write ψ(m) := φ_∞ − φ(m).
- m·ψ(m) is BOUNDED OSCILLATING in [−0.85, +0.56] over m=25..51 — no decay of m·ψ(m).
  So ĥ(m) = ĥ_∞-shifted + O(1/m) with an irreducible 1/m·oscillation. The oscillation frequency is
  b·log 2 ≈ 1.983 rad/unit-m: this is the Fourier signature of the H′ JUMP at t = log 2 — the prime
  itself is what forbids faster coefficient decay. Second IBP cannot remove it (H″ carries w_p·δ_{log2});
  a pointwise |Q| ≤ C/|k−l|^{1+ε} far-band bound is therefore FALSE, not merely hard — consistent with
  the measured binned envelope: max|Q|·|k−l| is flat ~0.9–1.0 across distances 2–26 (no extra decay).

## 3. Consequence: absolute (Gershgorin) row sums log-diverge; the classical fix is Hilbert's inequality

Substituting ψ into the Loewner form: Q[k,l] = kl(ψ(k)−ψ(l))/(k²−l²). Decompose ψ(m) = γ/m + β·osc(m)/m + O(1/m²):
- γ/m term ⟹ EXACTLY −γ/(k+l): a Cauchy (Hilbert-matrix) kernel. Rows non-summable (log-divergent),
  but operator norm ≤ π·γ (Hilbert's inequality, classical).
- β·osc/m term ⟹ modulated discrete Hilbert kernels 1/(k∓l) with sin/cos(mb·log2) diagonal modulation:
  pointwise ≤ β/|k−l| (matches the measured 1.15/d envelope) — again non-summable rows, but operator
  norm ≤ π·β-type constants (diagonal unitary modulation preserves the Hilbert-transform norm).
- O(1/m²) remainder ⟹ absolutely summable, plain Gershgorin.

So: for row k, the far tail behaves like (k·ψ(k))/l; since k·ψ(k) is bounded-oscillating and does NOT
vanish, Σ_l |Q[k,l]| over the infinite tail grows like k·ψ(k)·½·log — **pure absolute-row-sum Gershgorin
cannot certify the infinite tail, at any cutoff.** (Empirical foreshadowing: margins shrank as the window
widened 27→51 — that was the log.) The lemma must be restructured:

    λmin(T_tail) ≥ min_k d₀(k) − ‖E_near‖_row-sum − ‖E_far‖_op

near band by the v2 identity (obligation 1, unchanged); far field by explicit two-term Fourier asymptotics
of ĥ (jump at log 2 → exact β; endpoint/log-singularity → γ) + Hilbert/Cauchy operator-norm constants;
1/m² remainder by trivial summation. Every ingredient is classical with explicit constants.

## 4. Measured answer key (updates the v2 key)

- Cocycle/Loewner form: exact (2.24e-13).
- m·ψ(m) ∈ [−0.85, +0.56] for m=25..51 (bounded, oscillating at b·log2).
- Diagonal (obligation 3 targets): Q[k,k] for k=15..51 ranges 1.775 (k=19) to 3.448 (k=49), non-monotone;
  min over k≥25 is 2.143 (k=25). The 0.31·log(kβ) model under-predicts ~2×.
- Binned envelope: max|Q|·d ≈ 0.65–0.97 flat over d=2–26 (all pairs); no 1/d^{1+ε} regime in-window.
- K0=23 low block: Arb-CERTIFIED PD this session (results/arb_certify_k0_K23_nl50000.log, 12 pivots,
  worst mode-1 [4e-8 ± 4.35e-9]).

## 5. What this does NOT do

No claim about the composed theorem. Obligations restated, still OPEN: (1) near-band C₁/|k−l| with
explicit C₁ via the identity; (2') two-term ĥ asymptotics + Hilbert operator bound for the far field
(replaces the unattainable pointwise 1+ε target); (3) diagonal lower bound d₀(k) ≥ ~1.9 for k ≥ 25;
(4) interval-Schur composition against the certified K0=23 block.

## 6. UPDATE (same session, later turn): u_c and β have exact closed forms — verified to m=2001

Computing ĥ(m) EXACTLY via IBP on the kernel's closed forms (uc_closed_form_check.py — no quadrature;
H(0) = 0 exactly since the Lerch sum telescopes to Φ₀/4), the two-term law for odd m is:

    u(m) := (2m/A)·ĥ(m) = −1 − (4w_p/π)·sin(m·b·log 2) − 2/(π²m) + 1.720/m + O(1/m²)

- **u_c = −1 EXACT** (log-singularity coefficient ½ × commensurability: −(2/π)Si(mπ) → −1; all
  sin(mπ) boundary terms vanish). **β = 4w_p/π = 0.624052 EXACT** (jump IBP). Verified: residual×m
  = +1.720 constant to FOUR DIGITS over m = 25..2001 (results/uc_closed_form_check_out.txt).
  The prefactor identity Q[k,l] = (2kl/(A(k²−l²)))(ĥ(l)−ĥ(k)) re-verified against the raw kernel
  at rel. err ~3e-12.
- Supersedes the fitted u_c = −0.907 / β = 0.683 (obl2prime lane): those were finite-window (m ≤ 51)
  fit bias from the +1.72/m term (mean +1.72/m over [15,51] ≈ +0.06..0.09 — exactly the 9% gap).
  Also corrects §2 of THIS note: the reported m·ψ range [−0.85, +0.56] was contaminated by using the
  window tail-mean for φ_∞ (linear-in-m bias); the true law gives m·ψ(m) = 1 + 0.624·sin(θm) ∈ [0.38, 1.62].
- **Margin consequence:** the Cauchy-channel Hilbert bound is |u_c|·π/2 = 1.5708 (not 1.425), so
  "d₀ − |u_c|π/2" = 2.057 − 1.571 = 0.486 (not 0.632) — IF the β-channel contributes nothing, which
  is measured-plausible (combined norm converges from below through 1.02 at M=12801) but UNPROVEN;
  crude triangle charges it π/2·β ≈ 0.980 and then the k=25 corner FAILS (2.551 > 2.057). The
  knife-edge is real.
- **De-risking structure (proposed obligation 4′):** d₀(k) = log(kb) − 1.8527 ± 0.3612 grows without
  bound while every norm constant is fixed. So: certify modes [25, M*] as a finite Arb block
  (mechanical, same class as tonight's K0=23 cert), and apply the CRUDE triangle norm bound only to
  modes > M*, coupling handled by the three-term test λ_A·(d_C − ‖Off_C‖) ≥ ‖B‖² with a rectangle
  Schur bound on B. Working point M* = 801: d_C ≥ 5.52 worst-phase, d_C − 2.551 = 2.97,
  λ_A·2.97 ≈ 6.1 > ‖B‖² ≈ 3.7. Closes with ONLY classical constants — no interference lemma needed;
  any future sharp bound on the β-channel just shrinks M* (less compute), never gates feasibility.
