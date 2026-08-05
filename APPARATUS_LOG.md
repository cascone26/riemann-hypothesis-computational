# RH Apparatus — Build Log (append-only)

## 2026-08-06 — build kickoff (plan: APPARATUS_BUILD_PLAN.md)
- ENV oriented: Lean4 (v4.18/v4.29) + RiemannLean project EXISTS with a ZERO-SORRY Handle
  (TailLemma.tail_kernel) + partial Robin lane (4 sorries); cvxpy+scs+sympy/mpmath ready;
  Ollama (no models); 8GB/8core/15GB-free (heavy work -> Olympus).
- ORGAN 1a (verifier restore): mathlib cache restoring (`lake exe cache get`, running).
- ORGAN 2 (SDP/Weil, self-certifying, parallel):
  - v1: WRONG (min-eig -1.41) — used product φ_j(log n)φ_k(log n) for prime term. My built-in
    "RH predicts PSD" check caught it. [Demonstrates the apparatus thesis: my Weil-form error
    rate is the binding constraint the verifier removes.]
  - v2: VALIDATED the explicit-formula identity (Σ over zeros = arch+pole−prime; terms ~0.83
    cancel to 1e-4) — arch/pole/prime/signs correct.
  - v3: CORRECT positivity matrix (cross-correlation prime term). Q PSD (min-eig ≈ −3e-13).
    Eigenvalues [≈0,≈0,≈0,0.033]: three NEAR-NULL — the SATURATION ("barely true") appearing
    in the SDP lane; the Gaussian-cos basis at a=0,4,8 is near-extremal (W→0). Q saved (weil_Q_v3.npy).
  - NEXT: (a) determine if the near-null eigenvalues are exactly 0 (rank deficiency) or >0
    (saturation from above) via exact/higher-precision Gram; (b) exact-rational PSD certificate
    on a full-rank sub-basis (first self-certifying Handle); (c) hierarchy: grow basis, track the
    margin decay = shape constraint on any certificate.

## 2026-08-05 — ORGAN 2 FIRST HANDLE (self-certifying, no Lean)
- Margin basis found: phi_a, a∈{15,25,35,45} → Weil form min-eig=+0.0481032, grid-CONVERGED
  to 6 digits (stable rmax 60→120, pmax 100k→300k). ("wide 20,40,60,80" min-eig=−0.228 was a
  grid-TRUNCATION artifact — healed to +1.105 with wider r-grid; internal "RH⇒PSD" check flagged
  it correctly, NOT an RH violation.)
- RIGOROUS certificate (Arb enclosure prec=300 → exact-rational LDL^T): all 4 pivots positive
  rationals, Sylvester minors all +, ρ≤3.45e-24 ≪ 0.048 margin. HANDLE minted (LEDGER §731,
  proof pointer weil_cert_15_25_35_45.json). Machinery proven end-to-end.
- Honest: clear-margin nearly-diagonal case. NEXT = hierarchy toward coupled/saturated bases.
- ORGAN 1a: mathlib cache re-downloading (slow, background, ~5%). Non-blocking.

## 2026-08-05 — ORGAN 1 OPERATIONAL (the multiplier)
- 1a: mathlib restored; RiemannLean.TailLemma builds clean (8252 jobs, exit 0) → zero-sorry
  Handle tail_kernel RE-VERIFIES. Verifier live.
- 1b/1c: verify_loop.py — propose→verify→mint. Gate VERIFIED (verify-before-ship on the gate):
  true lemma PASS+recorded, `sorry` REJECTED, FALSE stmt FAIL. First loop-Handle:
  two_sub_two_cos_nonneg. handles.jsonl live.
- PLAN CORRECTION: the 4 Robin sorries = full Robin/Lagarias↔RH equivalences (not tractable);
  organ 1 = many small Handles via the loop + toolkit formalization, not one-shot RH.
- BOTH organs 1 (live) + 2 (§731/§732) operational. LEDGER §733.
