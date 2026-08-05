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
