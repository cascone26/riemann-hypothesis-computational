# arXiv Submission Notes

## Target: math.NT (Number Theory)

## Cross-list categories
- math-ph (Mathematical Physics) — for the Hilbert-Polya/operator sections
- cs.LO (Logic in Computer Science) — for the Lean 4 formalization

## MSC 2020 Primary
- 11M26 — Nonreal zeros of zeta and L-functions; Riemann-type hypotheses

## MSC 2020 Secondary
- 11M06 — Analytic continuation; summation methods
- 11Y35 — Analytic computations and problems
- 03B35 — Mechanization of proofs and logical operations

## Keywords
Riemann hypothesis, Li criterion, Li coefficients, pair positivity, functional pair,
Hilbert-Polya operator, Nyman-Beurling criterion, zero-free regions,
Lean 4 formalization, Robin's inequality, Lagarias inequality, WKB quantization

## Submission Checklist
- [x] main.tex compiles cleanly (last compile: prior to paper/main.pdf exists)
- [x] Bibliography inline (no separate .bib file needed)
- [x] No external figures (no .eps/.pdf/.png files needed)
- [x] Repository URL updated: https://github.com/cascone26/riemann-hypothesis-computational
- [x] MSC codes added as comment in preamble
- [ ] Author email/affiliation — Jacob may want to add: Grand Canyon University
- [ ] arXiv submission form: https://arxiv.org/submit
  - File: arxiv_submission/main.tex
  - Primary: math.NT
  - Cross: math-ph, cs.LO (optional)
  - Comments field: "11 pages of prose + 4 tables + 1 algorithm. Lean 4 source
    at the repository. No figures."

## Title
"Computational Investigations of the Riemann Hypothesis:
Parameterized Operator Search, Spectral Geometry, and Formal Verification"

## Author
Jacob Cascone

## Key claims to double-check before submission
1. λ_n > 0 for n ≤ 85 — Theorem 3 (unconditional) ✓
2. λ_n > 0 for n ≤ 7.1×10^6 using Odlyzko 2M — Theorem 4 (conditional on LMFDB verification) ✓
3. 145 certifiable runs, 13,881 right endpoints in scan to 10^7 ✓
4. d_5000 = 0.0723, NB distance monotone decreasing, full rank ✓
5. WKB optimal α ≈ 0.9, 2.6× better than BK ✓
6. Connes ε_N ~ 0.606·N^{-0.226} ✓

## Known limitations
- Lean proofs have `sorry` (statements formalized, not proofs)
- Corollary covers n ≤ 2824 formally; beyond is computational certificate
- No claim of proving RH; finite-n partial results only
