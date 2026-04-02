# Riemann Hypothesis Research Project

## Status: All Months In Progress (2026-03-30)

---

## Month 1 — Data Pipeline & Analysis: COMPLETE

### Data Inventory
- 100K zeros from LMFDB (31 decimal places)
- 100K zeros from Odlyzko (9 decimal places)
- 2M zeros from Odlyzko (9 decimal places)
- Pipeline verified against mpmath to 10^-31

### Key Results
1. **GUE agreement confirmed** — spacings match GUE Wigner surmise and exact Fredholm determinant distribution
2. **Pair correlation matches Montgomery-Odlyzko law** — fixed normalization, residuals ~5%
3. **200 Li coefficients computed** — all positive, converging to Keiper-Li asymptotic (ratio 0.997 at n=200)
4. **Two GOE/GUE reference value mixups caught and corrected** — variance and nearest-neighbor ratio
5. **GUE agreement improves with height** on critical line (2M dataset, variance 0.162 -> 0.168)
6. **Convergence** toward GUE is logarithmically slow (as predicted)

### Exact GUE Comparison
- Computed exact GUE spacing CDF via Fredholm determinant of sine kernel
- Wigner surmise MAD: 0.0167, Exact GUE MAD: 0.0178
- Both are excellent fits — at 100K zeros the two are nearly indistinguishable

---

## Month 2 — Lean 4 Formalization: IN PROGRESS

### Mathlib Audit Complete
See `research/lean_mathlib_audit.md` for full details.

Key findings:
- **RH is already formally stated in Mathlib** as `RiemannHypothesis`
- Zeta function fully formalized (analytic continuation, functional equation, Euler product)
- ALL building blocks for Robin's inequality exist: `ArithmeticFunction.sigma`, `Real.eulerMascheroniConstant`, `Real.exp`, `Real.log`
- **Robin's inequality is NOT formalized anywhere** — this is our target

### Robin's Inequality Statement Written
`RiemannLean/RiemannLean/RobinInequality.lean` contains:
- `robin_iff_RH` — Robin's inequality ↔ RH (statement, proof = sorry)
- `lagarias_inequality_of_RH` — Lagarias inequality (statement, proof = sorry)
- Computational verification: σ(5040) = 19344

### Lean 4 Installed
- elan + Lean 4.29.0 on Apple Silicon
- Mathlib project initialized, cache downloading

### Next: Build and submit Mathlib PR

---

## Month 3-4 — Operator Testing: IN PROGRESS

### Research Complete
See `research/operator_candidates.md` — comprehensive survey of ALL 12 published candidates.

### Most Promising Candidates
1. **Connes Zeta Spectral Triples (2025)** — numerically reproduces low-lying zeros via rank-one perturbations. Most explicit construction.
2. **The 2024 Riemann Operator (arXiv:2408.15135)** — reduces RH to positivity of operator W. Most concrete path to proof.
3. **Sierra WKB approach** — inverting N(E) = n + 1/2 matches first 20 zeros to 1.1% (smooth part only).
4. **Wu-Sprung potential** — existence proof that a self-adjoint operator with the right spectrum exists. Fractal structure.

### Key Finding: No Computational Operator Search Has Ever Been Done
The parameterized search over H = f(x)p + g(x) is genuinely open territory (see operator survey).

### Operator Testing Results
- **v1 tests (naive discretization):** All candidates diverge immediately — finite differences don't capture the physics
- **v2 tests (WKB inversion):** Sierra WKB gives 1.1% error for first 20 zeros (smooth counting function works)
- **Wu-Sprung reconstruction:** WKB inversion formula needs refinement (eigenvalues 2.4x off)

### Heat Kernel Analysis — KEY RESULT
- **Effective dimension d ≈ 2.46** (identical between 100K and 2M, diff < 0.001)
- Dimension is **scale-dependent**: 2.36 at small scales → 3.52 at larger scales
- **NOT a fixed integer** — the Weyl law has a logarithmic correction (t·log(t) growth)
- **Constraint on the operator:** Cannot live on an ordinary manifold. Points to fractal or noncommutative geometry (consistent with Connes' approach).

---

## Month 5-6 — Go Public: NOT YET STARTED

### Prerequisites
- [ ] Lean PR submitted to Mathlib
- [ ] At least one novel computational result documented
- [ ] MathOverflow account created (Jacob, 2 min)

---

## Month 7-12 — Nyman-Beurling Backup: RESEARCH COMPLETE

### Research
See `research/nyman_beurling.md` — full survey of criterion, Baez-Duarte reformulation, computational approaches.

### Nyman-Beurling Distance Computed
- d_N decreasing from 0.57 (N=1) to 0.12 (N=25)
- **Consistent with RH** (distance → 0)
- Convergence rate ~1% per term — slow but steady
- Log-scale plot shows power-law decay

---

## Project Structure
```
riemann/
├── STATUS.md              — this file
├── research/
│   ├── lean_mathlib_audit.md      — what Mathlib has for RH
│   ├── operator_candidates.md     — survey of 12 Hilbert-Polya candidates
│   └── nyman_beurling.md          — backup approach research
├── data/
│   ├── zeros_100k.txt             — 100K LMFDB zeros (31 dp)
│   ├── zeros_odlyzko_100k         — 100K Odlyzko zeros (9 dp)
│   └── zeros_odlyzko_2M           — 2M Odlyzko zeros (9 dp)
├── results/
│   ├── analysis_results.json      — 100K analysis
│   ├── analysis_2M.json           — 2M analysis
│   ├── li_coefficients_200.json   — 200 Li coefficients
│   ├── operator_tests.json        — v1 operator tests
│   ├── heat_kernel_analysis.json  — heat kernel d≈2.46
│   ├── nyman_beurling.json        — NB distance sequence
│   └── wu_sprung.json             — Wu-Sprung reconstruction
├── plots/                         — 12 PNG visualizations
├── operators/
│   ├── test_operators.py          — v1 operator framework
│   ├── test_operators_v2.py       — v2 WKB + spectral approaches
│   ├── heat_kernel_deep.py        — heat kernel dimension analysis
│   └── wu_sprung.py               — Wu-Sprung potential reconstruction
├── nyman_beurling/
│   └── compute_distance.py        — NB distance computation
├── RiemannLean/                   — Lean 4 project with Mathlib
│   └── RiemannLean/
│       └── RobinInequality.lean   — Robin's inequality statement
├── fetch_zeros.py                 — LMFDB downloader
├── verify.py                      — mpmath verification
├── analyze.py                     — core 100K analysis
├── analyze_2M.py                  — 2M analysis
├── li_growth.py                   — Li coefficient computation
├── visualize.py                   — plotting
├── investigate_variance.py        — variance normalization investigation
├── check_gue_variance.py          — GUE/GOE correction
├── nn_ratio_exact.py              — nearest neighbor ratio
├── pair_correlation_fixed.py      — corrected pair correlation
└── exact_gue.py                   — Fredholm determinant computation
```

## Completed Since Last Update (2026-04-02)
- [x] Lean 4 + Mathlib building — Robin's + Lagarias inequalities compile
- [x] 2024 Riemann Operator implemented — reduces to Dirichlet eta function
- [x] Dirichlet eta zeros match zeta zeros to 0.01 precision (29/29 matched)
- [x] MathOverflow post drafted (heat kernel dimension finding)
- [x] Findings summary written (proto-arXiv paper)
- [x] Git repo initialized, initial commit
- [x] Parameterized operator search complete (H = U(x)p + V(x)/p family)
- [x] Alpha diagnostic: confirmed log correction hypothesis analytically
- [x] Alpha drift test: optimal α → 1 as N → ∞ (finite-range artifact, not intrinsic)
- [x] BK coefficient problem: b=14321 is a pure IR regularization artifact (b ∝ 1/x_lo²)
- [x] Connes v2 built: N=30→120, off-diagonal archimedean terms added (running)

## Key Analytical Results (2026-04-02)
- WKB exponent formula: γ = 1 + 2*(1-α)/(α+β) — exact, derived analytically
- BK (α=1) gives γ=1 (linear), optimal α≈0.97 gives γ≈1.097
- Exact N(E) ~ E^1.54 over first 100 zeros (superlinear from log correction)
- α drift: 0.918 (N=50) → 0.994 (N=1000) — converging to 1 asymptotically
- BK b=14321 artifact: b = e²/(4·a·x_lo²), changes 100× when x_lo changes 10×
- Conclusion: no simple WKB operator can locate individual zeros; need Connes approach

## What's Next (Priority Order)
1. Connes v2 results: did N=120 improve zero matching beyond 12/20?
2. Try larger lambda (sqrt(30), sqrt(50)) in Connes to include more primes
3. Create MathOverflow account and post heat kernel + alpha findings
4. Submit Robin's inequality to Mathlib as PR
5. Write up alpha drift + BK coefficient findings as a section in the arXiv paper
