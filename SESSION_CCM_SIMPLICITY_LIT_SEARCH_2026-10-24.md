# CCM Exact Operator Simplicity — Literature Search, 2026-10-24

**Goal:** Find a genuinely NEW angle for analytically proving simplicity of the smallest eigenvalue of tau_N for all N in the Connes-Consani-Moscovici exact operator framework.

**Prior search angles (NOT retried):**
1. "simple ground state trace class perturbation large"
2. "non-degenerate eigenvalue random matrix perturbation additive"
3. "Toeplitz compact operator simple spectrum eigenvalue multiplicity"
4. "generic simplicity symmetric operator perturbation transversality"

---

## Angle A: Citation Search (Direct & Recent)

**Search queries:**
- "arXiv:2511.22755 cited by Connes Consani Moscovici exact operator"
- "Connes Consani Moscovici 2025 2026 follow-up simplicity eigenvalue"

**Findings:**

### Direct Citations Found
- **arXiv:2607.02828** — "A finite Guinand–Weil dictionary and archimedean tail order for the truncated Weil quadratic form" (Akiva Groskin, 2026-07)
  - Cites and extends CCM's truncated Weil form Q_W^N_λ
  - Establishes two-sided certification rules with explicit budget B_T ~ (2N+1)ρ log(T)/(π²T)
  - Provides finite-cutoff to cutoff-free positivity certification
  - **On simplicity**: The paper ASSUMES simplicity of the smallest eigenvalue and works forward from that assumption. Does NOT offer a proof method.

- **GitHub implementations** (akivag613/connes-cvs-, nrensen/connes-cvs-): Active numerical/computational projects; no new analytical simplicity proofs in codebases.

### Related Recent Work
- Connes's own blog post (2026-07) on Zeta Spectral Triples restates the CCM construction but does not sketch a proof of the simplicity assumption.
- "Weil positivity in compact windows: a finite reduction, certified two-sided bounds, and a Landau-Widom decay law" (arXiv:2608.24827) — focuses on positivity certification, not eigenvalue multiplicity.

### Assessment
**Not directly applicable.** Follow-up work acknowledges the simplicity assumption as an open structural requirement, but none of the citations published since 2511.22755 (late 2025) propose or prove a simplicity theorem for this operator. The problem remains open at the state-of-the-art.

---

## Angle B: Perron-Frobenius / Krein-Rutman Route

**Search queries:**
- "Krein-Rutman theorem simple eigenvalue positive operator"
- "Perron-Frobenius archimedean operator number theory"

**Findings:**

### Krein-Rutman Theorem
**Canonical statement** (Generalized Krein-Rutman, arXiv:1606.04377):
- For a **strongly positive** (or positivity-improving) **compact linear operator** T on a Banach space with a cone, the spectral radius ρ(T) is:
  1. A **positive real number**
  2. An **algebraically simple eigenvalue** (multiplicity exactly 1)
  3. Associated with a **strictly positive eigenvector**
  4. The unique eigenvalue of maximal absolute value

### Key Strength
**This directly guarantees simplicity automatically** — no perturbation theory, no transversality argument, no generic assumption. If the extremal eigenvalue is the spectral radius of a strongly positive operator, multiplicity = 1 is automatic.

### Applicability Check: Can τ or a Transformation Be Shown to Be Strongly Positive?

The "exact operator" τ is a real symmetric matrix of dimension 2N+1 with:
- Archimedean (digamma integral) base term
- Rank-finite prime/pole corrections (rank ≤ const as N → ∞)

**Obstacles to direct application:**
1. **Sign structure unclear**: The problem statement does not confirm entrywise non-negativity or positivity-improving property of τ itself. If τ has mixed signs (likely, given it's a perturbation of a spectral-triple operator), Krein-Rutman does not apply directly.
2. **Transformations possible but unexamined**: One might construct e^{-tτ} for suitable t, or -τ if τ's spectrum structure allows it, to make the relevant operator positive. This has NOT been explored in prior sessions or the literature.
3. **Compact operator requirement**: τ_N is finite-dimensional (dimension 2N+1), hence all operators on it are "compact" in the finite-rank sense, so this condition is trivially met.

### Perron-Frobenius Connections
Search for "Perron-Frobenius archimedean" yielded no specific results connecting the Perron-Frobenius theorem to "archimedean" operators in number theory. The term "archimedean" in the Connes-Consani context refers to the archimedean place (real/complex completions) in arithmetic geometry, not to ordered Banach spaces or Perron-Frobenius matrices. **Not directly applicable**.

### Assessment
**HIGH POTENTIAL IF POSITIVITY CAN BE ESTABLISHED.** Krein-Rutman is a genuinely powerful tool that, if applicable, would solve the simplicity problem completely for all N. The open question is: can τ_N (or a canonical transformation of it) be proven to be a positive/positivity-improving operator? This has NOT been attempted in prior work and represents a novel angle.

**Recommendation:** Investigate whether the digamma-plus-prime-correction structure of τ admits a positivity-improving interpretation or transformation. If yes, Krein-Rutman closes the problem.

---

## Angle C: Sturm-Liouville / Oscillation Theory Route

**Search queries:**
- "Sturm oscillation ground state simplicity integral operator"

**Findings:**

### Sturm Oscillation Theorem
**Classical statement**: For a Sturm-Liouville problem on [a, b], the nth eigenfunction has exactly n-1 zeros (nodes) in the interior of [a, b].

**Corollary for ground state (n=1)**: The principal eigenfunction (ground state) has **zero nodes** — it is sign-definite (all positive or all negative).

### Connection to Multiplicity
In some contexts (e.g., Schrödinger operators on ℝ), the fact that the ground state is sign-definite can be used to argue that its eigenvalue is simple. However:
1. This typically requires additional structure (e.g., the operator is self-adjoint and semi-bounded below).
2. The argument is qualitative ("if the ground state has no nodes and is strictly positive, then it is unique up to scaling").
3. **Quantitative multiplicity one requires further argument**, not just oscillation.

### Integral Operator Setting
The Sturm-Liouville framework is naturally differential-equation based (d²/dx²). For integral operators with digamma/kernel structure:
- Classical oscillation theory gives the nodal structure of eigenfunctions.
- Whether this immediately implies multiplicity = 1 is context-dependent and is NOT stated as a general theorem in the search results.

### Assessment
**Adjacent but not directly conclusive.** Oscillation theory is relevant to understanding the structure of the ground state (likely sign-definite due to nodal count = 0), but the search did not return a theorem of the form "if ground state has no nodes, eigenvalue is simple" for integral operators. This remains an open path to explore but is not a known standalone result.

---

## Angle D: Numerical Range / Birman-Schwinger Route

**Search queries:**
- "Birman-Schwinger simple eigenvalue trace class perturbation"
- "multiplicity one smallest eigenvalue Schrödinger operator integral kernel"

**Findings:**

### Birman-Schwinger Principle
**Statement**: Let H₀ be a self-adjoint operator and V a rank-k perturbation (k finite). Then z₀ is an eigenvalue of H = H₀ + V if and only if 1 is an eigenvalue of the Birman-Schwinger operator K(z₀), with **geometric multiplicities matched** (and algebraic multiplicities matched under trace-class conditions).

### Application to τ
The "exact operator" τ = base + rank-finite prime correction fits this framework:
- H₀ = archimedean base (spectral triple operator)
- V = prime/pole corrections (rank ≤ const)

**However:** The Birman-Schwinger principle tells you *how eigenvalues arise under perturbation* — i.e., given H₀ with spectrum S, the spectrum of H = H₀ + V can be tracked via the Birman-Schwinger operator. **It does not directly prove that the result of the perturbation has simple spectrum.**

### Multiplicity via Integral Kernels
For **positive integral kernels** K(x, y), the Perron-Frobenius theorem guarantees the largest eigenvalue is simple. If τ_N's action can be interpreted as a positive integral operator (or a power/transform of one), this would apply. **But the sign structure is not confirmed** in the CCM paper or follow-ups.

### Assessment
**Applicable to analyzing spectrum propagation, but not a direct simplicity proof.** Birman-Schwinger is a tool for tracking how the unperturbed spectrum changes under perturbation, not for proving the perturbed spectrum is non-degenerate. Useful for understanding the mechanism of spectral change, but does not address the core open problem directly.

---

## Summary: State of the Literature

| Angle | Finding | Applicability | Status |
|-------|---------|----------------|--------|
| **A: Citations** | arXiv:2607.02828 + GitHub implementations; all ASSUME simplicity, none PROVE it | Not applicable — problem still open | Settled |
| **B: Krein-Rutman** | Powerful theorem: strongly positive operator → automatic simple extremal eigenvalue | HIGH IF positivity of τ or transformation can be shown | **Open, high-value path** |
| **C: Sturm-Liouville** | Oscillation theory gives nodal count of ground state; alone does not imply multiplicity = 1 | Adjacent; would need supplementary argument | Inconclusive |
| **D: Birman-Schwinger** | Framework for tracking spectral perturbation; does not directly prove simplicity of result | Mechanism-understanding, not simplicity proof | Not applicable |

---

## Conclusion

**Genuinely new and applicable?** 

**Angle B (Krein-Rutman) is the most promising novel angle.** If the exact operator τ_N (or a canonical transformation of it, e.g., e^{-tτ_N} with the right sign, or -τ_N reframed as a largest-eigenvalue problem) can be proven to be a **strongly positive or positivity-improving operator**, then Krein-Rutman's theorem **automatically guarantees simplicity of the smallest (or largest) eigenvalue for all N**, bypassing the need for perturbation transversality or other generic arguments.

**The open question is structural:** Does the digamma-plus-prime-correction operator admit a positivity interpretation?

Angles A, C, D do not offer new applicable theorems; A shows the problem is still open in the literature; C is adjacent theory requiring supplementary work; D is mechanism-understanding, not a direct simplicity proof.

**Recommendation for next session:** Investigate the sign structure and positivity-improving properties of τ_N (or explicit transformations of it). If even a candidate positive operator can be identified, Krein-Rutman becomes the path to a full analytic proof.

---

## Sources

- [Zeta Spectral Triples (arXiv:2511.22755)](https://arxiv.org/pdf/2511.22755)
- [A finite Guinand–Weil dictionary (arXiv:2607.02828)](https://arxiv.org/pdf/2607.02828)
- [A Generalized Krein-Rutman Theorem (arXiv:1606.04377)](https://arxiv.org/pdf/1606.04377)
- [Weil positivity in compact windows (arXiv:2608.24827)](https://arxiv.org/pdf/2608.24827)
- [High-Precision Approximation of Riemann Zeros (arXiv:2605.20224)](https://arxiv.org/pdf/2605.20224)
- [Simplicity of extremal eigenvalues of the Klein-Gordon equation (arXiv:1008.2469)](https://arxiv.org/pdf/1008.2469)
- [The Birman–Schwinger principle applications (arXiv:2507.02251)](https://arxiv.org/pdf/2507.02251)
- [Eigenvalue Properties of Schrödinger Operators](https://projecteuclid.org/ebooks/advanced-studies-in-pure-mathematics/Spectral-and-Scattering-Theory-and-Applications/chapter/Eigenvalue-Properties-of-Schr%C3%B6dinger-Operators/10.2969/aspm/02310027.pdf)
- [akivag613/connes-cvs (GitHub)](https://github.com/akivag613/connes-cvs-)
- [nrensen/connes-cvs (GitHub)](https://github.com/nrensen/connes-cvs-)
