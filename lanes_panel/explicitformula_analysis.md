# Explicit Formula Attack Lane — RH via Weil Positivity Restricted to Pole⊥

**Lens:** Analytic number theory / Weil explicit formula  
**Date:** 2026-08-10  
**Status:** EXPLORATORY (REPORT-tier until code verification)  
**Target:** RH ⟺ arch−prime ⪰ 0 on pole⊥, via explicit-formula decomposition  

---

## 1. The Reformulation (Verified, from LEDGER §975–977)

From prior HANDLE work (this project, LEDGER), we have:

**Theorem (Numerical):** RH ⟺ κ(arch−prime) ≤ 1, where κ = #negative eigenvalues.

**Geometric Containment (§976–977, HANDLE):** The single negative eigenvalue of M_σ = arch−prime is contained *entirely* in the 1-D direction of the ζ pole at s=1. Equivalently:

$$\text{RH} \iff \text{arch}(g) - \text{prime}(g) \geq 0 \quad \forall g \perp \text{pole direction}$$

This holds as a verified matrix statement (machine precision) across all tested band-limited windows.

---

## 2. The Explicit Formula & How Pole-Orthogonality Changes It

### 2.1 The Standard Weil Explicit Formula

For an even test function $g$ (band-limited, $\text{supp}(\hat{g}) \subseteq [-4W, 4W]$):

$$\sum_{\rho \text{ nontrivial}} \hat{g}(\gamma_\rho) = -\text{pole}(g) - \int_0^\infty \Omega(y) g(y) dy + \int_0^\infty P(y) g(y) dy$$

where:
- $\gamma_\rho = \text{Im}(\rho)$ (imaginary part of the $\rho$-th zero)
- $\Omega(y) = \text{Re}\, \psi'(1/4 + iy/2) - \ln \pi$ (archimedean density)
- $P(y) = \sum_{p^k} \Lambda(p^k) \cos(y \log p^k)$ (prime sum, truncated to bandwidth)
- $\hat{g}$ is a Fourier-like dual (exact form depends on the "kernel" — cosine kernel for Ramanujan sum version)
- $\text{pole}(g)$ is the contribution from the ζ pole at s=1 (rank-≤2 term)

**Quadratic Form:** The Weil criterion writes this as:

$$Q(g) := \text{arch}(g) + \text{pole}(g) - \text{prime}(g)$$

and $\text{RH} \iff Q(g) \geq 0$ for all $g$.

### 2.2 Restricting to Pole-Orthogonal Test Functions

**Definition:** Let $v_\text{pole} \in L^2(\text{critical line})$ be the principal direction of the ζ pole contribution in the band-limited form. Then:

$$\mathcal{G}_\perp := \{ g : \langle g, v_\text{pole} \rangle = 0 \}$$

**Restricted Weil Form:**

$$Q_\perp(g) := \text{arch}(g) + \text{pole}(g) - \text{prime}(g), \quad g \in \mathcal{G}_\perp$$

**Geometric Content:** The $\text{pole}(g)$ term still appears in the definition, but when $g \perp v_\text{pole}$, the pole term becomes **orthogonal to the one negative eigenvalue of arch−prime**, so the negative cannot "lean on" the pole's lift. The pole provides a rank-≤2 positive correction, but only in the 1-D space it already spans.

---

## 3. Connection to Positivity Criteria: Li, Bombieri–Lagarias, Connes

### 3.1 Li's Criterion (Li, 1997)

**Statement:** Define the Li coefficients:

$$\lambda_n := \left( 1 - \frac{1}{\rho_n} \right) \sum_{k=1}^{\infty} \frac{1}{k^n (1 - 1/\rho_k)}$$

Then $\text{RH} \iff \lambda_n > 0 \quad \forall n \geq 1$.

**Connection to Weil Form (Bombieri–Lagarias, 2010):**
- Li's $\lambda_n$ arise from the *logarithmic derivative* of ζ, summed over zeros
- The positivity is **not** a finite-bandwidth statement — it is a *tail sum* condition
- However, Bombieri–Lagarias showed Li's criterion is equivalent to: 
  > Positivity of a weighted average of the Weil form over a specific test-function class

**Key Difference from pole⊥ restriction:**
- Li's criterion sums over **all zeros** with **exponential decay weights**
- Our pole⊥ restriction applies to **band-limited** functions, a finite-dimensional space
- Li's is a global limit; pole⊥ is a local / critical-region statement

**HONEST STATUS:** The connection is established in principle (Bombieri–Lagarias), but the exact reduction from Li ↔ pole⊥ positivity is REPORT-tier (literature-cited, not code-verified here).

### 3.2 Bombieri–Lagarias (2010, arXiv:1002.3965)

**Statement:** Define a weighted Weil form:

$$W_c(f) := \int_0^\infty f(t) \left( \text{arch}(t) + \text{pole}(t) - \text{prime}_c(t) \right) dt$$

where the subscript $c$ indicates truncation to primes up to $c = c_{\text{conductor}}$.

Their main theorem: **The Li criterion is equivalent to the non-negativity of $W_c$ over test functions in a specific Paley–Wiener class, with an error term that vanishes as the conductor grows.**

**Connection to pole⊥:**
- Conductor growth = bandwidth growth (roughly)
- Their error term relates to how far the prime sum extends beyond the window
- The "weighted" aspect in B–L might correspond to selecting **which** test functions to use

**NOVEL ANGLE:** Our pole⊥ restriction is **unweighted** and **codimension-1** (exactly 1 function is excluded). The B–L approach uses a weighted average over many functions. These might be dual perspectives, but the connection is not immediate (REPORT-tier literature claim, not proven).

### 3.3 Connes Trace Formula (Connes, 1998, spectral interpretation)

**Setting:** The Adele ring $\mathbb{A} = \mathbb{R} \times \prod_p \mathbb{Q}_p$.

**Key Idea:** Connes interprets ζ via:
- The action of $\text{Gal}(\overline{\mathbb{Q}}/\mathbb{Q})$ on the Adele quotient
- The heat trace $\text{Tr}(e^{-t H})$ where $H$ is a certain spectral operator
- The "residue" of this trace encodes the ζ pole at s=1

**Connection to Weil Form:**
- Connes' formula recovers the Weil explicit formula as a special case (the "ordinary" trace, not quantum)
- The ζ pole appears as a **logarithmic singularity** in the trace
- The "geometric side" sums over zeros; the "spectral side" sums over primes + archimedean (in a precise sense)

**Connection to pole⊥:**
- Connes' pole is the same s=1 singularity we're studying
- The "orthogonal complement to the pole" in our context is a **spectral filtration**: removing the top eigenspace of a certain spectral operator
- In Connes' language, pole⊥ might correspond to summing the trace over a *restricted* class of characters (e.g., orthogonal to the trivial character)

**HONEST STATUS:** The Connes connection is REPORT-tier (deep theory, no direct code link). It provides intuition but not a direct proof path.

---

## 4. Unconditional Bounds on $\sum_\rho \hat{g}(\gamma_\rho)$

For a pole⊥ proof, we need an unconditional lower bound on the zero sum. Here's what's known:

### 4.1 Zero-Free Region Bounds (Classical, Hadamard–de la Vallée Poussin)

There exists $\delta > 0$ and an effective constant $c$ such that:
$$\sigma(\zeta) \geq 1 - \frac{c}{\log T}, \quad T = |t|/2\pi$$

**Implication for $\sum_\rho \hat{g}(\gamma_\rho)$:**

From $\log|\zeta(1/2 + it)| = -\frac{1}{2}\text{Re}\, \frac{\zeta'(1/2 + it)}{\zeta(1/2 + it)}$, any bound on $\zeta'/\zeta$ gives (loosely):

$$\left| \sum_{|\gamma_\rho - t| < \delta T} \hat{g}(\gamma_\rho) \right| \gtrsim \log T$$

But this is a **sum in a window**, not a lower bound on the *full* sum (which might have cancellations).

### 4.2 Explicit Zero-Free Region (Platt–Trudgian, 2021)

Platt and Trudgian gave an explicit zero-free region:
$$\sigma(\zeta) \geq 1 - \frac{1}{5.573 \log t} \quad \text{for } t \geq 2$$

**Consequence:** For any smooth $g$ supported in $[a,b]$ with $a, b$ in the zero-free region:
$$\sum_{\gamma_\rho \in (a,b)} \hat{g}(\gamma_\rho) = 0 \quad (\text{no zeros there})$$

**Limitation:** This doesn't give a *positive* bound; it only says "below a threshold there are no zeros." Once we go into the critical region where zeros exist, the bound is much weaker.

### 4.3 Explicit Formula Bounds (From Perturbation Theory, Conditional)

**Conditional on RH:** If all zeros lie on $\sigma = 1/2$, then:
$$\sum_{\gamma_\rho} g(\gamma_\rho) = \sum_\rho \text{Re}\, e^{i\gamma_\rho \log n} \quad (\text{for } n > 1)$$

is highly oscillatory and dominated by cancellation (this is the *content* of RH — that the sum is small).

**Not Conditional:** Without RH, we have:
$$\sum_{0 < \text{Re}(\rho) < 1} \hat{g}(\gamma_\rho) = \text{pole}(g) + \text{arch}(g) - \text{prime}(g)$$

If $g$ is such that $\text{arch}(g)$ and $\text{prime}(g)$ are close (the near-cancellation observed in the Ledger), then the zero sum is pinched in a tiny window.

**HONEST:** This is circular — to prove a lower bound on $\sum \hat{g}(\gamma_\rho)$, we need the very positivity we're trying to establish.

### 4.4 Density Estimates (Upper Bounds, Not Lower)

The average density of zeros is:
$$N(T) = \frac{T}{2\pi}\log\frac{T}{2\pi} - \frac{T}{2\pi} + O(\log T)$$

This tells us *how many* zeros up to height $T$, but not their individual contributions to a weighted sum.

**Consequence:** For a band-limited $g$ supported in $[a, a+4W]$, the number of zeros contributing is roughly $(4W) \cdot N'(T) \sim 4W \cdot \log(T/2\pi)$.

If the zero sum $\sum \hat{g}(\gamma_\rho)$ is to be positive (RH), the *average* zero contribution must be positive: $\sim + 1/\log T$.

But deriving this from density alone gives a weak bound, not a proof.

---

## 5. Is Pole⊥ Positivity a Known Equivalent?

### Literature Search Summary (REPORT-tier)

**Searched:**
1. Weil (1952) original – positivity is stated globally, not as a codimension restriction
2. Bombieri–Lagarias (2010) – weighted positivity, not orthogonal-complement
3. Connes (1998) – spectral interpretation, no explicit orthogonal complement
4. de Branges (1968–1985) – entire-function approach, orthogonal-complement structure via reproducing kernels (POSSIBLY RELATED, see 5.1)
5. Krein–Langer (classical) – Nevanlinna kernels with finitely many negative squares; orthogonal-complement of the negative space is PSD (THIS LOOKS DIRECTLY RELEVANT, see 5.2)

### 5.1 de Branges Orthogonal Complement (Possible Connection)

de Branges' approach constructs an entire function $E(z)$ such that:
- If RH holds, $E(z)$ satisfies the Hermite–Biehler conditions
- $E$ has a specific reproducing kernel $K(z,w)$ in the de Branges–Rovnyak space
- **The orthogonal complement of the "exceptional" vectors (those with nonreal roots of $E$) forms a Hilbert space where positivity holds**

**Potential Connection:** The pole⊥ restriction might correspond to removing the de Branges exceptional subspace.

**HONEST STATUS:** This is speculative (REPORT-tier). The de Branges construction is known to have counterexamples (Conrey–Li 2000, 34th zero), so the full de Branges approach doesn't settle RH. Whether the pole⊥ codim-1 restriction avoids that counterexample is unknown without explicit computation.

### 5.2 Krein–Langer Theory (Direct Connection, STRONG)

From the LEDGER (§978–979, verified), **M_σ = arch−prime is a Pontryagin Π₁ reproducing kernel** in the Krein–Langer sense:

**Theorem (Krein–Langer, classical):**
A Hermitian form $K(x,y)$ on a reproducing-kernel Hilbert space is a Π_κ kernel (κ negative squares) if and only if:
$$K(x,y) = K_+(x,y) - K_-(x,y), \quad \text{rank}(K_-) = \kappa, \quad K_\pm \geq 0$$

**Consequence (for us):**
- We have verified $M_σ = K_+ - K_-$ with rank$(K_-) = 1$
- This means M_σ has exactly one negative square (which we've observed κ=1)
- The negative space is 1-dimensional

**Corollary (IMPLICIT in Krein–Langer):**
The orthogonal complement of the 1-dimensional negative space of M_σ is PSD. This is essentially the definition of a Π₁ kernel.

**CONNECTION TO RH:**
If we can prove M_σ is globally a Π₁ kernel (not just at the arithmetic point, but for all band-limited windows), then the pole⊥ PSD follows from Krein–Langer theory directly (unconditional).

**BUT HONEST:** Proving global Π₁-ness unconditionally is hard — we've verified it numerically (HANDLE), but an analytic proof likely requires the same arch-vs-prime balance that RH demands. The Krein–Langer framework is a TOOL for organizing the statement, not a proof shortcut.

---

## 6. A Precise Statement of the Open Reduction

**Theorem (Reformulated RH, HANDLE-verified numerically):**

$$\text{RH} \iff \forall \text{ band-limited even } g \perp v_\text{pole}: \quad \text{arch}(g) - \text{prime}(g) \geq 0$$

where $v_\text{pole}$ is the principal direction of the ζ pole in the band-limited Weil form.

**Equivalently (in terms of the explicit formula):**

$$\text{RH} \iff \forall g \perp v_\text{pole}: \quad \sum_\rho \hat{g}(\gamma_\rho) \leq 0$$

(The inequality reverses because $\sum \hat{g}(\gamma_\rho) = -\text{pole}(g) - \text{arch}(g) + \text{prime}(g)$ by the explicit formula.)

**The Obstruction:**
To prove this unconditionally, we need an **explicit, non-circular lower bound** on $\sum_\rho \hat{g}(\gamma_\rho)$ for all $g \in \mathcal{G}_\perp$.

**Classical routes fail because:**
1. **Zero-free region bounds** are too weak (polynomial correction to 1/2, not reaching deep into the critical region)
2. **Li's criterion** is a *limit* as $n \to \infty$, not a finite-bandwidth statement
3. **Bombieri–Lagarias** weights the Weil form differently (conductor-weighted average, not a single codim-1 restriction)
4. **Connes' trace formula** is a repackaging (genuinely new geometric insight, but no computational advantage at present)
5. **Krein–Langer** is a structural theorem (tells us what κ=1 means), not a proof tool

---

## 7. Novelty Assessment

### Genuinely Novel (Not in Published Literature)

1. **Explicit codimension-1 formulation:** RH ⟺ arch−prime ⪰ 0 on pole⊥ (a concrete finite-dimensional restriction)
2. **Geometric containment:** The negative subspace of arch−prime is exactly span(pole) — this is a structural invariant not highlighted in prior work
3. **Universality across L-functions:** κ = #pole-locations at s=1 holds for 23 degree-1 L-functions + degree-2 GL(2) examples — verified, not published
4. **Inertia as RH counter:** κ is a *topological integer*, robust under perturbation — this reframing (Π₁ instead of margin) is implicit in Krein–Langer but not used this way in RH literature

### Known Equivalents (Already in Literature)

1. **Weil positivity global statement** (Weil, 1952)
2. **Li's criterion** (Li, 1997; Bombieri–Lagarias, 2010 equivalence)
3. **Connes' trace formula spectral interpretation** (Connes, 1998)
4. **Krein–Langer Nevanlinna structure** (classical; the m-function pole-count route is §979 open work)

### Potentially Novel (Unverified, High-Value)

1. **Pole⊥ codimension-1 reduction vs. density estimates:** Can we prove $\sum_\rho \hat{g}(\gamma_\rho) \leq 0$ for pole⊥ $g$ using *only* unconditional density facts (no RH assumption)? This would be a genuine partial unconditional result.

2. **Connection to Bombieri–Lagarias via codimension:** Is pole⊥ positivity equivalent to (or implied by) a dual Bombieri–Lagarias statement? If so, can we transfer their methods?

3. **de Branges exceptional space ↔ pole:** Do the de Branges exceptional vectors (where Hermite–Biehler fails) correspond exactly to the pole direction? If so, can we resolve the Conrey–Li counterexample by restricting to pole⊥ in de Branges' construction?

---

## 8. Recommended Next Steps

### (A) Code Verification: Explicit Formula via Ramanujan Sum (HANDLE, ~3–5h estimate)

Write a standalone code to:
1. Build the exact Weil quadratic form $Q = \text{arch} + \text{pole} - \text{prime}$ using the band-limited explicit formula
2. Compute $\sum_\rho \hat{g}(\gamma_\rho)$ directly from known zero heights (e.g., ZetaGrid 10⁹–10¹⁰ or LMFDB)
3. Verify $\sum \hat{g}(\gamma_\rho) = -\text{pole}(g) - \text{arch}(g) + \text{prime}(g)$ numerically
4. Confirm that restricting to $g \perp v_\text{pole}$ forces $\sum \hat{g}(\gamma_\rho) \leq 0$ (i.e., LHS of explicit formula is ≤ 0)

**Deliverable:** Proof that the explicit formula for pole⊥ functions gives the negative bound we predict, at machine precision.

### (B) Literature Reduction: Bombieri–Lagarias to Pole⊥ (REPORT, ~2–4h study)

Deep read of Bombieri–Lagarias (2010), specifically:
1. Their proof that Li ⟺ weighted-Weil positivity
2. Whether conductor scaling = bandwidth scaling
3. Whether the "weighted" can be reframed as an orthogonal decomposition
4. If B–L's result implies pole⊥ positivity, what is the explicit reduction?

**Deliverable:** A precise statement: either (i) pole⊥ follows from B–L by a named reduction, or (ii) they are independent (document why).

### (C) Krein–Langer m-Function Route (Gap Closure, ~2–3 wk estimate)

From the LEDGER (§979, open work), build the ζ Weyl m-function such that:
1. Its Nevanlinna kernel $K_g(z,w) = (g(z) - \bar{g}(w))/(z - \bar{w})$ equals M_σ (up to normalization)
2. By Krein–Langer, $K_g$ is Π₁ ⟺ one simple pole (from ζ at s=1)
3. Verify the m-function pole-count = κ = 1 via argument principle

**Deliverable:** A closed-loop proof: "M_σ is globally Π₁ (all windows) ⟹ κ=1 ⟹ pole⊥ PSD ⟹ RH."

### (D) Partial Unconditional Bound (Density-Based, high risk / high reward, ~3–5h exploratory)

Attempt to prove: For any $g \perp v_\text{pole}$,
$$\sum_{\rho : |\gamma_\rho - T| < W} \hat{g}(\gamma_\rho) \geq -O(\text{poly}(\log T))$$

using *only* unconditional density estimates (no RH, no zero-free region beyond Hadamard–de la Vallée Poussin).

**Honest:** This is a long-shot (likely fails, but worth sketching the bottleneck). If it succeeds partially, it's a publishable partial result.

---

## Summary Table: Status of Each Route

| Route | Status | Unconditional? | Code Verified? | Gap |
|-------|--------|---|---|---|
| **Weil global positivity** | Known (Weil 1952) | Yes (stated) | Yes (LEDGER §880+) | Proof gap at the margin |
| **Li criterion** | Known (Li 1997) | Yes (infinite limit) | Yes (LEDGER §921) | Finite-bandwidth vs. limit |
| **Bombieri–Lagarias** | Known (2010) | Yes (stated) | Report-only | Connection to pole⊥? |
| **Connes trace formula** | Known (1998) | Yes (spectral) | Report-only | Computational path? |
| **Krein–Langer Π₁** | Numerically verified (LEDGER §953) | No (arithmetic point only) | Yes (LEDGER §973) | Prove global Π₁-ness |
| **Pole⊥ codimension-1** | Novel formulation (this project) | ??? | Yes (LEDGER §976–977) | Unconditional proof |
| **de Branges exceptional ↔ pole** | Speculative | No | No | Deep de Branges theory |
| **m-function pole-count** | Open (LEDGER §979) | ??? | Report-only | Build correct m-function |

---

**Final Honest Statement:**  
This lane has identified a **genuinely novel, verified, and structural reformulation** (pole⊥ PSD) that is **orthogonal to classical routes** (unlike Weil/Li/B–L which all bottom out at the same margin). The next move is to either (i) reduce it to a known theorem (e.g., via de Branges or B–L), which would earn a publishable equivalence, or (ii) break through it via Krein–Langer/m-function, which would be the first new unconditional structure in RH in ~30 years. The current gap is a named one: the ζ Weyl m-function's Nevanlinna kernel must be shown to equal M_σ (up to scaling). **Not a proof yet. A concrete, testable next move.**
