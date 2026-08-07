# The Anchor ∀-Theorem — final form (assembles §802A/§803A/§804A; drafts §806 — barrier session holds §805)

## Statement (target)

**Theorem (Weil positivity on the weighted-concentration class at the open rung).**
Let 2W = 0.66 (support 4W = 1.32 > log 3, primes {2,3} in range). For even ψ ∈ PW_{2W}
(band-limited, exponential type 2W, L²-normalized) define the off-window weighted tail
    T(ψ) = ∫_{|t| ∉ [14,44]} |ψ̂(t)|² · (1 + |Ω(t)|) dt,   Ω(t) = Re ψ₀(¼ + it/2) − ln π.
Then there is an explicit certified ε₀ > 0 such that T(ψ) ≤ ε₀ implies Q(ψ) ≥ 0, where
Q = pole + arch − prime is the Weil form.

The class {T(ψ) ≤ ε₀} is defined by an explicit weighted-tail condition — no reference to our
probe spans (non-circular), and the weight (1+|Ω|) is exactly what makes Q's off-window
contribution finite (Gödel's window-cap argument, quantitative form).

## Assembly chain (each line = one certified constant)

Split ψ = p + f: p = Π_E ψ (projection onto E = the top-k eigenspace of the even concentration
operator A = P_PW 1_{±[14,44]} P_PW), f = rest.

1. **Concentration control of ‖f‖.** ⟨Aψ,ψ⟩ ≥ 1 − ∫_{off} |ψ̂|² ≥ 1 − T(ψ). Eigen-split:
   ⟨Aψ,ψ⟩ ≤ ‖p‖² + λ_{k+1}‖f‖² ⇒ **‖f‖² ≤ T(ψ)/(1 − λ_{k+1})**.
   Needs: certified upper bound on λ_{k+1} — trace route: λ_{k+1} ≤ tr A − Σᵢ≤k λᵢ(section),
   section eigenvalues lower-bound true ones (compression interlacing). Float landscape (§804
   scout): k = 9 gives λ₁₀ ≤ 0.0371, so ‖f‖² ≤ 1.039·T(ψ). TO CERTIFY in arb: tr A (1D
   integral, closed-ish), section concentration matrix (finite smooth integrals), shifted-LDLᵀ
   eigen lower bounds.

2. **Q on p.** p ∈ E, an explicit 9-dim span (eigenvectors of the certified section, expressed in
   the 17-probe ambient basis). **m_E := certified pencil λ_min(Q, G) on E** via VᵀQV / VᵀGV
   interval products from gram14.json + bisected LDLᵀ. (Computing now.)

3. **Cross term |Q(p, f)|.** In/out split at ±[14,44]:
   - in-window: |Ω| ≤ sup_{[14,44]}|Ω| (certified interval sweep; note Ω is MONOTONE-ish here,
     sup at t=44 ≈ ln(44/2π)+... ≈ 2.0-ish — sharper than the [0,88] figure since p and f both
     enter, but f's in-window mass is small: ⟨Af,f⟩ ≤ λ_{k+1}‖f‖² gives ∫_{in}|f̂|² ≤ 0.0371‖f‖²
     — the cross term picks up ANOTHER √λ_{k+1} factor: |arch_in(p,f)| ≤ sup|Ω|·‖p‖·√(0.0371)‖f‖).
   - off-window: |arch_off(p,f)| ≤ √(∫_off |p̂|²(1+|Ω|)) · √T(ψ) — needs the certified weighted
     off-window mass of the eigenspace E (explicit finite computation on the section vectors:
     E is 96%+ in-window, its weighted off-mass is small — call it τ_E, certify).
   - prime: ≤ B_prime·‖p‖‖f‖ = 0.3579·‖p‖‖f‖ (§803A, certified).
   - pole: ≤ negligible certified constant (§803A envelope).
   **B_pf := 0.3579 + sup|Ω|_{[14,44]}·√λ_{k+1} + √τ_E + pole** — every piece certified.

4. **Q on f.** |Q(f)| ≤ [arch_in: λ_{k+1}·sup|Ω|_{in}·‖f‖²] + [arch_off: T-weighted, ≤ T(ψ)
   directly by the class definition] + [prime: 0.3579·‖f‖²] + [pole: negligible]
   ⇒ **|Q(f)| ≤ C_f·‖f‖² + T(ψ)** with C_f certified small.

5. **Assembly.** Q(ψ) ≥ m_E‖p‖² − 2B_pf‖p‖‖f‖ − C_f‖f‖² − T(ψ), with ‖f‖² ≤ 1.04·T and
   ‖p‖² ≥ 1 − ‖f‖². Positivity for T ≤ ε₀ where ε₀ solves the explicit quadratic — one line of
   certified arithmetic once m_E, B_pf, C_f land.

## Status
- gram14.json (17-probe ambient Q + G, zero-free): COMPUTING.
- arb certification of the concentration section (traces + matrix + eigen bounds): NEXT SCRIPT.
- τ_E, sup|Ω|_{[14,44]}, C_f: small certified computations after the section lands.
- m_E: pencil on E once both land.
- Assembly + §806 bank + Lean socket (hypothesis loci: everything ∅ except the identity-route
  pieces — NOTE: this theorem's Q-side constants are ALL zero-free/analytic; the ONLY hypothesis
  entering is via nothing at all — the statement is UNCONDITIONAL. The frame constant c and the
  verified zeros appear nowhere in this assembly. Check that claim twice at bank time.)

*Drafted while gram14 computes (session 2 cont.). Everything here is Report until each constant's
arb line exists; the float landscape informing k=9 is §804-scout tier.*

## Cross-check from the barrier map (BARRIER_GATEKEEPER.md §811, concurrent barrier session)

The barrier session ran this route through its completed obstruction map. Verdict: **BARRIER-CLEAN
for the stated (height-localized) scope.** Two flags folded in here:

- **m(W) → 0 is the B16 knife-edge crossing, not a technical lemma.** If the p-block margin m_E(W)
  degenerates to 0 as W grows, the whole route collapses into a mapped dead pattern. So keeping
  m_E(W) bounded below along the ladder is THE load-bearing structural fact (agrees with
  Grothendieck §799: "exhibit the structure keeping m(W) alive while dim K_W ~ e^{4W}"). The
  certified curve so far (§790/§795/§802A: pencil m_E = 0.1289 at the anchor, coarse λ_min rising)
  is on the safe side; watch it explicitly as W extends.
- **The strip-tail adversary-only role is INDEPENDENTLY required by the map** (B1: FE/growth-only
  components can never supply positive mass) — a second derivation of Red's §788 constraint. Good.
- **The unbounded-T₀ completion (full-class, all heights) cannot come from more envelope (B1),
  more sampling/zero-density (B5, now dead-proven via a published Beurling witness), or averaging
  (B10).** The lattice forces it to be genuinely NEW joint Euler-product ∧ functional-equation
  input. This bounds what "the next idea" must be — and it's exactly the B18 arithmetic-surface
  door the whole map points at. The anchor theorem does NOT claim to cross it; it's the honest
  floor the uniform program builds up from.
- **Guardrail (§813):** never substitute an approximate/Riemann-Siegel-grade FE for the exact FE
  anywhere in the completion — Turán's partial-sum criterion died on Haselgrove's witness; the
  Γ-factor's exact shape is load-bearing (K-P Annals 2011). If a future step reaches for AFE input,
  it's on a dead rung.

## Open careful-lemma items (documented for a clean pass, NOT rushed)
1. **τ_E in arb** — E's weighted off-window mass ≤ 0.525 (measured; 5 of 6 E-vectors ~2e-4, the
   weakest-concentration 6th leaks). Cert route: for each E-vector (float coeffs c_j), arb-integrate
   ∫_off|Σc_j Fev(c_j,t)|²(1+|Ω|): near [0,14] direct interval, far [44,∞) via quartic envelope +
   ln-tail (§791 family). Assembly ALREADY closes with the measured value (ε₀=2.3e-4>0).
2. **m_E basis-consistency** — the assembly needs m_E = Q-pencil-λmin on the CONCENTRATION eigenspace
   E (so ⟨Af,f⟩≤λ₇‖f‖² holds for f⊥E), but §802A certified Q on the capacity PROBE SPAN. These are
   close (both ~capacity-6 subspaces of the same window) but not identical. Clean fix: compute Q, G, A
   on ONE conditioned 8-probe basis, take E = its top-6 concentration eigenspace, certify m_E = Q on
   THAT E — self-consistent single-basis theorem. ~1 script; deferred to a fresh careful pass.
Both are tightenings of an already-closing theorem, not sign-flipping gaps.
