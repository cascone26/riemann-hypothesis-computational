# The razor edge, unified: our a=1 is the Weil-form face of Λ = 0 ("RH is barely true")

**Session 2026-07-30, Track 1 of the pivot brainstorm.** Consolidates §470–§474c into one statement and
places it honestly in the literature.

## The finding (ours, this session)

Across five independent computational/structural lanes, the same zero-margin criticality:

| # | Lane | Object | Result | Ledger |
|---|------|--------|--------|--------|
| 1 | KLMN form bound | a = sup ⟨ψ,𝒱ψ⟩/⟨ψ,L₀ψ⟩ | a = 1.00000 exactly (not <1, not >1) | §470 |
| 2 | Weil quadratic form | min-eig(arch+pole−prime) | ≈ 0 to machine precision (PSD-singular) | §470 |
| 3 | Keiper–Li | λ_n | > 0 Arb-certified n≤150, but a criterion with no slack (λ_n>0 ⟺ RH, boundary case) | §472 |
| 4 | Dirichlet decomposition | 𝒱 = 2C_N·I − E_jump | positivity = exact balance of two positive forms vs a divergent constant | §473 |
| 5 | Off-line detector | min-eig vs δ = Re ρ − ½ | PSD at δ=0 exactly; indefinite ~ −δ² off-line | §474b |

Plus the witness: the saturating direction ψ* is the **zero-node function** (nodes exactly at the γ's, §474a),
localized (PR → 6.3, §473b), with detection stiffness decaying steeply in height (§474c).

**One sentence:** the positivity that is equivalent to RH holds with *exactly zero margin*, saturated by the
zero-interpolating test function, and fails quadratically the instant a zero leaves the line.

## The honest literature placement (lit-gate applied)

This razor-edge phenomenon is **not novel to us**. It is a theorem of Rodgers–Tao (2018): the de
Bruijn–Newman constant satisfies **Λ ≥ 0**, hence RH ⟺ Λ = 0 — their own gloss: *"the Riemann Hypothesis,
if true, is only barely so."* The heat-flow encoding (H_t with all-real zeros iff t ≥ Λ) sits at its critical
point exactly at RH.

What our session adds is **the explicit map of this criticality into the Weil/operator encoding**, with the
critical direction exhibited concretely:

- Λ = 0 boundary (heat flow) ⟷ a = 1 boundary (KLMN/Weil form): both say "the RH-equivalent positivity has
  zero slack."
- Rodgers–Tao's mechanism (zeros attract under backward heat flow; any Λ < 0 would relax to GUE and
  contradict) ⟷ our §474c stiffness decay (the restoring force against off-line displacement vanishes at
  height ∞ — no uniform spectral gap protects the line).
- The Lehmer-pair phenomenon (near-degenerate zero pairs forcing Λ ≥ −ε pre-Rodgers–Tao) ⟷ our localized
  saturator ψ* (near-null directions of the Weil form supported near close zero pairs).

So: **independent confirmation of the field's sharpest known "razor edge" theorem, in a different encoding,
with the saturating witness made explicit.** Claim level: characterization/confirmation (Report), not new
mathematics. The five-lane table is a useful map; the equivalences Weil ⟺ Li are Bombieri–Lagarias 1999;
a=1-as-identity is the explicit formula itself.

## Why this matters strategically (the methodology theorem)

Five re-encodings, one criticality ⇒ **any positivity re-encoding of RH terminates at a zero-margin edge,
because the positivity IS the theorem.** There is no slack anywhere to exploit, in any encoding — Λ ≥ 0 says
this is not an artifact of our instruments but a proven property of the problem. Therefore the filter for
future attacks:

> A viable approach needs a source of rigidity EXTERNAL to the positivity being proved
> (reciprocity, functional equation structure, geometry, a theorem ladder with partial credit).

The reduced statement this line leaves behind (sharpest form of the open core in our encoding):

    E_jump(ψ) + ⟨ψ, L₀ψ⟩ ≥ 2 κ_crit C_N ‖ψ‖²   in the renormalized N→∞ limit
    (positive prime jump form + confining prolate ⊵ divergent PNT constant) ⟺ a ≤ 1 ⟺ RH.

## Where the machine points next (Track 2 candidates, feasibility-scouted separately)

1. **Semilocal support ladder** — Weil positivity for restricted-support test functions; archimedean rung
   proven (Connes–Consani); first prime rung(s) the live increment. Our §473 jump-form + §469 semilocal
   instruments apply directly.
2. **Vasyunin / Bettin–Conrey** — the Nyman–Beurling Gram matrix has closed-form arithmetic entries with a
   *reciprocity law* = external rigidity of exactly the kind the filter demands.
3. **Λ upper bound (Polymath15 lane)** — the one lane with a genuine theorem ladder (Λ ≤ 0.2 → lower);
   real unconditional partial progress, heavy engineering, our Arb toolkit applies.

*(Scout reports banked in the ledger as they land; this doc is the stable synthesis.)*
