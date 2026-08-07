# Barrier Gatekeeper — route audit protocol + audits (LEDGER §811)

Purpose: turn BARRIER_MAP.md from a museum into a gate. Any proof route gets audited: list the
ζ-properties it consumes, check each barrier row, name where the open completion points sit. A route
that consumes only properties from a witness-saturated lattice row is dead on arrival; a clean route
must show it consumes EP ∧ FE jointly, globally, individually.

## Protocol (run on any route before investing)
1. List the properties consumed: which parts of {Euler product, FE/Γ-factor, pole, counting
   precision, verified zeros ≤ T₀, positivity structure, family averaging} does each component use?
2. Check the lattice: does any COMPONENT claim standalone force? If a component uses a
   witness-saturated property set alone, it may only play a bounded/defensive role.
3. Check B8–B12: global (not per-prime)? class-level (not single test function)? individual (not
   averaged)? non-circular (what exactly is proven unconditionally vs equivalent-to-RH)?
4. Name the completion points: where does the route's finite/restricted scope end, and which
   barrier waits there? Any claim of the FULL theorem must name the new property input that crosses it.

## Audit 1 — Riemann8 height-localized rung architecture (REDUCTION.md, §788–§806, 2026-08-07)

Route: Q(ψ) = pole + arch − prime on Fejér-probe classes (support 4W, centers ≤ C_max), corner /
frame-window / strip-tail decomposition, p ⊕ r composition with pencil-normalized joint certificates,
verified zeros ≤ T₀.

**Property consumption — the fatal test PASSES:** Q is the Guinand–Weil form: pole (EP/ζ-pole), arch
(FE's Γ-factor), prime (EP). The route consumes EP ∧ FE jointly — it lives in the lattice's one
non-witness-saturated locus. Not disqualified.

**Row-by-row:**
- B4/B8 (local/pointwise): PASS — matrix certificates over spans, not pointwise comparisons.
- B9 (single test function): PASS — a certified λ_min over a CLASS is exactly the repair of the
  single-h identity trap.
- B10 (averaging): PASS — ζ's own verified zeros, no family average anywhere.
- B13 (sharp extremal witness): PASS — pencil normalization works with the ~23-dim near-null instead
  of fighting it; no magic function sought.
- B1-as-component: CONFIRMING INSTANCE — the strip tail (§791, cosh⁴ envelope) uses only FE/growth
  data, and the architecture correctly assigns it a purely DEFENSIVE role ("can never supply positive
  mass," Red's constraint). The map says FE/growth-only components can't force zeros on-line; the
  route already obeys. Independent agreement between the map and Red's §788 attack.
- B5-as-component: PASS with note — the frame inequality is L²/sampling information, but it is fed
  by ζ's actual certified zeros (individual data), not generic mean values, and composes with the
  corner + tail rather than claiming standalone force.
- B16 (knife-edge vs effective bounds): THE CRUX, correctly named — m(W) bounded below while
  dim K_W ~ e^{4W} is exactly where a knife-edge could reappear. §795's measured rising margins are
  the evidence it doesn't; the uniform-in-W statement is the open analytic work. The map adds: if
  m(W) → 0, the route degenerates into the B16 ceiling — so the m(W) lower bound is not a technical
  lemma, it IS the barrier crossing.
- B12 (circularity) + verification height (OBSTRUCTION_MAP #6): the two completion points, correctly
  acknowledged in the architecture itself (Red's §788: frame inequality FALSE on any class of dim > #zeros; "needs
  something genuinely new" for unbounded height). The map sharpens what "genuinely new" must mean:
  crossing T₀-finiteness cannot come from more envelope (B1-saturated), more sampling (B5-saturated),
  or more averaging (B10) — it must be a NEW joint EP∧FE consumption (e.g. the FE forcing zero
  supply above T₀ the way the integrality class forces rigidity, or B18's arithmetic surface
  supplying the polarization).

**Verdict: BARRIER-CLEAN for its stated scope** (height-localized rungs, unconditional). The full
theorem's two completion points (uniform-in-W margin, unbounded-T₀ supply) sit exactly on B16 and
B12/#6, and the lattice constrains any crossing to new joint EP∧FE input. No component overdraws its
property budget; the architecture's own honesty devices (Red's constraints, adversary-only tail)
match the map independently.

---
*Add an audit section per route. A route that fails the fatal test gets one line and a citation to
the lattice row that kills it.*
