# Superprism — the next lane after the Krein picture (2026-08-01, Direction 4)

The session built a complete Krein/de Branges picture: the |ζ|²-string (zeros in coefficients, S(T) in
phases, Möbius in the phase-twist), the zero-drum (primes in coefficients, §545/556), Suzuki's screw
(Weil form), all linked by Krein duality (§531), with the E-function showing our string ≠ Suzuki's space
(§555). Deep constraint from §500: RH is global-cancellation-only — no local/pointwise certificate exists.
Everything measured is δ²-soft ("barely true"). Where does this picture point next?

## Five ranked directions

**1. (TOP) The Weil form in the string basis — concrete positivity.**
Suzuki's Q_W (= our A+B−P) is positive ⟺ RH, but his realization is abstract. We now have an EXPLICIT
orthonormal basis: the |ζ|²-string's polynomials p_n. Express Q_W as a matrix in {p_n} and study its
positivity/structure. Unlike the abstract de Branges space, this is a computable infinite matrix whose
positivity is RH. Pairs the two things the session made explicit (the Weil form via screw §532, the string
basis §555). CONCRETE, and connects our machinery to the actual RH-equivalent.

**2. The Krein intertwiner — the map between the two de Branges spaces.**
§555: our quadrature space ≠ Suzuki's ξ-space, but they are Krein-dual and share the zeros. The explicit
INTERTWINING operator (which sends one canonical system to the other) is the bridge from the computable side
to the RH-equivalent side. Construct it from the two E-functions.

**3. Hilbert–Pólya forcing — what pins σ=½.**
§552/553: the exponent reads Re ρ; the reconstruction is best at ½, but nothing FORCES it — Suzuki assumes
RH. The gap: what extra structure makes the string coefficients real/positive in the way that pins zeros to
the line? The a=1 razor (§474) is the saturated form of this. Candidate: the reality of the Jacobi
coefficients {a_n=0, b_n>0} of the |ζ|² measure IS automatic; the question is whether the DUAL (drum)
coefficients' reality forces the zeros real. A reality/self-adjointness argument, not positivity.

**4. The zeros↔primes fixed point.**
Zeros determine primes (drum coefficients, §545) and primes determine zeros (explicit formula / string).
This is a fixed-point/self-consistency structure. RH may be a uniqueness/stability statement of that fixed
point. Concrete scout: iterate zeros → drum-coeffs → primes → explicit-formula → zeros′ and study whether
the true (on-line) zeros are the attracting fixed point and whether off-line configurations are unstable.

**5. Quantitative "how barely true" — the softness law → Λ_dBN.**
The δ²-stiffness grows as τ^1.6 (§525), the Möbius gap as 0.009·ln n (§519b), detection decays with height
(§474c). These quantify "barely true." Connect the rate to the de Bruijn–Newman constant Λ=0 boundary
(Rodgers–Tao) — a quantitative RH-margin law from the operator softness.

## Pick: Direction 1 (Weil form in string basis).
Most concrete, most directly RH-equivalent, reuses the session's two explicit assets. Scout below.

# Superprism round 3 (2026-08-03, post-§584, Fable) — after the crown
BOARD: the phase-arithmetic principle is complete (full Selberg class, observed+verified+derived, theorem-
candidate with one analytic gap). Krein dictionary, drum, Weil saturation all closed. Everything δ²-soft;
§500 metatheorem stands: RH content lives in global/asymptotic structure only.

KEY REALIZATION: the phase field is a FUNCTION of the positions (arg L′(ρ)=πS_L(γ)−π/2, and S = counting
fluctuation = determined by the zero POSITIONS). So the crown result says: the zero positions, through
their OWN counting function, encode the inverse-Euler arithmetic that the SAME positions encode via
Landau/explicit-formula. That is a SELF-CONSISTENCY (bootstrap) constraint on the zero configuration:
  positions ⟹ S-field ⟹ exp(iπS) twist ⟹ must reproduce Π_p(Euler)⁻¹ ⟸ primes ⟸ positions (Landau).
RH-adjacent question: is the true (on-line) configuration the unique/stiff solution of this bootstrap?

## Ranked directions
1. (TOP) BOOTSTRAP STIFFNESS: perturb ONE zero on-line by ε (no pair-symmetry protection!) and measure the
   degradation of position⟺phase self-consistency (twist-from-primes vs actual-phase agreement, comb
   quality). If LINEAR in ε — the bootstrap is a stiff variational characterization of the true zero set
   (contrast: every off-line channel is δ²-soft). Cheap, decisive, genuinely new angle.
2. COMB-AS-RH-WITNESS: off-line pair (½±δ+iγ) in the readout: predict deficit ∝ δ²log²n (pair-symmetry
   soft but with growing coefficient) — measure the coefficient, complete the softness map for the new channel.
3. ENVELOPE BOUND (the theorem gap): analytic, Gonek–Ng discrete moments — the path to a real theorem. Multi-session.
4. DEGREE-3 (sym²Δ): further universality; incremental; needs careful Γ-factors.
5. Writeup/consolidation of MOBIUS_PHASE_RESULT toward shareable form (Jacob's call, "no papers" rule).
PICK: 1 then 2, bank, then 3 as the standing frontier.
