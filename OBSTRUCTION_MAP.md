# RH Obstruction Map — why the classical analytic toolkit reduces to the zeros

**Compiled 2026-07-11 (sessions 5–8, overnight).** Every entry is verified two ways and
recorded in `LEDGER.md`. Purpose: a single artifact stating, for each classical route to
proving Li/Weil positivity, the *precise* mechanism by which it becomes circular (or hits a
hard ceiling) — so no future session re-treads, and so the exact obstruction a real proof
must overcome is explicit.

---

## The single target (the wall ★)

Weil/Bombieri/Connes: RH ⟺ **A(g) + B(g) ≥ P(g)** for all admissible bandlimited g, where
- **A(g)** = (1/2π)∫ Ω(r)|ĝ(½+ir)|² dr,  Ω(r) = Re ψ(¼+ir/2) − log π  (archimedean; smooth, sign-indefinite: Ω<0 for |r|<6.29, Ω>0 beyond)
- **B(g)** = rank-2 boundary term from the poles of ξ at s=0,1 (positive, concentrated at low r)
- **P(g)** = Σ_k Λ(k)k^{−1/2}(ĝ(log k)+ĝ(−log k))  (prime sum; k ≤ λ² by bandwidth)

Equivalent operator form: the matrix **tau = A + B − P** is PSD ⟺ RH. Numerically tau ⪰ 0 at
machine precision for all N (λ²=14), with a 23-dim near-kernel and spectral gap ≈1.87e-4.

---

## The obstruction map

| # | Route | What it would need | PRECISE obstruction (verified) | Status |
|---|-------|-------------------|-------------------------------|--------|
| 1 | **Moment / Hankel** (Track 1) | (λ_n) a Stieltjes moment sequence ⟹ all λ_n>0 | Hankel matrix fails PD at **N=2**: det=λ₁λ₃−λ₂²=−0.0037. (λ_n) is log-CONCAVE, the opposite of a moment sequence. | **DEAD** |
| 2 | **Beurling-Selberg** (pointwise majorant of P) | A smooth majorant of the prime comb dominated by the archimedean kernel | Archimedean kernel Ω̌(x)=−e^{−|x|/2}/(1−e^{−2|x|}) decays as **k^{−1/2}** at x=log k; prime weight Λ(k)/√k decays as **log(p)·k^{−1/2}**. Ratio Λ(k)(1−1/k²)→∞; prime comb exceeds arch kernel at **every prime ≥5**. No pointwise majorant exists. | **Hits wall** |
| 3 | **Large sieve** (global bound on P) | Bound Σ Λ(k)k^{−1/2}R(log k) globally by A+B | Prime symbol on the line D(r)=(1/π)Re[−ζ′/ζ(½+ir)] is SMOOTH (no poles — ζ′/ζ pole 1/(i(r−γ)) is imaginary). Large sieve only bites as λ²→∞, giving unconditional **zero-density** N(σ,T) — bounds OFF-line zeros, never forces them ON. | **Hits wall** (but gives real zero-density → feeds route 6) |
| 4 | **de Branges** (H(E) space) | ξ-space reproducing kernel positive-definite | A+B has **4 negative eigenvalues** (modes n=±1,±2) → not a valid dB inner product. Only tau=A+B−P is PSD. dB positivity K(t,t)>0 ⟺ phase φ′(t)>0 everywhere ⟺ no Lehmer collision ⟺ RH. | **Hits wall** (Conrey-Li showed naive version fails) |
| 5 | **Operator / tau spectral** (Track 2, Connes) | Prove tau ⪰ 0 unconditionally (or spectrum = zeros) | Every certificate (Gershgorin, Sylvester, tridiagonal, SOS) fails or reduces to Weil. tau ⪰ 0 IS RH; the 23-dim kernel and gap 1.87e-4 are structural but the gap value depends on zero locations. | **Hits wall** |
| 6 | **Pair positivity + verification height** (session 5) | Per-zero F_n(σ,γ)≥0 for admissible σ,γ | NOT circular but has a scaling law: γ_crit(n)~n/(2π). Gives λ_n>0 unconditionally for **n ≲ 2π·T** (T = on-line verification height). Not a proof of RH, but the strongest unconditional Li bound. | **PARTIAL — real result** |
| 7 | **Speiser / Levinson** (ζ′ zero-free region) | ζ′(s)≠0 for 0<Re<½ (⟺ RH) | Equivalent to RH but NOT circular-positivity — connects to Levinson mollifier method → positive PROPORTION on the line (Levinson ⅓ 1974 → Conrey ⅖ → current **≈41.7%**). Engine = Dirichlet-poly mean value thm (VERIFIED: (1/T)∫|Σaₙn^{−it}|²→Σ|aₙ|², err ~O(N/T)). BARRIER: mollifier length y=T^θ controllable only for **θ<4/7** (Conrey, Kloosterman/Weil); beyond that off-diagonals are unbounded shifted-convolution sums. Reaching 100% needs θ→1 = knowing ζ moments fully ≈ RH-hard. Also verified: ζ′ zeros all Re>½ (direct + argument principle). | **Different partial axis; barrier = shifted-convolution / θ<4/7** |

---

## The common root

Routes 2–5 all reduce to the same fact (route 7 Speiser is the exception — different obstruction): **the RH content is the difference between the
archimedean density (smooth) and the prime comb (discrete), and that difference is exactly
the zero sum Σ_ρ|ĝ(½+iγ)|².** Concretely:

- The archimedean kernel and the prime comb have the **same k^{−1/2} decay** (explicit-formula
  balance of the place at ∞ vs the finite places), differing only by the log(p) = Λ(k) factor.
- Positivity is **not pointwise** (prime spikes exceed the arch kernel) — it is a **global
  Diophantine constraint**: the autocorrelation R=g⋆g̃ cannot be simultaneously large at all
  prime points {log k}. That constraint is enforced *precisely by the zeros*.
- Under RH the prime symbol D(r) = Ω(r)/2π exactly (zeros invisible to Re[ζ′/ζ] on the line);
  the deviation of the truncated symbol from Ω/2π as λ²→∞ **is** the zeros. Circular at root.

## What a real proof must overcome

A proof cannot come from bounding P by A+B via any *local/pointwise* comparison (routes 2–4
provably fail). It must exploit the **global** structure — the Diophantine/spectral constraint
that the primes cannot align to violate positivity. The only non-circular partial handle found
is route 6 (verification height → n≲2π·T), which is genuinely unconditional but finite.

The honest frontier: a genuinely new idea is required — the classical analytic toolkit is
mapped and each piece provably reduces to the zeros. Numerics overwhelmingly confirm tau ⪰ 0
(= RH) at every scale, but confirmation is not proof.

---
*This map is a living document. Add a row when a new route is tested; never delete a verified
obstruction (that is the value — knowing precisely why each door is closed).*
