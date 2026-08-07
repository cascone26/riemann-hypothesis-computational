# The Strip-Envelope Tail Lemma (Riemann8, loaded by roundtable §788; banks as §791)

**Purpose.** Close the last zero-hypothesis in the open-rung pair. After §788.1/.2 the PD half is
zero-free (unconditional). The indefiniteness half (§476c) still bounds its zero tail one-sidedly,
valid only for on-line configurations above G=800. This lemma bounds the tail for **arbitrary**
zero configurations in the critical strip, with explicit constants, using the closed form of the
Fejér probes — no Phragmén–Lindelöf, no L¹ time-side loss.

## Statement

Let ψ = Σⱼ vⱼ F_{cⱼ} with F_c(r) = fej(W(r−c)) + fej(W(r+c)), fej(y) = (sin y / y)², centers
cⱼ ≤ B₁ = 44, and let the zero term of the Weil form be Σ_ρ ĝ(t_ρ − iσ_ρ) over nontrivial zeros
ρ = ½ + σ_ρ + i t_ρ (so σ_ρ ∈ [−½, ½], off-line pairs ρ ↔ 1−ρ̄ sharing t_ρ), with ĝ = ψ̂².

**Lemma.** For every zero configuration in the strip (no on-line-ness assumed anywhere above G):

  | Σ_{|t_ρ| > G} ĝ(t_ρ − iσ_ρ) |  ≤  S_strip(G) := cosh⁴(W/2) · S_online(G),

where S_online(G) = |v|₁² · (4/W⁴) · 2(ln G + ½)/(3(G−B₁)³) is exactly the §476c on-line tail
bound (density N(t+1) − N(t) ≤ 2 ln t, which counts ALL strip zeros — von Mangoldt/Backlund — not
just on-line ones).

## Proof (three lines, closed form)

1. **Pointwise envelope, complex argument.** For z = t + iσ, |σ| ≤ ½:
   |fej(W(z−c))| = |sin(W(z−c))|² / |W(z−c)|² ≤ cosh²(Wσ) / (W(t−c))²,
   since |sin(x+iy)|² = sin²x + sinh²y ≤ cosh²(y) and |W(z−c)|² ≥ W²(t−c)².
   Hence |F_c(z)| ≤ 2 cosh²(Wσ) / (W(t−B₁))² for t > 2B₁ (both mirrors dominated by the near one),
   and |ψ̂(z)|² ≤ |v|₁² · 4 cosh⁴(Wσ) / (W⁴ (t−B₁)⁴).

2. **Off-line pairs.** A pair (ρ, 1−ρ̄) contributes ĝ(t−iσ) + ĝ(t+iσ) = 2 Re[ψ̂(t+iσ)²], of
   magnitude ≤ 2|ψ̂(t+iσ)|² — i.e. per zero, the same 2·envelope as the on-line case. On-line
   zeros above G contribute 2ψ̂(t)² ∈ [0, 2·envelope]. So every zero above G, wherever it sits in
   the strip, is covered by ± 2·(envelope with cosh⁴).

3. **Sum by density.** Identical to §476c's integral bound with the constant carried through. ∎

## The numbers (Arb-certified in operators/strip_tail_cert.py)

- W = 0.33 (4W = 1.32), σ_max = ½: **cosh⁴(0.165) = 1.0557** — the strip costs 5.6%, not the
  factor ~2 estimated via the e^{4W|σ|} L¹ route at the table (that route is valid but lossy; the
  closed form is the right instrument).
- §476c restated: components prime₃ = −2.020e-4, zero-sum(≤800) = +1.176e-4, and now the tail
  **two-sided** ±cosh⁴·7.62e-5 ≈ ±8.04e-5 ⇒ Q_no3(v) ∈ **[−1.6486e-4, −3.9887e-6]** — entirely
  negative, **Arb-CERTIFIED** (operators/strip_tail_cert.out, banked §791).
  **The indefiniteness certificate survives with NO on-line hypothesis above 800 at all** — margin
  thinner (−0.40e-5 vs −0.80e-5) but strictly negative.
- Citing Platt–Trudgian to 3×10¹² instead: S_strip(3e12) ~ 1e-36 — negligible; margin returns to
  the full −0.8e-5. Either hypothesis set is now completely enumerated; the T₀=800 version needs
  only: zeros to 800 on-line + counting bound. Nothing else.

## What this lemma is NOT (Red's constraint, §788 round)

It bounds strip zeros **adversarially** — enemies only. It cannot serve as the frame/sampling
inequality's supply of positive mass (that would need unverified zeros as friends = on-line-ness =
the RH-fragment). The continuum rung theorem's positivity must come from the height-localized
class + verified zeros; this lemma only guards the boundary. Both halves of the open-rung pair now
carry fully enumerated hypotheses:
- **PD half:** none (zero-free direct route, §788).
- **Indefiniteness half — FOUR hypotheses (Gauss + Red, §791-addenda):** (1) the zero list is
  **complete** to height G (Gauss: N(800)≈490.7→491 ✓; N(2000)≈1517 ✓), (2) all on-line
  (Platt–Trudgian), (3) the counting bound N(t+1)−N(t) ≤ 2 ln t, (4) **ordinate enclosures** —
  mpmath dps-30 values correct within the declared 1e-9 balls (Red: a distinct trust point; PT does
  not supply ordinates. Hardening path: interval-Newton re-enclosure on ξ). (§791, G=2000 margin
  −7.9e-5 per addendum-2.)

**Adversarial status:** verified by Gauss (roundtable, §791-addendum) — envelope, pair bookkeeping,
and enclosure independently recomputed; STANDS. Class-descriptor correction absorbed: the probe
class is pure Fejér with center spacing s=2.0 (quartic height decay), NOT Gaussian-windowed.
Margin discipline: build nothing on the thin G=800 margin −3.99e-6; the G=2000 raise
(strip_tail_g2000.py) restores ~−8e-5.
