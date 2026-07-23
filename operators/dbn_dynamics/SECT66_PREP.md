# §178d analytic prep — reading of Connes 2602.04022 §"approximation k_λ of θ_x" (banked 2026-07-23)

Purpose: have the analytic frame ready the instant the §178d decay-law compute lands, so the empirical
rate can be interpreted (not just reported). Source: `membrane_pulls/connes_2602/rhready.tex`
§sectprolate (l.1094), §poisson (l.1144), §sectconv (l.1168), §miss (l.1174). Report-tier (reading).

## What Connes' construction is (confirms our tool is faithful — Handle-adjacent)
- Prolate wave op `pw_λ = -∂_x((λ²-x²)∂_x) + (2πλx)²` on L²[-λ,λ]; eigenfns = prolate spheroidal h_{n,λ}.
- Angle operator α_λ between projections P_λ (time-limit) and P̂_λ (band-limit): P_λ cos²α_λ = P_λ P̂_λ P_λ.
  Even prolate eigenvalues χ_m, with F(h_{2m,λ}) = χ_m h_{2m,λ}, χ_m² = ν_m, sign (-1)^m.
- **k_λ(u) = E(h_λ)(u)** on [λ⁻¹,λ], where h_λ = the UNIQUE (up to scalar) combo of h_{0,λ}, h_{4,λ} with
  vanishing integral. Computation of E(h_λ)(u) "only involves the sum over integers ≤ λ²".
  ⇒ our `bridge_overlap.k_lambda_on` (even eigvecs i0,i4; b0=-V[0,i0]/V[0,i4] to zero the integral;
  nmax=floor(λ²)+1) IS exactly this. Faithful.
- Near-null eigenvalue of A_λ (QW_λ(f,f)=⟨A_λ f|f⟩): ε(λ) ~ 1-χ_2, **exponential-of-exponential** (Fuchs Thm 1):
  `1-χ_2 ~ (2^14/3)√2 π⁵ e^{-4π e^L + (9/2)L}`, L=2log λ, e^L=λ²=c ⇒ e^{-4πc}. = our §171 envelope, SOURCE-CONFIRMED.

## THE CRUCIAL DISTINCTION (keep Report/Handle clean when reporting §178d)
Two DIFFERENT convergences in the paper — do not conflate:
1. **k̂_λ → Ξ** (Fourier/zeta side, §sectconv Fact): ALREADY CONTROLLED by Connes — polynomial rate,
   diff on line Im(z)=α bounded by `c·λ^{-1/2-α}(1-2α)⁻¹`. This is DONE analytically.
2. **k_λ → θ_x** (§miss "Remaining steps"): the OPEN piece — "it still remains to show that k_λ is a
   sufficiently good approximation of θ_x, λ=√x." θ_x defined earlier (l.1045), θ̂_x all-real-zeros ⟺ RH.

**What §178d actually measures:** dist(k_λ, ψ_c) = how well the prolate guess k_λ matches the TRUE
finite-c near-null eigenvector ψ_c of A_λ (Connes: k_λ is "an approximation of the eigenvector associated
to the smallest eigenvalue ε(λ)"). This is Handle-earning (real mpmath computation), and it directly
quantifies convergence #(the eigenvector approximation). It is RELATED to but NOT identical with the open
limit estimate #2 (k_λ→θ_x as λ→∞). Report the measured law as the eigenvector-approximation rate; the
λ→∞ limit estimate remains the open analytic target. Don't overclaim §178d = "closed §6.6".

## ANALYTIC LEAD to test against the decay law (the reason to have this ready)
k_λ uses ONLY the h_0,h_4 combo. The true ψ_c has small components along HIGHER even prolates (h_8,h_12,…).
Hypothesis: **dist(k_λ, ψ_c) ≈ contamination from the next even prolate mode**, i.e. governed by the ratio
of the next angle-operator eigenvalue (1-χ_4 vs 1-χ_2) or the next Fuchs term. If so, dist(c) should decay
at the Fuchs *gap between consecutive even prolate modes*, NOT at e^{-4πc} itself. Concrete test once the
decay law is in hand: fit log dist(c) vs c and compare the slope to (a) e^{-4πc} [same as ε], (b) the
1-χ_4 envelope, (c) the poly λ^{-1/2} rate of convergence #1. Whichever matches identifies WHAT the
correction is, which is the analytic handle.

## Early empirical (from the running §178d, c=3):
N=40→dist 6.957e-4, N=60→6.894e-4 (PLATEAU, <1%), log10λmin=-7.29 stable, 1-ov~2.4e-7. So dist(k_λ,ψ_3)≈6.9e-4
is the TRUE finite-c correction at c=3. Awaiting c=4,5,6 to fit dist(c). (c=13 from §178c never plateaued
at N≤80 — deep near-null needs more N; small-c plateaus are the clean law source.)

NEXT after law lands: (1) fit slope, (2) test the 3 hypotheses above, (3) LEDGER §178d + PROCESS append
with the number as proof pointer. Model: analytic estimate attempt = Opus.
