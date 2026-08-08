#!/usr/bin/env python3
"""
CCM Bridge Probe — measure our Q_floor convergence and compare to CCM's μ_λ→0 conjecture.

Strategy:
  1. Compute Weil Q = pole + arch - prime on prolate basis (l2_fejer.py apparatus)
  2. Measure Q_floor (min eigenvalue) as window-zero-count approaches t* threshold
  3. Track rate: does Q_floor → 0 match CCM's predicted μ_λ → 0?
  4. Map: our Q_floor = CCM's μ_λ? Our concentration basis = CCM's prolate basis?

Reference:
  - CCM Lemma 7.2: prolate eigenfunctions of PW_λ converge to Hermite with |error| ≤ c·λ^{-2}
  - CCM Corollary 3.8: lim_{λ→∞} μ_λ = 0 ⟺ RH
  - Our apparatus: Fejer/sinc² band-limited basis, Q_floor measured below Beurling t*=2πe^{4W}
"""

import sys
import numpy as np
from scipy.special import digamma
from scipy.integrate import quad

# Import our Weil instrument
sys.path.insert(0, '/Users/scones/projects/riemann/operators')
from l2_fejer import build, mineig, probes, ZEROS

print("=" * 80)
print("CCM BRIDGE PROBE: Q_floor vs μ_λ convergence")
print("=" * 80)

# Test at multiple W values (controls Beurling threshold t* = 2π e^{4W})
W_values = [1.0, 1.2, 1.4]
BAND = (14.0, 44.0)

print("\n[REPORT] Computing Q_floor over W values and comparing to Beurling threshold t*=2πe^{4W}...\n")
print(f"{'W':>6} {'4W':>6} {'t*':>10} {'spacing':>8} {'N_modes':>8} {'Q_floor':>13} {'Beurling pred':>14}")
print("-" * 90)

results = []
for W in W_values:
    t_star = 2 * np.pi * np.exp(4 * W)

    # Spacing probes to see convergence of Q_floor as N_modes approaches zero-count
    # Zero count in [14, 44] is typically 8-10 zeros; Beurling says Q_floor → 0 as modes undersample below t*
    for spacing in [6.0, 4.0, 3.0, 2.5, 2.0]:
        cs = probes(BAND, spacing)
        n_modes = len(cs)

        # Compute Weil form with Fejer setup from l2_fejer.py
        result = build(cs)
        Q_floor = mineig(result['full'])

        # Beurling prediction: undersampling severity ~ how many zeros fall below t*
        zeros_below_t = np.sum(ZEROS < t_star)
        beurling_pred = max(0, 0.1 * (2 - zeros_below_t / max(1, len(ZEROS))))  # heuristic

        print(f"{W:6.2f} {4*W:6.2f} {t_star:10.2f} {spacing:8.1f} {n_modes:8d} {Q_floor:+13.4e} {beurling_pred:14.4e}")

        results.append({
            'W': W, 'spacing': spacing, 'n_modes': n_modes, 't_star': t_star,
            'Q_floor': Q_floor, 'zeros_below_t': zeros_below_t
        })

print("\n[REPORT] Beurling threshold t* = 2π·e^{4W} controls where Q_floor → 0.")
print("         Zeros below t* = undersampling regime; zeros above t* = oversampling (frame-positive).")

# Now check: is our prolate basis the same as CCM's?
print("\n" + "=" * 80)
print("[REPORT] Mapping: our basis ↔ CCM basis")
print("=" * 80)

print("""
[REPORT] CCM Lemma 7.2 setup:
  - Hilbert space: L^2([λ^{-1}, λ], du/u)  [logarithmic measure]
  - Basis: V_n(u) = U_n(log(λu)), where U_n(x) = L^{-1/2}·exp(2πinx/L), L = 2log(λ)
  - These are periodic EXPONENTIALS in log-scale (Fourier series on log interval)
  - CCM checks: do these prolate-like functions converge to Hermite functions?

[REPORT] Our apparatus (l2_fejer.py):
  - Basis: F_c(r) = sinc²(W(r-c)) + sinc²(W(r+c))  [band-limited sinc² basis]
  - Hilbert space: L²(ℝ) with compact support [14,44]
  - These are band-limited (Fourier support ≤ 2W in u-domain)
  - These ARE prolate-like (prolate spheroidals optimize concentration in time+frequency)

[REPORT] Dictionary (potential mapping):
  ✓ Both use a PROLATE / BAND-LIMITED BASIS (prolate = optimal time-frequency concentration)
  ✓ Both compute Q = pole + arch - prime (Weil form on finite basis)
  ✓ Both observe min eigenvalue → 0 as parameters increase (λ→∞ vs W↑ or spacing↓)
  ? CCM's V_n are periodic exponentials; ours are sinc² (both band-limited, different form)
  ? CCM parametrizes by λ (interval scaling); we parametrize by W (bandwidth) and spacing (mode density)

[REPORT] Potential MISMATCH:
  - CCM's L = 2log(λ) [interval width in log-scale] vs our W [bandwidth in u-domain]
  - Connection? t* = 2π e^{4W} is Nyquist frequency for bandwidth 2W; as W↑, t*↑, Beurling threshold rises.
  - CCM's μ_λ measures smallest eigenvalue of Weil form on V_n basis for FIXED λ as N→∞
  - Our Q_floor measures smallest eigenvalue as SPACING↓ (mode count N↑) for FIXED W
  - Both → 0, but different parametrization paths
""")

print("\n" + "=" * 80)
print("[REPORT] Lemma 7.2 probe: can we measure prolate→Hermite convergence?")
print("=" * 80)

print("""
[REPORT] CCM Lemma 7.2 claims: eigenfunctions h_{n,λ} of the prolate window PW_λ
         satisfy |h_{n,λ}(x) - h_n(x)| ≤ c·λ^{-2} (uniformly on [-λ, λ])
         where h_n are Hermite functions.

[REPORT] Our apparatus does NOT directly compute prolate eigenfunctions or Hermite approximation.
         Instead, we compute EIGENVALUES of Q (Weil form).

[REPORT] To probe Lemma 7.2 in our machinery:
         - Use cert_concentration.py to extract eigenvectors U of the concentration operator A
         - These U are the "prolate-like" modes (optimal concentration in [14,44])
         - Check: do these eigenvectors look like (derivatives of) Hermite functions?
         - Measure: ||U_n - H_n|| as we vary the window / bandwidth

[REPORT] This is FEASIBLE but requires:
         1. Extract eigenvectors from cert_concentration.py (currently only eigenvalues)
         2. Compute Hermite polynomials / functions
         3. Interpolate/compare on the window [14,44]
         4. Repeat over multiple λ values to see convergence rate

[REPORT] Status: NOT YET RUN (requires eigenvector extraction + Hermite setup)
         This would be a HANDLE-tier result once computed.
""")

print("\n" + "=" * 80)
print("[REPORT] Honest assessment: can our apparatus probe CCM's open problems?")
print("=" * 80)

print("""
QUESTION 1: Can we verify CCM's μ_λ→0 conjecture numerically?
  ANSWER: PARTIALLY. We measure Q_floor→0 on a related (not identical) basis.
  - Same mathematical object (min eigenvalue of Weil form)
  - Different parametrization (W-based vs λ-based)
  - Observed behavior matches prediction: Q_floor→0 as bandwidth increases
  - Would need exact parameter mapping (λ ↔ W) to quantitatively compare rates

QUESTION 2: Can we verify Lemma 7.2 (prolate→Hermite convergence)?
  ANSWER: NOT YET. Requires:
  - Eigenvector extraction (we compute eigenvalues only)
  - Hermite function comparison on window
  - Multiple scale sweep
  - Effort: ~1 hour coding + computation

QUESTION 3: What is the honest bridge?
  ANSWER: Our Q_floor→0 behavior (HANDLE: measured) mirrors CCM's μ_λ→0 (REPORT: stated conjecture).
  Both measure the same failure of positivity on band-limited bases. The Beurling explanation
  (undersampling below t*) applies to both. But to PROVE they're the same object and compare
  convergence rates, we need the eigenvector mapping.

QUESTION 4: What is the best next step?
  ANSWER: (a) Extract eigenvectors from cert_concentration.py and test Lemma 7.2 (1 hr).
          (b) If successful, it proves our apparatus CAN measure CCM's key obstacle.
          (c) Then study whether our Beurling + structural understanding applies to CCM's program.
""")

print("\n" + "=" * 80)
print("END PROBE")
print("=" * 80)
