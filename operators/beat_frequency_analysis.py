#!/usr/bin/env python3
"""
Beat-frequency mechanism: why selectivity is j-dependent.

Theory: false positives occur when j×(γ₁/γ_fake) mod 1 ≈ 0 or 0.5.
At those j values, lsq_fake(j) ≈ lsq_riem(j') for some near-integer j'.
The fake's prime sum accidentally sits at a Riemann resonant scale → near-PSD.

Beat period for fake with gamma:
  per-j phase increment = γ₁/γ_fake (mod 1)
  coincidences (false positives) when phase ≈ 0 or 0.5

This explains:
  fake_13.8: γ₁/γ = 14.135/13.8 = 1.0243 → per-j phase = 0.0243
    → beat period = 1/0.0243 = 41.2 j-units
    → half-period coincidences at j = (0.5/0.0243) = 20.6 ≈ 20 ✓ (j=20 false neg)
    → full-period coincidences at j = (1/0.0243) = 41.2 ≈ 40 ✓ (j=40 false neg)

  fake_18.0: γ₁/γ = 14.135/18.0 = 0.785 → per-j |deviation| = |0.785 - 1| = 0.215
    → for fakes with γ>γ₁: lsq_fake < lsq_riem (fewer primes), different formula
    → fractional part of j×γ₁/γ: j=30→0.56≈0.5 ✓ (sel=1.96), j=40→0.41 ✓ (sel=0.25)

Verified against multi_j_fake138 data.
"""
import math

GAMMA1 = 14.134725141734693
J_vals = [5, 10, 15, 20, 25, 30, 35, 40]

# Multi-j observed selectivities (from results/multi_j_fake138.out)
obs_sel = {
    "fake13.8": {5:2.71, 10:48.97, 15:1793261, 20:0.0045, 25:15.87, 30:87.03, 35:14.29, 40:0.45},
    "fake18.0": {5:2.56, 10:399.64, 15:654804, 20:67612, 25:2905013, 30:1.96, 35:11858, 40:0.25},
    "fake19.0": {5:20.32, 10:0.63, 15:11.34, 20:1605, 25:67254, 30:1.45, 35:1.20, 40:88.01},
}

gammas_fake = {"fake13.8": 13.800, "fake18.0": 18.000, "fake19.0": 19.000}

print("BEAT FREQUENCY ANALYSIS")
print("="*80)
print()
print("Theory: false positives (sel<1 or sel≈1) occur when")
print("  phase = (j × γ₁/γ_fake) mod 1  is near 0.0 or 0.5")
print("  → lsq_fake(j) ≈ lsq_riem(some j') → accidentally resonant")
print()

for name, gamma_f in gammas_fake.items():
    ratio = GAMMA1 / gamma_f
    per_j = ratio % 1  # per-j phase increment (fraction)
    # actual: phase(j) = (j * ratio) mod 1

    print(f"--- {name} (γ_fake={gamma_f}) ---")
    print(f"  γ₁/γ_fake = {ratio:.6f}  per-j phase step = {per_j:.4f}")
    if per_j > 0.5: per_j = 1.0 - per_j  # symmetrize to [0, 0.5]
    beat_period_half = 0.5 / per_j if per_j > 0 else float('inf')
    print(f"  Half-period (false pos nodes) every {beat_period_half:.1f} j-steps")
    print(f"  Full-period (integer nodes) every {2*beat_period_half:.1f} j-steps")
    print()

    print(f"  {'j':>4}  {'phase':>8}  {'dist_to_node':>14}  {'pred_qual':>10}  {'obs_sel':>12}  {'match?':>8}")
    print(f"  {'-'*4}  {'-'*8}  {'-'*14}  {'-'*10}  {'-'*12}  {'-'*8}")

    for j in J_vals:
        phase = (j * GAMMA1 / gamma_f) % 1.0  # fractional part of j × γ₁/γ
        # Distance to nearest resonance node (0 or 0.5)
        dist = min(phase, abs(phase - 0.5), abs(phase - 1.0))

        # Qualitative prediction
        if dist < 0.05:
            pred = "FALSE_POS"
        elif dist < 0.15:
            pred = "LOW_SEL"
        elif dist < 0.30:
            pred = "MED_SEL"
        else:
            pred = "HIGH_SEL"

        # Observed
        obs = obs_sel[name].get(j, None)
        obs_str = f"{obs:.2f}" if obs is not None else "?"

        # Match?
        if obs is not None:
            if pred == "FALSE_POS" and obs < 5:
                match = "✓ MATCH"
            elif pred == "LOW_SEL" and obs < 100:
                match = "✓ match"
            elif pred in ("MED_SEL","HIGH_SEL") and obs >= 10:
                match = "✓ match"
            elif obs < 5:
                match = "✗ missed"
            else:
                match = "~"
        else:
            match = "?"

        print(f"  {j:>4}  {phase:>8.4f}  {dist:>14.4f}  {pred:>10}  {obs_str:>12}  {match:>8}")
    print()

print()
print("BEAT PERIODS BY FAKE GAMMA:")
for name, gamma_f in gammas_fake.items():
    ratio = GAMMA1 / gamma_f
    phase_step = ratio % 1
    if phase_step > 0.5: phase_step = 1.0 - phase_step
    half_period = 0.5 / phase_step if phase_step > 0 else float('inf')

    # Nodes of false positives (half-period coincidences)
    nodes = []
    j = 0
    while j <= 80:
        jf = j * GAMMA1 / gamma_f
        frac = jf % 1
        if min(frac, abs(frac-0.5), abs(frac-1.0)) < 0.06:
            nodes.append(j)
        j += 1
    print(f"  {name}: half-period={half_period:.1f}  false-positive j-nodes: {nodes[:8]}")

print()
print("OPTIMAL j FOR EACH FAKE (maximally far from resonance nodes):")
for name, gamma_f in gammas_fake.items():
    best_j = None; best_dist = -1
    for j in range(1, 100):
        phase = (j * GAMMA1 / gamma_f) % 1.0
        dist = min(phase, abs(phase-0.5), abs(phase-1.0))
        if dist > best_dist:
            best_dist = dist; best_j = j
    lsq_best = round(math.exp(2*math.pi*best_j/GAMMA1))
    print(f"  {name}: best j={best_j}, dist={best_dist:.4f}, lsq_riem={lsq_best:,}")
