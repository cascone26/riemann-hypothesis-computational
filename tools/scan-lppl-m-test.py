#!/usr/bin/env python3
"""
LPPL power-law exponent m test for RH.

Theory: If λ₁(Λ) ~ (tc-Λ)^(σ+i·ω) near a resonance tc,
and RH holds (all zeros on σ=½ line), then m = σ = ½ for ALL blobs.

The LPPL model fits:
  log|λ₁| = A + B*(tc-Λ)^m * cos(ω*ln(tc-Λ) + φ)

The m exponent (power-law envelope) should ≈ ½ = 0.5 if RH is true.
Observed m systematically different from ½ would suggest either:
  (a) RH violations at some zeros (σ ≠ ½)
  (b) Finite-N truncation effects on the effective exponent
  (c) The LPPL model is fitting log|λ₁| rather than λ₁ itself

Usage: python3 scan-lppl-m-test.py /tmp/lppl_results_all.json
"""
import json, math, sys
from collections import defaultdict

ZEROS = [14.134725, 21.022040, 25.010858]
ZLBLS = ['rho1', 'rho2', 'rho3']

RH_PREDICTION = 0.5

data = json.load(open(sys.argv[1] if len(sys.argv) > 1 else '/tmp/lppl_results_all.json'))
hits = [h for h in data if h.get('hit')]
all_blobs = [h for h in data if 'm' in h]

print(f"=== LPPL m-exponent test ===")
print(f"RH prediction: m = {RH_PREDICTION}")
print()

# Global m distribution
m_vals = [h['m'] for h in all_blobs]
m_mean = sum(m_vals)/len(m_vals)
m_sd = math.sqrt(sum((x-m_mean)**2 for x in m_vals)/(len(m_vals)-1))
print(f"All blobs: n={len(all_blobs)} m_mean={m_mean:.3f} m_sd={m_sd:.3f}")

# m histogram
m_counts = {}
for m in m_vals:
    k = round(m, 1)
    m_counts[k] = m_counts.get(k, 0) + 1
print("m histogram (all blobs):")
for k in sorted(m_counts):
    print(f"  m={k:.1f}: {'#'*m_counts[k]:<20} ({m_counts[k]})")

print()
print("=== m per zero cluster ===")
for zi, lbl in zip(ZEROS, ZLBLS):
    blobs = [h for h in hits if h['nearest_zero'] == zi]
    if not blobs: continue
    ms = [h['m'] for h in blobs]
    mean_m = sum(ms)/len(ms)
    delta_from_rh = mean_m - RH_PREDICTION
    print(f"{lbl}: n={len(blobs)} m_mean={mean_m:.3f} delta_from_RH={delta_from_rh:+.3f}")
    mc = {}
    for m in ms:
        k = round(m, 1)
        mc[k] = mc.get(k, 0) + 1
    print(f"  m dist: {dict(sorted(mc.items()))}")
    print(f"  hits: {['%s@%s m=%.1f w=%.2f' % (h['blob_start'],h.get('label','?'),h['m'],h['omega']) for h in blobs]}")
    print()

print("=== m=0.5 blobs only (closest to RH prediction) ===")
m_half_blobs = [h for h in hits if h['m'] == 0.5]
print(f"Hits with m=0.5: {len(m_half_blobs)}")
by_zero_mhalf = defaultdict(list)
for h in m_half_blobs:
    by_zero_mhalf[h['nearest_zero']].append(h['omega'])
for zi, lbl in zip(ZEROS, ZLBLS):
    ws = by_zero_mhalf.get(zi, [])
    if ws:
        mean_w = sum(ws)/len(ws)
        print(f"  {lbl}: n={len(ws)} mean_w={mean_w:.4f} theory={zi:.4f} delta={mean_w-zi:+.4f}")
