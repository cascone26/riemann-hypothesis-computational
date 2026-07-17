#!/usr/bin/env python3
"""
Post-processor for margin_j40.json.

The main script had a bug: margin was computed only when ev0_riem > 0,
but ev0_riem = -14462 (deeply negative), so all margins were NaN.

The CORRECT quantity to measure selectivity:

  depth_ratio(γ̃) = |ev[0](γ₁)| / |ev[0](γ̃)|

When ev[0](γ₁) is the DEEPEST (most negative), depth_ratio ≥ 1 for all γ̃.
depth_ratio = 1 at γ̃=γ₁ (reference); depth_ratio → large for γ̃ far from γ₁.

The 5782× claim from abc_credence means: at γ̃=13.8 (|Δγ|=0.33):
  depth_ratio = |ev0(γ₁)| / |ev0(13.8)| = 5782
  i.e. |ev0(13.8)| = 14462 / 5782 ≈ 2.5

For the power-law fit, we want:
  log(depth_ratio) ~ α × log(|Δγ|) + const
i.e. depth_ratio ~ A × |Δγ|^α (power law in offset from γ₁)

We expect α > 0 and depth_ratio → 1 as Δγ → 0.
"""
import json, math, numpy as np, sys

data_file = "results/margin_j40.json"
try:
    with open(data_file) as f:
        data = json.load(f)
except FileNotFoundError:
    print(f"Not done yet: {data_file} not found.")
    sys.exit(0)

ev0_riem = data["ev0_riem"]
print(f"Margin j40 post-processor")
print(f"ev[0]_riem = {ev0_riem:.6e}  (at γ₁=14.134725, lsq=52738111, N=4200)")
print()

if not data.get("points"):
    print("No points computed yet. Re-run after margin_j40.py completes.")
    sys.exit(0)

print(f"{'γ̃':>14}  {'Δγ':>+8}  {'lsq_γ̃':>14}  {'ev[0]':>14}  {'depth_ratio':>14}  label")
print("-"*100)

points_with_ratio = []
for p in data["points"]:
    gamma = p["gamma"]
    ev0 = p["ev0"]
    delta = p["delta"]
    lsq = p["lsq"]
    label = p["label"]

    # Depth ratio: |ev0_riem| / |ev0| (Riemann dip deeper → ratio > 1 for fakes)
    if ev0 != 0 and ev0_riem != 0:
        depth_ratio = abs(ev0_riem) / abs(ev0)
    else:
        depth_ratio = float("nan")

    print(f"{gamma:>14.6f}  {delta:>+8.4f}  {lsq:>14}  {ev0:>14.4e}  {depth_ratio:>14.4f}  {label}")
    p["depth_ratio"] = float(depth_ratio) if not math.isnan(depth_ratio) else None
    points_with_ratio.append(p)

# Power-law fit on fake points (not other Riemann zeros)
fakes = [p for p in points_with_ratio
         if "fake" in p["label"] and abs(p["delta"]) > 0
         and p["depth_ratio"] is not None and p["depth_ratio"] > 0]

print()
print(f"=== Power-law fit: depth_ratio(|Δγ|) ~ A × |Δγ|^α ===")
print(f"  (depth_ratio = |ev0_riem|/|ev0(γ̃)|, tells how much deeper Riemann is)")
print()

if len(fakes) >= 3:
    dg = np.array([abs(p["delta"]) for p in fakes])
    dr = np.array([p["depth_ratio"] for p in fakes])
    valid = (dg > 0) & (dr > 0) & np.isfinite(dr) & np.isfinite(dg)
    if valid.sum() >= 2:
        coeffs = np.polyfit(np.log(dg[valid]), np.log(dr[valid]), 1)
        alpha = coeffs[0]
        A = math.exp(coeffs[1])
        print(f"  depth_ratio(|Δγ|) ≈ {A:.4f} × |Δγ|^{alpha:.4f}")
        print(f"  α = {alpha:.4f}")
        if abs(alpha - 2.0) < 0.5:
            print(f"  *** α ≈ 2 — QUADRATIC selectivity ***")
            print(f"      depth_ratio ≥ C|Δγ|² → min_ev lower-bounded quadratically in distance to zero")
        elif abs(alpha - 1.0) < 0.3:
            print(f"  *** α ≈ 1 — LINEAR selectivity ***")
        elif alpha < 0:
            print(f"  *** α < 0 — depth_ratio DECREASES with |Δγ| (inverse: closer to γ₁ = deeper) ***")
        print(f"  Coefficient A = {A:.4e}")
        print()

        # Verify 5782× claim at |Δγ|=0.33 (γ̃=13.8)
        predicted_5782 = A * (0.3347 ** alpha)
        print(f"  At |Δγ|=0.3347 (γ̃=13.8): predicted depth_ratio = {predicted_5782:.1f}")
        print(f"  Claimed from abc_credence: 5782×")
        if abs(math.log(predicted_5782) - math.log(5782)) < 0.5:
            print(f"  ✓ Consistent with 5782× claim")
        else:
            print(f"  ✗ INCONSISTENT with 5782× claim")

        # Check: at γ̃=13.8 directly
        fake_13_8 = next((p for p in fakes if "13.8" in p["label"]), None)
        if fake_13_8:
            print(f"\n  Direct measurement at γ̃=13.8: depth_ratio = {fake_13_8['depth_ratio']:.2f}")
            print(f"  |ev0(13.8)| = {abs(fake_13_8['ev0']):.4e}")
        data["power_law_alpha_corrected"] = float(alpha)
        data["power_law_A_corrected"] = float(A)
    else:
        print("  Not enough valid points for fit")
else:
    print(f"  Only {len(fakes)} fake points computed so far — fit needs ≥ 3")

# Other Riemann zeros
zeros = [p for p in points_with_ratio if "zero_g" in p["label"]]
if zeros:
    print()
    print("=== Other Riemann zeros (at N=4200 — not their own Dip1) ===")
    for z in zeros:
        dr = z.get("depth_ratio")
        print(f"  {z['label']}: γ̃={z['gamma']:.4f}  |Δγ|={abs(z['delta']):.4f}  "
              f"ev[0]={z['ev0']:.4e}  depth_ratio={dr:.4f}" if dr else
              f"  {z['label']}: γ̃={z['gamma']:.4f}  ev[0]={z['ev0']:.4e}  depth_ratio=N/A")
    print("  (Expect depth_ratio ≈ 1 only at their OWN resonant Dip1, not at N=4200)")

# Save updated JSON
data["points"] = points_with_ratio
out_file = "results/margin_j40_corrected.json"
with open(out_file, "w") as f:
    json.dump(data, f, indent=2)
print(f"\nSaved corrected margins → {out_file}")
