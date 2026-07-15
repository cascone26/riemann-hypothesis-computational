#!/usr/bin/env python3
"""
density_analysis.py — Analyze resonance density as a function of lsq.

Key question: How does density of deep minima in λ_min(τ) scale with lsq?

Expected from theory:
- N(γ,T) ~ (T/2π) ln(T/2π) — zero density grows as T ln T
- Modes n_k = γ_k L / 2π — at larger L, more modes exist within [0,N]
- Multi-zero alignment probability grows with more zeros in the window
- Therefore: density of resonances should INCREASE with lsq

Combine data from:
1. gamma4_scan (lsq=700-980, step=5)
2. resonance_hunt_1000_2000 (lsq=1000-1178, step=2, score-filtered)

Estimate density in windows and fit scaling law.
"""
import numpy as np

# Data from gamma4_scan (step=5, all computed):
DEEP_700_980 = {715, 780, 805, 845, 870, 895, 900, 960}
SPAN_700_980 = (700, 980, 5)  # start, end, step → 56 points total

# From earlier scans (step=5):
DEEP_560_700 = {625, 635, 665, 695}
SPAN_560_700 = (560, 700, 5)  # 28 points

# From validate_predictor (step=2, focused regions only):
DEEP_800_CLUSTER = {802, 804}  # lsq=790-810 region
DEEP_1237_CLUSTER = {1230, 1236, 1238, 1248, 1250}
DEEP_1382_CLUSTER = {1380, 1382, 1398}

# From resonance_hunt_1000_2000 (step=2, score-filtered):
DEEP_1000_1178 = {
    1020, 1032, 1038, 1062, 1070, 1092, 1102, 1104,
    1108, 1110, 1116, 1118, 1122, 1124, 1128, 1130,
    1138, 1142, 1152, 1162, 1164, 1170, 1172
}
SPAN_1000_1178 = (1000, 1178, 2)

def density_in_window(deep_set, start, end, step):
    """Deep minima per 100 lsq units, corrected for step size."""
    # Total points in range
    n_total = len(range(start, end+1, step))
    n_deep = sum(1 for x in deep_set if start <= x <= end)
    span = end - start
    # Density in 'deep per step' then extrapolated
    density_per_100 = (n_deep / n_total) * (100 / step)
    density_raw = n_deep / (span / 100)
    return n_deep, n_total, span, density_raw, density_per_100

print("Resonance Density Analysis")
print("="*70)
print()
print("Window analysis (deep minima per 100 lsq units):")
print(f"{'Range':>15} {'N_deep':>7} {'N_total':>9} {'step':>5} {'density/100':>12}")
print("-"*60)

windows = [
    ("560-700", DEEP_560_700, 560, 700, 5),
    ("700-980", DEEP_700_980, 700, 980, 5),
    ("1000-1178", DEEP_1000_1178, 1000, 1178, 2),
]

centers = []
densities = []
for label, deep, start, end, step in windows:
    n_deep, n_total, span, d_raw, d_corr = density_in_window(deep, start, end, step)
    center = (start + end) / 2
    print(f"{label:>15} {n_deep:>7} {n_total:>9} {step:>5} {d_raw:>12.2f}")
    centers.append(center)
    densities.append(d_raw)

print()
print("Power law fit: density ~ C × lsq^α")
if len(centers) >= 2:
    log_c = np.array([np.log(c) for c in centers])
    log_d = np.array([np.log(d) for d in densities])
    coeffs = np.polyfit(log_c, log_d, 1)
    alpha = coeffs[0]
    C = np.exp(coeffs[1])
    print(f"  α = {alpha:.4f}  (expected ~1 from T log T counting)")
    print(f"  C = {C:.4e}")
    for c, d in zip(centers, densities):
        pred = C * c**alpha
        print(f"  lsq={c:.0f}: observed={d:.2f}, predicted={pred:.2f}")

print()
print("Cluster structure:")
print("  lsq=1100-1130: 8 deep in 30 units → dense cluster")
print("  lsq=1160-1172: 4 deep in 12 units → dense cluster")
print("  Gaps: 1080-1090 (0 deep), 1132-1136 (0 deep)")
print()
print("Score vs depth correlation:")
# From lsq=1100-1178 data
# High score but NOT deep: 1100(0.461), 1106(0.557), 1112(0.588), 1114(0.602), 1120(0.707)
# Low score but deep: 1092(0.365), 1138(0.499)
print("  Score>0.50 NOT deep: 1106(sc=0.557), 1112(0.588), 1114(0.602), 1120(0.707)")
print("  Score<0.50 ARE deep: 1092(sc=0.365), 1138(0.499)")
print("  → Score is necessary but not sufficient; mode DIVERSITY matters")
print()
print("Rate of false positives (score>0.35 but NOT deep):")
all_computed_1100_1130 = [1092,1094,1096,1098,1100,1102,1104,1106,1108,1110,1112,1114,1116,1118,1120,1122,1124,1126,1128,1130]
deep_1100_1130 = {1092,1102,1104,1108,1110,1116,1118,1122,1124,1128,1130}
not_deep = [x for x in all_computed_1100_1130 if x not in deep_1100_1130]
print(f"  Computed: {len(all_computed_1100_1130)}, Deep: {len(deep_1100_1130)}, Not deep: {len(not_deep)}")
print(f"  Hit rate: {len(deep_1100_1130)/len(all_computed_1100_1130):.1%}")
print(f"  Not deep: {not_deep}")
