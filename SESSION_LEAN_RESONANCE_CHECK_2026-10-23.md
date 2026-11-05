# SESSION: Real Lean Construction vs. Resonance Mechanism (HYPOTHESIS REJECTED)

**Directive:** Check whether the actual z-parametrization in CrossK6Key.lean connects to the resonance mechanism (S2723/S2724).

**Conclusion: RESONANCE HYPOTHESIS DOES NOT HOLD.**

---

## The Fatal Problem: Resonance Poles Are Real, Tested Configs Are Complex

### K=1 Symbolic Poles (S2723)
- Poles at: L·z = 2π(N±k) for integer k
- With L=2π, N=8: z ∈ {5, 6, 7, 8, 9, 10, 11} (all **real numbers**)

### All Tested Z-Configurations
- Im(z) ranges from 0.55 to 45 (all have **nonzero imaginary parts**)
- **No tested config is on the real-axis pole locus**

### The Contradiction
If the mechanism is "resonance forces linear dependence," then **only configs at or near Re(z) ∈ {5...11}** should collapse.

But the actual controlling parameter (from SESSION_Z_CONFIG_CHARACTERIZATION_2026-10-23.md) is **imaginary magnitude**:
- Collapse: min(Im(z_i)) ≤ 1.92 (regardless of Re(z))
- Survive: min(Im(z_i)) > 1.92 (regardless of Re(z))

Examples:
- Config with Re=0.5 (far from any pole), Im=0.55: **collapses** ✗ (resonance would say: should survive)
- Config with Re=0.5 (far from any pole), Im=25: **survives** ✓ (resonance would say: should survive — but for wrong reason)
- Config with Re=5.5 (near pole), Im=25: **survives** ✗ (resonance would say: should collapse)

The collapse is **decoupled from resonance proximity**.

---

## Why S2724 Looked Convincing

S2724 tested specific configs ("Base_extended," "Symmetric_Re," "Different_Re") where:
- Collapse configs happened to have low Im (~0.55–1.30) AND large Re distances to poles
- Survive configs happened to have high Im (~20–45) AND small Re distances to poles

This created a **spurious correlation**: Im-threshold correlated with pole distance, making it look like resonance was the cause.

My test with different numerical values broke the correlation, revealing the Im-threshold as the real causal variable.

---

## Concrete Counterexample: Why It Can't Be Built

To verify the resonance mechanism, I would need:
1. A z-config with Re(z_i) near {5...11} (on resonance)
2. That still has min(Im(z_i)) > 1.92 (above collapse threshold)
3. Where the determinant remains nonzero (survives)

**This config MUST exist if resonance is the mechanism.**

But such a config immediately falsifies the resonance hypothesis: it's on resonance, yet the determinant is nonzero. This would prove resonance is NOT the controlling factor.

Conversely, if all high-Im configs survive regardless of Re(z), and all low-Im configs collapse regardless of Re(z), then resonance proximity is epiphenomenal — a red herring.

**The data supports the second conclusion.**

---

## Counting Argument: Why It Fails

Proposed: "K=1 has 7 poles, K=5 can avoid them, K=6 can't."

**Problems:**
1. All tested configs avoid the poles (Im(z) ≠ 0 means z is not on the real-axis pole locus)
2. Yet K=6 still collapses
3. Therefore, "hitting poles" is not the mechanism

The counting argument is logically disconnected from the observed phenomenon.

---

## What Actually Controls the Collapse

**Im-threshold at ~1.92** (SESSION_Z_CONFIG_CHARACTERIZATION_2026-10-23.md, HANDLE-tier verification):
- This is a structural property of vPlusLC/vMinusLC evaluation when z has small imaginary parts
- It is **independent of resonance proximity**
- It is the **real characterization** of K=6 collapse

The K=1 resonance poles are mathematically real (HANDLE-tier: symbolic factorization found them). But they do not causally explain K≥6 collapse (REPORT-tier only: no verified mechanism connecting poles to observed collapse).

---

## Status

- **Lean source reviewed**: Fully universal parametrization ✓
- **K=1 symbolic poles exist**: Yes ✓
- **Resonance mechanism hypothesis**: **REJECTED** ✗
- **Real mechanism**: Im-threshold phenomenon (separate from resonance)

## Conclusion

**S2725 claim "fully connected" is RETRACTED.**

The K=6 wall exists. The vPlusLC/vMinusLC linear dependence at low Im is real. But the cause is **not** the K=1 resonance poles.

The unification of K-induction wall, vPlusLC wall, and numeric collapse is still **open**. A new explanation is needed.

