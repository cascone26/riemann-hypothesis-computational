#!/usr/bin/env python3
"""
D1 Canonical System: Extract H(t) Metric and Jacobi Coefficients from ξ(1/2+z)

The de Branges-Sakhnovich theory provides a formula relating a Hermite-Biehler
function E(z) to a canonical system. The key is the log-derivative.

For our case: E(z) = ξ(1/2+z) (after spectral coordinate change).

The canonical system is: dY/dt = z J H(t) Y,  where H(t) ⪰ 0.

One way to extract H from E: use the fact that the "characteristic function" W(z)
of the canonical system satisfies W(z) = E(z) up to a phase.

Concretely, for a Hermite-Biehler function E(z), we have:
  E(z) = (P(z) + i Q(z)) / (P(z) - i Q(z))   where P, Q are real entire functions
OR
  E(z) = B(z) / B̄(z̄)  where B entire with no UHP zeros

The metric H(t) can be extracted via:
  log E(iy) = ∫_0^∞ (yt / (1 + y²t²)) · f(t) dt   [Herglotz representation]
where f(t) is related to H(t).

Alternatively, a direct route: if E(z) = e^{iθ(z)}, then Im(log E) = θ,
and the spectral density dμ(λ) is related to dθ/dλ.

Let's implement this numerically.
"""

import mpmath as mp
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.linalg import eigh
import sys

mp.mp.dps = 35  # High precision

def xi(s):
    """Riemann xi function"""
    return mp.mpf('0.5') * s * (s - 1) * mp.power(mp.pi, -s/2) * mp.gamma(s/2) * mp.zeta(s)

def xi_log_deriv(s):
    """d/ds log ξ(s) = ξ'(s) / ξ(s)"""
    h = mp.mpf('1e-8')
    return (xi(s + h) - xi(s - h)) / (2 * h * xi(s))

# ============================================================================
# APPROACH: Extract H(t) from log-derivative phase and amplitude
# ============================================================================

print("=" * 80)
print("D1 CANONICAL SYSTEM: H-Metric Extraction from ξ(1/2+z)")
print("=" * 80)
print()

print("""
Key fact: For a Hermite-Biehler function E(z) = ξ(1/2+z), the metric H(t)
in the canonical system dY/dt = z J H(t) Y can be related to:

  (a) The phase: arg ξ(1/2+iy)  on the real axis y
  (b) The amplitude: |ξ(1/2+iy)| on the real axis y

The Herglotz representation gives:
  Im(log E(iy)) = ∫_0^∞ (2yt / (t² + y²)) ρ(t) dt   [ρ(t) = spectral density]

And the metric is H(t) = ρ(t) (after normalization).

Strategy:
1. Compute ξ(1/2 + iy) for y ∈ [0, 100] (real axis in spectral variable)
2. Extract phase and amplitude
3. Use inverse Stieltjes / Herglotz transform to recover ρ(t) from Im(log E)
4. Set H(t) = ρ(t) · I or a more sophisticated symmetric matrix
5. Solve the canonical ODE and extract Jacobi coefficients
6. Compare to Lanczos ground truth
""")

print()
print("-" * 80)
print("STEP 1: Sample ξ(1/2+iy) on the real axis")
print("-" * 80)
print()

# Sample on the real axis (imaginary spectral variable)
y_vals = np.logspace(-2, 2, 41)  # 0.01 to 100
xi_vals = []
phase_vals = []
amp_vals = []

for y in y_vals:
    try:
        z = mp.mpc(mp.mpf('0.5'), mp.mpf(str(y)))
        xi_z = xi(z)
        xi_vals.append(complex(xi_z))
        phase = mp.arg(xi_z)
        amp = abs(xi_z)
        phase_vals.append(float(phase))
        amp_vals.append(float(amp))
    except:
        pass

if xi_vals:
    print(f"Sampled ξ(1/2+iy) at {len(xi_vals)} points:")
    for i in [0, len(y_vals)//4, len(y_vals)//2, 3*len(y_vals)//4, -1]:
        if i < len(y_vals):
            y = y_vals[i]
            print(f"  y={y:.4f}:  ξ(1/2+{y:.4f}i) = {xi_vals[i]:.6f}")
    print()

print("-" * 80)
print("STEP 2: Phase and amplitude structure")
print("-" * 80)
print()

print("Phase arg(ξ(1/2+iy)):")
for i in [0, len(y_vals)//2, -1]:
    if i < len(phase_vals):
        print(f"  y={y_vals[i]:.4f}:  arg = {phase_vals[i]:.6f} rad = {np.degrees(phase_vals[i]):.2f}°")

print("\nAmplitude |ξ(1/2+iy)|:")
for i in [0, len(y_vals)//2, -1]:
    if i < len(amp_vals):
        print(f"  y={y_vals[i]:.4f}:  |ξ| = {amp_vals[i]:.8f}")

print()

# ============================================================================
# APPROACH 2: Direct numerical metric construction via spectral interpolation
# ============================================================================

print("-" * 80)
print("STEP 3: Construct H(t) via spectral density estimation")
print("-" * 80)
print()

print("""
Instead of inverting the Herglotz transform (which is ill-posed), use a
constructive approach:

1. Assume H(t) = ρ(t) I where ρ(t) is smooth and supported on [0, T_max]
2. Generate a CANDIDATE canonical system from a smooth ρ
3. Compute the resulting characteristic function W_ρ(z)
4. Check: does W_ρ match ξ(1/2+z)?  (Or at least match Jacobi coefficients?)

Smooth candidate: ρ(t) = A e^{-γt}  or  ρ(t) = 1/(1 + (t/σ)²)
""")

# Load ground truth b_n
try:
    G = np.load('/Users/scones/projects/riemann/operators/zeros2000.npy')
except:
    print("WARNING: zeros2000.npy not found, using placeholder")
    G = np.array([0.5j + np.sqrt(1j*k) for k in range(2000)])

x = 1.0/(0.25 + G.astype(float)**2)
w = x.copy()
w = w / w.sum()

def lanczos_jacobi(nodes, wt, K):
    """Return (a_n, b_n) via Lanczos"""
    q = np.sqrt(wt)
    qp = np.zeros_like(q)
    a = np.zeros(K)
    b = np.zeros(K)
    bprev = 0.0

    for k in range(K):
        v = nodes * q
        aa = np.dot(q, v)
        a[k] = aa
        v = v - aa * q - bprev * qp
        v = v - q * np.dot(q, v) - qp * np.dot(qp, v)
        nb = np.linalg.norm(v)
        if k + 1 < K:
            b[k+1] = nb
        qp = q
        q = v / nb if nb > 0 else v
        bprev = nb

    return a, b

a_truth, b_truth = lanczos_jacobi(x, w, 20)

print("Ground truth Jacobi (Lanczos):")
print("  n    a_n              b_n")
for n in range(min(10, len(b_truth))):
    print(f"  {n}    {a_truth[n]:.8e}   {b_truth[n]:.8e}")

print("\n" + "-" * 80)
print("STEP 4: Canonical system ODE solver")
print("-" * 80)
print()

print("""
Canonical system: dY/dt = z J H(t) Y,  Y(0) = I

Where Y(t) = [Y11(t) Y12(t); Y21(t) Y22(t)]  and  J = [0 -1; 1 0]

Expanding: dY_{11}/dt = -z Y_{21}(t)
           dY_{21}/dt = +z Y_{11}(t)  [if H is scalar ρ(t) I]

This is a Hamiltonian system. For a given H(t) (the metric), we integrate
to find Y(t), and extract the Weyl m-function and its limit κ → ∞.

The Jacobi coefficients are then extracted from the m-function's limit.

For our case, we'd need to:
1. Guess or estimate H(t) from ξ(1/2+z) [HARD part]
2. Solve the ODE numerically
3. Extract b_n from the spectral data
4. Check if b_n > 0 and match Lanczos
""")

print("""
HONEST ASSESSMENT: This route requires either
  (A) An explicit formula for H(t) in terms of ξ (from literature or derivation)
  (B) Numerically inverting the Herglotz transform (ill-posed, unstable)
  (C) A phenomenological fit of H(t) to match Lanczos data (backwards-engineering)

The canonical-system framework is THEORETICALLY sound: H ⪰ 0 ⟹ b_n > 0 automatic.
But EXTRACTING H from ξ is the HARD part — it's where the de Branges-Sakhnovich
machinery gets technical.
""")

print()
print("=" * 80)
print("REAL NEXT STEP")
print("=" * 80)
print("""
The canonical system / de Branges route is VALID but TECHNICAL.

Before diving into the ODE solver / Herglotz inversion, check:
1. Does Remling's book "Spectral Theory of 1D Schrödinger Operators" (2022)
   treat the canonical system for ξ(s)?
2. Does de Branges' original "Hilbert Spaces of Entire Functions" have an
   explicit metric for this family?
3. Can we extract H(t) from the known functional equation of ξ?

If (3), we might be able to prove H(t) ⪰ 0 unconditionally (without RH),
which would give b_n > 0 automatically.

CANDIDATE BREAKTHROUGH: ξ(s) satisfies ξ(s) = ξ(1-s) (functional equation).
This symmetry might directly constrain H(t) to be PSD.
""")

print()
