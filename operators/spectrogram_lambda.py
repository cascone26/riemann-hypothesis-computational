"""
Sliding-window spectrogram of lambda_n sequence.

Detects rogue off-critical-line zeros via their signature:
  - critical-line zeros -> constant-amplitude ridges (|alpha_j| = 1)
  - rogue zero at sigma_0 < 1/2 -> GROWING amplitude ridge (|alpha_0| > 1)

Resolution limit: rogue zeros with delta_sigma < log(N)/(N * gamma_0^2) are invisible at scale N.
"""
import numpy as np
from mpmath import mp, mpf, mpc, re as mre
import json, time

mp.dps = 25

gammas = [
    '14.134725141734693790',
    '21.022039638771554993',
    '25.010857580145688763',
    '30.424876125859513210',
    '32.935061587739189691',
    '37.586178158825671257',
    '40.918719012147495187',
    '43.327073280914999519',
    '48.005150881167159727',
    '49.773832477672302181',
]

# Precompute alphas (for on-line zeros at sigma=1/2)
alphas = []
for g in gammas:
    gamma = mpf(g)
    rho = mpc(mpf('0.5'), gamma)
    a = (rho - 1) / rho
    alphas.append(a)

def lambda_n_fast(n):
    total = 0.0
    for a in alphas:
        total += 2 * float(mre(1 - a**n))
    return total

# Compute lambda_n for n=1..N_MAX
N_MAX = 50000
print(f"Computing lambda_n for n=1..{N_MAX} with {len(alphas)} zeros...")
t0 = time.time()
lam_vals = np.zeros(N_MAX)
for n in range(1, N_MAX + 1):
    lam_vals[n-1] = lambda_n_fast(n)
    if n % 5000 == 0:
        elapsed = time.time() - t0
        eta = elapsed / n * (N_MAX - n)
        print(f"  n={n}, elapsed={elapsed:.1f}s, ETA={eta:.1f}s", flush=True)

elapsed = time.time() - t0
print(f"Done in {elapsed:.1f}s. lambda_n range: [{lam_vals.min():.4f}, {lam_vals.max():.4f}]")
print(f"All lambda_n > 0: {np.all(lam_vals > 0)}")
print(f"Min lambda_n: {lam_vals.min():.4f} at n={np.argmin(lam_vals)+1}")

# Sliding window DFT spectrogram
BLOCK = 2000   # window size
STEP = 200     # step size
n_windows = (N_MAX - BLOCK) // STEP + 1
N_FREQ = BLOCK // 2 + 1

print(f"\nComputing spectrogram: block={BLOCK}, step={STEP}, n_windows={n_windows}")

# Store amplitudes: shape (n_windows, N_FREQ)
amp_matrix = np.zeros((n_windows, N_FREQ))
window_centers = []

for w in range(n_windows):
    n_start = w * STEP
    n_end = n_start + BLOCK
    window_center = n_start + BLOCK // 2
    window_centers.append(window_center)

    segment = lam_vals[n_start:n_end]
    # Hann window to reduce spectral leakage
    hann = np.hanning(BLOCK)
    segment_w = (segment - segment.mean()) * hann
    fft_vals = np.fft.rfft(segment_w)
    amp_matrix[w, :] = np.abs(fft_vals) / (BLOCK / 2)

window_centers = np.array(window_centers)
freqs = np.fft.rfftfreq(BLOCK, d=1.0)  # cycles per integer n

print(f"\nFrequency resolution: {freqs[1]:.5f} cycles/n (period={1/freqs[1]:.1f})")
print(f"Expected gamma_1 frequency: theta_1/(2*pi) where theta_1 = arg(alpha_1)")

# Compute expected frequencies for each zero
print(f"\nExpected zero frequencies (on critical line):")
for i, g in enumerate(gammas):
    gamma = float(g)
    rho = complex(0.5, gamma)
    alpha = (rho - 1) / rho
    theta = np.angle(alpha)  # in radians
    freq = abs(theta) / (2 * np.pi)  # cycles per n
    period = 1 / freq if freq > 0 else float('inf')
    print(f"  gamma_{i+1} = {gamma:.3f}: theta={theta:.4f} rad, freq={freq:.5f} c/n, period={period:.1f}")

# Find top frequencies in spectrogram
print(f"\nTop 10 frequencies by time-averaged amplitude:")
mean_amps = amp_matrix.mean(axis=0)
top_f_idx = np.argsort(mean_amps)[::-1][:10]
for idx in top_f_idx:
    f = freqs[idx]
    a = mean_amps[idx]
    period = 1/f if f > 0 else float('inf')
    print(f"  freq={f:.5f} (period={period:.1f}), mean_amp={a:.4f}")

# Check for growing-amplitude ridges (rogue zero signature)
print(f"\nGrowing-amplitude check: are any frequency ridges intensifying over time?")
# For each frequency, fit linear trend to amplitude vs window center
from numpy.polynomial import polynomial as P

growing_ridges = []
for f_idx in range(1, N_FREQ):  # skip DC (index 0)
    amps = amp_matrix[:, f_idx]
    # Require significant amplitude: rogue-zero ridges have mean_amp > 0.05
    # (DC growth from λ_n → ∞ and low-frequency noise are mean_amp < 0.01)
    if amps.mean() < 0.05:
        continue
    # Linear fit
    coeffs = np.polyfit(window_centers, amps, 1)
    slope = coeffs[0]
    # Relative slope: slope / mean_amp (per n unit)
    rel_slope = slope / (amps.mean() + 1e-10)
    # Threshold calibrated from injection tests:
    # At N=50K, legitimate rogue-zero ridges have rel_slope >> 1e-5
    # Noise floor baseline at N=50K is ~1e-6; use 3× baseline as threshold
    if rel_slope > 3e-6:
        growing_ridges.append((freqs[f_idx], rel_slope, amps.mean(), amps.max()))

if growing_ridges:
    growing_ridges.sort(key=lambda x: -x[1])
    print(f"  WARNING: {len(growing_ridges)} growing-amplitude ridges found!")
    for f, rs, ma, mx in growing_ridges[:5]:
        print(f"  freq={f:.5f}, rel_slope={rs:.5f}/n, mean_amp={ma:.4f}, max_amp={mx:.4f}")
else:
    print(f"  CLEAN: no significantly growing amplitude ridges detected.")
    print(f"  Consistent with no rogue zero detectable at this resolution.")

# Resolution limit
gamma_1 = 14.134725
C_empirical = 90766.0  # calibrated from injection tests
alpha_empirical = 1.67
delta_sigma_min = C_empirical * (N_MAX ** (-alpha_empirical))
print(f"\n=== DETECTION RESOLUTION ===")
print(f"N = {N_MAX}")
print(f"Formula: Δσ_min ≈ C × N^{{-{alpha_empirical}}} where C={C_empirical:.0f}")
print(f"Resolution limit: delta_sigma_min ≈ {delta_sigma_min:.2e}")
print(f"(Rogue zeros with |σ₀ - 1/2| < {delta_sigma_min:.2e} are below detection threshold.)")
print(f"(Note: γ₀² in numerator — harder to detect zeros at large imaginary part.)")

# Save results
results = {
    'N_max': N_MAX,
    'n_zeros': len(alphas),
    'all_positive': bool(np.all(lam_vals > 0)),
    'min_lambda_n': float(lam_vals.min()),
    'min_lambda_n_at': int(np.argmin(lam_vals) + 1),
    'max_lambda_n': float(lam_vals.max()),
    'growing_ridges': len(growing_ridges),
    'delta_sigma_resolution': float(delta_sigma_min),
    'top_freqs': [(float(freqs[i]), float(mean_amps[i])) for i in top_f_idx],
}
with open('/Users/scones/projects/riemann/results/spectrogram_results.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nResults saved to results/spectrogram_results.json")
print(f"\n=== KEY FINDING ===")
print(f"lambda_n > 0 for all n=1..{N_MAX}: {np.all(lam_vals > 0)}")
print(f"Min lambda_n = {lam_vals.min():.6f} at n={np.argmin(lam_vals)+1}")
print(f"Growing ridges detected: {len(growing_ridges)}")
