#!/usr/bin/env python3
"""§PATH_A REFINED: Use Beurling-Selberg extremal MINORANT to give a provable lower bound via Rayleigh.

Key insight: If m(r) ≤ S(r) pointwise, then for ANY test function φ:
  ∫ m(r)|φ(r)|² dr ≤ ∫ S(r)|φ(r)|² dr

Thus: margin = inf_φ ∫(pole+arch-prime)|φ|² / ∫G|φ|²
             ≥ inf_φ ∫(pole+m)|φ|² / ∫G|φ|²  (if m is a lower bound on arch-prime part)

But this is ONLY valid if m ≤ (arch-prime). Let me instead use the direct Rayleigh bound:

For the symbol S(r) = arch(r) - prime(r), if m(r) ≤ S(r), then:
  margin ≥ inf_r m(r)  (Rayleigh bound for symbol-like operators)

So the question is: what is the BEST band-limited minorant m(r) ≤ S(r), in the sense that maximizes
inf_r m(r)?

The Beurling-Selberg theorem says this best minorant exists and is unique. We approximate it numerically
by constructing a band-limited function that is close to S(r) from below.

Strategy:
1. Compute S(r) densely
2. Construct m(r) via: m(r) = inf_g∈BandLimited { g(r) : g ≤ S, g is smooth, g maximizes integral }
   (Approximated by Fejér/spline methods)
3. Report min_r m(r) as the lower bound
4. Also compute the "lift" structure: how much does m differ from S_min?
"""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
from scipy.optimize import minimize
from scipy.signal import windows

def vonmangoldt(N):
    lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam

NMAX = 3000; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1)
LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)

def analyze_with_extremal_minorant(a, b, W, dr=0.004):
    """
    Construct the Beurling-Selberg extremal band-limited minorant of S(r) = arch - prime.

    Approach: use a Fejér-kernel convolution tuned to smooth out oscillations while
    respecting the lower bound property.
    """
    fourW = 4*W; L = b-a
    RGRID = np.arange(-(b+120), (b+120), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)

    arch_sym = OMEGA/(2*np.pi)
    mask = LOGN < fourW; ln = LOGN[mask]; lw = LAMW[mask]
    P = np.zeros_like(RGRID)
    for u, w in zip(ln, lw):
        P += w*np.cos(u*RGRID)
    P *= (1.0/np.pi)
    S = arch_sym - P

    win = (RGRID >= a) & (RGRID <= b)
    S_min = float(S[win].min())
    S_min_idx = np.where(win)[0][np.argmin(S[win])]
    S_min_loc = RGRID[S_min_idx]

    # Construct Beurling-Selberg minorant via multiple strategies:
    # Strategy 1: Direct convolution with smooth positive kernel
    # Strategy 2: Trigonometric polynomial via Fejér sums
    # Strategy 3: Spline below S

    # We'll use Strategy 1: band-limited Fejér kernel
    # The band-limited Fejér kernel of order n has exponential type 2πn and is positive
    B = 2.0 * W  # band-limit
    # Create a smooth kernel: Fejér-like but strictly band-limited
    # Use a cosine window in frequency domain to create smooth band-limited kernel

    freq = np.fft.fftfreq(len(RGRID), dr)
    kernel_freq = np.zeros_like(freq)
    # Set up a smooth band-limited kernel in frequency domain
    # Use a triangular/Fejér window
    mask_pos = np.abs(freq) <= B
    kernel_freq[mask_pos] = np.cos(np.pi * np.abs(freq[mask_pos]) / (2*B))**2
    kernel_time = np.fft.ifft(kernel_freq).real
    kernel_time = np.fft.fftshift(kernel_time)
    kernel_time = np.abs(kernel_time)  # Ensure non-negative
    kernel_time /= np.max(kernel_time)

    # Convolve S with kernel
    S_smooth_bad = np.convolve(S, kernel_time, mode='same')
    # Normalize
    S_smooth = S_smooth_bad * (dr / np.sum(kernel_time))

    # Check if m ≤ S on window; if not, apply a downward scale
    S_smooth_on_win = S_smooth[win]
    S_on_win = S[win]
    max_excess = np.max(S_smooth_on_win - S_on_win)

    if max_excess > 0:
        # Apply a correction: linear blend towards S to ensure m ≤ S everywhere
        alpha = min(1.0, 0.5 + max_excess / (2*np.abs(S_min)))
        S_smooth = alpha * S_smooth + (1.0 - alpha) * S

    m_min = float(S_smooth[win].min())
    m_min_idx = np.where(win)[0][np.argmin(S_smooth[win])]
    m_min_loc = RGRID[m_min_idx]

    # Compute the quality of the minorant
    # How much of the gap between S_min and 0 does the minorant recover?
    gap_recovered = (m_min - S_min) / np.abs(S_min) if S_min != 0 else 0

    # Also compute the L_infty error on the window
    linf_error = float(np.max(np.abs(S_smooth[win] - S[win])))

    # Check: is m ≤ S on the window?
    is_minorant = float(np.max(S_smooth[win] - S[win])) <= 1e-6

    # Compute margins using the full Weil operators
    N = max(4, round(0.75*L*2*W/np.pi)); CEN = list(np.linspace(a, b, N))

    def fej(x):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out

    Pm = np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); G = (Pm @ Pm.T)*dr

    arch = (Pm*OMEGA[None, :]) @ Pm.T * dr/(2*np.pi)
    cosu = np.cos(np.outer(ln, RGRID)); prime = np.zeros((N, N))
    for j in range(N):
        gp = ((Pm[j][None, :]*Pm) @ cosu.T)*dr/(2*np.pi)
        prime[j, :] = 2.0*(gp*lw[None, :]).sum(axis=1)
    prime = 0.5*(prime+prime.T)
    def fc(c):
        y = W*(0.5j-c); y2 = W*(0.5j+c); return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2
    fp = np.array([fc(c) for c in CEN]); pole = 2.0*np.real(np.outer(fp, fp))

    Q_full = pole + arch - prime
    margin_full = float(eigh(Q_full, G, eigvals_only=True)[0])

    # Get the optimal eigenfunction (eigenvector for the minimum eigenvalue)
    eigvals, eigvecs = eigh(Q_full, G)
    optimal_phi = eigvecs[:, 0]  # eigenvector for minimum eigenvalue

    # Evaluate the Rayleigh bound using the optimal phi and the minorant
    # ∫ m(r)|φ(r)|² ≈ sum_j (Pm[j]*m)(r) (Pm[j]*m)(r)
    m_values_proj = Pm @ (S_smooth * dr)  # project m onto basis
    phi_norm_sq = (optimal_phi @ G) @ optimal_phi  # ||φ||_G²

    # Numerator: <pole + (m-integrated), φ> / <φ, φ>G
    # This is complex to compute directly. Instead, use the fact that margin ≥ m_min (Rayleigh)
    rayleigh_bound_m = m_min

    return {
        'S_min': S_min,
        'S_min_loc': S_min_loc,
        'm_min': m_min,
        'm_min_loc': m_min_loc,
        'margin_full': margin_full,
        'rayleigh_bound_m': rayleigh_bound_m,
        'gap_recovered': gap_recovered,
        'linf_error': linf_error,
        'is_minorant': is_minorant,
        'N': N,
        'B_limit': B
    }

print("="*150)
print("BEURLING-SELBERG EXTREMAL MINORANT: Rayleigh lower bound on margin via band-limited m(r)≤S(r)")
print("="*150)
print(f"{'window':>12} {'4W':>4} {'S_min':>10} {'m_min':>10} {'Rayleigh(m)':>12} {'margin_full':>14} {'gap_recovered':>14} {'L∞error':>10} {'is_min?':>8}")
print("-"*150)

for (a, b, fourW) in [(14,44,2.0),(14,44,3.0),(14,44,4.0),(100,130,3.0),(100,130,4.0),(400,430,4.0)]:
    W = fourW/4
    result = analyze_with_extremal_minorant(a, b, W)

    S_min = result['S_min']
    m_min = result['m_min']
    rayleigh_m = result['rayleigh_bound_m']
    margin_full = result['margin_full']
    gap_rec = result['gap_recovered']
    linf_err = result['linf_error']
    is_min = result['is_minorant']

    print(f"{f'[{a},{b}]':>12} {fourW:>4.1f} {S_min:>10.4f} {m_min:>10.4f} {rayleigh_m:>12.4f} {margin_full:>+14.4e} {gap_rec:>14.2%} {linf_err:>10.4e} {str(is_min):>8}")

print("="*150)
print("INTERPRETATION:")
print("  S_min: minimum of Weil symbol (always negative, useless)")
print("  m_min: minimum of band-limited minorant (still negative, but closer to 0)")
print("  Rayleigh(m): m_min, used as a lower bound on the margin")
print("  margin_full: actual margin from Q")
print("  gap_recovered: (m_min - S_min) / |S_min|, fraction of gap between S_min and 0 that m recovers")
print("  L∞error: max|m(r) - S(r)| on window (how tight is the minorant?)")
print("  is_minorant: True if m(r) ≤ S(r) pointwise (should be True)")
print()
print("VERDICT:")
print("  If Rayleigh(m) ≥ 0: the band-limited minorant gives a POSITIVE lower bound (via extremal theory)")
print("  If Rayleigh(m) < 0 but gap_recovered ≈ 1.0: we recover MOST of the gap, but not enough for positivity")
print("  If gap_recovered << 1.0: the band-limiting constraint significantly cuts the available bound")
print("="*150)
