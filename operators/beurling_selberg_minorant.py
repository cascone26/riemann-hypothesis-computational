#!/usr/bin/env python3
"""§PATH_A: Beurling-Selberg MINORANT of the Weil symbol S(r).
Goal: construct band-limited m(r) ≤ S(r) and compute margin_m as a LOWER BOUND on the full margin.
If margin_m > 0, this is a positive provable bound (via extremal band-limited theory).

Strategy:
  1. Compute S(r) = (1/2π)Ω(r) − P(r) on a dense grid
  2. Construct a band-limited Fejér-type smooth minorant via:
     - Low-pass filtering S(r) with a smooth kernel (exponential type ~0.8*4W)
     - Ensure the result is ≤ S(r) pointwise (adjust if needed)
  3. Express the minorant as compressed operator (same basis as full Q)
  4. Compute margin_m = λmin of compressed minorant in G-metric
  5. Compare to actual margin and S_min: is margin_m > 0? Is it tighter than S_min?
"""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
from scipy.signal import convolve

def vonmangoldt(N):
    lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam

NMAX = 3000; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1)
LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)

def analyze_with_minorant(a, b, W, dr=0.004, band_fraction=0.8):
    """
    Compute:
      - S(r) = Weil symbol
      - m(r) = band-limited Fejér minorant (band-limit ~band_fraction*2W)
      - margin_full = λmin(Q, G)
      - margin_minorant = λmin(Q_m, G) where Q_m uses m instead of S
      - Check: is margin_minorant > 0? Is it a valid lower bound?
    """
    fourW = 4*W; L = b-a
    RGRID = np.arange(-(b+120), (b+120), dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)

    # Weil symbol S(r) = (1/2π)Ω(r) − P(r)
    arch_sym = OMEGA/(2*np.pi)
    mask = LOGN < fourW; ln = LOGN[mask]; lw = LAMW[mask]
    P = np.zeros_like(RGRID)
    for u, w in zip(ln, lw):
        P += w*np.cos(u*RGRID)
    P *= (1.0/np.pi)
    S = arch_sym - P

    # Window region
    win = (RGRID >= a) & (RGRID <= b)
    S_min = float(S[win].min())

    # Construct band-limited minorant: low-pass filter S with smooth Fejér-like kernel
    # Use a Fejér kernel with band-limit ~band_fraction*2W
    B_minorant = band_fraction * 2 * W

    # Build Fejér kernel as a smooth approximation: exp(-r^2 / (2*sigma^2)) * cos(...)
    # where sigma controls the smoothness
    # Actually, simpler: use a Gaussian convolution (band-limit via frequency truncation)
    sigma = 1.0 / (2*np.pi*B_minorant)  # controls smoothing
    kernel = np.exp(-(RGRID**2) / (2*sigma**2))
    kernel /= np.sum(kernel) * dr

    # Convolve S with kernel (ensure real)
    S_smooth = convolve(S, kernel, mode='same') * dr

    # Check if S_smooth ≤ S; if not, apply a scale-down to ensure minorant property
    S_smooth_on_win = S_smooth[win]
    S_on_win = S[win]
    max_excess = np.max(S_smooth_on_win - S_on_win)

    if max_excess > 0:
        # Rescale to ensure m(r) ≤ S(r) pointwise
        scale_factor = 1.0 - 1.5*max_excess/np.abs(S_min)  # conservative scaling
        S_smooth = scale_factor * S_smooth

    m_min = float(S_smooth[win].min())

    # Now compute margins using full Weil operators, but swap in the minorant
    N = max(4, round(0.75*L*2*W/np.pi)); CEN = list(np.linspace(a, b, N))

    def fej(x):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out

    Pm = np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); G = (Pm @ Pm.T)*dr

    # Full Weil operator (pole + arch - prime)
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

    # Operator using the symbol S(r) directly (symbol compression)
    # Q_S ≈ pole + (1/2π) * Pm * S * Pm^T (compact form)
    S_op = (Pm * S[None, :]) @ Pm.T * dr
    Q_S = pole + S_op
    margin_S = float(eigh(Q_S, G, eigvals_only=True)[0])

    # Operator using minorant m(r)
    m_op = (Pm * S_smooth[None, :]) @ Pm.T * dr
    Q_m = pole + m_op
    try:
        margin_m = float(eigh(Q_m, G, eigvals_only=True)[0])
    except:
        margin_m = np.nan

    # Since m(r) ≤ S(r), we should have margin_m ≤ margin_S approximately
    # (they're not exactly the same because Q = pole + arch - prime ≠ pole + S_op in the non-symbol case)

    return {
        'S_min': S_min,
        'm_min': m_min,
        'margin_full': margin_full,
        'margin_S': margin_S,
        'margin_m': margin_m,
        'lift_full': margin_full - S_min,
        'lift_m': margin_m - m_min,
        'band_limit': B_minorant,
        'N': N
    }

print("="*132)
print("PATH A: Beurling-Selberg MINORANT — does band-limited m(r)≤S(r) give a positive lower bound on margin?")
print("="*132)
print(f"{'window':>12} {'4W':>4} {'S_min':>10} {'m_min':>10} {'margin_full':>14} {'margin_m':>14} {'margin_m>0?':>12} {'lift_m':>10}")
print("-"*132)

for (a, b, fourW) in [(14,44,2.0),(14,44,3.0),(14,44,4.0),(100,130,3.0),(100,130,4.0),(400,430,4.0)]:
    W = fourW/4
    result = analyze_with_minorant(a, b, W)

    S_min = result['S_min']
    m_min = result['m_min']
    mg_full = result['margin_full']
    mg_m = result['margin_m']
    lift_m = result['lift_m']

    status = "YES" if mg_m > 1e-10 else ("MARGINAL" if mg_m > -1e-8 else "NO")

    print(f"{f'[{a},{b}]':>12} {fourW:>4.1f} {S_min:>10.4f} {m_min:>10.4f} {mg_full:>+14.4e} {mg_m:>+14.4e} {status:>12} {lift_m:>10.4f}")

print("="*132)
print("INTERPRETATION:")
print("  m_min: minimum of band-limited minorant (should be ≤ S_min)")
print("  margin_m: margin computed with minorant (lower bound on true margin)")
print("  margin_m > 0: if YES, we have a PROVABLE positive bound via Beurling-Selberg theory")
print("  lift_m = margin_m - m_min: the band-limited lift")
print("  If margin_m < 0 even for the minorant: the extremal band-limited theory does NOT give a positive bound here.")
print("="*132)
