#!/usr/bin/env python3
"""Berry-Keating periodic-orbit gate: reconstruct the ACTUAL zero positions from PRIMES ONLY.

The explicit-formula density identity (a distributional identity, smoothed by a Gaussian of width σ):
  Σ_γ g_σ(T-γ)  =  d_smooth(T)  +  d_osc(T)
with
  d_smooth(T) = (1/2π) log(T/2π)                                  [Weyl / prolate semiclassical term]
  d_osc(T)    = -(1/π) Σ_p Σ_{k≥1} (log p / p^{k/2}) cos(kT log p) e^{-(σ k log p)²/2}
                                                                   [prime periodic orbits, Gaussian-damped]
  g_σ(u)      = (1/(σ√2π)) exp(-u²/(2σ²))                          [smoothing kernel]

LHS uses ONLY the zeros; RHS uses ONLY primes (+ the smooth Weyl term, no zeros). If they agree,
primes reconstruct the zeros — the numerical content of "V_prime realizes the spectrum". NON-CIRCULAR:
no zero appears on the RHS. Gaussian smoothing makes the prime sum converge (only p^k ≲ e^{1/σ} matter).
FALSIFIABLE: report the pointwise match + peak-recovery (does RHS peak at each real γ?) + correlation."""
import numpy as np
from sympy import primerange

Z = np.load('zeros10k.npy')            # 10000 Riemann zeros, γ up to ~9877

def d_smooth(T):
    return np.log(T/(2*np.pi))/(2*np.pi)

def d_osc_from_primes(T, sigma, pmax):
    """oscillatory density from primes only. T is an array."""
    out = np.zeros_like(T)
    primes = list(primerange(2, pmax+1))
    for p in primes:
        lp = np.log(p)
        for k in range(1, 40):
            damp = np.exp(-(sigma*k*lp)**2/2)
            if damp < 1e-14:
                break
            amp = lp / p**(k/2)
            out += -(1/np.pi) * amp * damp * np.cos(k*T*lp)
    return out

def zero_side(T, sigma, Zwin):
    """smoothed delta-comb over the actual zeros in window."""
    out = np.zeros_like(T)
    c = 1.0/(sigma*np.sqrt(2*np.pi))
    for g in Zwin:
        out += c*np.exp(-(T-g)**2/(2*sigma**2))
    return out

# --- test window: zeros around T~1000 (spacing ~1.2, need sigma<spacing to resolve) ---
for (Tlo, Thi, sigma) in [(1000, 1040, 0.18), (5000, 5040, 0.14), (9000, 9040, 0.12)]:
    pmax = int(np.exp(6.5/sigma))          # include primes to where damping ~ e^{-(6.5)²/2}~2e-10
    pmax = min(pmax, 200000)
    Zwin = Z[(Z > Tlo-8) & (Z < Thi+8)]
    T = np.linspace(Tlo, Thi, 4000)
    lhs = zero_side(T, sigma, Zwin)
    rhs = d_smooth(T) + d_osc_from_primes(T, sigma, pmax)
    # correlation of the FLUCTUATION (subtract smooth mean)
    lf = lhs - d_smooth(T); rf = rhs - d_smooth(T)
    r = np.corrcoef(lf, rf)[0,1]
    # peak recovery: find local maxima of RHS, match to nearest real zero
    from scipy.signal import find_peaks
    pk,_ = find_peaks(rhs, height=d_smooth((Tlo+Thi)/2)*1.5)
    Tpk = T[pk]
    matched = 0; errs=[]
    zin = Zwin[(Zwin>Tlo)&(Zwin<Thi)]
    for g in zin:
        if len(Tpk):
            j = np.argmin(np.abs(Tpk-g))
            e = abs(Tpk[j]-g)
            if e < sigma:
                matched += 1; errs.append(e)
    rms = np.sqrt(np.mean(np.array(errs)**2)) if errs else np.nan
    print(f"T∈[{Tlo},{Thi}] σ={sigma} primes≤{pmax}: {len(zin)} zeros | "
          f"corr(fluct)={r:.4f} | peaks matched {matched}/{len(zin)} | RMS pos err={rms:.4f} "
          f"(spacing≈{2*np.pi/np.log((Tlo+Thi)/2/(2*np.pi)):.3f})")
