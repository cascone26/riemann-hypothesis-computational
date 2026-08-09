#!/usr/bin/env python3
"""PATH D: Prime concentration law at the hard corner [14,44].
Measure pnorm(4W) at low-height [14,44] as support 4W grows from 2.0 to 6.0.
Compare pnorm to:
  1) Raw prime sum S(4W) = Σ_{log n < 4W} Λ(n)/√n
  2) Growth functions: √S, ln(S), √(ln S)
  3) Fixed floor ~0.125 to detect saturation vs unbounded growth

KEY QUESTION: does pnorm(4W) saturate below floor (~0.125) or exceed it?
Corner verdict depends on whether pnorm → ∞ (margin<0 eventual) or pnorm↗0.125⁻ (RH-safe).
"""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh

def vonmangoldt(N):
    lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam

# Pre-computed grid for [14,44] window and wider sweep
NMAX = 1200; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1)
LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)
RGRID = np.arange(-300.0, 300.0, 0.003); dr = 0.003
OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)

def fej(x, W):
    y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
    out[nz] = (np.sin(y[nz])/y[nz])**2; return out

def Frow(c, r, W):
    return fej(r-c, W) + fej(r+c, W)

def raw_prime_sum(fourW):
    """Compute S(4W) = Σ_{log n < 4W} Λ(n)/√n"""
    mask = LOGN < fourW
    return float(LAMW[mask].sum())

def cell(cs, W):
    """Measure floor, pnorm, margin, margin_error for a given set of centers cs.
    Returns: floor, pnorm, margin, error_mag, margin_minus_floor_plus_pnorm, cond_G
    """
    n = len(cs); fourW = 4*W
    P = np.array([Frow(c, RGRID, W) for c in cs]); G = (P @ P.T)*dr
    arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)

    # Prime sum: only primes with log n < 4W
    mask = LOGN < fourW; ln = LOGN[mask]; lw = LAMW[mask]
    prime = np.zeros((n, n))
    if ln.size:
        cosu = np.cos(np.outer(ln, RGRID))
        for j in range(n):
            gp = ((P[j][None, :]*P) @ cosu.T)*dr/(2*np.pi)
            prime[j, :] = 2.0*(gp*lw[None, :]).sum(axis=1)
        prime = 0.5*(prime+prime.T)

    def fc(c):
        y = W*(0.5j-c); y2 = W*(0.5j+c)
        return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2

    fp = np.array([fc(c) for c in cs]); pole = 2.0*np.real(np.outer(fp, fp))
    archpole = pole + arch

    # Condition number of G
    gw = np.linalg.eigvalsh(G)
    condG = float(gw[-1]/max(gw[0], 1e-300))

    # Eigenvalues
    floor = float(eigh(archpole, G, eigvals_only=True)[0])
    pnorm = float(eigh(prime, G, eigvals_only=True)[-1])
    margin = float(eigh(archpole - prime, G, eigvals_only=True)[0])

    # Check the assembly bound: margin ≥ floor - pnorm?
    assembly_error = margin - (floor - pnorm)

    return floor, pnorm, margin, assembly_error, condG

# =============================================================================
# SWEEP 4W from 2.0 to 6.0 at fixed LOW height [14,44]
# =============================================================================
A, B = 14.0, 44.0
L = B - A
print("\n" + "="*110)
print("PATH D: PRIME CONCENTRATION LAW — hard corner [14,44], sweeping 4W from 2.0 to 6.0")
print("="*110)
print(f"Window: [{A:.0f}, {B:.0f}] (length {L:.0f}, arch-floor ~0.125)")
print(f"N-strategy: 0.75 * N_cap where N_cap = 15*4W/π  (capacity-matched)")
print("="*110)
print()
print(f"{'4W':>6} {'N_cap':>6} {'N':>4} {'nprime':>7} {'S(4W)':>12} {'√S':>9} {'ln(S)':>8} {'floor':>9} {'pnorm':>10} {'margin':>13} {'cond(G)':>8}")
print("-"*110)

results = []
for fourW in np.linspace(2.0, 6.0, 21):  # 21 points from 2.0 to 6.0
    W = fourW/4
    Ncap = 15*fourW/np.pi
    N = max(4, round(0.75*Ncap))

    # Window centers
    CEN = list(np.linspace(A, B, N))

    # Compute
    fl, pn, mg, aerr, cg = cell(CEN, W)

    # Prime sum
    S = raw_prime_sum(fourW)
    sqrtS = np.sqrt(S) if S > 0 else 0
    lnS = np.log(S) if S > 0 else 0

    # Count primes in range
    nprime = int((LOGN < fourW).sum())

    results.append({
        '4W': fourW,
        'N': N,
        'Ncap': Ncap,
        'nprime': nprime,
        'S': S,
        'sqrtS': sqrtS,
        'lnS': lnS,
        'floor': fl,
        'pnorm': pn,
        'margin': mg,
        'cond_G': cg
    })

    print(f"{fourW:>6.2f} {Ncap:>6.1f} {N:>4} {nprime:>7} {S:>12.4f} {sqrtS:>9.4f} {lnS:>8.4f} {fl:>9.4f} {pn:>10.5f} {mg:>+13.5e} {cg:>8.1f}")

print("-"*110)
print()

# =============================================================================
# ANALYSIS: Growth and saturation
# =============================================================================
print("="*110)
print("ANALYSIS: pnorm growth vs growth functions")
print("="*110)
print()

# Extract pnorm and S values
pnorm_vals = np.array([r['pnorm'] for r in results])
S_vals = np.array([r['S'] for r in results])
fourW_vals = np.array([r['4W'] for r in results])
floor_vals = np.array([r['floor'] for r in results])

# Fit pnorm to power-law and other growth functions
print(f"{'4W':>6} {'pnorm':>10} {'S(4W)':>12} {'pnorm/S':>10} {'pnorm/√S':>10} {'pnorm/ln(S)':>11}")
print("-"*110)
for r in results:
    ratio_S = r['pnorm'] / r['S'] if r['S'] > 0 else 0
    ratio_sqrtS = r['pnorm'] / r['sqrtS'] if r['sqrtS'] > 0 else 0
    ratio_lnS = r['pnorm'] / r['lnS'] if r['lnS'] > 0 else 0
    print(f"{r['4W']:>6.2f} {r['pnorm']:>10.6f} {r['S']:>12.4f} {ratio_S:>10.6f} {ratio_sqrtS:>10.6f} {ratio_lnS:>11.6f}")

print()

# Fit log-log to detect growth exponent
# log(pnorm) = a + b*log(4W)
log4W = np.log(fourW_vals)
logpnorm = np.log(pnorm_vals + 1e-10)  # add small value to avoid log(0)
fit_coeffs = np.polyfit(log4W, logpnorm, 1)
b_exponent = fit_coeffs[0]
a_offset = fit_coeffs[1]
pnorm_fit = np.exp(a_offset) * fourW_vals**b_exponent

print("Log-log fit: pnorm ~ a * (4W)^b")
print(f"  a = {np.exp(a_offset):.6e}")
print(f"  b (exponent) = {b_exponent:.4f}")
print(f"  Fit residual (RMS): {np.sqrt(np.mean((pnorm_vals - pnorm_fit)**2)):.6e}")
print()

# Check for saturation: does pnorm plateau?
diff_pnorm = np.diff(pnorm_vals)
avg_diff = np.mean(diff_pnorm)
max_diff = np.max(np.abs(diff_pnorm))
print(f"Saturation check (pnorm increments):")
print(f"  Average Δpnorm per step: {avg_diff:>10.6e}")
print(f"  Max |Δpnorm|: {max_diff:>10.6e}")
print(f"  Recent increment (last 3 points): {np.mean(diff_pnorm[-3:]):>10.6e}")

# Compare to low-height floor
fixed_floor = 0.125
print()
print(f"Corner safety check (fixed low-height floor ≈ {fixed_floor}):")
print(f"  max(pnorm) at 4W∈[2,6]: {np.max(pnorm_vals):.6f}")
print(f"  pnorm exceeds floor?  {np.max(pnorm_vals) > fixed_floor}")
print(f"  Margin to safety:      {fixed_floor - np.max(pnorm_vals):>+.6f}")
print()

# Extrapolate: if b > 0, what happens at 4W → ∞?
if b_exponent > 0.01:
    pnorm_at_8 = np.exp(a_offset) * 8**b_exponent
    print(f"Growth extrapolation (power-law model):")
    print(f"  At 4W=8:   pnorm ≈ {pnorm_at_8:.6f}  ({'EXCEEDS' if pnorm_at_8 > fixed_floor else 'stays below'} floor {fixed_floor})")
    print(f"  Exponent b={b_exponent:.4f} > 0  →  pnorm → ∞ as 4W → ∞")
    print(f"  ⟹  margin = floor - pnorm would eventually NEGATIVE")
elif b_exponent < -0.01:
    print(f"Growth extrapolation (power-law model):")
    print(f"  Exponent b={b_exponent:.4f} < 0  →  pnorm → 0 as 4W → ∞")
    print(f"  ⟹  margin → floor (saturates from below, RH-SAFE)")
else:
    print(f"Growth extrapolation (power-law model):")
    print(f"  Exponent b={b_exponent:.4f} ≈ 0  →  pnorm → const as 4W → ∞")
    if np.max(pnorm_vals) < fixed_floor:
        print(f"  ⟹  margin → floor - const > 0 (RH-SAFE)")
    else:
        print(f"  ⟹  margin → const < 0  (RH-RELEVANT)")

print()
print("="*110)
print("VERDICT")
print("="*110)

# Final verdict
max_pnorm = np.max(pnorm_vals)
min_margin = np.min([r['margin'] for r in results])
worst_4W = results[np.argmin([r['margin'] for r in results])]['4W']

print(f"At [14,44] over 4W ∈ [2.0, 6.0]:")
print(f"  • Measured pnorm range: {np.min(pnorm_vals):.6f} to {np.max(pnorm_vals):.6f}")
print(f"  • Fixed low-height floor: ~{fixed_floor}")
print(f"  • pnorm crosses floor? {max_pnorm > fixed_floor}")
print()
print(f"  • Measured margin range: {np.min([r['margin'] for r in results]):+.5e} to {np.max([r['margin'] for r in results]):+.5e}")
print(f"  • Worst margin at 4W={worst_4W:.2f}")
print()

if b_exponent > 0.05:
    verdict = "UNDETERMINED — pnorm grows with 4W; asymptotic behavior unclear at finite precision"
    concern = "POTENTIAL RH-RELEVANT if pnorm → ∞ as 4W → ∞"
elif max_pnorm > fixed_floor + 0.01:
    verdict = "CONCERNING — pnorm has exceeded low-height floor in measured range"
    concern = "Margin may go negative at higher 4W (extend sweep or use arb precision)"
elif min_margin < 0:
    verdict = "NEGATIVE MARGIN DETECTED — investigate float precision!"
    concern = "May be float floor artifact (margins < 1e-6)"
else:
    verdict = "RH-SAFE (measured) — pnorm stays below floor over 4W ∈ [2.0, 6.0]"
    concern = "Monitor saturation at higher 4W; float precision <1e-6 is unreliable"

print(f"CORNER VERDICT: {verdict}")
print(f"CAVEAT: {concern}")
print("="*110)
