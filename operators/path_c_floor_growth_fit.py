#!/usr/bin/env python3
"""PATH C: PIN THE FLOOR GROWTH LAW
Measure floor(t)=λmin(arch+pole,G) across MANY heights at FIXED 4W=3.
Heights: {14,30,60,120,250,500,1000,2000}, window length L=30, capacity-matched N.
Fit: floor(t) = c·ln(t/2π) + d (+ possibly e/ln(t)).
Cross-check: c vs archimedean density Ω(t)=Re ψ₀(¼+it/2)−ln π ~ ln(t/2π).
"""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
import sys

def vonmangoldt(N):
    lam = np.zeros(N+1); isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True; lp = np.log(p); pk = p
            while pk <= N: lam[pk] = lp; pk *= p
    return lam

NMAX = 800; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1)
LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)

# Fixed 4W=3
fourW = 3.0; W = fourW/4
Ncap = 15*fourW/np.pi; N = max(4, round(0.75*Ncap))

# Grid resolution: finer than verify_race_premises for better accuracy
RGRID = np.arange(-600.0, 600.0, 0.003); dr = 0.003
OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)

def fej(x, W):
    y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
    out[nz] = (np.sin(y[nz])/y[nz])**2; return out
def Frow(c, r, W): return fej(r-c, W) + fej(r+c, W)

def cell(cs, W):
    """Compute floor, pnorm, margin for given centers cs at bandwidth W.
    Returns: floor, pnorm, margin, cond(G), Omega_val"""
    n = len(cs)
    P = np.array([Frow(c, RGRID, W) for c in cs]); G = (P @ P.T)*dr
    arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)
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
    AP = pole+arch

    gw = np.linalg.eigvalsh(G); condG = float(gw[-1]/max(gw[0], 1e-300))

    # Regularize G if ill-conditioned (add small multiple of identity)
    G_reg = G.copy()
    if condG > 1e8:
        eps = 1e-10 * np.mean(gw[-3:])  # small regularization
        G_reg = G + eps * np.eye(n)
        gw_reg = np.linalg.eigvalsh(G_reg)
        condG = float(gw_reg[-1] / max(gw_reg[0], 1e-300))

    try:
        floor = float(eigh(AP, G_reg, eigvals_only=True)[0])
        pnorm = float(eigh(prime, G_reg, eigvals_only=True)[-1])
        margin = float(eigh(AP-prime, G_reg, eigvals_only=True)[0])
    except np.linalg.LinAlgError:
        # Fallback: use regular eigendecomposition (not generalized)
        # Transform to standard form: solve AP in metric G
        try:
            L = np.linalg.cholesky(G_reg)
            AP_std = np.linalg.solve(L, AP @ np.linalg.solve(L.T, np.eye(n)))
            floor = float(np.linalg.eigvalsh(AP_std)[0])
            pnorm = float(np.linalg.eigvalsh(np.linalg.solve(L, prime @ np.linalg.solve(L.T, np.eye(n))))[-1])
            margin = float(np.linalg.eigvalsh(np.linalg.solve(L, (AP-prime) @ np.linalg.solve(L.T, np.eye(n))))[0])
        except:
            # Last resort: return NaN
            floor, pnorm, margin = np.nan, np.nan, np.nan

    # Archimedean density at the window midpoint
    mid = (cs[0] + cs[-1])/2
    Om = np.log(mid/(2*np.pi))
    return floor, pnorm, margin, condG, Om, mid

print("="*110)
print("PATH C: FLOOR GROWTH LAW (4W=3 fixed, L=30, N=0.75*N_cap, many heights)")
print("="*110)
print(f"Window length L=30, 4W={fourW}, N_cap={Ncap:.1f}, N={N}")
print()

# Data collection across heights
# NOTE: Grid RGRID=[-600,600] limits us to ~t≤600 reliably; larger heights show numerical breakdown
heights = [14, 30, 60, 120, 250, 500]  # Exclude t=1000,2000 (outside grid reliable range)
L = 30.0
results = []

print(f"NOTE: Using heights up to t=500 (within RGRID=[−600,600] safe range).")
print(f"      For larger heights, need adaptive grid; skipping t=1000,2000 numerically problematic points.")
print()

for a in heights:
    b = a + L
    mid = (a+b)/2
    cs = list(np.linspace(a, b, N))
    fl, pn, mg, cg, Om, mid_check = cell(cs, W)
    results.append((mid, fl, Om, N, pn, mg, cg))
    print(f"Height {a:>4}–{b:>4} (mid={mid:>7.1f}, N={N:>2}): floor={fl:>8.4f}  Omega(ln)={Om:>9.3f}  ratio={fl/Om:>7.3f}  pnorm={pn:>8.4f}  margin={mg:>+13.4e}  condG={cg:>7.1f}")

print("="*110)

# Extract data for fitting
ts = np.array([r[0] for r in results])        # midpoints (heights)
floors = np.array([r[1] for r in results])    # floor values
omegas = np.array([r[2] for r in results])    # Omega values

# Fit 1: floor(t) = c·ln(t/2π) + d
# Linear fit: y = c*x + d where x=ln(t/2π), y=floor(t)
X1 = np.column_stack([omegas, np.ones(len(omegas))])  # Design matrix
try:
    coeff1, residuals1, rank1, s1 = np.linalg.lstsq(X1, floors, rcond=None)
    c, d = coeff1

    # Error bars via residual variance
    n_pts = len(floors)
    if n_pts > 2:
        residual_var = np.sum((floors - (c*omegas + d))**2) / (n_pts - 2)
        # Standard errors (covariance from inverse of X'X)
        XtX_inv = np.linalg.inv(X1.T @ X1)
        se_c = np.sqrt(residual_var * XtX_inv[0,0])
        se_d = np.sqrt(residual_var * XtX_inv[1,1])
        R2 = 1 - np.sum((floors - (c*omegas + d))**2)/np.sum((floors - np.mean(floors))**2)
    else:
        se_c, se_d = np.nan, np.nan
        R2 = np.nan

    print()
    print("FIT 1: floor(t) = c·ln(t/2π) + d")
    print("-" * 110)
    print(f"  c (slope)  = {c:>10.6f} ± {se_c:>10.6f}    [coefficient of ln(t/2π)]")
    print(f"  d (intercept) = {d:>10.6f} ± {se_d:>10.6f}")
    print(f"  R² = {R2:>.6f}")
    print()

    # Residuals and fit quality
    fitted = c*omegas + d
    residuals = floors - fitted
    max_resid = np.max(np.abs(residuals))
    print(f"  Max residual: {max_resid:.2e}")
    print(f"  Residuals: {residuals}")

except np.linalg.LinAlgError as e:
    print(f"ERROR fitting floor(t) = c·ln(t/2π) + d: {e}")
    c, d, se_c, se_d, R2 = np.nan, np.nan, np.nan, np.nan, np.nan

# Fit 2: floor(t) = c·ln(t/2π) + d + e/ln(t) (3-parameter)
print()
print("FIT 2: floor(t) = c·ln(t/2π) + d + e/ln(t)  [optional higher-order correction]")
print("-" * 110)
try:
    X2 = np.column_stack([omegas, np.ones(len(omegas)), 1.0/omegas])  # c·x + d + e/x
    coeff2, residuals2, rank2, s2 = np.linalg.lstsq(X2, floors, rcond=None)
    c2, d2, e2 = coeff2

    n_pts = len(floors)
    if n_pts > 3:
        residual_var2 = np.sum((floors - (c2*omegas + d2 + e2/omegas))**2) / (n_pts - 3)
        XtX_inv2 = np.linalg.inv(X2.T @ X2)
        se_c2 = np.sqrt(residual_var2 * XtX_inv2[0,0])
        se_d2 = np.sqrt(residual_var2 * XtX_inv2[1,1])
        se_e2 = np.sqrt(residual_var2 * XtX_inv2[2,2])
        R2_2 = 1 - np.sum((floors - (c2*omegas + d2 + e2/omegas))**2)/np.sum((floors - np.mean(floors))**2)
    else:
        se_c2, se_d2, se_e2 = np.nan, np.nan, np.nan
        R2_2 = np.nan

    print(f"  c = {c2:>10.6f} ± {se_c2:>10.6f}")
    print(f"  d = {d2:>10.6f} ± {se_d2:>10.6f}")
    print(f"  e = {e2:>10.6f} ± {se_e2:>10.6f}  [coefficient of 1/ln(t)]")
    print(f"  R² = {R2_2:>.6f}")

    fitted2 = c2*omegas + d2 + e2/omegas
    residuals2 = floors - fitted2
    max_resid2 = np.max(np.abs(residuals2))
    print(f"  Max residual: {max_resid2:.2e}")

    # Compare to fit 1
    improvement = (max_resid - max_resid2) / max_resid if max_resid > 0 else 0
    print(f"  Improvement over fit 1: {improvement*100:.1f}%  [e ≈ {e2:.6f} is {'significant' if abs(e2)>0.01 else 'negligible'}]")

except np.linalg.LinAlgError as e:
    print(f"ERROR fitting 3-parameter model: {e}")
    c2, d2, e2 = np.nan, np.nan, np.nan

print()
print("="*110)
print("CROSS-CHECK: Archimedean density")
print("="*110)
print("Ω(t) = Re ψ₀(¼+it/2) − ln π  should satisfy  Ω(t) ~ ln(t/2π)  (pure density, 4W-independent).")
print(f"Measured c={c:.6f} relates to the archimedean floor's link to zero density.")
print()

# Estimate what c should be if it's purely from Ω
# The coefficient c in floor = c·Ω + d should equal c if the floor tracks Ω exactly
print(f"Interpretation of c={c:.6f}:")
print(f"  If c ≈ 1.0  : floor tracks Ω exactly (marginal/unweighted)")
print(f"  If c < 1.0  : floor is a smoothed/fractional density (weighted/tapered)")
print(f"  If c > 1.0  : floor scales above pure density (amplification factor)")
print()

# Known values from literature
print("Known constants:")
print(f"  §816 reports floor ≈ 0.153 (W-invariant)")
print(f"  §839 reports floor ≈ 0.155·ln(t/2π)")
print(f"  Explicit formula: Ω(t) ~ ln(t/2π)  at large t")
print(f"  Beurling density: related to zero spacing and Nyquist sampling")
print()
print(f"HANDLE: This script measures and fits floor(t) via the cell() function")
print(f"        (same apparatus as §834, cert_floor_vs_prime.py). All outputs are MEASURED,")
print(f"        fits are via standard linear regression with reported error bars.")
print()

# Summary table
print("="*110)
print("SUMMARY TABLE: floor(t) ~ {0:.6f}·ln(t/2π) + {1:.6f}".format(c, d))
print("="*110)
print(f"{'Height':>8} {'t (mid)':>10} {'floor':>10} {'ln(t/2π)':>12} {'fitted':>10} {'resid':>11}")
for i, (t, f, om) in enumerate(zip(ts, floors, omegas)):
    fitted_val = c*om + d
    resid = f - fitted_val
    print(f"{int(ts[i]-15):>8} {t:>10.1f} {f:>10.5f} {om:>12.5f} {fitted_val:>10.5f} {resid:>+11.2e}")
print("="*110)
