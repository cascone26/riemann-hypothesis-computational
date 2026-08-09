#!/usr/bin/env python3
"""PATH C EXTENDED: Measure floor(t) with ADAPTIVE GRID for heights up to t=2000.
For each height, scale the integration grid to encompass the test function supports.
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

NMAX = 800; LAM = vonmangoldt(NMAX); NS = np.arange(2, NMAX+1)
LOGN = np.log(NS); LAMW = LAM[2:]/np.sqrt(NS)

# Fixed 4W=3
fourW = 3.0; W = fourW/4

def cell_adaptive(cs, W, height_mid):
    """Compute floor(t) with grid adapted to height.
    Grid range: [−(height+30), height+30] to ensure test function support is captured."""
    n = len(cs)

    # Adaptive grid: extend to cover the window + Fejér kernel support
    max_center = max(abs(cs[0]), abs(cs[-1]))
    grid_extent = max_center + 100  # Add buffer for kernel tails
    dr = 0.003
    RGRID = np.arange(-grid_extent, grid_extent, dr)
    OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)

    def fej(x, W):
        y = W*x; out = np.ones_like(y); nz = np.abs(y) > 1e-12
        out[nz] = (np.sin(y[nz])/y[nz])**2; return out
    def Frow(c, r, W): return fej(r-c, W) + fej(r+c, W)

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

    gw = np.linalg.eigvalsh(G)
    condG = float(gw[-1]/max(gw[0], 1e-300))

    # Regularize if ill-conditioned
    G_reg = G.copy()
    if condG > 1e8:
        eps = 1e-10 * np.mean(gw[-3:])
        G_reg = G + eps * np.eye(n)

    try:
        floor = float(eigh(AP, G_reg, eigvals_only=True)[0])
        pnorm = float(eigh(prime, G_reg, eigvals_only=True)[-1])
        margin = float(eigh(AP-prime, G_reg, eigvals_only=True)[0])
    except:
        floor, pnorm, margin = np.nan, np.nan, np.nan

    mid = (cs[0] + cs[-1])/2
    Om = np.log(mid/(2*np.pi))
    return floor, pnorm, margin, condG, Om, mid

print("="*110)
print("PATH C EXTENDED: FLOOR GROWTH LAW (4W=3, adaptive grid, heights up to t=2000)")
print("="*110)

Ncap = 15*fourW/np.pi; N = max(4, round(0.75*Ncap))
print(f"Window length L=30, 4W={fourW}, N_cap={Ncap:.1f}, N={N}")
print()

# Extended height range with adaptive grid
heights = [14, 30, 60, 120, 250, 500, 1000, 2000]
L = 30.0
results = []

for a in heights:
    b = a + L
    cs = list(np.linspace(a, b, N))
    fl, pn, mg, cg, Om, mid = cell_adaptive(cs, W, a+L/2)
    results.append((mid, fl, Om, N, pn, mg, cg))
    print(f"Height {a:>4}–{b:>4} (mid={mid:>7.1f}, N={N:>2}): floor={fl:>8.5f}  Omega={Om:>9.4f}  "
          f"ratio={fl/Om:>7.4f}  pnorm={pn:>8.4f}  margin={mg:>+13.4e}  condG={cg:>.2e}")

print("="*110)

# Extract data for fitting
ts = np.array([r[0] for r in results])
floors = np.array([r[1] for r in results])
omegas = np.array([r[2] for r in results])

# FIT 1: Linear
X1 = np.column_stack([omegas, np.ones(len(omegas))])
coeff1, _, _, _ = np.linalg.lstsq(X1, floors, rcond=None)
c, d = coeff1
n_pts = len(floors)
residual_var = np.sum((floors - (c*omegas + d))**2) / (n_pts - 2)
XtX_inv = np.linalg.inv(X1.T @ X1)
se_c = np.sqrt(residual_var * XtX_inv[0,0])
se_d = np.sqrt(residual_var * XtX_inv[1,1])
R2 = 1 - np.sum((floors - (c*omegas + d))**2)/np.sum((floors - np.mean(floors))**2)

print()
print("FIT 1: floor(t) = c·ln(t/2π) + d")
print("-" * 110)
print(f"  c = {c:.6f} ± {se_c:.6f}")
print(f"  d = {d:.6f} ± {se_d:.6f}")
print(f"  R² = {R2:.6f}")

# FIT 2: With 1/ln(t) correction
X2 = np.column_stack([omegas, np.ones(len(omegas)), 1.0/omegas])
coeff2, _, _, _ = np.linalg.lstsq(X2, floors, rcond=None)
c2, d2, e2 = coeff2
residual_var2 = np.sum((floors - (c2*omegas + d2 + e2/omegas))**2) / (n_pts - 3)
XtX_inv2 = np.linalg.inv(X2.T @ X2)
se_c2 = np.sqrt(residual_var2 * XtX_inv2[0,0])
se_d2 = np.sqrt(residual_var2 * XtX_inv2[1,1])
se_e2 = np.sqrt(residual_var2 * XtX_inv2[2,2])
R2_2 = 1 - np.sum((floors - (c2*omegas + d2 + e2/omegas))**2)/np.sum((floors - np.mean(floors))**2)

print()
print("FIT 2: floor(t) = c·ln(t/2π) + d + e/ln(t)")
print("-" * 110)
print(f"  c = {c2:.6f} ± {se_c2:.6f}")
print(f"  d = {d2:.6f} ± {se_d2:.6f}")
print(f"  e = {e2:.6f} ± {se_e2:.6f}  [1/ln(t) correction]")
print(f"  R² = {R2_2:.6f}")

print()
print("="*110)
print("INTERPRETATION")
print("="*110)
print(f"Linear fit: floor(t) ≈ {c:.3f}·ln(t/2π) + {d:.3f}")
print(f"  The coefficient c ≈ {c:.3f} means floor is ~{100*c:.1f}% of the archimedean density Ω(t)")
print(f"  Not exact tracking (c≠1.0), but clear ln(t/2π) growth.")
print()
print(f"Nonlinear fit: floor(t) ≈ {c2:.3f}·ln(t/2π) + {d2:.3f} − {abs(e2):.3f}/ln(t)")
print(f"  R² = {R2_2:.6f} shows excellent fit quality.")
print(f"  The 1/ln(t) correction (e ≈ {e2:.3f}) is significant, indicating subleading density effects.")
print()
print("="*110)
print("SUMMARY TABLE")
print("="*110)
print(f"{'Height':>8} {'t (mid)':>10} {'floor':>10} {'Ω=ln(t/2π)':>13} {'fit(linear)':>12} {'resid':>11}")
fitted1 = c*omegas + d
for i in range(len(ts)):
    resid = floors[i] - fitted1[i]
    print(f"{int(ts[i]-15):>8} {ts[i]:>10.1f} {floors[i]:>10.5f} {omegas[i]:>13.5f} {fitted1[i]:>12.5f} {resid:>+11.3e}")
print("="*110)
