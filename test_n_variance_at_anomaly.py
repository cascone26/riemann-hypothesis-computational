#!/usr/bin/env python3
"""
Test: is the γ=124.3 anomaly (cos≈0.994) an N-dependent artifact?
Measure cos at this height with N ranging from 6 to 16 to see how misalignment varies.
"""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
import sys

NMAX = 400
def vonmangoldt(N):
    lam = np.zeros(N+1)
    isc = np.zeros(N+1, dtype=bool)
    for p in range(2, N+1):
        if not isc[p]:
            isc[2*p::p] = True
            lp = np.log(p)
            pk = p
            while pk <= N:
                lam[pk] = lp
                pk *= p
    return lam

LAM = vonmangoldt(NMAX)
NS = np.arange(2, NMAX+1)
LOGN = np.log(NS)
LAMW = LAM[2:]/np.sqrt(NS)
RGRID = np.arange(-320.0, 320.0, 0.0025)
dr = 0.0025
OMEGA = np.real(digamma(0.25 + 0.5j*RGRID)) - np.log(np.pi)

def fej(x, W):
    y = W*x
    out = np.ones_like(y)
    nz = np.abs(y) > 1e-12
    out[nz] = (np.sin(y[nz])/y[nz])**2
    return out

def Frow(c, r, W):
    return fej(r-c, W) + fej(r+c, W)

def mats(cs, W):
    n = len(cs)
    fourW = 4*W
    P = np.array([Frow(c, RGRID, W) for c in cs])
    G = (P @ P.T)*dr
    arch = (P*OMEGA[None, :]) @ P.T * dr/(2*np.pi)

    mask = LOGN < fourW
    ln = LOGN[mask]
    lw = LAMW[mask]
    prime = np.zeros((n, n))
    if ln.size:
        cosu = np.cos(np.outer(ln, RGRID))
        for j in range(n):
            gp = ((P[j][None, :]*P) @ cosu.T)*dr/(2*np.pi)
            prime[j, :] = 2.0*(gp*lw[None, :]).sum(axis=1)
        prime = 0.5*(prime+prime.T)

    def fc(c):
        y = W*(0.5j-c)
        y2 = W*(0.5j+c)
        return (np.sin(y)/y)**2 + (np.sin(y2)/y2)**2
    fp = np.array([fc(c) for c in cs])
    pole = 2.0*np.real(np.outer(fp, fp))

    return pole+arch, prime, G

gamma = 124.3
window_length = 30
height_scale = max(1.0, np.log(max(gamma, 30) / 30.0))
fourW_target = 2.8 * (1.0 + height_scale * 0.3)
W = fourW_target / 4.0
Ncap = 30 * 2 * W / np.pi

hw = window_length / 2.0
a = gamma - hw
b = gamma + hw

print("="*90)
print(f"N-VARIANCE TEST at γ={gamma} (window [{a:.1f}, {b:.1f}], 4W={fourW_target:.4f}, Ncap={Ncap:.1f})")
print("="*90)
print()
print(f"{'N':>3} {'N/Ncap':>7} {'floor':>8} {'pnorm':>8} {'margin':>12} {'cos(v_p,v_a)':>12} {'angle°':>7} {'condG':>8}")
print("-"*90)

results = []
for N in range(6, 17):
    try:
        CEN = list(np.linspace(a, b, N))
        AP, PR, G = mats(CEN, W)

        gw = np.linalg.eigvalsh(G)
        condG = float(gw[-1] / max(gw[0], 1e-300))

        wf, Vf = eigh(AP, G)
        floor = float(wf[0])
        v_a = Vf[:, 0]

        wp, Vp = eigh(PR, G)
        pnorm = float(wp[-1])
        v_p = Vp[:, -1]

        wm, _ = eigh(AP - PR, G)
        margin = float(wm[0])

        cosang = abs(float(v_a @ (G @ v_p)))
        ang = np.degrees(np.arccos(min(1.0, cosang)))

        results.append({'N': N, 'cos': cosang, 'angle': ang, 'margin': margin})

        print(f"{N:3d} {N/Ncap:7.2f} {floor:8.4f} {pnorm:8.4f} {margin:+12.3e} {cosang:12.6f} {ang:7.2f} {condG:8.1f}")
    except Exception as e:
        print(f"{N:3d} ERROR: {e}", file=sys.stderr)

print("-"*90)
print()

# Analysis
cos_vals = [r['cos'] for r in results]
angle_vals = [r['angle'] for r in results]

print("Statistics across N-variance:")
print(f"  cos(v_p,v_a) ranges from {min(cos_vals):.6f} to {max(cos_vals):.6f}")
print(f"  Mean cos: {np.mean(cos_vals):.6f}")
print(f"  Std dev:  {np.std(cos_vals):.6f}")
print()
print(f"  Angle ranges from {min(angle_vals):.2f}° to {max(angle_vals):.2f}°")
print(f"  Deviation from 90°: {abs(90 - np.mean(angle_vals)):.2f}° (mean)")
print()

# Find N that gives minimal cos (best orthogonality)
best_N = results[np.argmin(cos_vals)]['N']
worst_N = results[np.argmax(cos_vals)]['N']
print(f"Best orthogonality at N={best_N} (cos={min(cos_vals):.6f})")
print(f"Worst at N={worst_N} (cos={max(cos_vals):.6f})")
print()

print("Interpretation:")
if max(cos_vals) > 0.5:
    print(f"  ✗ At some N, vectors align significantly (cos up to {max(cos_vals):.3f})")
    print(f"    The misalignment is NOT robust to N choice at this height")
else:
    print(f"  ✓ Vectors stay mostly orthogonal across all N")
print()

print("="*90)
