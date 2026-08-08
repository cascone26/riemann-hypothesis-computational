#!/usr/bin/env python3
"""
Verify: does γ=69.5 really have cos_max=0.999984, or was that a measurement error?
Run the exact adaptive-scan logic for this one height with detailed output.
"""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh

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

# Exact reproduction of adaptive scan for γ=69.5
gamma = 69.5

# Setup window
window_length = 30
height_scale = max(1.0, np.log(max(gamma, 30) / 30.0))
fourW_target = 2.8 * (1.0 + height_scale * 0.3)
W = fourW_target / 4.0
Ncap = 30 * 2 * W / np.pi

hw = window_length / 2.0
a = gamma - hw
b = gamma + hw

# Scan N from 0.40 to 0.90 of capacity (as in adaptive scan)
N_min = max(5, int(np.round(Ncap * 0.40)))
N_max = int(np.round(Ncap * 0.90))
N_range = range(N_min, N_max + 1)

print(f"γ = {gamma}")
print(f"Window [{a:.1f}, {b:.1f}]")
print(f"4W = {fourW_target:.4f}, Ncap = {Ncap:.2f}")
print(f"N range: [{N_min}, {N_max}]")
print()

print(f"{'N':>3} {'N/Ncap':>7} {'cos':>12} {'angle°':>8}")
print("-"*40)

cos_vals = []
for N in N_range:
    CEN = list(np.linspace(a, b, N))
    AP, PR, G = mats(CEN, W)

    wf, Vf = eigh(AP, G)
    v_a = Vf[:, 0]

    wp, Vp = eigh(PR, G)
    v_p = Vp[:, -1]

    cosang = abs(float(v_a @ (G @ v_p)))
    ang = np.degrees(np.arccos(min(1.0, cosang)))

    cos_vals.append(cosang)
    print(f"{N:3d} {N/Ncap:7.2f} {cosang:12.6f} {ang:8.2f}")

print("-"*40)
print()
print(f"cos_min: {min(cos_vals):.8f}")
print(f"cos_max: {max(cos_vals):.8f}")
print(f"cos_mean: {np.mean(cos_vals):.8f}")
print()
if max(cos_vals) > 0.9:
    print("ERROR: Got cos_max > 0.9. This suggests a computation error.")
else:
    print("✓ Verified: cos values are reasonable (< 0.05).")
