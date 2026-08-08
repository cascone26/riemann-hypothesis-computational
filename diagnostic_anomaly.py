#!/usr/bin/env python3
"""
Diagnose the anomaly at γ=124.3 where cos(v_p,v_a) ≈ 0.9957 (nearly aligned).
This violates §835's claim of height-stable misalignment.
"""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
import mpmath

# Setup (copied from main)
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

print("="*100)
print("DIAGNOSTIC: γ=124.3 anomaly — why cos(v_p,v_a) ≈ 0.9957?")
print("="*100)
print()

# Reproduce the anomalous case
gamma = 124.3
capacity_fraction = 0.60
window_length = 30
N = 11  # As observed in the run

height_scale = max(1.0, np.log(max(gamma, 30) / 30.0))
fourW_target = 2.8 * (1.0 + height_scale * 0.3)
W = fourW_target / 4.0
Ncap = 30 * 2 * W / np.pi

hw = window_length / 2.0
a = gamma - hw
b = gamma + hw

print(f"Measurement parameters:")
print(f"  Height γ:          {gamma}")
print(f"  Window [a,b]:      [{a:.1f}, {b:.1f}]")
print(f"  Window width:      {window_length}")
print(f"  4W (target):       {fourW_target:.4f}")
print(f"  W:                 {W:.4f}")
print(f"  N_cap:             {Ncap:.2f}")
print(f"  N (fixed):         {N}")
print(f"  N/N_cap:           {N/Ncap:.3f}")
print()

CEN = list(np.linspace(a, b, N))
print(f"Window centers:")
for i, c in enumerate(CEN):
    print(f"  {i}: {c:.2f}")
print()

AP, PR, G = mats(CEN, W)

print(f"Matrix sizes and conditioning:")
gw = np.linalg.eigvalsh(G)
condG = float(gw[-1] / max(gw[0], 1e-300))
print(f"  Gram eigenvalues: {gw}")
print(f"  Gram cond(G):     {condG:.1f}")
print()

# Arch+pole
wf, Vf = eigh(AP, G)
floor = float(wf[0])
v_a = Vf[:, 0]
print(f"Archimedean+pole spectrum (G-generalized eigvalues):")
print(f"  Full spectrum: {wf}")
print(f"  Floor (min):   {floor:.6f}")
print(f"  v_a (bottom eigvec): {v_a}")
print()

# Prime
wp, Vp = eigh(PR, G)
pnorm = float(wp[-1])
v_p = Vp[:, -1]
print(f"Prime spectrum:")
print(f"  Full spectrum: {wp}")
print(f"  pnorm (max):   {pnorm:.6f}")
print(f"  v_p (top eigvec): {v_p}")
print()

# Margin
wm, _ = eigh(AP - PR, G)
margin = float(wm[0])
print(f"Margin Q (arch+pole-prime):")
print(f"  Full spectrum: {wm}")
print(f"  Margin (min):  {margin:.6e}")
print()

# Misalignment (THE ANOMALY)
cosang = abs(float(v_a @ (G @ v_p)))
ang = np.degrees(np.arccos(min(1.0, cosang)))
print(f"MISALIGNMENT (THE ANOMALY):")
print(f"  cos(v_p, v_a) = v_a^T · G · v_p = {cosang:.8f}")
print(f"  Angle:                            {ang:.4f}°")
print()

# Debug: check if v_a and v_p are actually close
dot_without_G = float(v_a @ v_p)
print(f"Debug: inner products:")
print(f"  v_a · v_p (no G):                {dot_without_G:.8f}")
print(f"  v_a · (G @ v_p):                 {cosang:.8f}")
print(f"  ‖v_a‖_G (G-norm of v_a):        {float(np.sqrt(v_a @ (G @ v_a))):.8f}")
print(f"  ‖v_p‖_G (G-norm of v_p):        {float(np.sqrt(v_p @ (G @ v_p))):.8f}")
print()

# Is this a numerical artifact? Check if v_a ≈ v_p
v_diff = v_a - v_p
print(f"Eigenvector difference analysis:")
print(f"  ‖v_a - v_p‖:                    {np.linalg.norm(v_diff):.8f}")
print(f"  max|v_a - v_p|:                 {np.max(np.abs(v_diff)):.8f}")
print()

# Check cross-eigenvector products
print(f"Full (arch+pole) × prime eigenvector products:")
print(f"  v_a (arch bottom) vs all prime eigvecs:")
for i, wp_val in enumerate(wp):
    v_prime_i = Vp[:, i]
    prod = abs(float(v_a @ (G @ v_prime_i)))
    print(f"    prime eig#{i} (val={wp_val:.6f}): cos={prod:.6f}")
print()

print("Hypothesis checks:")
print(f"  Is v_a ≈ v_p? (eigenvector nearly identical)")
if np.linalg.norm(v_diff) < 0.1:
    print(f"    YES — they're very close (diff={np.linalg.norm(v_diff):.6f})")
else:
    print(f"    NO — they're distinct (diff={np.linalg.norm(v_diff):.6f})")
print()

if margin > 0:
    print(f"  ✓ Margin still positive (Q ≥ 0), so no RH violation at float precision")
    print(f"    But orthogonality claim (§835) is WEAKER at this height")
else:
    print(f"  ✗ Margin negative (Q < 0) — RH VIOLATED if real")
print()

print("="*100)
