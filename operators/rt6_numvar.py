#!/usr/bin/env python3
"""§454 (spectral rigidity — Hilbert-Pólya evidence): number variance Σ²(L) of the zeros vs GUE & Poisson.
GUE (a self-adjoint operator's spectrum): Σ²(L) ~ (1/π²)(ln(2πL)+γ_E+1) — LOG growth (rigid). Poisson (random):
Σ²(L)=L (linear). Unfold zeros via exact θ (unit mean spacing), compute Σ²(L). Rigidity (log, not linear) =
evidence the zeros ARE a self-adjoint spectrum (Θ exists)."""
import numpy as np, mpmath as mp
mp.mp.dps=20
Z=np.load('zeros150.npy')
x=np.array([float(mp.siegeltheta(g))/np.pi+1 for g in Z])   # unfolded positions (unit mean spacing)
x=x-x[0]
def numvar(L,nsamp=4000):
    starts=np.linspace(0,x[-1]-L,nsamp)
    counts=np.array([np.sum((x>=s)&(x<s+L)) for s in starts])
    return counts.var()
EG=0.5772156649
print("§454 number variance Σ²(L): zeros vs GUE(log) vs Poisson(linear). unfolded 150 zeros.")
print(f"{'L':>5} {'Σ²(zeros)':>11} {'GUE':>8} {'Poisson':>8}")
for L in [1,2,3,5,8,12,16,20]:
    sv=numvar(float(L)); gue=(1/np.pi**2)*(np.log(2*np.pi*L)+EG+1)
    print(f"{L:>5} {sv:11.4f} {gue:8.4f} {L:8.1f}")
print("""
reading: Σ²(zeros) tracking the GUE LOG curve (not the Poisson line) => the zeros have GUE spectral RIGIDITY =
strong evidence they are the spectrum of a self-adjoint operator (Hilbert-Pólya Θ exists). This is the statistical
face of the same conclusion: Θ = xp (§449) with GUE fluctuations; RH = Θ self-adjoint. (150 zeros is modest; the
rigidity should still be visibly sub-linear vs Poisson.)
""")
