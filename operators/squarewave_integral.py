"""Continuous form of the elementary readout: does the INTEGRAL
  I(u) = (1/K) int_{T0}^{T} sgn Z(t) e^{i theta(t)} e^{iut} (theta'(t)/pi) dt
reproduce the mu-comb that the zero-sum gives? (K = normalizing count = (theta(T)-theta(T0))/pi)

If yes: the mu-content comes from the absolutely continuous part of the counting measure,
weight = bounded unit-modulus square wave of Hardy Z -> classical Fourier analysis territory.
Z via vectorized RS main sum (error O(t^-1/4), only matters near |Z| tiny - measure ~0).
"""
import numpy as np
from scipy.special import loggamma as slg

T0, T, dt = 50.0, 5000.0, 0.004
t = np.arange(T0, T, dt)
th = np.imag(slg(0.25+0.5j*t)) - t/2*np.log(np.pi)
thp = np.gradient(th, dt)

# vectorized RS main sum for Z
M = np.sqrt(t/(2*np.pi)).astype(int)
Mmax = M.max()
Z = np.zeros_like(t)
for n in range(1, Mmax+1):
    mask = M >= n
    Z[mask] += 2*np.cos(th[mask] - t[mask]*np.log(n))/np.sqrt(n)
sq = np.sign(Z)

base = sq*np.exp(1j*th)*thp/np.pi
K = (th[-1]-th[0])/np.pi
print(f"T=[{T0},{T}]  K~{K:.0f} zeros' worth  grid {len(t)} pts")
res = []
for n, mu in [(2,-1),(3,-1),(4,0),(5,-1),(6,1),(7,-1),(8,0),(10,1),(11,-1),(13,-1),(14,1),(15,1),(30,-1)]:
    I = 1j*np.trapezoid(base*np.exp(1j*np.log(n)*t), dx=dt)/K
    ok = "OK " if (mu==0 and abs(I.real)<0.05) or (mu!=0 and np.sign(I.real)==mu) else "MISS"
    res.append(ok)
    print(f"  n={n:<3} mu={mu:+d}  Re I = {I.real:+.4f}  Im = {I.imag:+.4f}  {ok}")
# off-support control
for u in [np.log(2)+0.3, 2.2]:
    I = 1j*np.trapezoid(base*np.exp(1j*u*t), dx=dt)/K
    print(f"  off-support u={u:.3f}  |I| = {abs(I):.4f}")
print("verdict:", "CONTINUOUS FORM CARRIES THE COMB" if all(r=="OK " for r in res) else "partial: %s"%res)
print("DONE")
