#!/usr/bin/env python3
"""
Diagnose the discrepancy: compare Omega values from scipy vs mpmath.
"""

import numpy as np
from mpmath import mp, mpc, digamma as mp_digamma, log
from scipy.special import digamma as scipy_digamma

mp.dps = 50

LAMBDA_SQ = 14
L = 2.0 * mp.log(mp.sqrt(mp.mpf(LAMBDA_SQ)))
LOGPI_mp = mp.log(mp.pi)
LOGPI_np = float(np.log(np.pi))

print(f"L = {L}")
print(f"LOGPI_mp = {LOGPI_mp}")
print(f"LOGPI_np = {LOGPI_np}")

# Test digamma values at several r points
test_rs = [1.0, 10.0, 50.0, 100.0]

print("\n" + "="*80)
print("Comparing digamma values:")
print("="*80)

for r in test_rs:
    # scipy version (from ccm_gap_scaling.py)
    z_scipy = 0.25 + 0.5j * r
    psi_scipy = scipy_digamma(z_scipy)
    omega_scipy = np.real(psi_scipy) - LOGPI_np

    # mpmath version (from our lite script)
    z_mp = mp.mpf(0.25) + mp.mpc(0, mp.mpf(r) / 2.0)
    psi_mp = mp_digamma(z_mp)
    omega_mp = psi_mp.real - LOGPI_mp

    print(f"\nr = {r}")
    print(f"  scipy Omega:   {float(omega_scipy):.6e}")
    print(f"  mpmath Omega:  {float(omega_mp):.6e}")
    print(f"  Difference:    {float(omega_mp) - float(omega_scipy):.6e}")
    print(f"  Ratio (mp/sp): {float(omega_mp) / float(omega_scipy):.6f}")

# Also test if Omega is truly zero-free
print("\n" + "="*80)
print("Checking if Omega is zero-free (no sign changes in [1, 60]):")
print("="*80)

rs = np.linspace(1, 60, 100)
omegas_scipy = []
omegas_mp = []

for r in rs:
    z_scipy = 0.25 + 0.5j * r
    psi_scipy = scipy_digamma(z_scipy)
    omega_scipy = np.real(psi_scipy) - LOGPI_np
    omegas_scipy.append(omega_scipy)

    z_mp = mp.mpf(0.25) + mp.mpc(0, mp.mpf(r) / 2.0)
    psi_mp = mp_digamma(z_mp)
    omega_mp = psi_mp.real - LOGPI_mp
    omegas_mp.append(float(omega_mp))

omegas_scipy = np.array(omegas_scipy)
omegas_mp = np.array(omegas_mp)

print(f"Scipy Omega min: {omegas_scipy.min():.6e}, max: {omegas_scipy.max():.6e}")
print(f"mpmath Omega min: {omegas_mp.min():.6e}, max: {omegas_mp.max():.6e}")
print(f"Scipy has sign change? {np.any(omegas_scipy[:-1] * omegas_scipy[1:] < 0)}")
print(f"mpmath has sign change? {np.any(omegas_mp[:-1] * omegas_mp[1:] < 0)}")

# Check if the issue is with the grid resolution
print("\n" + "="*80)
print("Analysis: Grid resolution impact (2001 vs 60000)")
print("="*80)

print("With G=2001 points, we integrate over the full domain with sparse sampling.")
print("The arch matrix entries depend on the integral:")
print("  arch_{n,m} = (1/2π) ∫ Omega(r) e_hat_n(r) e_hat_m*(r) dr")
print("\nWith coarse grid, fine oscillations in e_hat or Omega may be missed.")
print("For N=20, the Fourier modes go up to s_n = 2π*20/L ≈", float(2*np.pi*20/L))
print("With grid spacing h = 2π*20*2/L / 2000 ≈", float(2*np.pi*20*2/L / 2000))
