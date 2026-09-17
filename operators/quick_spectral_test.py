"""Quick test of direct spectral method with minimal truncation."""
import numpy as np
from scipy.integrate import solve_ivp
import sys, os

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from canonical_weyl_solver import (
    explicit_formula_symbol,
    solve_canonical_system_finite_a,
    extract_weyl_m_function
)

print("="*80)
print("QUICK SPECTRAL κ TEST")
print("="*80)

# Test just ONE z-value and ONE U-value with minimal a-range
z_test = complex(1.0, 1.0)
U_test = 1
a_values = [1, 2, 3, 4, 5]

print(f"\nTest setup:")
print(f"  z = {z_test}")
print(f"  U = {U_test}")
print(f"  a values: {a_values}")
print()

m_values = []
print("  a    |m(z)|       Re(m)          Im(m)")
print("-" * 60)

for a in a_values:
    try:
        X_final = solve_canonical_system_finite_a(z_test, a)
        if X_final is not None:
            m_a = extract_weyl_m_function(X_final, z_test, a)
            if not np.isnan(m_a):
                m_values.append(m_a)
                print(f"  {a}    {abs(m_a):.6e}  {m_a.real:+.6e}  {m_a.imag:+.6e}")
            else:
                print(f"  {a}    NaN")
                m_values.append(m_a)
        else:
            print(f"  {a}    Solve failed")
            m_values.append(np.nan + 1j*np.nan)
    except Exception as e:
        print(f"  {a}    Error: {str(e)[:40]}")
        m_values.append(np.nan + 1j*np.nan)

m_values = np.array(m_values)

print()
print("Analysis:")
print(f"  Valid points: {np.sum(~np.isnan(m_values.real))}")
print(f"  |m| range: {np.nanmin(np.abs(m_values)):.3e} to {np.nanmax(np.abs(m_values)):.3e}")
print(f"  Final |m|: {abs(m_values[-1]):.3e}")

# Look for growth pattern
valid_m = m_values[~np.isnan(m_values.real)]
if len(valid_m) > 1:
    magnitudes = np.abs(valid_m)
    diffs = np.diff(magnitudes)
    print(f"\n  Magnitude progression:")
    for i, (a, mag) in enumerate(zip(a_values[:len(magnitudes)], magnitudes)):
        if i > 0:
            print(f"    a={a}: |m|={mag:.3e} (Δ={diffs[i-1]:+.3e})")
        else:
            print(f"    a={a}: |m|={mag:.3e}")

    max_jump = np.max(np.abs(diffs))
    if max_jump > 1.0:
        print(f"\n  POSSIBLE POLE SIGNATURE: large jump of {max_jump:.3e} detected")
    else:
        print(f"\n  NO STRONG POLE SIGNATURE: max jump {max_jump:.3e}")

print("\n" + "="*80)
print("DONE")
print("="*80)
