#!/usr/bin/env python3
"""
toy_lambda_flow.py — LEDGER §131 (Lens B / Λ≤0 functional prep): validate the dBN zero-flow dynamics on a TOY system
WITH a complex conjugate pair (the Λ>0 scenario), to (i) confirm the collision-onto-the-line mechanism at a toy "Λ",
(ii) test which functionals are monotone, (iii) demonstrate the CIRCULARITY of pure-flow functionals for the upper bound.
Flow (RT): dx_j/dt = 2 Σ_{k≠j} 1/(x_j − x_k), extended to complex x_j. Forward t ⟹ zeros more real (Im→0).
"""
import numpy as np
from scipy.integrate import solve_ivp

def flow_rhs(t, y):
    x = y[:len(y)//2] + 1j*y[len(y)//2:]
    dx = np.zeros_like(x)
    for j in range(len(x)):
        s = 0j
        for k in range(len(x)):
            if k!=j: s += 1.0/(x[j]-x[k])
        dx[j] = 2*s
    return np.concatenate([dx.real, dx.imag])

# toy: 8 real zeros ~arithmetic progression + one complex pair a±ib (RH-violating config)
real0 = np.array([-7.,-5.,-3.,-1.,1.,3.,5.,7.])
cpair = np.array([0.3+0.8j, 0.3-0.8j])          # a complex conjugate pair (off the real axis)
x0 = np.concatenate([real0, cpair])
y0 = np.concatenate([x0.real, x0.imag])

def functionals(x):
    E = sum(1.0/abs(x[j]-x[k])**2 for j in range(len(x)) for k in range(len(x)) if k!=j)   # RT energy
    Im2 = np.sum(x.imag**2)                                                                  # complex-zero measure
    # a candidate ARITHMETIC-flavored functional: weighted by 1/|x|² (mimics the 2/γ² Weil weight)
    Wq = sum(1.0/(abs(x[j])**2+1) * abs(x[j]-x[k])**-2 for j in range(len(x)) for k in range(len(x)) if k!=j)
    return E, Im2, Wq

print("Toy dBN flow: complex pair 0.3±0.8i + 8 real zeros. Flow FORWARD (t↑); does Im→0 (collision at toy-Λ)?")
print(f"{'t':>7} {'max|Im|':>10} {'ΣIm²':>12} {'E(RT)':>12} {'Wq(arith?)':>12}")
sol = solve_ivp(flow_rhs, [0, 4.0], y0, t_eval=np.linspace(0,4,17), rtol=1e-9, atol=1e-12)
prev=None
for i,t in enumerate(sol.t):
    x = sol.y[:len(x0),i] + 1j*sol.y[len(x0):,i]
    E,Im2,Wq = functionals(x)
    mono_E = '' if prev is None else ('↓' if E<prev[0] else '↑')
    mono_I = '' if prev is None else ('↓' if Im2<prev[1] else '↑')
    print(f"{t:>7.2f} {np.max(np.abs(x.imag)):>10.4e} {Im2:>12.4e}{mono_I} {E:>12.4e}{mono_E} {Wq:>12.4e}",flush=True)
    prev=(E,Im2)
print("\nREAD: Im→0 = the complex pair collides onto the real axis at toy-Λ (forward flow makes real). ΣIm² & E monotone?")
print("CIRCULARITY: ΣIm²(t=0) is exactly the RH-violation measure — monotone but its t=0 value = the answer. Pure-flow")
print("functionals only give the LOWER bound (RT). The Λ≤0 functional must couple to the ARITHMETIC (the primes/ξ),")
print("not just the zero-positions — testing whether any position-only functional can distinguish RH-true from -false.")
