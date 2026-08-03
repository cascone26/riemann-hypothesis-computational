"""Quantitative check of the 625 proof skeleton:
J(u) =? sum_n mu(n) n^{-1/2} M_n(u)/K with M_n(u) = int |zeta| e^{i(2 theta + (u - log n) t)} dt/(2pi...)
Compute M_n resonance profile numerically: is M peaked at u=log n, and does the skeleton
resum to the measured J?"""
import numpy as np
from scipy.special import loggamma as slg
T0, T, dt = 50.0, 3000.0, 0.004
t = np.arange(T0, T, dt)
th = np.imag(slg(0.25+0.5j*t)) - t/2*np.log(np.pi)
M = np.sqrt(t/(2*np.pi)).astype(int)
Z = np.zeros_like(t)
for n in range(1, M.max()+1):
    m = M >= n
    Z[m] += 2*np.cos(th[m] - t[m]*np.log(n))/np.sqrt(n)
absZ = np.abs(Z)   # |zeta| = |Z|
K = (th[-1]-th[0])/np.pi
# resonance kernel: M(v) = (1/K) int |Z| e^{i(2 theta(t) + v t)} dt  -- peak expected at stationary pts
def Mker(v):
    return np.trapezoid(absZ*np.exp(1j*(2*th + v*t)), dx=dt)/K
print("resonance kernel |M(v)| around v = -log n + u ... scan v:")
for v in np.arange(-3.2, 0.3, 0.25):
    print(f"  v={v:+.2f}  |M|={abs(Mker(v)):.4f}")
# skeleton resum at u = log 2 and log 6: J_skel = sum_n mu(n) n^-1/2 M(u - log n)
from sympy import mobius
w_t = np.sign(Z)*np.exp(1j*th)
for u in (np.log(2), np.log(6), np.log(5)):
    Jdir = np.trapezoid(w_t*np.exp(1j*u*t), dx=dt)/K
    Jskel = sum(int(mobius(n))*n**-0.5*Mker(u - np.log(n)) for n in range(1, 40) if mobius(n) != 0)
    print(f"u={u:.3f}: J_direct={Jdir.real:+.4f}{Jdir.imag:+.4f}i  J_skeleton={Jskel.real:+.4f}{Jskel.imag:+.4f}i")
print("DONE")
