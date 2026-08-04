"""Linked-cluster resummation test (item 4):
A_zeros(u,X) = (1/K)[ int e^{iPsi}(theta'/pi)dt + [S e^{iPsi}] - i int S(t) Psi'(t) e^{iPsi} dt ]
with Psi = Phi_X + ut, Phi_X = -pi S_X (prime oscillator field), S ~ S_Y (Y >> X, zero-free).
Predict the §632 dip-recovery: measured -0.3035/-0.2959/-0.3199/-0.3255 at X=300/1e3/3e3/1e4."""
import numpy as np
from scipy.special import loggamma as slg
from sympy import primerange
T0, T = 15.0, 9877.0
dt = 0.005
t = np.arange(T0, T, dt)
th = np.imag(slg(0.25+0.5j*t)) - t/2*np.log(np.pi)
thp = np.gradient(th, dt)
K = (th[-1]-th[0])/np.pi
def build_S(X, taper=True):
    S = np.zeros_like(t)
    dS = np.zeros_like(t)   # derivative
    for p in primerange(2, X):
        lp = np.log(p)
        m = 1
        while p**m < X:
            n = p**m; u_n = m*lp
            if taper:
                z = (np.log(n) - (np.log(X)-np.log(2)))/(2*np.log(2))
                w = 1.0 if z <= 0 else (0.0 if z >= 1 else 0.5*(1+np.cos(np.pi*z)))
            else:
                w = 1.0
            if w > 0:
                amp = w/(np.pi*np.sqrt(n)*m)
                S += -amp*np.sin(t*u_n)
                dS += -amp*u_n*np.cos(t*u_n)
            m += 1
    return S, dS
print("building S_Y (Y=30000)...", flush=True)
SY, _ = build_S(30000)
import sys
u = np.log(2)
for X in (300, 1000, 3000, 10000):
    SX, dSX = build_S(X, taper=False)
    Phi = -np.pi*SX
    Psi = Phi + u*t
    Psip = -np.pi*dSX + u
    eiP = np.exp(1j*Psi)
    A_smooth = np.trapezoid(eiP*thp/np.pi, dx=dt)/K
    bnd = (SY[-1]*eiP[-1] - SY[0]*eiP[0])/K
    A_dS = (bnd - 1j*np.trapezoid(SY*Psip*eiP, dx=dt))/K
    A = A_smooth + A_dS
    print(f"X={X:>6}: smooth {A_smooth.real:+.4f}  dS-part {A_dS.real:+.4f}  TOTAL {A.real:+.4f}", flush=True)
print("measured (632): -0.3035 / -0.2959 / -0.3199 / -0.3255")
print("DONE")
