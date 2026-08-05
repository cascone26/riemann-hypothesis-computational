"""§650 finite-T error decay at 7.5x: identity split (smooth + tail + finite-T) on a high window."""
import numpy as np
from scipy.special import loggamma as slg
from sympy import primerange
z = np.load("zeros100k.npy")
T0, T1, dt = 60000.0, 74000.0, 0.004
t = np.arange(T0, T1, dt)
th = np.imag(slg(0.25+0.5j*t)) - t/2*np.log(np.pi)
thp = np.gradient(th, dt)
K = (th[-1]-th[0])/np.pi
S_true = np.searchsorted(z, t) - th/np.pi - 1.0
def build(X, taper_at=None):
    S = np.zeros_like(t)
    lim = taper_at*2 if taper_at else X
    for p in primerange(2, lim):
        lp = np.log(p); m = 1
        while p**m < lim:
            n = p**m
            if taper_at:
                zz = (np.log(n) - (np.log(taper_at)-np.log(2)))/(2*np.log(2))
                w = 1.0 if zz <= 0 else (0.0 if zz >= 1 else 0.5*(1+np.cos(np.pi*zz)))
            else:
                w = 1.0
            if w > 0:
                S += -(w/np.pi)*np.sin(t*lp*m)/(n**0.5*m)
            m += 1
    return S
def buildd(X):
    S = np.zeros_like(t); dS = np.zeros_like(t)
    for p in primerange(2, X):
        lp = np.log(p); m = 1
        while p**m < X:
            n = p**m
            amp = 1/(np.pi*np.sqrt(n)*m)
            S += -amp*np.sin(t*lp*m); dS += -amp*lp*m*np.cos(t*lp*m)
            m += 1
    return S, dS
print("building fields...", flush=True)
SY = build(None, 30000)
SY2 = build(None, 100000)
u = np.log(2)
X = 3000
SX, dSX = buildd(X)
Psi = -np.pi*SX + u*t
Psip = -np.pi*dSX + u
eiP = np.exp(1j*Psi)
def piece(R):
    return (((R[-1]*eiP[-1]-R[0]*eiP[0]) - 1j*np.trapezoid(R*Psip*eiP, dx=dt))/K).real
d_true = piece(S_true - SY)
d_tail = piece(SY2 - SY)
print(f"HIGH window T=[60k,74k], X=3000: deficit {d_true:+.4f} = tail {d_tail:+.4f} + finite-T {d_true-d_tail:+.4f}")
print("(at T~1e4 was: deficit -0.0269 = tail -0.0088 + finite-T -0.0182; half-window low/high: -0.0207/-0.0162)")
print("DONE")
