"""Box (iii) assembly: deficit(X) = -i int R_Y Psi' e^{iPsi} dt / K.
R_Y = S_true - S_Y(30k,taper). Test: R_tail = S_{Y2=100k,taper} - S_Y (pure prime tail,
zero-free) vs the remainder (finite-T Landau corrections + continuous spectrum)."""
import numpy as np
from scipy.special import loggamma as slg
from sympy import primerange
g = np.load("zeros10k.npy")
T0, T, dt = 15.0, 9877.0, 0.01
t = np.arange(T0, T, dt)
th = np.imag(slg(0.25+0.5j*t)) - t/2*np.log(np.pi)
thp = np.gradient(th, dt)
K = (th[-1]-th[0])/np.pi
S_true = np.searchsorted(g, t) - th/np.pi - 1.0
def build(X, taper_at=None):
    S = np.zeros_like(t); dS = np.zeros_like(t)
    lim = taper_at*2 if taper_at else X
    for p in primerange(2, lim):
        lp = np.log(p); m = 1
        while p**m < lim:
            n = p**m
            if taper_at:
                z = (np.log(n) - (np.log(taper_at)-np.log(2)))/(2*np.log(2))
                w = 1.0 if z <= 0 else (0.0 if z >= 1 else 0.5*(1+np.cos(np.pi*z)))
            else:
                w = 1.0
            if w > 0:
                amp = w/(np.pi*np.sqrt(n)*m)
                S += -amp*np.sin(t*lp*m)
                if not taper_at: dS += -amp*lp*m*np.cos(t*lp*m)
            m += 1
    return S, dS
print("building fields...", flush=True)
SY, _ = build(None, taper_at=30000)
SY2, _ = build(None, taper_at=100000)
R_true_base = S_true - SY
R_tail = SY2 - SY
u = np.log(2)
print(" X     deficit(true R)   tail part   finite-T remainder", flush=True)
for X in (300, 1000, 3000, 10000):
    SX, dSX = build(X)
    Psi = -np.pi*SX + u*t
    Psip = -np.pi*dSX + u
    eiP = np.exp(1j*Psi)
    def piece(R):
        return (((R[-1]*eiP[-1]-R[0]*eiP[0]) - 1j*np.trapezoid(R*Psip*eiP, dx=dt))/K).real
    d_true = piece(R_true_base)
    d_tail = piece(R_tail)
    print(f"{X:>6}   {d_true:+.4f}         {d_tail:+.4f}     {d_true-d_tail:+.4f}", flush=True)
print("DONE")
