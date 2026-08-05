"""Matched-edge deficit scan: hold X/(T/2pi) FIXED at 0.3 (below edge) across T windows;
does the finite-T deficit decay with T at fixed edge-ratio? (settles 688b)"""
import numpy as np
from scipy.special import loggamma as slg
from sympy import primerange
z = np.load("zeros100k.npy")
def build_S(t, X, taper_at=None):
    S = np.zeros_like(t)
    lim = taper_at*2 if taper_at else X
    for p in primerange(2, lim):
        lp = np.log(p); m = 1
        while p**m < lim:
            n = p**m
            if taper_at:
                zz = (np.log(n)-(np.log(taper_at)-np.log(2)))/(2*np.log(2))
                w = 1.0 if zz<=0 else (0.0 if zz>=1 else 0.5*(1+np.cos(np.pi*zz)))
            else: w = 1.0
            if w>0: S += -(w/np.pi)*np.sin(t*lp*m)/(n**0.5*m)
            m += 1
    return S
def buildd(t, X):
    S=np.zeros_like(t); dS=np.zeros_like(t)
    for p in primerange(2,X):
        lp=np.log(p); m=1
        while p**m<X:
            n=p**m; amp=1/(np.pi*np.sqrt(n)*m)
            S+=-amp*np.sin(t*lp*m); dS+=-amp*lp*m*np.cos(t*lp*m); m+=1
    return S,dS
u=np.log(2)
for T0,T1 in ((3000,5000),(9000,11000),(30000,34000),(68000,72000)):
    dt=0.006
    t=np.arange(T0,T1,dt)
    th=np.imag(slg(0.25+0.5j*t))-t/2*np.log(np.pi)
    K=(th[-1]-th[0])/np.pi
    Tmid=0.5*(T0+T1); edge=Tmid/(2*np.pi)
    X=int(0.3*edge)   # matched: X = 0.3 * edge
    Xtap=max(20000,3*X)
    S_true=np.searchsorted(z,t)-th/np.pi-1.0
    SY=build_S(t,None,min(Xtap,30000))
    SX,dSX=buildd(t,X)
    Psi=-np.pi*SX+u*t; Psip=-np.pi*dSX+u; eiP=np.exp(1j*Psi)
    d=(((S_true[-1]*eiP[-1]-S_true[0]*eiP[0])-1j*np.trapezoid(S_true*SY*0+S_true*Psip*eiP,dx=dt))/K).real
    print(f"T~{Tmid:.0f} edge={edge:.0f} X={X}(0.3edge): deficit {d:+.5f}", flush=True)
print("DONE")
