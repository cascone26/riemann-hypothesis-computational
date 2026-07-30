#!/usr/bin/env python3
"""Weyl, §427: does F_op(T)=N_op(T)-<N>(T) correlate with S(T)=N_true(T)-<N>(T)?
Uses the naive bare B(E)=0 extraction at large X (X=80,100) -- per the room's own referee finding
(von Neumann/Arb), bare and phase-corrected converge to the SAME asymptotic spectrum, so bare-at-large-X
is a clean, bug-free stand-in for 'the X-stable corrected basis' without inheriting my earlier WKB
root-dropping bug.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import pearsonr

Z = np.load('zeros500.npy')
FP=8*np.pi**2; w=2*np.pi*np.sqrt(2); sp=np.sqrt(2.0)

def find_op_roots(X, gmax=95.0, dg=0.02):
    def rhs(x,y,alpha):
        xi,xip=y
        return [xip, -(2*x*xip+(FP*x*x-alpha)*xi)/(x*x-2)]
    def Bval(gamma):
        alpha=-(gamma/2.0)**2; dxi0=-(2*FP-alpha)/(2*sp)
        xs=np.linspace(X-8,X,300)
        s=solve_ivp(rhs,[sp+1e-6,X],[1.0,dxi0],args=(alpha,),method='DOP853',rtol=1e-11,atol=1e-13,t_eval=xs)
        u=xs*s.y[0]; M=np.column_stack([np.sin(w*xs),np.cos(w*xs)])
        (A,B),*_=np.linalg.lstsq(M,u,rcond=None); return B
    gs=np.arange(3.0,gmax,dg)
    Bs=np.array([Bval(g) for g in gs])
    roots=[]
    for i in range(len(gs)-1):
        if Bs[i]*Bs[i+1]<0:
            a,b=gs[i],gs[i+1]
            for _ in range(26):
                m=(a+b)/2
                if Bval(a)*Bval(m)<0: b=m
                else: a=m
            roots.append((a+b)/2)
    return np.array(roots)

def Nsmooth(T):
    return (T/(2*np.pi))*np.log(T/(2*np.pi*np.e)) + 7.0/8.0

def N_count(roots, T):
    return np.searchsorted(roots, T, side='right')

gmax = 95.0
grid = np.linspace(14.0, gmax, 300)
Ntrue = np.array([np.searchsorted(Z, t, side='right') for t in grid])
Ssm = Nsmooth(grid)
S_true = Ntrue - Ssm

print("Finding op roots at X=80 and X=100 (this is the heavy part)...")
for X in [80.0, 100.0]:
    roots = find_op_roots(X, gmax=gmax)
    Nop = np.array([N_count(roots, t) for t in grid])
    Fop = Nop - Ssm
    rho, pval = pearsonr(Fop, S_true)
    rms = np.sqrt(np.mean((Fop - S_true)**2))
    mean_spacing = 2*np.pi/np.log(gmax/(2*np.pi))
    print(f"\nX={X}: n_op_roots(<{gmax})={len(roots)}  n_true_zeros(<{gmax})={int(np.sum(Z<gmax))}")
    print(f"  Pearson rho(F_op, S_true) = {rho:.4f}  (p={pval:.2e})")
    print(f"  RMS(F_op - S_true) = {rms:.4f}   mean spacing at T={gmax}: {mean_spacing:.4f}")
    print(f"  RMS / spacing = {rms/mean_spacing:.4f}")
    print(f"  max|F_op|={np.max(np.abs(Fop)):.3f}  max|S_true|={np.max(np.abs(S_true)):.3f}")
