#!/usr/bin/env python3
"""Weyl, §427/428 fast pass: F_op(T) vs S(T) correlation, lighter grid, single X, lower rtol
(Report-tier scipy, not trying for precision -- trying to get a number out before the room needs it)."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import pearsonr

Z = np.load('zeros500.npy')
FP=8*np.pi**2; w=2*np.pi*np.sqrt(2); sp=np.sqrt(2.0)

def find_op_roots(X, gmax=95.0, dg=0.1):
    def rhs(x,y,alpha):
        xi,xip=y
        return [xip, -(2*x*xip+(FP*x*x-alpha)*xi)/(x*x-2)]
    def Bval(gamma):
        alpha=-(gamma/2.0)**2; dxi0=-(2*FP-alpha)/(2*sp)
        xs=np.linspace(X-8,X,150)
        s=solve_ivp(rhs,[sp+1e-6,X],[1.0,dxi0],args=(alpha,),method='RK45',rtol=1e-9,atol=1e-11,t_eval=xs)
        u=xs*s.y[0]; M=np.column_stack([np.sin(w*xs),np.cos(w*xs)])
        (A,B),*_=np.linalg.lstsq(M,u,rcond=None); return B
    gs=np.arange(3.0,gmax,dg)
    Bs=np.array([Bval(g) for g in gs])
    roots=[]
    for i in range(len(gs)-1):
        if Bs[i]*Bs[i+1]<0:
            a,b=gs[i],gs[i+1]
            for _ in range(20):
                m=(a+b)/2
                if Bval(a)*Bval(m)<0: b=m
                else: a=m
            roots.append((a+b)/2)
    return np.array(roots)

def Nsmooth(T):
    return (T/(2*np.pi))*np.log(T/(2*np.pi*np.e)) + 7.0/8.0

gmax = 95.0
X = 80.0
print("Finding op roots at X=80, dg=0.1 (fast pass)...")
roots = find_op_roots(X, gmax=gmax)
print(f"found {len(roots)} roots < {gmax}: {np.round(roots,3)}")

grid = np.linspace(14.0, gmax, 200)
Ntrue = np.array([np.searchsorted(Z, t, side='right') for t in grid])
Nop   = np.array([np.searchsorted(roots, t, side='right') for t in grid])
Ssm   = Nsmooth(grid)
S_true = Ntrue - Ssm
F_op   = Nop - Ssm

rho, pval = pearsonr(F_op, S_true)
rms = np.sqrt(np.mean((F_op - S_true)**2))
sign_agree = np.mean(np.sign(F_op) == np.sign(S_true))
mean_spacing = 2*np.pi/np.log(gmax/(2*np.pi))
print(f"\nPearson rho(F_op,S_true) = {rho:.4f}  (p={pval:.2e})")
print(f"sign-agreement = {sign_agree*100:.1f}%")
print(f"RMS(F_op-S_true) = {rms:.4f}   mean spacing@{gmax} = {mean_spacing:.4f}  ratio={rms/mean_spacing:.4f}")
print(f"max|F_op|={np.max(np.abs(F_op)):.3f}  max|S_true|={np.max(np.abs(S_true)):.3f}")
print(f"F_op range: [{F_op.min():.3f},{F_op.max():.3f}]  monotone-check (first vs last 10%): "
      f"{np.mean(F_op[:20]):.3f} -> {np.mean(F_op[-20:]):.3f}")
