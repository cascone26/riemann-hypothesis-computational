#!/usr/bin/env python3
"""Cycle 12: Temple's inequality rigorous LOWER bound on lambda_full with trial vector v_arch.
H=full op (Qfull, mass a*I). Trial psi=c_arch (arch ground eigvec). Temple:
  lambda_full >= m - sigma^2/(eta - m),   m=<psi|H|psi>, sigma^2=<H^2>-m^2, eta<=lambda_2(full), eta>m.
If bound>0 on a range of a => RIGOROUS lambda_full>0 there (beyond Suzuki's small-a Thm 1.4).
Exploratory eta = true lambda_2(full); a real proof needs a rigorous eta lower bound (next step)."""
import numpy as np, suzuki_screw_lambda_a as S
def temple(a,M=3000,K=60):
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=S.vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; diffs=np.arange(-M,M+1)*h; idx=np.arange(M+1)
    gf=np.array([S.screw_g(d,Lam,Nmax) for d in diffs]); ga=np.array([S.screw_g_no_primes(d,Lam,Nmax) for d in diffs])
    Gf=gf[(idx[:,None]-idx[None,:])+M]; Ga=ga[(idx[:,None]-idx[None,:])+M]
    kk=np.arange(1,K+1); arg=np.pi*np.outer((x+a),kk)/(2*a); Vp=(kk*np.pi/(2*a))*np.cos(arg)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h; Wp=Vp*w[:,None]
    Qf=Wp.T@Gf@Wp; Qf=(Qf+Qf.T)/2; Qa=Wp.T@Ga@Wp; Qa=(Qa+Qa.T)/2
    ea,Va=np.linalg.eigh(Qa); psi=Va[:,0]
    evf=np.linalg.eigvalsh(Qf); lam_full=evf[0]/a; lam2=evf[1]/a
    nn=float(psi@psi)
    m=float(psi@Qf@psi)/(a*nn)                      # Rayleigh of full op at arch vector
    H2=float(psi@Qf@Qf@psi)/(a*a*nn)                # <H^2>
    sig2=H2-m*m
    eta=lam2
    bound=m - sig2/(eta-m) if eta>m else float('nan')
    return lam_full,lam2,m,sig2,eta,bound
print(f"{'a':>5} {'lam_full':>9} {'lam2':>8} {'m(trial)':>9} {'sigma^2':>9} {'eta':>7} {'Temple LB':>10}")
for a in [0.5,0.8,1.0,1.3,1.6,2.0,2.5]:
    lf,l2,m,s2,eta,b=temple(a)
    print(f"{a:5.1f} {lf:+9.4f} {l2:+8.4f} {m:+9.4f} {s2:9.4f} {eta:+7.3f} {b:+10.4f}")
print("\nTemple LB>0 => rigorous lambda_full>0 (given eta a valid lower bd on lambda_2). Sign+size of LB is the result.")
