#!/usr/bin/env python3
"""
subgoalH_dbn_trunc.py — LEDGER §158 (verify §157's dBN monotone signal): is dλ_min/ds truncation-STABLE, or an artifact
of the un-renormalized velocity field v_j=2Σ_{k≠j}1/(γ_j−γ_k) (conditionally convergent, biased by the cutoff)? Compute
the slope at increasing zero counts N=1500,3000,6000,10000. STABLE ⟹ real flow-monotonicity (§157 candidate holds).
DRIFTS with N ⟹ the monotone signal is a truncation/mean-field-drift artifact, not the true dBN flow.
Also: compare the collective-flow λ_min change to what a naive §129 single-zero picture would give (is it different?).
"""
import numpy as np

ALL = np.loadtxt("/Users/scones/projects/riemann/data/zeros_100k.txt", max_rows=10000)

def lam_min_flowed(a, s, Nz, M=1400):
    G = ALL[:Nz]
    D = G[:,None]-G[None,:]; np.fill_diagonal(D, np.inf); V=2.0*np.sum(1.0/D,axis=1)
    g = G + s*V
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    dif=np.arange(-M,M+1)*h
    C0=-2.0*np.sum(1.0/g**2); gv=C0+2.0*np.sum(np.cos(np.outer(dif,g))/g**2,axis=1)
    idx=np.arange(M+1); Gm_=gv[(idx[:,None]-idx[None,:])+M]
    X=x[:,None]; Y=x[None,:]; N=(X**2+Y**2)/(4*a)-np.abs(X-Y)/2+a/6
    sw=np.sqrt(w); Gm=sw[:,None]*Gm_*sw[None,:]; Gm=(Gm+Gm.T)/2; Nm=sw[:,None]*N*sw[None,:]; Nm=(Nm+Nm.T)/2
    c=sw/np.linalg.norm(sw); P=np.eye(M+1)-np.outer(c,c); Gp=P@Gm@P; Np=P@Nm@P
    ev_N,V_N=np.linalg.eigh((Np+Np.T)/2); keep=ev_N>1e-12*ev_N.max()
    S=V_N[:,keep]/np.sqrt(ev_N[keep]); Ar=S.T@((Gp+Gp.T)/2)@S
    return np.linalg.eigvalsh((Ar+Ar.T)/2)[0]

a=1.5; ds=0.02
print(f"§158 truncation test of dλ_min/ds under the dBN flow (a={a}):\n")
print(f"{'N zeros':>8} {'λ_min(s=0)':>13} {'λ_min(-ds)':>13} {'λ_min(+ds)':>13} {'slope dλ/ds':>13}")
for Nz in [1500, 3000, 6000, 10000]:
    l0=lam_min_flowed(a,0.0,Nz); lm=lam_min_flowed(a,-ds,Nz); lp=lam_min_flowed(a,+ds,Nz)
    slope=(lp-lm)/(2*ds)
    print(f"{Nz:>8} {l0:>13.5e} {lm:>13.5e} {lp:>13.5e} {slope:>13.4e}",flush=True)
print("\nREAD: if 'slope' converges as N↑ ⟹ REAL dBN flow-monotonicity of the arithmetic λ_min (§157 candidate solid).")
print("If 'slope' keeps drifting with N ⟹ the monotone signal is the un-renormalized velocity's truncation artifact —")
print("the true dBN velocity needs mean-field renormalization, and the naive collective flow signal is not trustworthy.")
