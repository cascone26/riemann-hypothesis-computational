#!/usr/bin/env python3
"""
subgoalH_groundstate_zeros.py — LEDGER §161 (C→A: the CORRECT CvS object): the proven foothold (Connes–van Suijlekom
Thm 6.1) is about the zeros of η̂_c = Fourier transform of the GROUND-STATE eigenvector (even sector) — a SPARSE set —
proven on the critical line (⟺ real z) at every finite c via Carathéodory–Fejér (PSD form, one near-null direction). The
OPEN crux: do these ground-state zeros converge to the true γ as the cutoff→∞? (No rate analysis exists.)
§160/g_a_W computed the WRONG object (whole D_{a,θ} spectrum = dense comb). HERE: extract the ground state ψ_1 of the
compact Weil operator, form ψ̂_1(z)=∫ψ_1(x)e^{izx}dx, find its REAL zeros, (1) verify they are REAL (foothold on our
tooling), (2) compare to γ_1,γ_2,... and track convergence vs a — the concrete CvS convergence, correctly specified.
"""
import numpy as np, sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

GAM=np.array([14.134725,21.022040,25.010858,30.424876,32.935062,37.586178,40.918719,43.327073,48.005151,49.773832])

def ground_state(a, M=1800):
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    dif=np.arange(-M,M+1)*h; gv=np.array([screw_g(d,Lam,Nmax) for d in dif]); idx=np.arange(M+1)
    G=gv[(idx[:,None]-idx[None,:])+M]
    X=x[:,None]; Y=x[None,:]; N=(X**2+Y**2)/(4*a)-np.abs(X-Y)/2+a/6
    sw=np.sqrt(w); Gm=sw[:,None]*G*sw[None,:]; Gm=(Gm+Gm.T)/2; Nm=sw[:,None]*N*sw[None,:]; Nm=(Nm+Nm.T)/2
    c=sw/np.linalg.norm(sw); P=np.eye(M+1)-np.outer(c,c); Gp=P@Gm@P; Np=P@Nm@P
    ev_N,V_N=np.linalg.eigh((Np+Np.T)/2); keep=ev_N>1e-12*ev_N.max()
    S=V_N[:,keep]/np.sqrt(ev_N[keep]); Ar=S.T@((Gp+Gp.T)/2)@S
    wv,U=np.linalg.eigh((Ar+Ar.T)/2); f=(S@U[:,0]); f=f/np.where(sw>0,sw,1.0)
    return x, f, w, wv[0]

def psihat(x, f, w, z):
    # ∫ f(x) e^{izx} dx ; f is real → split; find zeros of the ANALYTIC-signal magnitude on real z after de-phasing
    return np.array([np.sum(f*np.exp(1j*zz*x)*w) for zz in z])

def find_zeros(z, vals):
    # zeros of a real-ish function: use the real part if f even (ψ̂ real) or |·| minima
    a=np.abs(vals); med=np.median(a); zs=[]
    for i in range(1,len(z)-1):
        if a[i]<a[i-1] and a[i]<a[i+1] and a[i]<0.1*med:
            d=(a[i-1]-a[i+1])/(2*(a[i-1]-2*a[i]+a[i+1])+1e-30); zs.append(z[i]+d*(z[1]-z[0]))
    return np.array(zs)

print("§161 CvS ground-state zeros: ψ̂_1(z) zeros — REAL (foothold)? converging to γ? (correct sparse object)\n")
zg=np.linspace(0.3,50,6000)
print(f"{'a':>5} {'λ_min':>11} {'parity':>7} {'#zeros<50':>9} {'ψ̂_1 zeros (first 6)':>34}  {'|z-γ| (γ1..4)':>16}")
for a in [1.5,2.0,2.5,3.0]:
    x,f,w,lm=ground_state(a)
    par = 'even' if np.sum(f*f[::-1]*w)>0 else 'odd'
    vals=psihat(x,f,w,zg); zs=find_zeros(zg,vals)
    zs=zs[zs>1.0]
    match=[float(np.min(np.abs(zs-g))) if len(zs) else np.nan for g in GAM[:4]]
    z6=[round(z,3) for z in zs[:6]]
    print(f"{a:>5.1f} {lm:>11.3e} {par:>7} {len(zs):>9} {str(z6):>34}  {str([round(m,2) for m in match]):>16}",flush=True)
print("\nREAD: if ψ̂_1 zeros are sparse & real & the low ones sit near γ_1,γ_2,... (|z−γ| small, TIGHTENING with a) ⟹ we")
print("reproduce the CvS convergence on our tooling with the CORRECT object. Then our tangency(§154)/collectivity(§156)")
print("findings can target WHY/whether the zeros converge — the open crux with no rate analysis. If NOT near γ ⟹ the")
print("compact-G_a ground state's transform differs from the CvS Fourier-Mellin one; refine the transform.")
