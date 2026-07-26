#!/usr/bin/env python3
"""§255 falsifier prong (PurpleBot): epsilon-jitter threshold sweep on the telescoping cascade.
Jitter ONE prime position (log 3 -> log 3 + eps), keep weights fixed, and measure at what eps the
telescope (ratio~1.0 each step, deep lam_min) breaks. Threshold curve eps* vs lam_min depth =
the quantitative "arithmetic budget" the identity consumes (§253b shifted comb was a +0.2 jolt;
this shrinks it). Scaffold = beurling_cascade_cvs.py (the RIGHT CvS/Groskin operator, §253b).
Run: python3 epsilon_jitter_sweep.py [c=9] [N=30] [dps=100] [T=250]
"""
import sys, os, mpmath as mp
VENV=os.path.expanduser("~/projects/riemann/cvs_venv/lib/python3.14/site-packages")
sys.path.insert(0, VENV); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import connes_cvs.operator as O
from cvs_fast import precompute_components

c   = int(sys.argv[1]) if len(sys.argv)>1 else 9
N   = int(sys.argv[2]) if len(sys.argv)>2 else 30
dps = int(sys.argv[3]) if len(sys.argv)>3 else 100
T   = int(sys.argv[4]) if len(sys.argv)>4 else 250
mp.mp.dps=dps
comp = precompute_components(c, N, T, dps)
idx=comp["idx"]; L=comp["L"]; DIM=2*N+1
Vp=mp.matrix(DIM,N+1); Vp[N,0]=mp.mpf(1); r2=1/mp.sqrt(2)
for k in range(1,N+1): Vp[N+k,k]=r2; Vp[N-k,k]=r2
arc={n:comp["arcpol"][n] for n in idx}; arcd={n:comp["arcpol_d"][n] for n in idx}

def assemble(valmap,dermap):
    Q=mp.matrix(DIM,DIM)
    for i in range(DIM):
        m=i-N
        for j in range(DIM):
            nn=j-N
            Q[i,j]=dermap[nn] if m==nn else (valmap[m]-valmap[nn])/(m-nn)
    for i in range(DIM):
        for j in range(i+1,DIM):
            a=(Q[i,j]+Q[j,i])/2; Q[i,j]=a; Q[j,i]=a
    return Vp.T*Q*Vp

def cascade(label, steps):
    """steps = ordered [(name,[(n,logpos,w)])]. Returns (log10 lam_min, max|ratio-1|, ratios)."""
    allpd=[pd[0] for _,pd in steps]
    fv={}; fd={}
    for n in idx:
        fv[n]=arc[n]+sum(O.psi_prime(mp.mpf(n),L,[t]) for t in allpd)
        fd[n]=arcd[n]+sum(O.psi_prime_deriv(mp.mpf(n),L,[t]) for t in allpd)
    Qe=assemble(fv,fd); eigs,vecs=mp.eigsy(Qe)
    k0=min(range(N+1),key=lambda i:eigs[i])
    psi=mp.matrix([vecs[r,k0] for r in range(N+1)]); lam0=eigs[k0]
    Qa=assemble(arc,arcd); Earch=(psi.T*Qa*psi)[0]
    hs=[]
    for name,pd in steps:
        val={n:O.psi_prime(mp.mpf(n),L,pd) for n in idx}
        der={n:O.psi_prime_deriv(mp.mpf(n),L,pd) for n in idx}
        h=(psi.T*assemble(val,der)*psi)[0]; hs.append((name,h))
    R=Earch; ratios=[]
    for i,(name,h) in enumerate(hs):
        R=R+h
        nh=abs(hs[i+1][1]) if i+1<len(hs) else None
        if nh and nh>0: ratios.append((name,float(abs(R)/nh)))
    llam=float(mp.log10(abs(lam0))) if lam0!=0 else -99.0
    dev=max((abs(r-1.0) for _,r in ratios), default=float('nan'))
    print(f"{label:>12}  lam_min log10={llam:8.2f}  max|ratio-1|={dev:10.4f}  ratios="+
          " ".join(f"{n}:{r:.3f}" for n,r in ratios), flush=True)
    return llam,dev,ratios

pdata,_=O.prime_powers_up_to(c)
def steps_with_jitter(eps):
    out=[]
    for (n,logn,w) in pdata:
        lp = logn + eps if n==3 else logn   # jitter ONLY log 3; weights untouched
        out.append((f"n{n}",[(n,lp,w)]))
    return out

print(f"# epsilon-jitter sweep on log 3 | c={c} N={N} dps={dps} L={mp.nstr(L,6)} | primes: "+
      ",".join(str(n) for n,_,_ in pdata), flush=True)
EPS = sys.argv[5:] or ["0","0.0001","0.0003","0.001","0.003","0.01","0.03","0.1","0.2"]
for eps in EPS:
    cascade(f"eps={eps}", steps_with_jitter(mp.mpf(eps)))
