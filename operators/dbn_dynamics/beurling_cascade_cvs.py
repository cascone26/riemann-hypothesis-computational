#!/usr/bin/env python3
"""§253 (correct): near-null TELESCOPING CASCADE for REAL vs BEURLING primes, on the RIGHT CvS/Groskin operator.
Reuses connes_cvs (the §176 machinery). Beurling = arbitrary (n, logpos, weight) triples into psi_prime.
Q: does R_after(step) ~ |h(next)| telescope (real, §176b) or HALT (Beurling, §L2 prediction)?
Run: python3 beurling_cascade_cvs.py [c=7] [N=30] [dps=80] [T=250]
"""
import sys, os, mpmath as mp
VENV=os.path.expanduser("~/projects/riemann/cvs_venv/lib/python3.14/site-packages")
sys.path.insert(0, VENV); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import connes_cvs.operator as O
from cvs_fast import precompute_components

c   = int(sys.argv[1]) if len(sys.argv)>1 else 7
N   = int(sys.argv[2]) if len(sys.argv)>2 else 30
dps = int(sys.argv[3]) if len(sys.argv)>3 else 80
T   = int(sys.argv[4]) if len(sys.argv)>4 else 250
mp.mp.dps=dps
comp = precompute_components(c, N, T, dps)
idx=comp["idx"]; L=comp["L"]; DIM=2*N+1
Vp=mp.matrix(DIM,N+1); Vp[N,0]=mp.mpf(1); r2=1/mp.sqrt(2)
for k in range(1,N+1): Vp[N+k,k]=r2; Vp[N-k,k]=r2
arc={n:comp["arcpol"][n] for n in idx}; arcd={n:comp["arcpol_d"][n] for n in idx}

def prime_maps(prime_data):
    val={n:O.psi_prime(mp.mpf(n),L,prime_data) for n in idx}
    der={n:O.psi_prime_deriv(mp.mpf(n),L,prime_data) for n in idx}
    return val,der

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

def run(label, steps):
    """steps = ordered list of (name, [(n,logpos,w)]) single-'prime' data. Build full op, cascade."""
    allpd=[pd[0] for _,pd in steps]  # flatten single triples
    fv={}; fd={}
    for n in idx:
        fv[n]=arc[n]+sum(O.psi_prime(mp.mpf(n),L,[t]) for t in allpd)
        fd[n]=arcd[n]+sum(O.psi_prime_deriv(mp.mpf(n),L,[t]) for t in allpd)
    Qe=assemble(fv,fd); eigs,vecs=mp.eigsy(Qe)
    k0=min(range(N+1),key=lambda i:eigs[i])
    psi=mp.matrix([vecs[r,k0] for r in range(N+1)]); lam0=eigs[k0]
    # arch-only expectation
    av,ad=arc,arcd; Qa=assemble(av,ad); Earch=(psi.T*Qa*psi)[0]
    print(f"\n== {label} ==  lam_min log10={float(mp.log10(abs(lam0))):.2f}  arch<psi|.|psi>={mp.nstr(Earch,5)}")
    print(f"   {'step':>10} {'h=<psi|Vn|psi>':>16} {'log10|R_after|':>15} {'|R|/|h_next|':>13}")
    R=Earch; hs=[]
    for name,pd in steps:
        Vv,Vd=prime_maps(pd); Qn=assemble(Vv,Vd); h=(psi.T*Qn*psi)[0]; hs.append((name,h))
    R=Earch
    print(f"   {'(arch)':>10} {mp.nstr(Earch,5):>16} {float(mp.log10(abs(Earch))):15.2f}")
    for i,(name,h) in enumerate(hs):
        R=R+h  # lam_min = Earch + sum_n h(n); each h peels off the residual if it telescopes
        nh=abs(hs[i+1][1]) if i+1<len(hs) else None
        ratio=float(abs(R)/nh) if nh and nh>0 else float('nan')
        lr=float(mp.log10(abs(R))) if R!=0 else -99
        lh=float(mp.log10(abs(h))) if h!=0 else -99
        print(f"   {name:>10} {mp.nstr(h,5):>16} {lr:15.2f} {ratio:13.3f}   log10|h|={lh:.2f}")
    return lam0

# REAL prime powers up to c
pdata,_=O.prime_powers_up_to(c)
real=[(f"n{n}",[(n,logn,w)]) for (n,logn,w) in pdata]
run("REAL primes (expect telescope: ratio~1.0 each step)", real)

# BEURLING B1: powers of 2 only, same log-weight form (RH-false; zeta_B zeros off critical line)
b1=[]
k=1
while k*mp.log(2)<L*0.999:
    n2=2**k; b1.append((f"2^{k}",[(n2, k*mp.log(2), mp.log(2)/mp.sqrt(n2))])); k+=1
run("BEURLING {2^k} (RH-false; expect HALT if cascade=integer-fingerprint)", b1)

# BEURLING B2: real weights, positions shifted +0.2 (breaks additive/PNT coherence)
b2=[(f"s{n}",[(n, logn+mp.mpf('0.2'), w)]) for (n,logn,w) in pdata if float(logn)+0.2<float(L)*0.999]
run("BEURLING shifted +0.2 (breaks additive coherence)", b2)
print("\nREAD: REAL ratios ~1.0 (telescope, §176b). If BEURLING ratios scatter/blow up => cascade HALTS =")
print("arithmetic-coherence fingerprint of the integers (RH-true) = concrete 'what a proof must consume' (§182).")
