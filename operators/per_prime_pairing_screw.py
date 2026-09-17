#!/usr/bin/env python3
"""Cycle 11: per-prime decomposition of the screw-frame pairing pair(a)=<v_arch|Q_prime|v_arch>.
g_prime(t)=Sum_n Lambda(n)/sqrt(n)*(|t|-log n)_+ ; per prime-power n contributes
  c_n = Lambda(n)/sqrt(n) * <v_arch| K_n |v_arch>/(a||v_arch||^2),  K_n from kernel (|x-y|-log n)_+.
Shows which n dominate, the tail decay, and total growth vs |lambda_arch(screw)|. Reuses v_arch from
the arch-only ground state (screw_g_no_primes)."""
import numpy as np, suzuki_screw_lambda_a as S
def analyze(a,M=3000,K=60):
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=S.vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; diffs=np.arange(-M,M+1)*h; idx=np.arange(M+1)
    ga=np.array([S.screw_g_no_primes(d,Lam,Nmax) for d in diffs]); Ga=ga[(idx[:,None]-idx[None,:])+M]
    kk=np.arange(1,K+1); arg=np.pi*np.outer((x+a),kk)/(2*a); Vp=(kk*np.pi/(2*a))*np.cos(arg)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h; Wp=Vp*w[:,None]
    Qa=Wp.T@Ga@Wp; Qa=(Qa+Qa.T)/2
    ea,Va=np.linalg.eigh(Qa); c=Va[:,0]; lam_arch=ea[0]/a; nrm=float(c@c)
    # per-prime-power kernels
    contribs={}
    for n in range(2,min(Nmax,int(np.exp(2*a)))+1):
        if Lam[n]>0 and np.log(n)<2*a:
            kn=np.maximum(np.abs(diffs)-np.log(n),0.0)   # (|t|-log n)_+
            Kn=kn[(idx[:,None]-idx[None,:])+M]; Qn=Wp.T@Kn@Wp
            cn=(Lam[n]/np.sqrt(n))*float(c@((Qn+Qn.T)/2)@c)/nrm/a
            contribs[n]=cn
    return lam_arch,contribs
for a in [1.5,2.0,2.5]:
    lam_arch,cc=analyze(a)
    tot=sum(cc.values()); items=sorted(cc.items(),key=lambda kv:-abs(kv[1]))
    print(f"\na={a}: lam_arch={lam_arch:+.4f}  |lam_arch|={abs(lam_arch):.4f}  pair_total={tot:+.4f}  (should ~ -lam_arch)")
    print(f"  top-6 |contrib| by n:  "+"  ".join(f"n={n}:{v:+.4f}" for n,v in items[:6]))
    cum=0; 
    for n,v in sorted(cc.items()):
        cum+=v
    print(f"  #prime-powers={len(cc)}  largest n={max(cc)}  sum(pos)={sum(v for v in cc.values() if v>0):+.3f} sum(neg)={sum(v for v in cc.values() if v<0):+.3f}")
