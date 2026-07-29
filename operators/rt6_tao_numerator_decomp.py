"""Tao: decompose sigma_off(N) = num_rms / diag  across N, to adjudicate the
'numerator is N-flat' claim the whole power-vs-log fork hinges on.
kappa_eff(sigma_off) = kappa(numerator) - kappa(diagonal). Riemann/Hilbert:
kappa(num)=0, kappa(diag)~3.5 => sigma_off slope ~ -0.25. Tao measured -0.083.
Direct decomposition settles which. s=42 sector, chunked Gram, no fit-extrapolation."""
import json, numpy as np
from math import log
def sieve_Lambda(M):
    spf=np.zeros(M+1,dtype=np.int64)
    for i in range(2,M+1):
        if spf[i]==0:
            for j in range(i,M+1,i):
                if spf[j]==0: spf[j]=i
    L=np.zeros(M+1)
    for n in range(2,M+1):
        p=int(spf[n]); m=n
        while m%p==0: m//=p
        if m==1: L[n]=log(p)
    return L
def triples_at_sum(s,lo=1):
    return [(a,b,s-a-b) for a in range(lo,s) for b in range(a+1,s) if s-a-b>b]
def admissible(t):
    for p in (3,5,7):
        if len(set(x%p for x in t))==p: return False
    return True
s=42
T=[t for t in triples_at_sum(s) if admissible(t)]; n=len(T)
NMAX=6_400_000; MAXO=60
Lam=sieve_Lambda(NMAX+MAXO+5); dL=Lam-1.0
offmask=~np.eye(n,dtype=bool)
def gram(N,chunk=250_000):
    G=np.zeros((n,n)); done=0
    while done<N:
        c=min(chunk,N-done); base=np.arange(done,done+c)
        Vc=np.array([dL[base+t[0]]*dL[base+t[1]]*dL[base+t[2]] for t in T])
        G+=Vc@Vc.T; done+=c
    return G/N
Ns=[100_000,200_000,400_000,800_000,1_600_000,3_200_000,6_400_000]
rows=[]
print(f"{'N':>9} {'num_rms':>9} {'diag_mean':>10} {'sig_off':>8}")
for N in Ns:
    G=gram(N); diag=np.diag(G); num=np.sqrt((G[offmask]**2).mean())
    d=np.sqrt(diag); R=G/np.outer(d,d); so=np.sqrt((R[offmask]**2).mean())
    rows.append((N,num,diag.mean(),so)); print(f"{N:>9} {num:>9.3f} {diag.mean():>10.2f} {so:>8.5f}")
def kappa_log(xs,ys):  # fit y = c*(log x)^kappa -> slope of log y vs log(log x)
    lx=np.log(np.log(np.array(xs))); ly=np.log(np.array(ys))
    return np.polyfit(lx,ly,1)[0]
def alpha_pow(xs,ys):  # fit y = c*x^-alpha -> -slope of log y vs log x
    return -np.polyfit(np.log(xs),np.log(ys),1)[0]
Ns_a=[r[0] for r in rows]; num=[r[1] for r in rows]; diag=[r[2] for r in rows]; so=[r[3] for r in rows]
# use last 5 rungs (past floor crossover)
sl=slice(2,None)
kn=kappa_log(Ns_a[sl],num[sl]); kd=kappa_log(Ns_a[sl],diag[sl]); ks=kappa_log(Ns_a[sl],so[sl])
print(f"\nlog-family exponents (last 5 rungs, y~(logN)^k):")
print(f"  kappa(numerator) = {kn:+.2f}   kappa(diagonal) = {kd:+.2f}   kappa(sigma_off)= {ks:+.2f}")
print(f"  check: k(num)-k(diag) = {kn-kd:+.2f}  vs measured k(sig_off) = {ks:+.2f}")
print(f"power-law: alpha(sigma_off) = {alpha_pow(Ns_a[sl],so[sl]):.3f}")
json.dump({"object":"Tao numerator/diagonal decomposition of sigma_off(N)",
 "version":"rt6-tao-numerator-decomp-2026-07-28","sector_s":s,"n":n,
 "ladder":[{"N":r[0],"num_rms":r[1],"diag_mean":r[2],"sig_off":r[3]} for r in rows],
 "kappa_numerator":kn,"kappa_diagonal":kd,"kappa_sigma_off":ks,
 "note":"if kappa(num)~0 and kappa(diag)~3.5 the sigma_off decay is steep-log; if numerator ALSO grows the effective decay is shallow"},
 open("results/rt6_tao_numerator_decomp.json","w"),indent=2)
print("wrote results/rt6_tao_numerator_decomp.json")
