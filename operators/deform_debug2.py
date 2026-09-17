import sys,os,math
sys.path.insert(0, os.path.join('..','external','groskin_2607.02828','code'))
from flint import arb, ctx
from arb_ldlt_certify import prime_powers_up_to
from dh_certified_tau import build_tau_general, margin
def vonmangoldt(N):
    lam=[0.0]*(N+1); comp=[False]*(N+1)
    for p in range(2,N+1):
        if not comp[p]:
            for m in range(2*p,N+1,p): comp[m]=True
            lp=math.log(p);pk=p
            while pk<=N: lam[pk]=lp;pk*=p
    return lam
c,N,prec=50,64,800; ctx.prec=prec
VM=vonmangoldt(c)
# how many composites have nonzero VM? (should be 0)
print("composites with VM!=0:", [n for n in range(2,c+1) if VM[n]!=0 and any(n%d==0 for d in range(2,n) if d*d<=n and (n//d)*(1)!=n) ]) 
# proper: n is prime power iff single prime factor
def is_pp(n):
    m=n;p=2
    while p*p<=m:
        if m%p==0:
            while m%p==0:m//=p
            return m==1
        p+=1
    return True
print("nonzero VM on non-prime-powers:", [n for n in range(2,c+1) if VM[n]!=0 and not is_pp(n)])
# Path A: prime-power positions only
pdA=prime_powers_up_to(c); wA=[arb(p).log()*(arb(q)**arb('-0.5')) for q,p in pdA]; posA=[arb(q).log() for q,p in pdA]
nA=margin(*build_tau_general(c,N,prec,nu=0.25,include_pole=True,prime_weights=wA,prime_positions=posA))[0]
# Path A' : same prime-power positions but built from VM (all n filtered to VM>0)
nzn=[n for n in range(2,c+1) if VM[n]>0]
wC=[arb(VM[n])*(arb(n)**arb('-0.5')) for n in nzn]; posC=[arb(n).log() for n in nzn]
nC=margin(*build_tau_general(c,N,prec,nu=0.25,include_pole=True,prime_weights=wC,prime_positions=posC))[0]
# Path B: all-n positions, zero on composites
wB=[arb(VM[n])*(arb(n)**arb('-0.5')) for n in range(2,c+1)]; posB=[arb(n).log() for n in range(2,c+1)]
nB=margin(*build_tau_general(c,N,prec,nu=0.25,include_pole=True,prime_weights=wB,prime_positions=posB))[0]
print(f"A prime-powers Groskin-weights: n_neg={nA}")
print(f"C prime-powers VM-weights(VM>0): n_neg={nC}")
print(f"B all-n VM-weights(0 on comp):   n_neg={nB}")
print(f"count: A={len(pdA)} C={len(nzn)} B={len(posB)} positions")
