import numpy as np
def build_dL(N):
    LIM=N+120; spf=np.zeros(LIM,dtype=np.int32)
    for i in range(2,LIM):
        if spf[i]==0: spf[i:LIM:i]=np.where(spf[i:LIM:i]==0,i,spf[i:LIM:i])
    Lam=np.zeros(LIM)
    for n in range(2,LIM):
        p=int(spf[n]); m=n
        while m%p==0: m//=p
        if m==1: Lam[n]=np.log(p)
    return Lam-1.0
def admissible(t): return all(len({x%p for x in t})<p for p in (3,5,7))
def triples(s): return [(a,b,c) for a in range(1,s) for b in range(a+1,s)
                        for c in (s-a-b,) if c>b and admissible((a,b,c))]
s=36; T=triples(s)
# how many off-diagonal pairs share exactly 1 offset value vs 0?
n=len(T); n0=n1=0
for i in range(n):
    for j in range(i+1,n):
        sh=len(set(T[i])&set(T[j]))
        if sh==0: n0+=1
        elif sh==1: n1+=1
        else: print("SHARED>=2!",T[i],T[j])
print(f"s={s} n={n}: off-diag pairs  0-shared={n0}  1-shared={n1}  ({100*n1/(n0+n1):.0f}% share an offset)")
print("\nC (=E[u_t u_t']) and R by shared-count, vs N -- is the 1-shared group the slow rate-setter?")
Ns=[100000,1000000,10000000,40000000]
for N in Ns:
    dL=build_dL(N); lo,hi=2,N-60; idx=np.arange(lo,hi)
    U=[dL[idx+t[0]]*dL[idx+t[1]]*dL[idx+t[2]] for t in T]
    D=np.array([ (U[i]*U[i]).mean() for i in range(n)])
    C0=[];C1=[];R0=[];R1=[]
    for i in range(n):
        for j in range(i+1,n):
            c=(U[i]*U[j]).mean(); r=c/np.sqrt(D[i]*D[j]); sh=len(set(T[i])&set(T[j]))
            (C0 if sh==0 else C1).append(c); (R0 if sh==0 else R1).append(r)
    rc0=np.sqrt(np.mean(np.array(C0)**2)); rc1=np.sqrt(np.mean(np.array(C1)**2))
    rr0=np.sqrt(np.mean(np.array(R0)**2)); rr1=np.sqrt(np.mean(np.array(R1)**2))
    print(f"  N={N:>9}  logN={np.log(N):5.2f}  D_min={D.min():7.2f}  |  rmsC0={rc0:6.3f} rmsR0={rr0:.4f}  |  rmsC1={rc1:7.2f} rmsR1={rr1:.4f}")
print("\n  Lemma2 covers 0-shared (C0 should be O(1)). 1-shared: C1 carries a dLambda^2 -> ~logN?  R1 is the (logN)^-1.3 rate-setter.")
print("  proof bound for 1-shared: |C|~logN over D>=(1/4)logN  =>  R=O(1), NOT o(1). Does R1 empirically ->0 anyway?")
