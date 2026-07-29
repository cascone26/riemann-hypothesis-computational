import numpy as np
def build_dL(N):
    LIM=N+80
    spf=np.zeros(LIM,dtype=np.int32)
    for i in range(2,LIM):
        if spf[i]==0:
            spf[i:LIM:i]=np.where(spf[i:LIM:i]==0,i,spf[i:LIM:i])
    Lam=np.zeros(LIM)
    # prime powers: n with n/spf^k ==1
    for n in range(2,LIM):
        p=int(spf[n]); m=n
        while m%p==0: m//=p
        if m==1: Lam[n]=np.log(p)
    return Lam-1.0

def admissible(t):
    return all(len({x%p for x in t})<p for p in (3,5,7))
def triples(s):
    return [(a,b,c) for a in range(1,s) for b in range(a+1,s)
            for c in (s-a-b,) if c>b and admissible((a,b,c))]

# fixed sector s=42; measure numerator E[u_t u_t'] and diagonal E[u_t^2] vs N
S=42; T=triples(S)
pairs=[(0,5),(1,7),(2,9),(3,11),(4,13)]  # a spread of shape-mismatched off-diag pairs
Ns=[100000,1000000,10000000]
print(f"s={S}, n_triples={len(T)}")
print(f"{'N':>10} | {'E[delta^2]':>10} {'diag E[u^2]':>11} | " + " ".join(f"num{p}" for p in range(len(pairs))))
diag_by_N=[]; num_by_N=[[] for _ in pairs]
for N in Ns:
    dL=build_dL(N); n=np.arange(2,N-80)
    U=[dL[n+t[0]]*dL[n+t[1]]*dL[n+t[2]] for t in T]
    d2=(dL[2:N-80]**2).mean()
    diag=np.mean([ (U[i]*U[i]).mean() for i in range(len(T))])
    diag_by_N.append(diag)
    nums=[]
    for k,(i,j) in enumerate(pairs):
        v=(U[i]*U[j]).mean(); nums.append(v); num_by_N[k].append(v)
    print(f"{N:>10} | {d2:>10.3f} {diag:>11.2f} | " + " ".join(f"{x:>6.2f}" for x in nums))

lnN=np.log(Ns)
kappa_diag=np.polyfit(np.log(lnN),np.log(diag_by_N),1)[0]
print(f"\ndiagonal E[u^2] ~ (log N)^kappa, kappa_diag = {kappa_diag:.2f}")
for k,(i,j) in enumerate(pairs):
    arr=np.array(num_by_N[k])
    trend=np.polyfit(lnN,arr,1)[0]
    print(f"  numerator pair{k} {T[i]}x{T[j]}: {arr.round(3)}  linear-in-lnN slope={trend:+.3f} (flat=>N-indep)")
print("\nVERDICT: if diag grows ~(logN)^2-3 and numerators are flat/O(1) with no N-trend,")
print("then sigma_off = num/diag ~ 1/(logN)^kappa = LOG-FAMILY (Riemann's structure), not power.")
