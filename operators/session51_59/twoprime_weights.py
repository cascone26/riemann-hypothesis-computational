# TWO-PRIME weight-criticality: lam^2 in [3,5), primes 2 & 3 active. Sweep w2, w3 INDEPENDENTLY
# around their arithmetic values (w2=log2/sqrt2, w3=log3/sqrt3). Q: does the true (w2,w3) sit at the
# positivity EDGE in each direction (min(tau_EVEN) pinned ~0 at arithmetic, negative just above)?
import mpmath as mp
mp.mp.dps=28
def min_even(LSQ,f2,f3,N=12):
    # f2,f3 = multipliers on arithmetic weights for p=2, p=3 (and p-powers use their own Lambda)
    L=2*mp.log(mp.sqrt(LSQ)); TWO_PI=2*mp.pi; LOGPI=mp.log(mp.pi)
    idx=list(range(-N,N+1)); D=len(idx)
    Rq=TWO_PI*N/L+55; G=900
    grid=[(-Rq+2*Rq*mp.mpf(i)/G) for i in range(G+1)]; hh=2*Rq/G
    wq=[hh]*(G+1); wq[0]=hh/2; wq[-1]=hh/2
    Om=[mp.re(mp.digamma(mp.mpf('0.25')+0.5j*r))-LOGPI for r in grid]
    def ehat(r,n):
        s=TWO_PI*n/L; d=s-r
        if abs(d)<mp.mpf('1e-22'): return mp.sqrt(L)
        return (1/mp.sqrt(L))*mp.e**(0.5j*d*L)*2*mp.sin(0.5*d*L)/d
    E=[[ehat(r,n) for r in grid] for n in idx]
    def arch(i,j): return mp.re(mp.fsum(E[i][g]*mp.conj(E[j][g])*wq[g]*Om[g] for g in range(G+1)))/(2*mp.pi)
    sh=mp.sinh(L/4)**2
    def Bnm(n,m):
        num=L**2-16*mp.pi**2*m*n; den=(L**2+16*mp.pi**2*m**2)*(L**2+16*mp.pi**2*n**2); return 32*L*sh*num/den
    def vm(k):
        n=k;p=0;f=2;ok=True
        while f*f<=n:
            if n%f==0:
                if p==0:p=f
                elif p!=f:ok=False;break
                while n%f==0:n//=f
            f+=1
        if n>1:
            if p==0:p=n
            elif p!=n:ok=False
        return (mp.log(p),p) if(ok and p>0) else (mp.mpf(0),0)
    # prime-power list with base prime, apply multiplier per base
    VM=[]
    for k in range(2,int(mp.floor(LSQ))+1):
        lam,base=vm(k)
        if lam>0:
            mult=f2 if base==2 else (f3 if base==3 else mp.mpf(1))
            VM.append((k,lam*mult))
    def qk(n,m,y):
        if n==m: return 2*(1-abs(y)/L)*mp.cos(TWO_PI*n*y/L)
        return (mp.sin(TWO_PI*m*y/L)-mp.sin(TWO_PI*n*y/L))/(mp.pi*(n-m))
    def Pnm(n,m): return mp.fsum(lp*mp.mpf(k)**mp.mpf('-0.5')*qk(n,m,mp.log(k)) for k,lp in VM)
    T=mp.matrix(D,D)
    for i,n in enumerate(idx):
        for j,m in enumerate(idx): T[i,j]=arch(i,j)+Bnm(n,m)-Pnm(n,m)
    def col_e(n):
        v=mp.matrix(D,1)
        if n==0: v[idx.index(0)]=1
        else: v[idx.index(n)]=1/mp.sqrt(2); v[idx.index(-n)]=1/mp.sqrt(2)
        return v
    Ec=[col_e(0)]+[col_e(n) for n in range(1,N+1)]; k=len(Ec)
    BE=mp.matrix(k,k)
    for a in range(k):
        for b in range(k): BE[a,b]=mp.fsum(Ec[a][i]*T[i,j]*Ec[b][j] for i in range(D) for j in range(D))
    return float(min(mp.eigsy(BE,eigvals_only=True)))
LSQ=mp.mpf('4.0')  # primes 2,3 (and 4=2^2) active
print(f"lam^2=4.0 (primes 2,3,4=2^2). arithmetic = (f2,f3)=(1,1). min(tau_EVEN):",flush=True)
print(f"  vary f2 (f3=1):   " + "  ".join(f"f2={f}:{min_even(LSQ,mp.mpf(str(f)),mp.mpf('1')):+.3e}" for f in [0.9,1.0,1.1]),flush=True)
print(f"  vary f3 (f2=1):   " + "  ".join(f"f3={f}:{min_even(LSQ,mp.mpf('1'),mp.mpf(str(f))):+.3e}" for f in [0.9,1.0,1.1]),flush=True)
print(f"  both true (1,1):  {min_even(LSQ,mp.mpf('1'),mp.mpf('1')):+.4e}",flush=True)
