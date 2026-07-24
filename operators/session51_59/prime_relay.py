# PRIME RELAY test: does min(tau_EVEN) dip toward 0 as lam^2 -> 3- (one prime), then get RESCUED
# when p=3 enters at lam^2=3+? If yes: RH-positivity = successive prime rescues of the even ground state.
import mpmath as mp
mp.mp.dps=30
def min_tau_even(LSQ,N=10):
    L=2*mp.log(mp.sqrt(LSQ)); TWO_PI=2*mp.pi; LOGPI=mp.log(mp.pi)
    idx=list(range(-N,N+1)); D=len(idx)
    Rq=TWO_PI*N/L+50; G=1000
    grid=[(-Rq+2*Rq*mp.mpf(i)/G) for i in range(G+1)]; hh=2*Rq/G
    wq=[hh]*(G+1); wq[0]=hh/2; wq[-1]=hh/2
    Om=[mp.re(mp.digamma(mp.mpf('0.25')+0.5j*r))-LOGPI for r in grid]
    def ehat(r,n):
        s=TWO_PI*n/L; d=s-r
        if abs(d)<mp.mpf('1e-24'): return mp.sqrt(L)
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
        return mp.log(p) if(ok and p>0) else mp.mpf(0)
    VM=[(k,vm(k)) for k in range(2,int(mp.floor(LSQ))+1) if vm(k)>0]
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
    return min(mp.eigsy(BE,eigvals_only=True))
print("lam^2   nprimes  min(tau_EVEN)   <- does it dip to 0 near 3-, jump up at 3+ (p=3 enters)?",flush=True)
for LSQ in ['2.80','2.90','2.97','2.99','3.01','3.10','3.30']:
    L=mp.mpf(LSQ); npr=len([k for k in range(2,int(mp.floor(L))+1) if k in (2,3,5,7,11,13) or True]) 
    npr=int(mp.floor(L))-1
    mte=min_tau_even(L)
    flag=" <== p=3 ENTERS" if LSQ=='3.01' else ""
    print(f"{LSQ}    {npr}       {mp.nstr(mte,5):>12}{flag}",flush=True)
