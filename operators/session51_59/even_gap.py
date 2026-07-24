# EVEN-BLOCK SPECTRAL GAP above the pinned null mode k_lambda. lowest 4 even eigenvalues of tau.
# Q: is the gap (2nd even eig) bounded below by a PRIME-built positive quantity? And does the prime
# term OPEN the gap (arch-only gap vs full gap)? That would be the analytic lever for even-positivity.
import mpmath as mp
mp.mp.dps=32
def even_spectrum(LSQ,N=12,use_primes=True):
    L=2*mp.log(mp.sqrt(LSQ)); TWO_PI=2*mp.pi; LOGPI=mp.log(mp.pi)
    idx=list(range(-N,N+1)); D=len(idx)
    Rq=TWO_PI*N/L+50; G=1100
    grid=[(-Rq+2*Rq*mp.mpf(i)/G) for i in range(G+1)]; hh=2*Rq/G
    wq=[hh]*(G+1); wq[0]=hh/2; wq[-1]=hh/2
    Om=[mp.re(mp.digamma(mp.mpf('0.25')+0.5j*r))-LOGPI for r in grid]
    def ehat(r,n):
        s=TWO_PI*n/L; d=s-r
        if abs(d)<mp.mpf('1e-26'): return mp.sqrt(L)
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
    VM=[(k,vm(k)) for k in range(2,int(mp.floor(LSQ))+1) if vm(k)>0] if use_primes else []
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
    ev=sorted(mp.eigsy(BE,eigvals_only=True))
    return [float(ev[t]) for t in range(4)]
print(f"{'lam^2':>6} | {'FULL: e0':>10} {'e1(gap)':>10} {'e2':>9} | {'ARCH: e0':>10} {'e1':>9}",flush=True)
for LSQ in [mp.mpf('2.2'),mp.mpf('2.5'),mp.mpf('2.72')]:
    ef=even_spectrum(LSQ,use_primes=True); ea=even_spectrum(LSQ,use_primes=False)
    print(f"{float(LSQ):6.2f} | {ef[0]:10.2e} {ef[1]:10.4f} {ef[2]:9.4f} | {ea[0]:10.4f} {ea[1]:9.4f}",flush=True)
print("\nDoes the prime OPEN the even gap e1 (full e1 >> arch e1)? Is e0~0 (=k_lambda null mode)?",flush=True)
