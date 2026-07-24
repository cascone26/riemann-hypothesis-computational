# PIVOT viability probe: does the # of near-null eigenvalues of tau match the # of zeta zeros within
# spectral reach? If YES -> eigenvalue-side realization (trace-formula) is live. If NO -> that door
# is a different structure than a naive zero-counter.
import mpmath as mp
mp.mp.dps=26
GAM=[14.134725,21.022040,25.010858,30.424876,32.935062,37.586178,40.918719,43.327073,48.005151,49.773832]
def probe(LSQ,N=12,thresh=1e-5):
    L=2*mp.log(mp.sqrt(LSQ)); TWO_PI=2*mp.pi; LOGPI=mp.log(mp.pi)
    idx=list(range(-N,N+1)); D=len(idx)
    reach=float(TWO_PI*N/L)  # spectral reach in r (=gamma) units
    Rq=TWO_PI*N/L+55; G=800
    grid=[(-Rq+2*Rq*mp.mpf(i)/G) for i in range(G+1)]; hh=2*Rq/G
    wq=[hh]*(G+1); wq[0]=hh/2; wq[-1]=hh/2
    Om=[mp.re(mp.digamma(mp.mpf('0.25')+0.5j*r))-LOGPI for r in grid]
    def ehat(r,n):
        s=TWO_PI*n/L; d=s-r
        if abs(d)<mp.mpf('1e-20'): return mp.sqrt(L)
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
    ev=sorted(float(e) for e in mp.eigsy(T,eigvals_only=True))
    nnull=sum(1 for e in ev if abs(e)<thresh)
    nzero=sum(1 for g in GAM if g<reach)
    return reach,nnull,nzero,ev[:6]
for LSQ in [mp.mpf('14'),mp.mpf('30'),mp.mpf('60')]:
    reach,nn,nz,lo=probe(LSQ)
    print(f"lam^2={float(LSQ):.0f}: reach={reach:.1f} | #near-null(|e|<1e-5)={nn}  #zeta-zeros<reach={nz} | lowest eig: {['%.1e'%e for e in lo]}",flush=True)
print("\nmatch (nn==nz) => eigenvalue-side realization live; mismatch => cascade != naive zero-count.",flush=True)
