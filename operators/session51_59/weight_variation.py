# VARIATIONAL LEVER: treat the prime-2 weight w as FREE. Sweep w around the arithmetic value
# w2 = Lambda(2)/sqrt2 = log2/sqrt2. Compute min(tau_EVEN)(w). Is w2 the ARGMAX (optimal rescue)?
# Where does min cross 0 (positivity boundary)? Does the arithmetic weight sit at a special point?
import mpmath as mp
mp.mp.dps=30
def min_even_of_w(LSQ,wfac,N=12):
    # wfac = w / w2  (1.0 = true arithmetic weight)
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
    lg2=mp.log(2); w2=lg2/mp.sqrt(2); w=wfac*w2
    def qk(n,m,y):
        if n==m: return 2*(1-abs(y)/L)*mp.cos(TWO_PI*n*y/L)
        return (mp.sin(TWO_PI*m*y/L)-mp.sin(TWO_PI*n*y/L))/(mp.pi*(n-m))
    T=mp.matrix(D,D)
    for i,n in enumerate(idx):
        for j,m in enumerate(idx): T[i,j]=arch(i,j)+Bnm(n,m)-w*qk(n,m,lg2)
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
for LSQ in [mp.mpf('2.5'),mp.mpf('2.72')]:
    print(f"\nlam^2={float(LSQ):.2f}  (w2=log2/sqrt2={float(mp.log(2)/mp.sqrt(2)):.4f}); sweep w/w2:",flush=True)
    print(f"  {'w/w2':>6} {'min(tau_EVEN)':>14}",flush=True)
    for wf in [0.0,0.5,0.8,0.9,1.0,1.1,1.2,1.5,2.0]:
        m=min_even_of_w(LSQ,mp.mpf(str(wf)))
        star=" <== arithmetic" if wf==1.0 else ""
        print(f"  {wf:6.2f} {m:14.6e}{star}",flush=True)
