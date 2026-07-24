# CASCADE HYPOTHESIS TEST: do the low near-null eigenvectors of tau=QW_lambda ALL have
# spectral function v_hat(r)=sum c_n ehat(r,n) vanishing at the zeta zeros? (hi-prec, dps=35)
import mpmath as mp
mp.mp.dps=35
LSQ=14; L=2*mp.log(mp.sqrt(LSQ)); TWO_PI=2*mp.pi; LOGPI=mp.log(mp.pi)
N=10; idx=list(range(-N,N+1)); D=len(idx)
Rq=TWO_PI*N/L+50; G=1400
grid=[(-Rq+2*Rq*mp.mpf(i)/G) for i in range(G+1)]; hh=2*Rq/G
wq=[hh]*(G+1); wq[0]=hh/2; wq[-1]=hh/2
Om=[mp.re(mp.digamma(mp.mpf('0.25')+0.5j*r))-LOGPI for r in grid]
def ehat(r,n):
    s=TWO_PI*n/L; d=s-r
    if abs(d)<mp.mpf('1e-28'): return mp.sqrt(L)
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
VM=[(k,vm(k)) for k in range(2,LSQ+1) if vm(k)>0]
def qk(n,m,y):
    if n==m: return 2*(1-abs(y)/L)*mp.cos(TWO_PI*n*y/L)
    return (mp.sin(TWO_PI*m*y/L)-mp.sin(TWO_PI*n*y/L))/(mp.pi*(n-m))
def Pnm(n,m): return mp.fsum(lp*mp.mpf(k)**mp.mpf('-0.5')*qk(n,m,mp.log(k)) for k,lp in VM)
print("building tau hi-prec...",flush=True)
T=mp.matrix(D,D)
for i,n in enumerate(idx):
    for j in range(i,D):
        m=idx[j]; v=arch(i,j)+Bnm(n,m)-Pnm(n,m); T[i,j]=v; T[j,i]=v
Ev,Q=mp.eigsy(T)  # eigenvalues Ev, eigenvectors columns of Q
order=sorted(range(D),key=lambda k: Ev[k])
gam=[mp.mpf(g) for g in [14.134725,21.022040,25.010858,30.424876,32.935062]]
def vhat(cvec,r): return abs(mp.fsum(cvec[k]*ehat(r,idx[k]) for k in range(D)))
# median scale of vhat off zeros
def medscale(cvec): 
    return mp.fsum(vhat(cvec,mp.mpf(z)) for z in [7,11,17,27,35])/5
print(f"CASCADE TEST: |vhat_k(gamma_j)|/scale for the lowest 8 eigenvectors",flush=True)
print("rank  eigenvalue    " + "  ".join(f"g{i+1}" for i in range(5)),flush=True)
for rank in range(8):
    k=order[rank]; c=[Q[i,k] for i in range(D)]; sc=medscale(c)
    ratios=[float(vhat(c,g)/sc) for g in gam]
    print(f"{rank:3d}  {mp.nstr(Ev[k],4):>12}  " + "  ".join(f"{x:.3f}" for x in ratios),flush=True)
print("(all rows small at every gamma => cascade hypothesis CONFIRMED: whole near-null tower vanishes at zeros)",flush=True)
