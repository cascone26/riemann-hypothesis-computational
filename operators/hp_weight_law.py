"""S2869: weight law. Quadrature nodes t_j (zeros of R_xi) vs nearest zeta zero, and alpha_j/w(t_j)."""
import sys, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
lam2=float(sys.argv[1]); dps=int(sys.argv[2]); N=int(sys.argv[3])
W02,WR,WP,L,_=parts(lam2,N,dps=dps,as_mp=True); S=W02-WR-WP; M=S.rows
nu=[2*mp.pi*n/L for n in range(-N,N+1)]
E,V=mp.eigsy(S); i0=min(range(M),key=lambda i:E[i]); eps=E[i0]; xi=[V[r,i0] for r in range(M)]
def polymul(a,b):
    r=[mp.mpf(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): r[i+j]+=x*y
    return r
P=[mp.mpf(0)]*M
for n in range(M):
    q=[mp.mpf(1)]
    for m in range(M):
        if m!=n: q=polymul(q,[-nu[m],mp.mpf(1)])
    for i in range(len(q)): P[i]+=xi[n]*q[i]
while abs(P[-1])<mp.mpf(10)**(-dps//2): P.pop()
roots=mp.polyroots(P[::-1],maxsteps=800,extraprec=6*dps)
t=sorted([mp.re(r) for r in roots]); mi=max(abs(mp.im(r)) for r in roots)
C=mp.matrix(M,len(t))
for i in range(M):
    for j in range(len(t)): C[i,j]=1/(nu[i]-t[j])
T=S-eps*mp.eye(M); Cp=mp.inverse(C.T*C)*C.T; A=Cp*T*Cp.T
off=max(abs(A[i,j]) for i in range(A.rows) for j in range(A.cols) if i!=j)
zz=[]; k=1
while True:
    g=mp.zetazero(k).imag
    if g>nu[-1]*1.6: break
    zz.append(g); k+=1
print(f"lam2={lam2} N={N} nu_N={mp.nstr(nu[-1],5)} eps={mp.nstr(eps,4)} max|Im t|={mp.nstr(mi,2)} max|offdiag A|={mp.nstr(off,2)} #alpha>0={sum(1 for j in range(len(t)) if A[j,j]>0)}/{len(t)}")
for j,x in enumerate(t):
    if x<=0: continue
    g=min(zz,key=lambda y:abs(y-x)) if zz else None
    w=4*mp.sin(x*L/2)**2/L
    print(f"  t={mp.nstr(x,12):>16}  nearest gamma={mp.nstr(g,12):>16}  |t-gamma|={mp.nstr(abs(x-g),3):>9}  alpha/w={mp.nstr(A[j,j]/w,8)}")
