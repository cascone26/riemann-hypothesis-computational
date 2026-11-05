"""S2868: finite Gauss-quadrature test. t_j = zeros of R_xi(z)=sum_n xi_n/(z-nu_n) (ground vector xi);
c_t(n)=1/(nu_n+t) (note sign convention: zero side used 1/(nu+gamma) at gamma=-t ... we test both).
Solve tau - eps I = C A C^T for A (2N x 2N); report off-diagonal size of A and signs of diag."""
import sys, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
lam2=float(sys.argv[1]); dps=int(sys.argv[2]); N=int(sys.argv[3])
W02,WR,WP,L,_=parts(lam2,N,dps=dps,as_mp=True); S=W02-WR-WP; M=S.rows
nu=[2*mp.pi*n/L for n in range(-N,N+1)]
E,V=mp.eigsy(S); i0=min(range(M),key=lambda i:E[i]); eps=E[i0]; xi=[V[r,i0] for r in range(M)]
# zeros of R(z)=sum xi_n/(z-nu_n): numerator polynomial P(z)=sum_n xi_n prod_{m!=n}(z-nu_m), degree<=2N
import itertools
P=[mp.mpf(0)]
def polymul(a,b):
    r=[mp.mpf(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): r[i+j]+=x*y
    return r
P=[mp.mpf(0)]*(M)
for n in range(M):
    q=[mp.mpf(1)]
    for m in range(M):
        if m!=n: q=polymul(q,[-nu[m],mp.mpf(1)])
    for i in range(len(q)): P[i]+=xi[n]*q[i]
while abs(P[-1])<mp.mpf(10)**(-dps//2): P.pop()
roots=mp.polyroots(P[::-1],maxsteps=400,extraprec=4*dps)
print(f"N={N}: eps={mp.nstr(eps,4)} deg P={len(P)-1}  max|Im root|={mp.nstr(max(abs(mp.im(r)) for r in roots),3)}")
t=sorted([mp.re(r) for r in roots])
print("  real zeros of xihat (t>0):",[mp.nstr(x,6) for x in t if x>0])
for sgn in [1,-1]:
    C=mp.matrix(M,len(t))
    for i in range(M):
        for j in range(len(t)): C[i,j]=1/(nu[i]+sgn*t[j])
    T=S-eps*mp.eye(M)
    # A = C^+ T C^+T
    Cp=mp.inverse(C.T*C)*C.T
    A=Cp*T*Cp.T
    resid=mp.mnorm(C*A*C.T-T,1)
    off=max(abs(A[i,j]) for i in range(A.rows) for j in range(A.cols) if i!=j); dg=[A[i,i] for i in range(A.rows)]
    w=[4*mp.sin(x*L/2)**2/L for x in t]
    print(f"  c_t=1/(nu{'+' if sgn>0 else '-'}t): resid={mp.nstr(resid,3)}  max|offdiag A|={mp.nstr(off,3)}  max|diag|={mp.nstr(max(abs(d) for d in dg),3)}  #diag>0={sum(1 for d in dg if d>0)}/{len(dg)}")
    print("     diag(A)/w(t):",[mp.nstr(dg[j]/w[j],5) for j in range(len(t))][:12])
