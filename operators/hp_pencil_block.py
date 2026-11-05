"""S2884: structure of B_{N+1} = [[B_N, x],[y^T, beta]] in nested G-orthonormal parity bases (Shi unshifted pencil)."""
import sys, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
lam2=float(sys.argv[1]); dps=int(sys.argv[2]); N1=int(sys.argv[3])  # build at N1, inspect nesting N1-1 inside N1
W02,WR,WP,L,_=parts(lam2,N1,dps=dps,as_mp=True); S=W02-WR-WP; M=S.rows; mp.mp.dps=dps
c=N1; nu=[2*mp.pi*n/L for n in range(-N1,N1+1)]
r2=1/mp.sqrt(2)
def ev(k):  # even combo of +-k (k>=1) and e0
    v=[mp.mpf(0)]*M
    if k==0: v[c]=1
    else: v[c+k]=r2; v[c-k]=r2
    return v
def od(k):
    v=[mp.mpf(0)]*M; v[c+k]=r2; v[c-k]=-r2; return v
# zero-mean even contrast vectors ordered by frequency: e_k - (const)e0 so mean zero: mean of ev(k) = sqrt2 (k>=1), ev(0) mean 1
evc=[]
for k in range(1,N1+1):
    v=ev(k); v[c]-=mp.sqrt(2); evc.append(v)   # zero-mean, nested in k
odc=[od(k) for k in range(1,N1+1)]             # automatically zero-mean, nested
def ip(a,b): return mp.fsum(a[i]*mp.fsum(S[i,j]*b[j] for j in range(M)) for i in range(M))
def gs(vs):
    out=[]
    for v in vs:
        w=v[:]
        for u in out:
            p=ip(u,w); w=[w[i]-p*u[i] for i in range(M)]
        nrm=mp.sqrt(ip(w,w)); out.append([x/nrm for x in w])
    return out
E=gs(evc); O=gs(odc)
SD=lambda v: [nu[i]*v[i] for i in range(M)]
B=mp.matrix(N1,N1)
for i in range(N1):
    for j in range(N1): B[i,j]=ip(E[i],SD(O[j]))
nB=mp.mnorm(B[0:N1-1,0:N1-1],'F'); x=[B[i,N1-1] for i in range(N1-1)]; y=[B[N1-1,j] for j in range(N1-1)]
print(f"N={N1}: ||B_N||_F={mp.nstr(nB,4)}  ||x||(old even -> new odd)={mp.nstr(mp.norm(x),4)}  ||y||(new even -> old odd)={mp.nstr(mp.norm(y),4)}  beta={mp.nstr(B[N1-1,N1-1],4)}")
print("   |B| row-by-row max over cols>row (upper part):",[mp.nstr(max([abs(B[i,j]) for j in range(i+1,N1)] or [0]),2) for i in range(N1)])
print("   |B| row-by-row max over cols<row (lower part):",[mp.nstr(max([abs(B[i,j]) for j in range(i)] or [0]),2) for i in range(N1)])
print("   diag:",[mp.nstr(B[i,i],5) for i in range(N1)])
sv=sorted([mp.sqrt(abs(e)) for e in mp.eigsy(B.T*B)[0]]); print("   singular values (=theta):",[mp.nstr(s,10) for s in sv[:6]])
