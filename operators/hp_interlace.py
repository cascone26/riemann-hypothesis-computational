"""S2882: do CCM node sets interlace between N and N+1?  t^{(N+1)}_k < t^{(N)}_k < t^{(N+1)}_{k+1}."""
import sys, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
lam2=float(sys.argv[1]); dps=int(sys.argv[2]); Ns=[int(a) for a in sys.argv[3:]]
def polymul(a,b):
    r=[mp.mpf(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): r[i+j]+=x*y
    return r
def nodes(N):
    W02,WR,WP,L,_=parts(lam2,N,dps=dps,as_mp=True); S=W02-WR-WP; M=S.rows; mp.mp.dps=dps
    nu=[2*mp.pi*n/L for n in range(-N,N+1)]
    E,V=mp.eigsy(S); i0=min(range(M),key=lambda i:E[i]); xi=[V[r,i0] for r in range(M)]
    P=[mp.mpf(0)]*M
    for n in range(M):
        q=[mp.mpf(1)]
        for m in range(M):
            if m!=n: q=polymul(q,[-nu[m],mp.mpf(1)])
        for i in range(len(q)): P[i]+=xi[n]*q[i]
    while abs(P[-1])<mp.mpf(10)**(-dps//2): P.pop()
    r=mp.polyroots(P[::-1],maxsteps=800,extraprec=6*dps)
    return sorted([mp.re(x) for x in r if mp.re(x)>0]), max(abs(mp.im(x)) for x in r)
prev=None
for N in Ns:
    t,mi=nodes(N)
    print(f"N={N}: {len(t)} positive nodes (max|Im|={mp.nstr(mi,2)}): {[mp.nstr(x,8) for x in t]}",flush=True)
    if prev is not None:
        a,b=prev,t  # a: N-1 nodes, b: N nodes
        ok=all(b[k]<a[k]<b[k+1] for k in range(len(a))) if len(b)==len(a)+1 else None
        dec=all(b[k]<a[k] for k in range(len(a)))
        print(f"   interlace(N-1 inside N) = {ok}   all t_k decreased = {dec}",flush=True)
    prev=t
