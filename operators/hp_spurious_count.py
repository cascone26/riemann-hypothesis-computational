import sys, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
lam2=float(sys.argv[1])
def polymul(a,b):
    r=[mp.mpf(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): r[i+j]+=x*y
    return r
for spec in sys.argv[2:]:
    N,dps=map(int,spec.split(':'))
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
    t=sorted([mp.re(x) for x in mp.polyroots(P[::-1],maxsteps=1500,extraprec=8*dps) if mp.re(x)>0])
    Z=[];k=1
    while True:
        g=mp.zetazero(k).imag; Z.append(g); k+=1
        if g>t[-1]+10: break
    def near(x):
        j=min(range(len(Z)),key=lambda i:abs(Z[i]-x)); gap=(Z[min(j+1,len(Z)-1)]-Z[max(j-1,0)])/2
        return abs(Z[j]-x)<0.25*gap
    res=[x for x in t if near(x)]; sp=[x for x in t if not near(x)]
    print(f"N={N}: {len(t)} positive nodes, resolved(near zeta)={len(res)}, spurious={len(sp)}: {[mp.nstr(x,5) for x in sp]}",flush=True)
