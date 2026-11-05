"""S2882b: interlacing in the PLANTED off-line world (zeros #a,#b -> gm ± i*f*h). Structural vs arithmetic test."""
import sys, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
lam2=float(sys.argv[1]); dps=int(sys.argv[2]); a,b=map(int,sys.argv[3].split(',')); f=mp.mpf(sys.argv[4]); Ns=[int(x) for x in sys.argv[5:]]
def polymul(p,q):
    r=[mp.mpf(0)]*(len(p)+len(q)-1)
    for i,x in enumerate(p):
        for j,y in enumerate(q): r[i+j]+=x*y
    return r
def nodes(N):
    W02,WR,WP,L,_=parts(lam2,N,dps=dps,as_mp=True); S=W02-WR-WP; M=S.rows; mp.mp.dps=dps
    nu=[2*mp.pi*n/L for n in range(-N,N+1)]
    def term(z):
        w=4*mp.sin(z*L/2)**2/L; return mp.matrix([[w/((nu[i]+z)*(nu[j]+z)) for j in range(M)] for i in range(M)])
    ga=mp.zetazero(a).imag; gb=mp.zetazero(b).imag; gm=(ga+gb)/2; h=(gb-ga)/2; z=gm+1j*f*h
    T=term(z)+term(mp.conj(z))+term(-z)+term(-mp.conj(z))
    S=S-(term(ga)+term(-ga)+term(gb)+term(-gb))+mp.matrix([[mp.re(T[i,j]) for j in range(M)] for i in range(M)])
    E,V=mp.eigsy(S); i0=min(range(M),key=lambda i:E[i]); xi=[V[r,i0] for r in range(M)]
    P=[mp.mpf(0)]*M
    for n in range(M):
        q=[mp.mpf(1)]
        for m in range(M):
            if m!=n: q=polymul(q,[-nu[m],mp.mpf(1)])
        for i in range(len(q)): P[i]+=xi[n]*q[i]
    while abs(P[-1])<mp.mpf(10)**(-dps//2): P.pop()
    r=mp.polyroots(P[::-1],maxsteps=800,extraprec=6*dps)
    return sorted([mp.re(x) for x in r if mp.re(x)>0]), max(abs(mp.im(x)) for x in r), E[i0]
prev=None
for N in Ns:
    t,mi,e=nodes(N)
    print(f"N={N}: eps={mp.nstr(e,3)} max|Im node|={mp.nstr(mi,2)} positive nodes: {[mp.nstr(x,6) for x in t]}",flush=True)
    if prev is not None and len(t)==len(prev)+1:
        print(f"   interlace={all(t[k]<prev[k]<t[k+1] for k in range(len(prev)))}  decreased={all(t[k]<prev[k] for k in range(len(prev)))}",flush=True)
    elif prev is not None: print(f"   (node count changed {len(prev)}->{len(t)})",flush=True)
    prev=t
