import sys, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
for lam2,N in [(36,24),(10,24)]:
    W02,WR,WP,L,pw=parts(lam2,N,dps=30,as_mp=True); mp.mp.dps=30
    M=2*N+1; nu=[2*mp.pi*n/L for n in range(-N,N+1)]; c=N
    def gen(X): return [mp.mpf(0) if j==c else mp.pi*(nu[c]-nu[j])*X[c,j] for j in range(M)]
    g02,gR,gP=gen(W02),gen(WR),gen(WP)
    V=nu[-1]; zs=[]; k=1
    while True:
        z=mp.im(mp.zetazero(k))
        if z>V+1: break
        zs.append(z); k+=1
    gz=[mp.pi*sum((4/L)*mp.sin(z*L/2)**2*(1/(nu[j]-z)+1/(nu[j]+z)) for z in zs) for j in range(M)]
    st=lambda g,j: -(g[j+1]-g[j])/mp.pi
    print(f"lam2={lam2} N={N}  step = s02 - sR - sP - szeros(removed)   [sR,sP enter with minus sign]")
    for j in [0,6,12,18,23,24,30,36,42,47]:
        a,b,cc,d=st(g02,j),st(gR,j),st(gP,j),st(gz,j)
        print(f"  j={j:2d} nu={float(nu[j]):8.3f}  s02={float(a):+.4e}  -sR={float(-b):+.4e}  -sP={float(-cc):+.4e}  -sZ={float(-d):+.4e}  total={float(a-b-cc-d):+.4e}")
