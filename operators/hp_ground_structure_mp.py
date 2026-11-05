import sys, numpy as np, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
lam2=float(sys.argv[1]); dps=int(sys.argv[2])
for N in [int(a) for a in sys.argv[3:]]:
    W02,WR,WP,L,_=parts(lam2,N,dps=dps,as_mp=True)
    S=W02-WR-WP; M=S.rows
    E,V=mp.eigsy(S); order=sorted(range(M),key=lambda i:E[i])
    def vec(k):
        v=[V[r,order[k]] for r in range(M)]; s=1 if v[N]>=0 else -1; return [s*x for x in v]
    g=vec(0)
    even=max(abs(g[r]-g[M-1-r]) for r in range(M)); odd=max(abs(g[r]+g[M-1-r]) for r in range(M))
    print(f"N={N} dps={dps}: eig0..3 ",[mp.nstr(E[order[k]],8) for k in range(4)])
    print(f"   e1/e0={mp.nstr(E[order[1]]/E[order[0]],8)}  ground: even-dev={mp.nstr(even,3)} odd-dev={mp.nstr(odd,3)}  sum={mp.nstr(mp.fsum(g),8)}")
    print("   signs:",''.join('+' if x>0 else '-' for x in g))
    off=[S[i,j] for i in range(M) for j in range(M) if i!=j]
    print(f"   off-diag >0 frac={sum(1 for x in off if x>0)/len(off):.3f}")
