"""Sign structure of S=W02-WR-WP (=CCM tau = Shi S_NL) and its ground state, for Perron-Frobenius simplicity."""
import sys, numpy as np, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
lam2=float(sys.argv[1]) if len(sys.argv)>1 else 14.0
for N in [int(a) for a in sys.argv[2:]] or [3,6,10]:
    W02,WR,WP,L,_=parts(lam2,N,dps=30)
    S=W02-WR-WP; M=len(S); idx=np.arange(-N,N+1)
    off=S[~np.eye(M,dtype=bool)]
    Dsg=np.diag((-1.0)**idx); C=Dsg@S@Dsg; offC=C[~np.eye(M,dtype=bool)]
    mpS=mp.matrix(S.tolist()); E,V=mp.eigsy(mpS)
    ev=sorted([(E[i],i) for i in range(M)],key=lambda t:t[0])
    g=[float(V[r,ev[0][1]]) for r in range(M)]
    g=np.array(g); g*=np.sign(g[N])
    print(f"N={N} L={L:.4f}: off>0 frac={np.mean(off>0):.3f}  checkerboard off<0 frac={np.mean(offC<0):.3f}")
    print("   eig0..3:",[mp.nstr(ev[k][0],6) for k in range(4)], " gap ratio e1/e0=",mp.nstr(ev[1][0]/ev[0][0],6))
    print("   ground vec signs:", ''.join('+' if x>0 else '-' for x in g), " sum=",round(g.sum(),4), " even?",np.allclose(g,g[::-1],atol=1e-8))
