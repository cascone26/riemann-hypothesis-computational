"""At L* = 2 pi k / gamma_j exactly (zeta zero ON node k): is the ground-vector component xi_k exactly 0?"""
import sys, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
N=int(sys.argv[1]); k=int(sys.argv[2]); j=int(sys.argv[3])
for dps in [int(a) for a in sys.argv[4:]]:
    mp.mp.dps=dps+20; g=mp.zetazero(j).imag; Ls=2*mp.pi*k/g
    W02,WR,WP,_,_=parts(mp.e**Ls,N,dps=dps,as_mp=True); S=W02-WR-WP; M=S.rows
    E,V=mp.eigsy(S); o=sorted(range(M),key=lambda i:E[i]); xi=[V[r,o[0]] for r in range(M)]
    print(f"dps={dps}: e0={mp.nstr(E[o[0]],5)}  |xi_k|={mp.nstr(abs(xi[N+k]),5)}  |xi_-k|={mp.nstr(abs(xi[N-k]),5)}  neighbors |xi_k-1|={mp.nstr(abs(xi[N+k-1]),4)} |xi_k+1|={mp.nstr(abs(xi[N+k+1]),4)}",flush=True)
