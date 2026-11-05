"""S2883: Shi's UNSHIFTED contrast pencil K_un y = theta G_un y (G_un=C^T S C, K_un=C^T S D C) — exactly nested in N.
Compare its positive eigenvalues with CCM nodes / zeta zeros, and test interlacing across N."""
import sys, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
lam2=float(sys.argv[1]); dps=int(sys.argv[2]); Ns=[int(a) for a in sys.argv[3:]]
prev=None
for N in Ns:
    W02,WR,WP,L,_=parts(lam2,N,dps=dps,as_mp=True); S=W02-WR-WP; M=S.rows; mp.mp.dps=dps
    D=mp.diag([2*mp.pi*n/L for n in range(-N,N+1)])
    # Helmert-type orthonormal basis of zero-mean vectors (any basis gives same spectrum): columns e_i - e_{i+1}, then orthonormalize
    C=mp.matrix(M,M-1)
    for j in range(M-1): C[j,j]=1; C[j+1,j]=-1
    Q,R=mp.qr(C); C=Q[:,0:M-1] if Q.cols>=M-1 else Q
    G=C.T*S*C; K=C.T*S*D*C
    ev=mp.eig(mp.inverse(G)*K,left=False,right=False)
    pos=sorted([mp.re(e) for e in ev if mp.re(e)>0]); mi=max(abs(mp.im(e)) for e in ev)
    print(f"N={N}: minEig(G_un)={mp.nstr(min(mp.eigsy(G)[0]),3)} max|Im|={mp.nstr(mi,2)} positive theta: {[mp.nstr(x,12) for x in pos[:6]]} ... top {mp.nstr(pos[-1],5)}",flush=True)
    if prev is not None and len(pos)==len(prev)+1:
        print(f"   interlace={all(pos[k]<prev[k]<pos[k+1] for k in range(len(prev)))}",flush=True)
    prev=pos
G=[mp.zetazero(k).imag for k in range(1,5)]; print("zeta:",[mp.nstr(g,12) for g in G])
