"""Gantmacher-Krein test: is the even (or odd) sector block of tau = W02-WR-WP (CCM = L*Shi S_NL) totally
positive, possibly after checkerboard signing?  Gasca-Pena: A is TP iff all INITIAL minors > 0
(contiguous minors using row 1 or column 1 block start)."""
import sys, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
lam2=float(sys.argv[1]); dps=int(sys.argv[2])
def sectors(S,N):
    M=S.rows; c=N  # index of n=0
    r2=1/mp.sqrt(2)
    # even basis: e0, (e_k+e_-k)/sqrt2 ; odd basis: (e_k-e_-k)/sqrt2, k=1..N
    def ev(k):
        v=[mp.mpf(0)]*M
        if k==0: v[c]=1
        else: v[c+k]=r2; v[c-k]=r2
        return v
    def od(k):
        v=[mp.mpf(0)]*M; v[c+k]=r2; v[c-k]=-r2; return v
    def block(vs):
        n=len(vs); B=mp.matrix(n,n)
        for i in range(n):
            Si=[mp.fsum(S[a,b]*vs[i][b] for b in range(M)) for a in range(M)]
            for j in range(n): B[j,i]=mp.fsum(vs[j][a]*Si[a] for a in range(M))
        return B
    return block([ev(k) for k in range(N+1)]), block([od(k) for k in range(1,N+1)])
def initial_minors(A):
    n=A.rows; out=[]
    for size in range(1,n+1):
        for start in range(0,n-size+1):
            # rows start..start+size-1, cols 0..size-1  and rows 0..size-1, cols start..
            for (rs,cs) in [(start,0),(0,start)]:
                sub=mp.matrix(size,size)
                for i in range(size):
                    for j in range(size): sub[i,j]=A[rs+i,cs+j]
                out.append((size,rs,cs,mp.det(sub)))
    return out
if __name__=="__main__":
 for N in [int(a) for a in sys.argv[3:]]:
     W02,WR,WP,L,_=parts(lam2,N,dps=dps,as_mp=True); S=W02-WR-WP
     Ev,Od=sectors(S,N)
     for name,B in [('even',Ev),('odd',Od)]:
         for sgn in ['plain','checker']:
             A=B.copy()
             if sgn=='checker':
                 for i in range(A.rows):
                     for j in range(A.cols): A[i,j]*=(-1)**(i+j)
             mins=initial_minors(A); neg=[m for m in mins if m[3]<=0]
             print(f"N={N} {name:4s} {sgn:7s}: initial minors={len(mins)}  nonpositive={len(neg)}", ("first bad: size=%d rs=%d cs=%d val=%s"%(neg[0][0],neg[0][1],neg[0][2],mp.nstr(neg[0][3],4))) if neg else "  ==> TP")
