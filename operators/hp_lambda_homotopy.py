"""lambda-homotopy (S2866): as L varies, node nu_k=2 pi k/L sweeps across zeta zero gamma_j at L*=2 pi k/gamma_j.
Bisect L where xihat_N(nu_k) (∝ (-1)^k xi_k of ground vector) changes sign; compare with L*."""
import sys, mpmath as mp
sys.path.insert(0,'operators'); from ccm_parts import parts
N=int(sys.argv[1]); k=int(sys.argv[2]); j=int(sys.argv[3]); dps=int(sys.argv[4])
mp.mp.dps=dps; g=mp.zetazero(j).imag; Ls=2*mp.pi*k/g
def f(L):
    W02,WR,WP,_,_=parts(mp.e**L,N,dps=dps,as_mp=True); S=W02-WR-WP; M=S.rows
    E,V=mp.eigsy(S); i0=min(range(M),key=lambda i:E[i]); xi=[V[r,i0] for r in range(M)]
    s=1 if xi[N]>0 else -1  # fix overall sign by central entry
    return s*(-1)**k*xi[N+k]
a,b=Ls-mp.mpf('0.02'),Ls+mp.mpf('0.02'); fa,fb=f(a),f(b)
print(f"N={N} node k={k} zero #{j} gamma={mp.nstr(g,12)} L*={mp.nstr(Ls,12)} (lam2*={mp.nstr(mp.e**Ls,8)})  f(L*-.02)={mp.nstr(fa,4)} f(L*+.02)={mp.nstr(fb,4)}",flush=True)
if fa*fb<0:
    for it in range(30):
        c=(a+b)/2; fc=f(c)
        if fa*fc<=0: b,fb=c,fc
        else: a,fa=c,fc
    print(f"  sign flip at L={mp.nstr((a+b)/2,12)}   L*={mp.nstr(Ls,12)}   diff={mp.nstr((a+b)/2-Ls,4)}  (node offset nu_k-gamma at flip = {mp.nstr(2*mp.pi*k/((a+b)/2)-g,4)})",flush=True)
else: print("  no sign flip in window")
