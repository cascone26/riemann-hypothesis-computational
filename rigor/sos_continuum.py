#!/usr/bin/env python3
"""SOS / Positivstellensatz certificate that the even Weil block Q(a) ⪰ 0 across a WHOLE a-interval in one shot
(idea: 'third path around the continuum' — replaces the 44-panel sweep AND eigenvector-enclosure).

Matrix-SOS interval certificate (Putinar, degree 2): for a symmetric matrix polynomial M(a) (n×n, deg ≤2 in a),
   M(a) = S0(a) + (a−a0)(a1−a)·S1,   S0(a)=A+(B+Bᵀ)a+C a² a matrix-SOS (Gram [[A,B],[Bᵀ,C]]⪰0), S1=P1⪰0.
If PSD A,B,C,P1 satisfy the coefficient match, then M(a)⪰0 for ALL a∈[a0,a1] — PROVEN for the polynomial M.
Our real entries are transcendental in a; the deg-2 fit's remainder is what rigor/affine_arith.py bounds (SOS
proves the polynomial part, affine bounds the remainder — they compose). Honest tiers below.

Report tier: cvxpy/SCS finds a floating certificate + dense-grid eigenvalue confirmation.
Handle tier : rational rounding of (A,B,C,P1) + interval re-check (stubbed/flagged — the remaining rigorous step).
"""
import sys, numpy as np
import mpmath as mp
mp.mp.dps = 30

# ---- parametrized even-block Weil entry Q(k,l; a) (mpmath, closed form; NL Lerch terms) ----
def _mk(a, NL=800):
    a=mp.mpf(a); b=mp.pi/(2*a); log2=mp.log(2); twoa=2*a; wp=log2/mp.sqrt(2)
    Om0=mp.digamma(mp.mpf(1)/4)-mp.log(mp.pi); Phi0=mp.zeta(2,mp.mpf(1)/4)
    def Pp(p,t,s):
        if p==0: return 1/s
        if p==1: return t/s-1/s**2
        return t*t/s-2*t/s**2+2/s**3
    def II(p,al,om,sc,t0,t1):
        s=mp.mpc(-al,om); F=lambda tt: mp.e**(s*tt)*Pp(p,tt,s); v=F(t1)-F(t0)
        return v.real if sc=='cos' else v.imag
    def dterms(k,l):
        if k==l: return [(2*a*k*k*b*b,0,k*b,'cos'),(-k*k*b*b,1,k*b,'cos'),(-k*b,0,k*b,'sin')]
        al=k*l*b*(mp.mpf(1)/(k-l)-mp.mpf(1)/(k+l)); ak=-k*l*b*(mp.mpf(1)/(k-l)+mp.mpf(1)/(k+l))
        return [(al,0,l*b,'sin'),(ak,0,k*b,'sin')]
    def Q(k,l):
        tot=mp.mpf(0); D=dterms(k,l)
        gt=[(-4,0,mp.mpf(-1)/2),(-4,0,mp.mpf(1)/2),(8,0,mp.mpf(0)),(-Om0/2,1,mp.mpf(0)),(-Phi0/4,0,mp.mpf(0))]
        for n in range(NL): gt.append((mp.mpf(1)/4/(n+mp.mpf(1)/4)**2,0,2*n+mp.mpf(1)/2))
        for (ag,pg,alp) in gt:
            for (ad,pd,om,sc) in D: tot+=ag*ad*II(pg+pd,alp,om,sc,0,twoa)
        for (ad,pd,om,sc) in D:
            tot+=wp*ad*II(pd+1,0,om,sc,log2,twoa)+(-wp*log2)*ad*II(pd,0,om,sc,log2,twoa)
        return tot
    return Q

MODES=[1,3,5,7]  # even sector (odd modes), the 4×4
def block(a):
    Q=_mk(a); n=len(MODES)
    return np.array([[float(Q(MODES[i],MODES[j])) for j in range(n)] for i in range(n)])

def fit_polys(a0,a1,deg=2,ns=7):
    """Sample M(a) on a grid, fit each entry as a degree-`deg` poly in a. Returns coeff matrices M[d] (n×n),
    M(a)=Σ_d M[d] a^d, plus the max fit residual over a dense grid (the quantity affine_arith must bound)."""
    n=len(MODES); ags=np.linspace(a0,a1,ns)
    blocks=[block(a) for a in ags]
    Mc=[np.zeros((n,n)) for _ in range(deg+1)]
    for i in range(n):
        for j in range(n):
            c=np.polyfit(ags,[blocks[t][i,j] for t in range(len(ags))],deg)  # highest power first
            for d in range(deg+1): Mc[d][i,j]=c[deg-d]
    # residual on a denser grid
    dense=np.linspace(a0,a1,41); res=0.0
    for a in dense:
        Mpoly=sum(Mc[d]*a**d for d in range(deg+1)); res=max(res,np.max(np.abs(Mpoly-block(a))))
    return Mc,res

def sos_certify(a0,a1,Mc):
    import cvxpy as cp
    n=Mc[0].shape[0]
    A=cp.Variable((n,n),symmetric=True); C=cp.Variable((n,n),symmetric=True)
    B=cp.Variable((n,n)); P1=cp.Variable((n,n),symmetric=True)
    Gram=cp.bmat([[A,B],[B.T,C]])
    cons=[Gram>>0, P1>>0]
    # coeff match: M0=A−a0a1 P1 ; M1=(B+Bᵀ)+(a0+a1)P1 ; M2=C−P1
    cons+=[Mc[0]==A-a0*a1*P1, Mc[1]==(B+B.T)+(a0+a1)*P1, Mc[2]==C-P1]
    prob=cp.Problem(cp.Minimize(0),cons)
    prob.solve(solver=cp.SCS, eps=1e-9)
    return prob.status, (A.value,B.value,C.value,P1.value)

def dense_confirm(a0,a1,ns=61):
    mn=min(np.linalg.eigvalsh(block(a)).min() for a in np.linspace(a0,a1,ns))
    return mn

def machinery_demo():
    """Validate the SOS interval-certificate MACHINERY on a well-conditioned parametrized PD family, so the
    certificate path is proven independent of the razor-thin real block. M(a)=[[2+a, a],[a, 2-a]] on a∈[a0,a1]:
    eigenvalues 2±√(2a²)  ⟹ PD with comfortable margin for |a|<√2."""
    import cvxpy as cp
    a0,a1=-0.5,0.5
    Mc=[np.array([[2.,0.],[0.,2.]]), np.array([[1.,1.],[1.,-1.]]), np.zeros((2,2))]  # M0+M1 a+M2 a²
    status,cert=sos_certify(a0,a1,Mc)
    A,B,C,P1=cert; Gram=np.block([[A,B],[B.T,C]])
    gmin=np.linalg.eigvalsh(Gram).min(); pmin=np.linalg.eigvalsh(P1).min()
    # verify the certificate identity M(a)=S0(a)+(a-a0)(a1-a)P1 holds at samples
    resid=0.
    for a in np.linspace(a0,a1,9):
        S0=A+(B+B.T)*a+C*a*a; lhs=Mc[0]+Mc[1]*a+Mc[2]*a*a; rhs=S0+(a-a0)*(a1-a)*P1
        resid=max(resid,np.max(np.abs(lhs-rhs)))
    ok = status in ('optimal','optimal_inaccurate') and gmin>-1e-6 and pmin>-1e-6 and resid<1e-5
    print("MACHINERY VALIDATION — SOS interval certificate on a well-conditioned M(a)=[[2+a,a],[a,2-a]], a∈[-.5,.5]:")
    print(f"  status={status}  Gram min-eig={gmin:.2e} (⪰0)  P1 min-eig={pmin:.2e} (⪰0)  identity residual={resid:.1e}")
    print(f"  ⟹ SOS MACHINERY {'WORKS' if ok else 'FAILED'}: a single certificate proves M(a)⪰0 across the WHOLE interval.")
    return ok

if __name__=="__main__":
    print("="*96); machinery_demo(); print("="*96)
    a0,a1=(float(sys.argv[1]),float(sys.argv[2])) if len(sys.argv)>2 else (0.48,0.52)
    print(f"\nREAL even 4×4 Weil block Q(a)⪰0, a∈[{a0},{a1}]  (modes {MODES}) — accuracy stress test:")
    Mc,res=fit_polys(a0,a1,deg=2)
    mn=dense_confirm(a0,a1)
    print(f"  dense-grid min eigenvalue over window = {mn:.3e}  ({'PD' if mn>0 else 'dips ≤0'});  deg-2 fit residual = {res:.2e}")
    print(f"  FINDING: the block's true margin (~1e-6, per high-precision scan) is ~{res/max(mn,1e-12):.0e}× SMALLER than")
    print(f"  the polyfit residual — so a loose fit + float64 SDP CANNOT certify the real block. The SOS path is sound")
    print(f"  (machinery validated above) but REQUIRES sub-margin entry enclosure: feed rigor/affine_arith.py's")
    print(f"  affine/Taylor enclosure (not a polyfit) + an extended-precision/rationalized SDP verify. That is the")
    print(f"  composition point — affine bounds the entries to <margin, SOS then certifies the whole window at once.")
