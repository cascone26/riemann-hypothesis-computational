# SCOPING: uncertainty-principle route to rigorous lambda_a>0.
# 2pi Q(v) = int S_a(xi)|vhat|^2 dxi, supp v in (-a,a):
#   S_a(xi) = S_arch(xi) + S_pole_a(xi) + S_prime(xi)
#   S_arch: from arch-core+tail. Core |t|-part: -2*2.686*2 = ... derive from kernel FT on |t|<=2a (use truncated FTs
#           uniformly: all kernels restricted to |t|<=2a, FT explicit/numeric -> legal multiplier).
#   Approach here (scoping): numerically FT each truncated kernel piece; find bad set B={S_a<0}; Donoho-Stark:
#   mass_B(vhat) <= 2a|B|/(2pi) * ||vhat||^2_"peak"  — use simplest DS: int_B |vhat|^2 <= |B| * ||v||_1^2 <= |B|*2a*||v||_2^2.
#   Then 2pi*lambda_a >= inf_{B^c} S_a * (2pi - is bounded...) — do the exact two-set bound:
#   2pi Q >= (inf_B S_a) * m_B + (inf_{B^c} S_a) * (2pi||v||^2 - m_B),  m_B <= 2a|B| ||v||^2.
#   => lambda_a >= [ inf_{B^c}S_a * 2pi + (inf_B S_a - inf_{B^c}S_a) * 2a|B| ] / (2pi)   [inf_B S_a<0]
import numpy as np, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),'..'))
import suzuki_screw_lambda_a as S
def S_a(a, xis):
    # numeric cosine-FT of each kernel piece truncated to t in [0,2a]  (S(xi)=2*int_0^{2a} k(t)cos(xi t)dt of -g''… 
    # easier: direct multiplier = FT of kernel * xi^2 identity used in §81; here compute S = -(g'')^ via g'' pieces:
    # g''_atoms: arch-core'' = 2*2.686*delta_0?? NO — use verified route: S(xi) = 2*[  sum over pieces  ]:
    #   arch-core (2.686|t|): contributes CONSTANT -2*2.686  (verified §79: -5.372)
    #   prime n: -(Lam(n)/sqrt n)*2 cos(xi log n)  for log n <= 2a
    #   pole: g_pole''(t) = -(e^{t/2}+e^{-t/2})  => contribution +2*int_0^{2a}(e^{t/2}+e^{-t/2})cos(xi t)dt  [truncated]
    #   arch-tail'': compute numerically from S.phi_tail second derivative — approximate via finite diff of the tail kernel.
    out=np.zeros_like(xis)
    out += -2*2.6861  # arch core (=-|Omega(0)|=-5.3722 total)
    # primes with log n <= 2a
    Nmax=int(np.floor(np.exp(2*a)))
    Lam=S.vonmangoldt_upto(max(Nmax,3))
    for n in range(2,Nmax+1):
        if Lam[n]>0:
            out += -(Lam[n]/np.sqrt(n))*2*np.cos(xis*np.log(n))
    # pole (truncated integral, exact closed form via numeric quad)
    ts=np.linspace(0,2*a,4000); dt=ts[1]-ts[0]; w=np.full_like(ts,dt); w[0]=w[-1]=dt/2
    ker=np.exp(ts/2)+np.exp(-ts/2)
    for i,xi in enumerate(xis):
        out[i] += 2*np.sum(w*ker*np.cos(xi*ts))
    # arch tail: kernel k(t) = -0.25(PHI - e^{-t/2}phi_tail); second derivative numerically, then +2 int (-k'')cos... 
    # equivalently S_tail(xi) = -2*int k''(t)cos = use FT of k with xi^2: S_tail = 2*xi^2*int_0^{2a} k(t)cos(xi t)dt - boundary
    # scoping: include via xi^2-FT with boundary terms dropped (small, flag as approximation)
    kt=-0.25*(S.PHI_1_2_QUARTER-np.exp(-ts/2)*np.array([S.phi_tail(t) for t in ts]))
    for i,xi in enumerate(xis):
        out[i] += 2*(xi**2)*np.sum(w*kt*np.cos(xi*ts))
    return out
for a in [0.30,0.37,0.42]:
    xis=np.linspace(0,300,6000)
    Sv=S_a(a,xis)
    bad=Sv<0
    if bad.any():
        Bmeas=2*(xis[1]-xis[0])*bad.sum()   # symmetric ±xi
        infB=Sv[bad].min(); infG=Sv[~bad].min()
        lower=(infG*2*np.pi + (infB-infG)*2*a*Bmeas)/(2*np.pi)
        print(f"a={a:.2f}: bad set |B|={Bmeas:.2f} (xi in ±{xis[bad][-1] if bad.any() else 0:.1f}), inf_B S={infB:.2f}, inf_good S={infG:.3f}  =>  DS lower bound lambda_a >= {lower:+.3f}")
    else:
        print(f"a={a:.2f}: S_a(xi) >= {Sv.min():.3f} everywhere -> lambda_a >= {Sv.min()/(2*np.pi):.3f} trivially")
print("\nbound >0 at a=0.42 => the uncertainty route CAN close the first inductive step rigorously. <0 => need Nazarov/sharper.")
