#!/usr/bin/env python3
"""
CYCLE 295 — arbitrary-precision port of suzuki_exact_Q.py's UNCONDITIONAL closed-form Suzuki form,
to certify λ_a>0 past the float64 eigenvalue-solve floor (~1e-13 for exact-Q, ~1e-16 zero-side).

Same closed forms as suzuki_exact_Q.py (validated there vs the old quadrature pipeline):
  Q = Qpole + Qcore + Qtail + Qprime,  λ_a = min eig(Q)/a  (unconditional: real Λ(n) + arch, NO zeros/RH).
Only special function: Lerch Φ(e^{-2t},1,¼) = mp.lerchphi. Certification: mp.cholesky succeeds ⟺ Q PD ⟺ λ_a>0.
VALIDATION: must reproduce the float64 exact-Q value at a=0.5 (λ_a≈1.02e-6) before trusting a≥0.9.
"""
import sys, mpmath as mp
mp.mp.dps = 50
PI = mp.pi
C_CORE = -(mp.digamma(mp.mpf(1)/4) - mp.log(PI))/2      # ≈ +2.6861

def vonmangoldt_upto(N):
    Lam = [mp.mpf(0)]*(N+1)
    comp = [False]*(N+1)
    for p in range(2, N+1):
        if not comp[p]:
            for m in range(2*p, N+1, p): comp[m]=True
            lp = mp.log(p); pk=p
            while pk<=N: Lam[pk]=lp; pk*=p
    return Lam

def E_p(p, a, k):
    om = k*PI/(2*a)
    return om*(mp.e**(-p*a) - ((-1)**k)*mp.e**(p*a))/(p*p+om*om)

def S_entry(s, a, k, l):
    L = 2*a; ok=k*PI/L; ol=l*PI/L
    if k==l:
        return (L-s)*mp.cos(ok*s) + mp.sin(ok*s)/ok
    if (k+l)%2==1: return mp.mpf(0)
    return (mp.sin(ol*s)-mp.sin(ok*s))/(ok-ol) + (mp.sin(ok*s)+mp.sin(ol*s))/(ok+ol)

def lerch1q(t): return mp.lerchphi(mp.e**(-2*t), 1, mp.mpf(1)/4)
EPS = mp.mpf(10)**(-mp.mp.dps-10)
GAMMA = mp.euler

def ker_reg(t):
    # ker(t)=½e^{-t/2}Φ(e^{-2t},1,¼); ker_reg = ker - ker_sing, bounded (log-singularity removed).
    # ker_reg = ½e^{-t/2}[Φ(e^{-2t},1,¼) + log(1-e^{-2t})];  ker_reg(0)=-½ψ(¼) finite.
    return mp.mpf('0.5')*mp.e**(-t/2)*(lerch1q(t) + mp.log(1-mp.e**(-2*t)))

def ker_sing(t):
    return -mp.mpf('0.5')*mp.e**(-t/2)*mp.log(1-mp.e**(-2*t))   # = ½Σ_{k≥1} e^{-(2k+½)t}/k

def sing_cos_moment(omega, L, Ksum=None):
    """∫_0^L ker_sing(t) cos(ωt) dt  EXACT: slow part via digamma, fast part = finite e^{-pL} sum.
       = ¼ Re[(1/b)(ψ(1+b)+γ)] , b=¼-iω/2   MINUS  ½Σ_{k≥1}(1/k)e^{-p_k L}(p_k cosωL-ω sinωL)/(p_k²+ω²)."""
    b = mp.mpf(1)/4 - 1j*omega/2
    slow = mp.mpf('0.25')*mp.re((mp.digamma(1+b)+GAMMA)/b)
    fast = mp.mpc(0)
    cwL = mp.cos(omega*L); swL = mp.sin(omega*L)
    if Ksum is None: Ksum = int(mp.mp.dps*2.5/ (2*float(L)) ) + 20   # e^{-2kL} decay -> #terms for dps
    for k in range(1, Ksum+1):
        p = 2*k + mp.mpf('0.5')
        fast += (mp.mpf(1)/k)*mp.e**(-p*L)*(p*cwL - omega*swL)/(p*p+omega*omega)
    return slow - mp.mpf('0.5')*mp.re(fast)

def _moments_closed(omega, L, Ksum):
    """Return (Mc, Msin, Mtsin) = ∫_0^L ker(t){cos ωt, sin ωt, t sin ωt} dt, ker(t)=½Σ_{k≥0}e^{-(2k+½)t}/(k+¼).
    Each = slow(digamma/trigamma) + finite e^{-p_kL} correction (k=0..Ksum). No Lerch."""
    w=omega; b=mp.mpf(1)/4 - 1j*w/2; d=b-mp.mpf(1)/4    # d=-iω/2
    psi_b=mp.digamma(b); psi_q=mp.digamma(mp.mpf(1)/4); dpsi_b=mp.polygamma(1,b)
    Mc_slow  = mp.mpf('0.25')*mp.re((psi_b-psi_q)/d)
    Msin_slow= mp.mpf('0.25')*mp.im((psi_b-psi_q)/d)
    Mtsin_slow=(mp.mpf(1)/8)*mp.im((psi_b-psi_q)/(d*d) - dpsi_b/d)
    cwL=mp.cos(w*L); swL=mp.sin(w*L)
    Mc_f=mp.mpf(0); Ms_f=mp.mpf(0); Mts_f=mp.mpf(0)
    for k in range(0,Ksum+1):
        p=2*k+mp.mpf('0.5'); den=p*p+w*w; e=mp.e**(-p*L); c=mp.mpf(1)/(k+mp.mpf(1)/4)
        Mc_f  += -c*mp.mpf('0.5')*e*(p*cwL - w*swL)/den
        Ms_f  += -c*mp.mpf('0.5')*e*(w*cwL + p*swL)/den
        # J1 fast part (e^{-pL} terms of ∫ t e^{-pt} sinωt dt):
        j1f = e*( -L*(w*cwL+p*swL)/den + swL/den - 2*p*(w*cwL+p*swL)/(den*den) )
        Mts_f += c*mp.mpf('0.5')*j1f
    return Mc_slow+Mc_f, Msin_slow+Ms_f, Mtsin_slow+Mts_f

def tail_moments(a, K):
    """Mc[j]=∫ker cos(ω_j t); Md[k]=∫ker ω_k(L-t)sin(ω_k t) = ω_k(L·Msin[k]-Mtsin[k]). All closed-form."""
    L=2*a; om=[k*PI/L for k in range(1,K+1)]
    Ksum=int(mp.mp.dps*2.4/(2*float(L)))+25
    Mc=[]; Md=[]
    for j in range(K):
        mc,ms,mts=_moments_closed(om[j],L,Ksum)
        Mc.append(mc); Md.append(om[j]*(L*ms - mts))
    return Mc, Md, om

def build_Q(a, K, NG=40):
    a=mp.mpf(a)
    Ehp=[E_p(mp.mpf('0.5'),a,k) for k in range(1,K+1)]
    Ehm=[E_p(mp.mpf('-0.5'),a,k) for k in range(1,K+1)]
    L=2*a; Nmax=int(mp.floor(mp.e**L))
    Lam=vonmangoldt_upto(max(Nmax,2))
    primes=[(n, Lam[n]) for n in range(2,Nmax+1) if Lam[n]>0 and mp.log(n)<L]
    logn=[mp.log(n) for (n,ln) in primes]; wpr=[ln/mp.sqrt(n) for (n,ln) in primes]
    Mc,Md,om=tail_moments(a,K)
    Q=[[mp.mpf(0)]*K for _ in range(K)]
    for k in range(1,K+1):
        for l in range(k,K+1):
            qpole = Ehp[k-1]*Ehm[l-1]+Ehm[k-1]*Ehp[l-1]
            qcore = (-2*C_CORE*a) if k==l else mp.mpf(0)
            if (k+l)%2==1: qtail=mp.mpf(0)
            elif k==l:     qtail=Md[k-1]
            else:          qtail=2*om[k-1]*om[l-1]*(Mc[k-1]-Mc[l-1])/(om[k-1]**2-om[l-1]**2)
            qprime = mp.mpf(0)
            for i in range(len(primes)):
                qprime += -wpr[i]*S_entry(logn[i],a,k,l)
            v=qpole+qcore+qtail+qprime
            Q[k-1][l-1]=v; Q[l-1][k-1]=v
    return mp.matrix(Q), a

def lambda_a_mp(a, K, NG=40):
    Q,a=build_Q(a,K,NG)
    E=mp.eigsy(Q, eigvals_only=True)
    lam_min=min([E[i] for i in range(K)])
    # certify PD via cholesky (raises if not PD)
    pd=True
    try: mp.cholesky(Q)
    except Exception: pd=False
    return lam_min/a, pd

if __name__=="__main__":
    print(f"mpmath exact-Q λ_a (dps={mp.mp.dps}), certify λ_a>0 unconditionally")
    print("VALIDATION @ a=0.5 (float64 exact-Q gave 1.019e-6):")
    lam,pd=lambda_a_mp(0.5, K=40)
    print(f"  a=0.50: λ_a={mp.nstr(lam,8)}  PD_cholesky={pd}")
    print("\nPUSH past float64 floor:")
    for a,K in [(0.8,50),(0.9,55),(1.0,60),(1.2,70)]:
        lam,pd=lambda_a_mp(a,K)
        print(f"  a={a}: λ_a={mp.nstr(lam,6)}  PD_cholesky={pd}  {'>0 CERTIFIED' if (lam>0 and pd) else 'CHECK'}")
