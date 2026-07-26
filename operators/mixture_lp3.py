#!/usr/bin/env python3
"""§324c: CORRECTED mixture-LP test (§324b part B was a complex-quadrature ARTIFACT — Φ=e^{-u} gives B=e^{-s²}
in closed form, zero-free; winding=2 was garbage). Clean statement:
  for ANY Stieltjes moment sequence m_k=∫Φ(u)u^{2k}du (Φ≥0), is F(s)=Σ_k (-1)^k m_k s^{2k}/((2k)! k!) ∈ LP?
Test via COEFFICIENTS (polyroots on the even polynomial truncation — reliable, no complex quadrature):
adversarial positive measures Φ (single/multi delta, uniform, gamma, adversarial bumps) → their m_k → truncated
F → all real roots?  Anchor: Φ=e^{-qu} ⇒ F=(1/q)e^{-(s/q)²} (Gaussian, zero-free) — analytic check of the method."""
import mpmath as mp
mp.mp.dps=40
K=26
def F_roots_from_measure(mfun):
    # m_k = ∫ Φ(u) u^{2k} du  via mfun(k); build even poly in s: sum (-1)^k m_k/((2k)! k!) s^{2k} => poly in w=s^2
    m=[mfun(k) for k in range(K)]
    cw=[(-1)**k*m[k]/(mp.factorial(2*k)*mp.factorial(k)) for k in range(K)]  # coeff of w^k, w=s^2
    # roots in w: want all real NEGATIVE (=> s real). polyroots on cw (low->high)
    rts=mp.polyroots(list(reversed(cw)),maxsteps=400,extraprec=250)
    return rts
def verdict(rts):
    # s real  <=>  w=s^2 real & <=0. Count w-roots that are real&<=0 vs complex.
    complexw=[r for r in rts if abs(mp.im(r))/max(abs(r),1)>mp.mpf('1e-10')]
    posw=[r for r in rts if mp.re(r)>1e-12 and abs(mp.im(r))/max(abs(r),1)<=1e-10]
    return len(complexw),len(posw),len(rts)
tests=[]
# (1) delta at u0: m_k=u0^{2k}
for u0 in [mp.mpf(1),mp.mpf(3)]:
    tests.append((f"delta u={u0}", lambda k,u0=u0: u0**(2*k)))
# (2) gamma Φ=u^p e^{-qu}: m_k=Γ(2k+p+1)/q^{2k+p+1}
for (p,q) in [(0,mp.mpf(1)),(2,mp.mpf(1)),(1,mp.mpf('0.5')),(5,mp.mpf(2))]:
    tests.append((f"gamma p={p} q={q}", lambda k,p=p,q=q: mp.gamma(2*k+p+1)/q**(2*k+p+1)))
# (3) uniform on [0,U]: m_k=U^{2k+1}/(2k+1)
for U in [mp.mpf(2),mp.mpf(5)]:
    tests.append((f"uniform[0,{U}]", lambda k,U=U: U**(2*k+1)/(2*k+1)))
# (4) adversarial two-bump: 0.9*delta(1)+0.1*delta(8) (very unequal, wide)
tests.append(("2delta 0.9@1+0.1@8", lambda k: mp.mpf('0.9')*1**(2*k)+mp.mpf('0.1')*mp.mpf(8)**(2*k)))
# (5) lognormal-ish Φ=exp(-(ln u)^2): m_k=∫ u^{2k} e^{-(ln u)^2} du/u ... use m_k = e^{k^2} (moments of lognormal) adversarial heavy-tail
tests.append(("lognormal m_k=e^{k^2}", lambda k: mp.e**(mp.mpf(k)**2)))
print("mixture-LP via coefficients: does F(s)=Σ(-1)^k m_k s^{2k}/((2k)!k!) have all-real zeros?",flush=True)
print(" measure Φ | #complex w-roots | #positive w-roots | total  => LP iff both 0",flush=True)
allLP=True
for name,mfun in tests:
    try:
        rts=F_roots_from_measure(mfun)
        nc,npos,tot=verdict(rts)
        lp = nc==0 and npos==0
        allLP=allLP and lp
        print(f"  {name:22s} | {nc:2d} | {npos:2d} | {tot} => {'LP' if lp else '*** NOT LP ***'}",flush=True)
    except Exception as e:
        print(f"  {name:22s} | polyroots failed: {e}",flush=True)
print("",flush=True)
print(f"ALL tested positive measures give LP mixture? {allLP}",flush=True)
print("anchor: gamma p=0 q=1 => F=e^{-s^2} exactly (Gaussian, zero-free) — method sanity.",flush=True)
print("DONE",flush=True)
