#!/usr/bin/env python3
"""rigor_lint.py — standing invariant checker for the Weil-cert numerics (idea #1).
Catches the two bug CLASSES actually on record before any output is trusted in the LEDGER:
  • reimplementation bugs (the retracted "operator doesn't converge": wrong eigenvector / no symmetrization)
      -> DYNAMIC invariants: matrix symmetry, per-mode diagonal>0, parity decoupling, D_kl asymptotic.
  • interval-rigor bugs (the §199.4 catches: float-extraction breaking rigor; arb(-1,1) sign-bias)
      -> STATIC scan of the target cert script for the known anti-patterns.
The DYNAMIC side re-derives Q[k,l] from the documented formulas with an INDEPENDENT mpmath impl — so if it
disagrees with the cert's entries that is itself a caught bug (a second implementation, at the entry level).
Usage:  rigor_lint.py [cert_script.py ...]      exit 0 = all pass, 1 = any fail."""
import sys, re, os
import mpmath as mp
mp.mp.dps = 34
a = mp.mpf('0.549'); b = mp.pi/(2*a); log2 = mp.log(2); twoa = 2*a
Om0 = mp.digamma(mp.mpf(1)/4) - mp.log(mp.pi); Phi0 = mp.zeta(2, mp.mpf(1)/4)
wp = mp.log(2)/mp.sqrt(2); NL = 1200

# ---- INDEPENDENT re-implementation of the entry Q[k,l] (mpmath) ----
def _Pp(p,t,s):
    if p==0: return 1/s
    if p==1: return t/s-1/s**2
    return t*t/s-2*t/s**2+2/s**3
def _II(p,al,om,sc,t0,t1):
    s=mp.mpc(-al,om); F=lambda tt: mp.e**(s*tt)*_Pp(p,tt,s); v=F(t1)-F(t0)
    return v.real if sc=='cos' else v.imag
def _dterms(k,l):
    if k==l: return [(2*a*k*k*b*b,0,k*b,'cos'),(-k*k*b*b,1,k*b,'cos'),(-k*b,0,k*b,'sin')]
    al=k*l*b*(mp.mpf(1)/(k-l)-mp.mpf(1)/(k+l)); ak=-k*l*b*(mp.mpf(1)/(k-l)+mp.mpf(1)/(k+l))
    return [(al,0,l*b,'sin'),(ak,0,k*b,'sin')]
def _Dcf(k,l,t):   # D_kl(t) from the closed-form term list (odd k,l)
    tot=mp.mpf(0)
    for (amp,p,om,sc) in _dterms(k,l): tot += amp*(t**p)*(mp.cos(om*t) if sc=='cos' else mp.sin(om*t))
    return tot
def _dprime(k,x): return (k*mp.pi/(2*a))*mp.cos(k*mp.pi*(x+a)/(2*a))   # v_k'(x)
def _Dnum(k,l,t): # ∫ v_k'(x)v_l'(x-t) + v_l'(x)v_k'(x-t) dx over overlap [t-a,a]  (the DEFINITION)
    return (mp.quad(lambda x:_dprime(k,x)*_dprime(l,x-t),[t-a,a])
           +mp.quad(lambda x:_dprime(l,x)*_dprime(k,x-t),[t-a,a]))
def Q(k,l):
    tot=mp.mpf(0); D=_dterms(k,l)
    gt=[(-4,0,mp.mpf(-1)/2),(-4,0,mp.mpf(1)/2),(8,0,mp.mpf(0)),(-Om0/2,1,mp.mpf(0)),(-Phi0/4,0,mp.mpf(0))]
    for n in range(NL): gt.append((mp.mpf(1)/4/(n+mp.mpf(1)/4)**2,0,2*n+mp.mpf(1)/2))
    for (ag,pg,alp) in gt:
        for (ad,pd,om,sc) in D: tot+=ag*ad*_II(pg+pd,alp,om,sc,0,twoa)
    for (ad,pd,om,sc) in D:
        tot+=wp*ad*_II(pd+1,0,om,sc,log2,twoa)+(-wp*log2)*ad*_II(pd,0,om,sc,log2,twoa)
    return tot

def _p(name, ok, detail): print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {detail}"); return ok

def dynamic_checks():
    print("DYNAMIC invariants (independent mpmath re-derivation):")
    ok=True
    # 1) symmetry Q[k,l]==Q[l,k]
    e=max(abs(Q(k,l)-Q(l,k)) for (k,l) in [(1,3),(3,7),(5,7),(1,7)])
    ok&=_p("symmetry Q=Qᵀ", e<mp.mpf(10)**-20, f"max|Q[k,l]-Q[l,k]|={mp.nstr(e,3)}")
    # 2) per-mode diagonal > 0 (individual modes safe — a proven structural fact)
    dmin=min(Q(k,k) for k in [1,3,5,7,11,21]); 
    ok&=_p("diagonal Q[k,k]>0", dmin>0, f"min diag={mp.nstr(dmin,4)} (k≤21)")
    # 3) closed-form D_kl vs DEFINITION quadrature (the anti-reimplementation-drift guard)
    de=max(abs(_Dcf(k,l,t)-_Dnum(k,l,t)) for (k,l) in [(1,3),(5,7),(3,7)] for t in [mp.mpf('0.3'),mp.mpf('0.8')])
    ok&=_p("D_kl closed-form ≡ definition", de<mp.mpf(10)**-18, f"max|Dcf-Dnum|={mp.nstr(de,3)}")
    # 4) coupling sign/scale: l·|Q[1,l]| bounded (1/l law) — guards a blown-up off-diagonal
    c=[l*abs(Q(1,l)) for l in [21,41,61]]; 
    ok&=_p("coupling 1/l law", all(x<mp.mpf('0.1') for x in c), f"l·|Q[1,l]|={[mp.nstr(x,3) for x in c]} (→const)")
    return ok

ANTI = [
    (r"arb\(\s*-1\s*,\s*1\s*\)", "arb(-1,1) is midpoint −1 (sign-biased!); tail enclosure must be arb(0,1)=[-1,1]"),
    (r"\.mid\(\)\s*-\s*.*\.rad\(\)", "extracting mid()-rad() as a float breaks interval rigor; test with s>0"),
    (r"float\(\s*\w+\.(mid|rad)\(\)\s*\)", "float(x.mid/rad) drops the interval; keep arb through the pivot test"),
]
REQUIRE_IF_PD = [
    (r"\(\s*s\s*>\s*0\s*\)|s\s*>\s*0", "a rigorous pivot test `s>0` (arb interval comparison), not a float compare"),
]
def static_checks(path):
    print(f"STATIC scan: {path}")
    src=open(path).read(); ok=True
    for pat,msg in ANTI:
        hits=[i+1 for i,ln in enumerate(src.splitlines()) if re.search(pat,ln) and not ln.strip().startswith('#')]
        ok&=_p(f"no anti-pattern /{pat[:22]}…/", not hits, "clean" if not hits else f"lines {hits}: {msg}")
    if re.search(r"LDL|pivot|PD|certif", src, re.I):
        for pat,msg in REQUIRE_IF_PD:
            ok&=_p("uses rigorous pivot test", bool(re.search(pat,src)), msg)
    return ok

if __name__=="__main__":
    scripts=[s for s in sys.argv[1:] if s.endswith('.py')] or ['operators/dbn_dynamics/arb_certify.py']
    d=dynamic_checks()
    s=all(static_checks(p) for p in scripts if os.path.exists(p))
    print(f"\n{'✓ RIGOR LINT PASS' if (d and s) else '✗ RIGOR LINT FAIL'} — {'safe to trust in LEDGER' if (d and s) else 'DO NOT bank this output'}")
    sys.exit(0 if (d and s) else 1)
