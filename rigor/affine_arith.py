#!/usr/bin/env python3
"""affine_arith.py — first-order AFFINE ARITHMETIC for the a-parameter, the cure for the dependency blowup
that defeats naive interval arithmetic on an a-panel (§199.6: plain Arb over an a-interval over-widens ~400×
because it treats every occurrence of `a` as independent; entries actually move together and ~linearly in a).

Each quantity is  x = x0 + x1·ε + [±r],  ε ∈ [−1,1] the SHARED normalized parameter deviation
(ε = (a − a_center)/(h/2)), x0 center, x1 first-order sensitivity, r a rigorous remainder radius. Because ε is
SHARED, correlated variation cancels instead of compounding — the standard dependency-problem cure.

Rebuilt clean (the first agent build was a facade: its entry returned ~1e-51 and its 'naive' baseline was a
synthetic proxy — kept at affine_arith.agent.py.bak). This version demonstrates on a REAL Weil-entry sub-integral
with a GENUINE flint-Arb naive baseline, showing the textbook a-cancellation kβ·2a = kπ.
"""
import mpmath as mp
mp.mp.dps = 40
mpf = mp.mpf

class AffineForm:
    __slots__ = ('x0','x1','r')
    def __init__(self, x0, x1=0, r=0):
        self.x0 = mpf(x0); self.x1 = mpf(x1); self.r = abs(mpf(r))
    def enclosure(self):
        rad = abs(self.x1) + self.r
        return (self.x0 - rad, self.x0 + rad)
    def width(self):
        lo,hi = self.enclosure(); return hi - lo
    def __add__(s, o):
        return AffineForm(s.x0+o.x0, s.x1+o.x1, s.r+o.r) if isinstance(o,AffineForm) else AffineForm(s.x0+o, s.x1, s.r)
    __radd__ = __add__
    def __neg__(s): return AffineForm(-s.x0, -s.x1, s.r)
    def __sub__(s, o): return s + (-o if isinstance(o,AffineForm) else -mpf(o))
    def __rsub__(s, o): return (-s) + o
    def __mul__(s, o):
        if not isinstance(o, AffineForm):
            o = mpf(o); return AffineForm(s.x0*o, s.x1*o, s.r*abs(o))
        x0,x1,r = s.x0,s.x1,s.r; y0,y1,q = o.x0,o.x1,o.r
        # ε² and all cross/remainder terms folded conservatively into the new remainder (rigorous over |ε|≤1)
        nr = abs(x1*y1) + abs(x0)*q + abs(y0)*r + abs(x1)*q + abs(y1)*r + r*q
        return AffineForm(x0*y0, x0*y1 + x1*y0, nr)
    __rmul__ = __mul__
    def _tayl(s, f, f1, f2max):
        d = abs(s.x1) + s.r
        nr = abs(f1(s.x0))*s.r + mpf('0.5')*f2max*d*d
        return AffineForm(f(s.x0), f1(s.x0)*s.x1, nr)
    def exp(s):
        lo,hi = s.enclosure(); return s._tayl(mp.exp, mp.exp, mp.exp(hi))         # f''=exp, max at hi
    def sin(s):
        return s._tayl(mp.sin, mp.cos, mpf(1))                                    # |sin''|≤1
    def cos(s):
        return s._tayl(mp.cos, lambda x:-mp.sin(x), mpf(1))                       # |cos''|≤1
    def recip(s):
        lo,hi = s.enclosure()
        if lo <= 0 <= hi: raise ZeroDivisionError("0 in enclosure — reciprocal not analytic")
        m = min(abs(lo),abs(hi)); return s._tayl(lambda x:1/x, lambda x:-1/(x*x), 2/m**3)  # |f''|=2/|x|³≤2/m³
    def __truediv__(s,o): return s * (o.recip() if isinstance(o,AffineForm) else AffineForm(mpf(1)/mpf(o)))
    def __repr__(s):
        lo,hi=s.enclosure(); return f"{float(s.x0):.6e}+{float(s.x1):.3e}ε±{float(s.r):.2e} =[{float(lo):.6e},{float(hi):.6e}]"

# ---------- REAL Weil-entry sub-integral, as an AffineForm in a ----------
# P_k(a) = ∫_0^{2a} (-4 e^{-t/2})·sin(kβ t) dt,  β = π/(2a)  — a genuine boundary×arch piece of Q[k,l].
# Closed form: -4·Im[(E-1)/s], s=-1/2+iω, ω=kβ, E=e^{(-1/2+iω)2a}=e^{-a}(cos(ω2a)+i sin(ω2a)); ω·2a=kπ (cancels!).
def a_form(a_center, h): return AffineForm(a_center, mpf(h)/2, 0)   # a = a_center + (h/2)ε
def Pk_affine(k, a_center, h):
    a = a_form(a_center, h); beta = AffineForm(mp.pi)/(a*2); w = beta*k; twoa = a*2
    ea = (a*(-1)).exp(); ang = w*twoa                # e^{-a}, and ω·2a (≈kπ; tiny linear coeff = the cancellation)
    Ere = ea*ang.cos(); Eim = ea*ang.sin()
    den = w*w + mpf('0.25')                          # |s|² = 1/4 + ω²
    im = (Eim*mpf('-0.5') + (Ere - 1)*(w*(-1))) / den   # Im[(E-1)·conj(s)]/|s|²
    return im*(-4)
def Pk_true(k, a):
    a=mpf(a); beta=mp.pi/(2*a)
    return -4*mp.quad(lambda t: mp.e**(-t/2)*mp.sin(k*beta*t), [0, 2*a])
def Pk_naive_arb(k, a0, a1):
    from flint import arb
    a = arb((a0+a1)/2, (a1-a0)/2)                    # ball [a0,a1] (arb is mid,rad) — every use independent (dependency problem)
    beta = arb.pi()/(2*a); w = k*beta; twoa = 2*a
    ea = (-a).exp(); ang = w*twoa
    Ere = ea*ang.cos(); Eim = ea*ang.sin()
    den = w*w + arb('0.25')
    return (-4)*((Eim*arb('-0.5') + (Ere-1)*(-w))/den)

if __name__=="__main__":
    print("=== AffineForm self-tests ===")
    e=AffineForm(0,1,0)
    q=(1+e)*(1+e); lo,hi=q.enclosure()               # (1+ε)² true range [0,4]
    print(f" (1+ε)² enclosure [{float(lo):.4f},{float(hi):.4f}] ⊇ true [0,4]?  {lo<=0 and hi>=4}")
    x=AffineForm(mpf('0.5'),mpf('0.1'),0)
    for nm,af,tru in [('exp',x.exp(),mp.exp),('sin',x.sin(),mp.sin),('cos',x.cos(),mp.cos),('recip',x.recip(),lambda v:1/v)]:
        lo,hi=af.enclosure(); ok=all(lo<=tru(mpf('0.5')+mpf('0.1')*eps)<=hi for eps in [-1,-0.5,0,0.5,1])
        print(f" {nm}: enclosure covers true over ε∈[-1,1]?  {ok}")
    print("\n=== REAL Weil sub-integral P_k(a)=∫_0^{2a}(-4e^{-t/2})sin(kβt)dt : naive-Arb vs affine (k=3) ===")
    print(f"{'panel h':>8} {'naive width':>14} {'affine width':>14} {'ratio':>9} {'affine⊇true?':>12}")
    for a0,a1 in [(0.48,0.56),(0.50,0.54),(0.51,0.53),(0.515,0.525)]:
        ac=(a0+a1)/2; h=a1-a0; k=3
        nb=Pk_naive_arb(k,a0,a1); nw=float(nb.rad())*2
        af=Pk_affine(k,ac,h); lo,hi=af.enclosure(); aw=float(hi-lo)
        cover=all(lo<=Pk_true(k,ac+(mpf(h)/2)*eps)<=hi for eps in [-1,-0.5,0,0.5,1])
        print(f"{h:>8.3f} {nw:>14.3e} {aw:>14.3e} {nw/aw:>8.1f}x {str(cover):>12}")
    print("\nREAD: naive-Arb width ~flat (dependency blowup); affine width ∝ h and far tighter, while rigorously")
    print("enclosing the true entry. This primitive unblocks the continuum; SOS (sos_continuum.py) then certifies")
    print("PD across a whole window from the affine enclosure.")
