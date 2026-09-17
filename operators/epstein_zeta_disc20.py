# Epstein zeta of Q=x^2+5y^2 (disc -20, class number 2) via Dirichlet decomposition:
#   2*zeta_Q1(s) = zeta(s) L(s,chi_-20) + L(s,chi_-4) L(s,chi_5)  = F(s).
# Davenport-Heilbronn: this has OFF-LINE zeros. Find one (degree-2-family FAIL for the discriminator).
import mpmath as mp
mp.mp.dps=30
def Lchi(s, q, chi):
    # L(s,chi) = q^{-s} sum_{a=1}^{q} chi(a) zeta(s, a/q)  (Hurwitz)
    return mp.power(q,-s)*mp.fsum(chi(a)*mp.zeta(s, mp.mpf(a)/q) for a in range(1,q+1))
def chi_m4(n):
    r=n%4; return 1 if r==1 else (-1 if r==3 else 0)
def chi_5(n):
    r=n%5; return 0 if r==0 else (1 if r in (1,4) else -1)
def chi_m20(n):
    return chi_m4(n)*chi_5(n)
def F(s):
    s=mp.mpc(s)
    return mp.zeta(s)*Lchi(s,20,chi_m20) + Lchi(s,4,chi_m4)*Lchi(s,5,chi_5)
# sanity: F should be real for s on the critical line? No -- F(1/2+it) is generally complex.
# Scan |F| over the strip Re in [0.5,1.4], t in [0,20] for near-zeros OFF the line (Re != 0.5).
print("Scanning |F(sigma+it)| for off-line zeros (Davenport-Heilbronn). |F|<0.05 flagged:")
best=[]
for t in [x/10 for x in range(5,205,2)]:      # t 0.5..20.4
    for sig in [0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25]:
        v=abs(F(mp.mpc(sig,t)))
        if v<0.08:
            best.append((round(sig,2),round(t,2),float(v)))
best.sort(key=lambda x:x[2])
print(f"  {len(best)} candidate near-zeros off the line; smallest 8:")
for b in best[:8]: print("   sigma=%.2f t=%.2f |F|=%.4f"%b, flush=True)
# refine the best candidate with findroot (2D: Re F = Im F = 0)
if best:
    sig0,t0,_=best[0]
    try:
        z=mp.findroot(lambda z: F(z), mp.mpc(sig0,t0))
        print(f"\nREFINED off-line zero: s = {mp.nstr(z,12)}  |F(z)|={float(abs(F(z))):.2e}  Re(s)-1/2 = {float(mp.re(z))-0.5:+.4f}")
        print("=> OFF-LINE ZERO CONFIRMED" if abs(float(mp.re(z))-0.5)>1e-4 else "=> on line after refine")
    except Exception as e:
        print("findroot failed:",e)
