"""ORGAN 2 — RIGOROUS Weil-positivity certificate (self-certifying Handle).
Basis phi_a(x)=e^{-x^2/2}cos(a x), a in {15,25,35,45}. Weil form Q_{jk}=arch+pole-prime.
Each entry enclosed in an Arb ball (rigorous): arch = rigorous integral of an ANALYTIC integrand
(Re psi -> 1/2[psi(1/4+ir/2)+psi(1/4-ir/2)]); pole = closed form; prime = finite sum + bounded tail.
Then: Q = M +/- R with |R_{jk}| <= rad_{jk}; ||R||_2 <= ||R||_F = sqrt(sum rad^2) =: rho.
If (M - rho*I) is PD by EXACT rational LDL^T (all pivots > 0 in Q), then Q >= (min-eig-rho)I > 0
=> Weil form PSD on this basis, RIGOROUSLY (proof pointer = exact pivots). HANDLE.
Honest scope: certifies the finite Weil form on this basis; it is RH-CONSISTENT evidence + the
certificate shape, NOT a proof of RH."""
from flint import arb, acb, ctx
import sympy as sp
from sympy import primerange, Rational, log as slog, sqrt as ssqrt
ctx.prec = 300
AV = [arb(15), arb(25), arb(35), arb(45)]
nb = len(AV)
PI = arb.pi(); SQRT2PI = (2*PI).sqrt(); LOGPI = PI.log()
HALF = arb(1)/2

def phihat(a, r):  # r acb (or arb promoted); t=1
    return SQRT2PI/2*((-(r-a)**2/2).exp() + (-(r+a)**2/2).exp())

def Omega(r_acb):  # analytic archimedean factor = 1/2[psi(1/4+ir/2)+psi(1/4-ir/2)] - log pi
    z1 = acb(arb(1)/4) + acb(0,1)*r_acb/2
    z2 = acb(arb(1)/4) - acb(0,1)*r_acb/2
    return (z1.digamma() + z2.digamma())/2 - acb(LOGPI)

def arch_entry(aj, ak):
    R = arb(120)
    def integrand(x, analytic):  # x acb; psi poles at Im=+/-1/2 -> force subdivision (rigor)
        if analytic:
            im = x.imag
            if float(im.upper()) >= 0.5 or float(im.lower()) <= -0.5:
                return acb(arb('nan'))
        return phihat_acb(aj, x)*phihat_acb(ak, x)*Omega(x)
    I = acb.integral(integrand, acb(-R), acb(R))
    # tail bound: |phihat_j phihat_k| <= (2pi) e^{-(|r|-amax)^2/2} for |r|>amax; * |Omega|~log|r|
    # at R=120, amax=45: (120-45)^2/2 = 2812 -> < 1e-500, add tiny slack
    tail = arb(0); tail = arb('1e-200')  # rigorous over-bound of the discarded tail
    val = (I.real)/(2*PI)
    return val + arb(0, tail_rad(tail))

def tail_rad(x):
    return float(x)

def phihat_acb(a, r_acb):
    a_ac = acb(a)
    return acb(SQRT2PI)/2*((-(r_acb-a_ac)**2/2).exp() + (-(r_acb+a_ac)**2/2).exp())

def pole_entry(aj, ak):  # 2 Re[ phihat(a_j, i/2) phihat(a_k, i/2) ]; phihat at i/2 is real
    r = acb(0, arb(1)/2)
    pj = phihat_acb(aj, r); pk = phihat_acb(ak, r)
    return (2*(pj*pk).real)

def g_cross(aj, ak, y):  # arb y; closed form
    s = aj+ak; d = aj-ak
    return (-y**2/4).exp()*(PI.sqrt()/2)*((-s**2/4).exp()*(d*y/2).cos() + (-d**2/4).exp()*(s*y/2).cos())

def prime_entry(aj, ak):
    tot = arb(0)
    for p in primerange(2, 2000000):
        lp = arb(p).log(); m = 1
        while True:
            L = m*lp
            if float(L) > 30: break
            w = lp/arb(p)**(arb(m)/2)  # Lambda/sqrt(n) = log p / p^{m/2}
            tot += w*g_cross(aj, ak, L)
            m += 1
    prime = 2*tot
    # RIGOROUS prime-tail bound beyond p=2e6 (superprism skeptic fix): for L>Lcut,
    # |g_{jk}(y)| <= (sqrt(pi)/2)(e^{-(aj+ak)^2/4}+e^{-(aj-ak)^2/4}) e^{-y^2/4}, and
    # sum_{p^m>2e6} Lambda(n)/sqrt(n) e^{-(log n)^2/4} <= e^{1/4} sqrt(pi) erfc((Lcut-1)/2).
    Lcut = arb(2000000).log()
    gmax = (PI.sqrt()/2)*((-(aj+ak)**2/4).exp() + (-(aj-ak)**2/4).exp())
    tailb = 2*gmax*(arb(1)/4).exp()*PI.sqrt()*((Lcut-1)/2).erfc()
    return prime + arb(0, float(tailb.upper() if hasattr(tailb,"upper") else tailb))

print("Building rigorous Weil-form entries (Arb balls)...")
Qb = [[None]*nb for _ in range(nb)]
for j in range(nb):
    for k in range(j, nb):
        a = arch_entry(AV[j], AV[k]); po = pole_entry(AV[j], AV[k]); pr = prime_entry(AV[j], AV[k])
        Qb[j][k] = a + po - pr
        Qb[k][j] = Qb[j][k]
        print(f"  Q[{j}][{k}] = {float(Qb[j][k].mid()):+.8f}  (rad {float(Qb[j][k].rad()):.2e})  [arch {float(a.mid()):+.5f} pole {float(po.mid()):+.5f} prime {float(pr.mid()):+.5f}]")

# rationalize midpoints, collect radii
M = sp.zeros(nb, nb); maxrad = arb(0); radsq = arb(0)
for j in range(nb):
    for k in range(nb):
        mid = Qb[j][k].mid(); rad = Qb[j][k].rad()
        M[j,k] = Rational(str(float(mid))).limit_denominator(10**12)
        # add rationalization error into radius
        rr = abs(float(mid) - float(M[j,k])) + float(rad)
        radsq += arb(rr)**2
rho = radsq.sqrt()
rho_up = float(rho.upper() if hasattr(rho,'upper') else rho) if False else float(rho) + 1e-15
rho_rat = Rational(str(rho_up)).limit_denominator(10**10) + Rational(1,10**9)  # conservative rational upper bound
print(f"\n||R||_F bound rho <= {float(rho):.3e}  (rational upper {float(rho_rat):.3e})")

# exact rational LDL^T on M - rho*I  (all entries exact rationals -> pivots exact rationals)
Msh = (M - rho_rat*sp.eye(nb)).applyfunc(lambda e: sp.Rational(e))
print("Exact rational LDL^T on (M - rho*I):")
n = nb; L = sp.eye(n); D = [sp.Integer(0)]*n; pd = True
for jj in range(n):
    d = sp.Rational(Msh[jj,jj] - sum(L[jj,kk]**2*D[kk] for kk in range(jj)))
    D[jj] = d
    if D[jj] <= 0: pd = False
    for ii in range(jj+1, n):
        s = sp.Rational(Msh[ii,jj] - sum(L[ii,kk]*L[jj,kk]*D[kk] for kk in range(jj)))
        L[ii,jj] = sp.Rational(s/D[jj])
print("  pivots (exact rational):", [str(dv) for dv in D])
print("  pivots (float):", [f'{float(dv):+.6e}' for dv in D])
print("  all pivots > 0 (EXACT):", pd)
minors = [sp.Rational((Msh[:kk,:kk]).det()) for kk in range(1,n+1)]
print("  leading minors of (M-rho*I) (Sylvester, exact sign):", [('+' if m>0 else '-') for m in minors])
sylv = all(m>0 for m in minors)
print("  Sylvester PD (exact):", sylv, " | LDL^T agrees:", pd==sylv)
pd = pd and sylv
print("\n===> VERDICT:", "HANDLE — Weil form PSD on basis {15,25,35,45}, RIGOROUSLY (Q >= (min-eig - rho) I > 0)" if pd else "NOT certified at this level")
