"""Yakaboylu HP operator — CORRECT test: the eigenvalue mechanism is the BOUNDARY CONDITION
Psi_s(0) = -f*Gamma(s)*eta(s), eigenfunction Psi_s(x)=f*int_0^inf J0(2sqrt(xt)) t^{s-1}/(1+e^t) dt (Eq 7-8).
Quantization Psi_s(0)=0 <=> eta(s)=0 <=> zeta(s)=0 (nontrivial). Eigenvalue i(1/2-s).
Test HARD: (A) verify Psi_s(0)=Gamma(s)eta(s) via direct quadrature vs mpmath. (B) hunt zeros of
Psi_s(0) on critical line, compare to real gamma_k. (C) expose the spurious 1-2^{1-s}=0 zeros."""
import mpmath as mp
mp.mp.dps = 30
G = [mp.mpf(x) for x in open("operators/_zeros150.txt").read().split()[:12]]

def Psi0_quad(s):   # direct: int_0^inf t^{s-1}/(1+e^t) dt  (J0(0)=1)
    return mp.quad(lambda t: t**(s-1)/(1+mp.e**t), [0, mp.inf])
def Gamma_eta(s):   # closed form: Gamma(s)*eta(s), eta(s)=(1-2^{1-s})zeta(s)
    return mp.gamma(s)*(1-2**(1-s))*mp.zeta(s)

print("=== (A) Psi_s(0) == Gamma(s)*eta(s) ?  [the construction's core identity] ===")
for s in [mp.mpf('2'), mp.mpc('0.5','14.134725'), mp.mpc('0.5','21.022040'), mp.mpc('0.7','10')]:
    q, c = Psi0_quad(s), Gamma_eta(s)
    print(f"  s={complex(s):>22}: quad={complex(q):.6g}  closed={complex(c):.6g}  |diff|={float(abs(q-c)):.2e}")

print("\n=== (B) zeros of Psi_s(0) on critical line s=1/2+i*gamma  ==  real zeta zeros? ===")
def bc_on_line(g):  # |Gamma(s)eta(s)| at s=1/2+ig; zero <=> zeta zero
    return abs(Gamma_eta(mp.mpc('0.5', g)))
for gk in G:
    g0 = mp.findroot(lambda g: (Gamma_eta(mp.mpc('0.5',g))).real if False else Gamma_eta(mp.mpc('0.5',g)),
                     gk, solver='muller')
    # findroot on complex-valued: hunt where eta=0. Use eta directly.
    print(f"  seed gamma={float(gk):.6f} -> BC-zero at {float(g0.real):.6f}  |resid eta|={float(abs((1-2**(1-mp.mpc('0.5',g0)))*mp.zeta(mp.mpc('0.5',g0)))):.1e}")

print("\n=== (C) spurious zeros from (1-2^{1-s})=0 : s=1+2*pi*i*k/ln2 (eta=0 but zeta!=0) ===")
for k in [1,2]:
    s = mp.mpc(1, 2*mp.pi*k/mp.log(2))
    print(f"  k={k}: s={complex(s):.4f}  1-2^(1-s)={complex(1-2**(1-s)):.2e}  zeta(s)={complex(mp.zeta(s)):.4f} (nonzero=>spurious)")
