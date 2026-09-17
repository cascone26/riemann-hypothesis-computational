"""
Suzuki 2301.00421 screw-line S_t(z) evaluator — clean (OCR-independent) build.

Definitions (all from the paper, no reliance on the garbled (1.6) OCR):
  xi(s)                                 completed zeta, real on the critical line
  E_xi(z)  = xi(1/2 - i z) + xi'(1/2 - i z)                         (1.3)
  E#_xi(z) = xi(1/2 + i z) + xi'(1/2 + i z)   [ = E_xi^sharp ]
  Theta_xi(z)  = E#_xi(z) / E_xi(z)                                 (1.4)
  Theta#_xi(z) = 1 / Theta_xi(z)   = E_xi(z) / E#_xi(z)
  P_t(z)  = sum_{gamma in Gamma} (1/m_gamma) * (e^{-i gamma t} - 1)/gamma * 1/(z - gamma)   (3.2)
            Gamma = { +/- gamma_n },  gamma_n = n-th zeta zero ordinate,  m_gamma = 1 (RH simple)
  S_t(z)  = (i/2) * (1 + Theta#_xi(z)) * P_t(z)                     (1.5)

VALIDATION GATE: under RH, Theta_xi is an inner function, so |Theta_xi(x)| = 1 for all real x.
If our E_xi / E#_xi transcription is correct this must hold numerically. Run first.
"""
import mpmath as mp

mp.mp.dps = 30

def xi(s):
    # completed zeta: xi(s) = 1/2 * s(s-1) * pi^{-s/2} * Gamma(s/2) * zeta(s)
    return mp.mpf('0.5') * s * (s - 1) * mp.pi**(-s/2) * mp.gamma(s/2) * mp.zeta(s)

def xip(s):
    return mp.diff(xi, s)

def E_xi(z):
    s = mp.mpf('0.5') - 1j*z
    return xi(s) + xip(s)

def Esharp_xi(z):
    s = mp.mpf('0.5') + 1j*z
    return xi(s) + xip(s)

def Theta_xi(z):
    return Esharp_xi(z) / E_xi(z)

def Thetasharp_xi(z):
    return E_xi(z) / Esharp_xi(z)

# cache zeta-zero ordinates
_ZC = []
def gammas(N):
    global _ZC
    while len(_ZC) < N:
        _ZC.append(mp.im(mp.zetazero(len(_ZC)+1)))
    return _ZC[:N]

def P_t(z, t, N=200):
    z = mp.mpc(z)
    acc = mp.mpc(0)
    for g in gammas(N):
        for gam in (g, -g):
            acc += (mp.e**(-1j*gam*t) - 1)/gam * 1/(z - gam)
    return acc

def S_t(z, t, N=200):
    return (1j/2) * (1 + Thetasharp_xi(z)) * P_t(z, t, N)


if __name__ == "__main__":
    print("=== VALIDATION GATE: |Theta_xi(x)| == 1 for real x (RH inner-function boundary) ===")
    ok = True
    for x in [0.5, 1.3, 3.0, 7.0, 14.5, 21.0]:
        val = Theta_xi(mp.mpf(x))
        m = abs(val)
        flag = "OK" if abs(m - 1) < 1e-12 else "**FAIL**"
        if abs(m-1) >= 1e-12: ok = False
        print(f"  x={x:6.2f}   |Theta_xi(x)| = {mp.nstr(m, 16)}   {flag}")
    print(f"GATE: {'PASSED' if ok else 'FAILED'} — E_xi/E#_xi transcription {'validated' if ok else 'WRONG'}")

    if ok:
        print("\n=== S_t(z) sample values (N=120 zero pairs) ===")
        for t in [0.5, 1.0]:
            for x in [0.3, 2.0, 5.0]:
                v = S_t(mp.mpf(x), mp.mpf(t), N=120)
                print(f"  t={t:.2f} x={x:.2f}:  S_t = {mp.nstr(v, 10)}   (Im/Re ratio {mp.nstr(mp.im(v)/mp.re(v) if mp.re(v)!=0 else mp.inf,4)})")
