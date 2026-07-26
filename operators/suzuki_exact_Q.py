#!/usr/bin/env python3
"""suzuki_exact_Q.py -- §242: QUADRATURE-FREE assembly of the Suzuki form Q_{kl}=∫∫g(x-y)v_k'(x)v_l'(y).

Fixes the a>1 wall of suzuki_screw_lambda_a.py (trapezoid over kinks + catastrophic e^{a/2}
cancellation). Method: integrate by parts twice: Q_{kl} = -∫∫ g''(x-y) v_k(x) v_l(y) (v_k(±a)=0),
then each piece of g (Suzuki eq 1.3) contributes in CLOSED FORM:
  g = boundary + prime + arch_core + arch_tail with
  boundary  = -4(e^{t/2}+e^{-t/2}-2):  g'' = -(e^{t/2}+e^{-t/2}); const drops EXACTLY.
      Q_pole = E_{1/2}(k)E_{-1/2}(l) + E_{-1/2}(k)E_{1/2}(l),
      E_p(k) = ∫e^{px}v_k = ω_k[e^{-pa} - (-1)^k e^{pa}]/(p²+ω_k²),  ω_k=kπ/(2a).
  arch_core = c|t|, c=-(ψ(¼)-logπ)/2:  g''=2cδ(t) ⟹ Q_core = -2c·a·δ_kl  (∫v_kv_l = a δ_kl).
  prime     = Σ Λ(n)/√n (|t|-log n)_+:  g'' = Σ Λ/√n [δ(t-log n)+δ(t+log n)] ⟹
      Q_prime = -Σ Λ(n)/√n · S_kl(log n),  S_kl(s)=∫v_k(y+s)v_l(y)dy + (k↔l)  (EXACT trig):
      S_kk(s) = (L-s)cos(ω_k s) + sin(ω_k s)/ω_k,   L=2a
      S_kl(s) = 0 for k+l odd (exact parity decoupling);
              = (sin(ω_l s)-sin(ω_k s))/(ω_k-ω_l) + (sin(ω_k s)+sin(ω_l s))/(ω_k+ω_l), k+l even.
  arch_tail = -¼[Φ(1,2,¼) - h(|t|)], h(u)=e^{-u/2}Φ(e^{-2u},2,¼):
      h'(u) = -2e^{-u/2}Φ(e^{-2u},1,¼)  (EXACT identity: (2n+½)/(n+¼)² = 2/(n+¼)).
      One by-parts: Q_tail = ½∫_0^L e^{-t/2}Φ(e^{-2t},1,¼) Ψ_kl(t) dt   (log-singular at 0 only),
      Ψ_kl(t) = Φ_kl(t)-Φ_kl(-t):  0 for k+l odd;
              = 2ω_kω_l(cos(ω_k t)-cos(ω_l t))/(ω_k²-ω_l²) for k≠l even sum;
              = ω_k(L-t)sin(ω_k t) for k=l.
      Quadrature: t = L·s³ Gauss-Legendre mapping (tames the log), NG≈320 nodes -> ~machine.
Sign checks: I=-∫∫K''vv verified on K=t²/2 (gives -∫v_k∫v_l = direct expansion). Parity decoupling
falls out of BOTH S and Ψ (matches §192's exact even/odd decoupling — internal consistency).
Also provides: exact RHS f_±(k)=E_{±1}(k), exact Fourier ŝ_k(z), W(a,θ;z) — full pipeline
quadrature-free.  VALIDATION gate in main: vs old lambda_a at a∈{0.3,0.5,0.8} + K/NG convergence +
Thm 1.5 reality check. Per-prime pieces exposed (S_matrices) for Hellmann-Feynman (§241 next phase).
"""
import numpy as np
from scipy.special import digamma
import mpmath as mp

PSI = float(digamma(0.25)) - float(np.log(np.pi))     # ψ(¼)-logπ ≈ -5.3722
C_CORE = -PSI / 2.0                                    # ≈ +2.6861

def vonmangoldt_upto(N):
    Lam = np.zeros(N + 1)
    is_comp = np.zeros(N + 1, dtype=bool)
    for p in range(2, N + 1):
        if not is_comp[p]:
            is_comp[2*p::p] = True
            lp = np.log(p); pk = p
            while pk <= N:
                Lam[pk] = lp; pk *= p
    return Lam

_EULER = 0.5772156649015329
def lerch1_quarter_t(t):
    """Φ(e^{-2t},1,¼) as a function of t>0. Series for t>=0.05; exact z->1 asymptotic
    Φ(z,1,s) = -ln(1-z) - ψ(s) - γ + O((1-z)ln(1-z)) for t<1e-3 (float64-safe: uses
    -ln(1-e^{-2t}) = -ln(2t) + t - t²/3·... via log1p/expm1); mpmath in between."""
    if t >= 0.05:
        z = np.exp(-2*t)
        tot, zn, n = 0.0, 1.0, 0
        while True:
            trm = zn/(n+0.25); tot += trm
            if trm < 1e-17*tot: break
            zn *= z; n += 1
            if n > 4000: break
        return tot
    with mp.workdps(25):
        return float(mp.lerchphi(mp.e**(-2*mp.mpf(float(t))), 1, mp.mpf(1)/4))

def E_p(p, a, kk):
    om = kk*np.pi/(2*a)
    return om*(np.exp(-p*a) - ((-1.0)**kk)*np.exp(p*a))/(p*p + om*om)

def S_mat(s, a, K):
    """S_kl(s) exact shift-overlap (symmetrized), 0<=s<L=2a."""
    L = 2*a
    kk = np.arange(1, K+1); om = kk*np.pi/L
    S = np.zeros((K, K))
    A = om[:, None]; B = om[None, :]
    par = ((kk[:, None] + kk[None, :]) % 2 == 0)
    with np.errstate(divide='ignore', invalid='ignore'):
        off = (np.sin(B*s) - np.sin(A*s))/(A - B) + (np.sin(A*s) + np.sin(B*s))/(A + B)
    dg = (L - s)*np.cos(om*s) + np.sin(om*s)/om
    S = np.where(par, off, 0.0)
    np.fill_diagonal(S, dg)
    return S

def tail_nodes(a, NG=320):
    """Gauss-Legendre on s∈(0,1], t=L s³; returns t nodes, weights (dt=3Ls²ds), kernel w(t)."""
    L = 2*a
    xs, ws = np.polynomial.legendre.leggauss(NG)
    s = 0.5*(xs + 1.0); wgt = 0.5*ws
    t = L*s**3; jac = 3*L*s**2
    ker = np.array([0.5*np.exp(-tt/2)*lerch1_quarter_t(tt) for tt in t])
    return t, wgt*jac*ker

def Q_tail(a, K, NG=320):
    L = 2*a
    kk = np.arange(1, K+1); om = kk*np.pi/L
    t, w = tail_nodes(a, NG)                              # (NG,)
    cos_t = np.cos(np.outer(om, t))                        # (K,NG)
    A = om[:, None, None]; B = om[None, :, None]
    par = ((kk[:, None] + kk[None, :]) % 2 == 0)
    with np.errstate(divide='ignore', invalid='ignore'):
        Psi = 2*A*B*(cos_t[:, None, :] - cos_t[None, :, :])/(A*A - B*B)   # (K,K,NG)
    Q = np.einsum('klt,t->kl', np.where(par[:, :, None], Psi, 0.0), w)
    # diagonal: Ψ_kk = ω_k (L-t) sin(ω_k t)
    sin_t = np.sin(np.outer(om, t))
    dg = np.einsum('kt,t->k', om[:, None]*(L - t)[None, :]*sin_t, w)
    np.fill_diagonal(Q, dg)
    return Q

def build_Q(a, K, NG=320, return_pieces=False):
    kk = np.arange(1, K+1)
    Ehalf_p = E_p(0.5, a, kk); Ehalf_m = E_p(-0.5, a, kk)
    Qpole = np.outer(Ehalf_p, Ehalf_m) + np.outer(Ehalf_m, Ehalf_p)
    Qcore = -2*C_CORE*a*np.eye(K)
    Qtail = Q_tail(a, K, NG)
    L = 2*a
    Nmax = int(np.floor(np.exp(L)))
    Lam = vonmangoldt_upto(max(Nmax, 2))
    Qprime = np.zeros((K, K)); pieces = {}
    for n in range(2, Nmax+1):
        if Lam[n] > 0 and np.log(n) < L:
            Sn = S_mat(np.log(n), a, K)
            Qn = -(Lam[n]/np.sqrt(n))*Sn
            Qprime += Qn
            if return_pieces: pieces[n] = Qn
    Q = Qpole + Qcore + Qtail + Qprime
    Q = (Q + Q.T)/2
    if return_pieces:
        return Q, {"pole": Qpole, "core": Qcore, "tail": Qtail, "prime": Qprime, "per_prime": pieces}
    return Q

def lambda_a_exact(a, K=80, NG=320):
    Q = build_Q(a, K, NG)
    return np.linalg.eigvalsh(Q)[0]/a

def shat(z, a, kk):
    """ŝ_k(z)=∫v_k e^{izx}dx = e^{-iza} ω_k(1-(-1)^k e^{2iza})/(ω_k²-z²)."""
    om = kk*np.pi/(2*a)
    zz = z if abs(abs(z)-om.min()) > 1e-9 else z + 1e-9
    den = om*om - zz*zz
    den = np.where(np.abs(den) < 1e-12, 1e-12, den)
    return np.exp(-1j*zz*a)*om*(1.0 - ((-1.0)**kk)*np.exp(2j*zz*a))/den

def W_family(a, K=80, NG=320, lam_shift=0.2):
    """Returns W(theta,z) callable + lam_a. Solve (Q-λ aI)c=f_± exactly like old pipeline."""
    from scipy.linalg import solve, eigvalsh
    kk = np.arange(1, K+1)
    Q = build_Q(a, K, NG)
    lam_a = np.linalg.eigvalsh(Q)[0]/a
    lam = lam_a - lam_shift
    Ta = Q - lam*a*np.eye(K)
    fplus = E_p(1.0, a, kk); fminus = E_p(-1.0, a, kk)
    cplus = solve(Ta, fplus, assume_a='sym')
    cminus = solve(Ta, fminus, assume_a='sym')
    npl = np.sqrt(abs(cplus @ fplus)); nmi = np.sqrt(abs(cminus @ fminus))
    cminus = cminus*(npl/nmi)
    def W(theta, z):
        vp = np.sum(cplus*shat(z, a, kk)); vm = np.sum(cminus*shat(z, a, kk))
        return (z - 1j)*vp + np.exp(1j*theta)*(z + 1j)*vm
    return W, lam_a

if __name__ == "__main__":
    import suzuki_screw_lambda_a as OLD
    print("VALIDATION GATE — exact-Q vs old quadrature pipeline")
    print(f"{'a':>5} {'λ_a exact':>13} {'λ_a old':>13} {'rel diff':>10}   K-conv(K=60 vs 90)  NG-conv(320 vs 640)")
    for a in [0.3, 0.5, 0.8]:
        le = lambda_a_exact(a, K=80, NG=320)
        lo = OLD.lambda_a(a, M=3000, K=40)
        le60 = lambda_a_exact(a, K=60, NG=320); le90 = lambda_a_exact(a, K=90, NG=320)
        leN = lambda_a_exact(a, K=80, NG=640)
        print(f"{a:5.2f} {le:13.6e} {lo:13.6e} {abs(le-lo)/abs(lo):10.2e}   "
              f"{abs(le60-le90):.2e}          {abs(le-leN):.2e}")
    print("\nλ_a at a>1 (OLD scheme broke here — §42 'spuriously negative'):")
    for a in [1.0, 1.5, 2.0, 2.5, 3.0]:
        le = lambda_a_exact(a, K=max(80, int(30*a)+40), NG=400)
        print(f"  a={a}: λ_a = {le:.6e}  ({'POSITIVE ✓ (Yoshida/RH-consistent)' if le > 0 else 'NEGATIVE — investigate'})")
    print("\nThm 1.5 reality check at a=1.5 (|W| on-axis dips vs off-axis 0.5i):")
    W, lam_a = W_family(1.5, K=90, NG=400)
    zg = np.linspace(5, 45, 4000)
    Wz = np.abs([W(0.0, z) for z in zg])
    med = np.median(Wz)
    mins = [(zg[i], Wz[i]) for i in range(1, len(zg)-1)
            if Wz[i] < Wz[i-1] and Wz[i] < Wz[i+1] and Wz[i] < 0.02*med]
    print(f"  λ_a={lam_a:.3e}; zeros(θ=0): {[round(z,3) for z,_ in mins]}")
    for z0, w0 in mins[:6]:
        print(f"    z={z0:7.3f}  |W|on={w0:.2e}  |W(z+0.5i)|off={abs(W(0.0, z0+0.5j)):.2e}  ratio={abs(W(0.0,z0+0.5j))/(w0+1e-30):.0f}x")
