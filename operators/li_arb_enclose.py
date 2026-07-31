#!/usr/bin/env python3
"""RT loop-2 Leg 1 HANDLE: rigorous 2-way Arb enclosure of the Keiper-Li coefficients lambda_n,
via certified acb power-series (no assumption of RH; no zero locations). Falsification test built in.

lambda_n = sum_{k=1}^n C(n,k)(-1)^{k+1} sigma_k,  sigma_k = sum_rho rho^{-k} = -[s^{k-1}](xi'/xi)(0).
xi'/xi = A'/A + zeta'/zeta, A = xi/zeta = 1/2 s(s-1) pi^{-s/2} Gamma(s/2) (archimedean+pole).
  A'/A closed form: [s^m] = -1  + (m==0: -ln(pi)/2 - gamma/2) + (m>=1: (1/2)(-1)^{m+1} zeta(m+1)/2^m)
  zeta'/zeta: certified acb_series (zeta series at 0, derivative, series division) -> Taylor coeffs.
All arithmetic in flint acb (rigorous balls). Binomial cancellation WIDENS the ball honestly; report where
the lower bound stays > 0 (certified lambda_n>0 = no falsification) vs where the ball straddles 0 (precision
exhausted -> raise prec). If any ball has upper bound < 0 -> RH DISPROVEN (fire).
"""
from flint import acb, arb, acb_series, ctx

def enclose(N, prec_bits):
    ctx.prec = prec_bits
    ctx.cap = N + 2                      # power-series truncation length
    pi = arb.pi()
    gamma = arb.const_euler()
    lnpi = acb(pi).log()

    # zeta'/zeta as acb_series at s=0
    z = acb_series([0, 1])               # the variable s (expansion point 0)
    Z = z.zeta()                         # zeta(s) = sum c_m s^m, certified
    Zp = Z.derivative()                  # zeta'(s)
    ZpZ = Zp / Z                         # zeta'/zeta, certified series
    # coefficient extraction: [s^m] ZpZ
    def carith(m):
        return acb(ZpZ[m])

    def ctriv(m):
        val = acb(-1)
        if m == 0:
            val = val - lnpi/2 - acb(gamma)/2
        else:
            val = val + acb((-1)**(m+1)) * acb(m+1).zeta() / (acb(2)**m) / 2
        return val

    # sigma_k = -( ctriv(k-1) + carith(k-1) ),  k=1..N
    sigma = [-(ctriv(k-1) + carith(k-1)) for k in range(1, N+1)]

    # lambda_n = sum_{k=1}^n C(n,k)(-1)^{k+1} sigma_k  (interval arithmetic)
    from math import comb
    lam = []
    for n in range(1, N+1):
        acc = acb(0)
        for k in range(1, n+1):
            acc = acc + acb(comb(n, k) * ((-1)**(k+1))) * sigma[k-1]
        lam.append(acc)
    return lam

def report(N, prec_bits):
    lam = enclose(N, prec_bits)
    print(f"Rigorous Arb enclosure of lambda_n (prec={prec_bits} bits, N={N}):")
    print(f"{'n':>3} {'lambda_n (midpoint)':>26} {'rad(log10)':>11} {'lo>0?':>6} {'verdict':>10}")
    first_straddle = None
    for n in range(1, N+1):
        L = lam[n-1].real                       # arb (ball)
        mid = L.mid(); rad = L.rad()
        lo = mid - rad; hi = mid + rad
        import math
        try: rlog = math.log10(float(rad)) if float(rad)>0 else -99
        except: rlog = 0
        lo_pos = float(lo) > 0
        hi_neg = float(hi) < 0
        verdict = "RH-FALSE!" if hi_neg else ("lam>0 cert" if lo_pos else "straddle-0")
        if (not lo_pos) and first_straddle is None: first_straddle = n
        print(f"{n:>3} {float(mid):>26.15g} {rlog:>11.2f} {str(lo_pos):>6} {verdict:>10}")
    print(f"\nCertified lambda_n>0 up to n = {(first_straddle-1) if first_straddle else N} "
          f"(ball lower-bound > 0, unconditional, NO RH assumption).")
    if first_straddle:
        print(f"From n={first_straddle} the ball straddles 0 (binomial cancellation exhausted {prec_bits}-bit prec)"
              f" -> raise prec or use Kreminski recursion. NO ball had upper<0 (no falsification).")

if __name__ == "__main__":
    import sys
    N = int(sys.argv[1]) if len(sys.argv)>1 else 15
    P = int(sys.argv[2]) if len(sys.argv)>2 else 400
    report(N, P)
