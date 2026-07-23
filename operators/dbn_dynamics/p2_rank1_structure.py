#!/usr/bin/env python3
"""§187: is p=2-only Weil positivity a RANK-1 UPDATE to an arch form with EXACTLY ONE negative eigenvalue?
The single p=2 place = w|v̂(log2)|², w=Λ(2)/√2 — a rank-1 (rank-2 with ±) form. §186: arch-alone min-eig ≈ −0.42
(a=0.549), full form +2.5e-5. Hypothesis: arch-only (A+B) has EXACTLY ONE negative eigenvalue, and the p=2 rank-1
term lifts precisely that mode to the +2.5e-5 floor ⟹ positivity ⟺ a scalar secular condition. TEST it:
 (1) full generalized spectrum of Qa (arch-only) vs Mass — count negatives, report λ_1(neg), λ_2.
 (2) rank + signs of the prime part Qp = Qf − Qa (rank-1? rank-2? PSD?).
 (3) interlacing check: does λ_min(Qf) sit between λ_min(Qa) and λ_2(Qa), consistent with a rank-1 lift?
Reuses weil_matrices from mpmath_weil_hiprec (additive)."""
import sys, math
import mpmath as mp
sys.path.insert(0, '/Users/scones/projects/riemann/operators/dbn_dynamics')
from mpmath_weil_hiprec import weil_matrices, gen_eig_sym

def spectrum(Q, Mass):
    E, V = gen_eig_sym(Q, Mass)
    return list(E), V

if __name__ == "__main__":
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    K   = int(sys.argv[2]) if len(sys.argv) > 2 else 24
    M   = int(sys.argv[3]) if len(sys.argv) > 3 else 600
    avals = [mp.mpf(x) for x in (sys.argv[4].split(',') if len(sys.argv) > 4 else ['0.45', '0.50', '0.549'])]
    mp.mp.dps = dps
    print(f"§187 p=2 rank-1 structure | dps={dps} K={K} M={M}")
    for a in avals:
        Qf, Qa, Mass, Vv, xs, w = weil_matrices(a, K, M)
        Ea, _ = spectrum(Qa, Mass)      # arch-only spectrum
        Ef, _ = spectrum(Qf, Mass)      # full spectrum
        # prime part as a matrix: Qp = Qa - Qf  (so full = arch - Qp; sign TBD — report both interpretations)
        Qp = mp.matrix(K, K)
        for i in range(K):
            for j in range(K): Qp[i, j] = Qf[i, j] - Qa[i, j]
        Ep, _ = spectrum(Qp, Mass)      # generalized eigenvalues of the prime part
        n_neg_arch = sum(1 for e in Ea if e < 0)
        n_nz_prime = sum(1 for e in Ep if abs(e) > mp.mpf(10) ** (-dps + 8))
        n_pos_prime = sum(1 for e in Ep if e > mp.mpf(10) ** (-dps + 8))
        print(f"\n=== a={mp.nstr(a,4)} ===")
        print(f"  arch-only spectrum (lowest 4): {[mp.nstr(e,5) for e in Ea[:4]]}")
        print(f"    #negative(arch) = {n_neg_arch}   λ_1={mp.nstr(Ea[0],6)}  λ_2={mp.nstr(Ea[1],6)}")
        print(f"  prime-part(Qf−Qa) spectrum (extremes): {[mp.nstr(e,5) for e in Ep[:2]]} ... {[mp.nstr(e,5) for e in Ep[-2:]]}")
        print(f"    rank(prime)≈{n_nz_prime}  (#pos={n_pos_prime}, #neg={n_nz_prime-n_pos_prime})  ⟹ rank-1/2 ⟺ single place")
        print(f"  full λ_min = {mp.nstr(Ef[0],6)}   (interlace: arch λ_1={mp.nstr(Ea[0],5)} ≤ full λ_min ≤ arch λ_2={mp.nstr(Ea[1],5)} ?)")
        ok = Ea[0] <= Ef[0] <= Ea[1] + mp.mpf(10)**(-6)
        print(f"    RANK-1-LIFT interlacing {'HOLDS' if ok else 'FAILS'}: {mp.nstr(Ea[0],4)} ≤ {mp.nstr(Ef[0],4)} ≤ {mp.nstr(Ea[1],4)}")
    print("\nREAD: if #negative(arch)=1 AND rank(prime)≤2 AND interlacing holds ⟹ p=2 positivity = a SCALAR secular")
    print("condition (rank-1 update lifting the single arch negative mode to ≥0). That is the analytically-attackable form.")
