#!/usr/bin/env python3
"""ROSETTA, completed & CORRECTED — the actual positive-definite form proving RH-for-curves.
(First attempt used End(H^1)=M_2 with Tr on H^1 -> came out INDEFINITE [-1,-1,1,1]: WRONG. The Rosati
positivity is on the endomorphism ALGEBRA End^0(E) with the REDUCED trace, not all of M_2.)
CORRECT: for an ordinary elliptic curve End^0(E)=Q(alpha), alpha=Frobenius (root of x^2-a_p x+p),
Rosati involution = conjugation. The positive form is the NORM: N(m+n*alpha)=m^2+a_p*mn+p*n^2,
Gram [[1,a_p/2],[a_p/2,p]], POSITIVE-DEFINITE <=> det=p-a_p^2/4>0 <=> a_p^2<4p <=> HASSE <=> |alpha|=sqrt p <=> RH-for-E.
Hasse (|a_p|<=2 sqrt p) is the GEOMETRIC input (intersection positivity / ample cone). Number-field GAP:
no finite endomorphism algebra with a norm form => no automatic positivity => the a=1 razor, unproven.
"""
import numpy as np
def norm_gram(a_p, p): return np.array([[1.0, a_p/2],[a_p/2, float(p)]])
if __name__ == "__main__":
    print("Rosati norm form on End^0(E): PD <=> disc a_p^2-4p < 0 <=> Hasse <=> RH-for-E\n")
    print(f"{'curve':>10} {'a_p':>5} {'disc':>6} {'eigs':>20} {'PD=RH?':>7}")
    for a_p,p in [(0,5),(-1,7),(3,11),(-2,13),(2,17),(5,23)]:
        ev=np.linalg.eigvalsh(norm_gram(a_p,p))
        print(f"{'E/F'+str(p):>10} {a_p:>5} {a_p*a_p-4*p:>6} {str(np.round(ev,3)):>20} {str(np.all(ev>0)):>7}")
    print("\nCounterfactual (a_p^2>4p, RH-violating; Hasse forbids on real curves):")
    for a_p,p in [(7,11),(9,17)]:
        ev=np.linalg.eigvalsh(norm_gram(a_p,p))
        print(f"  a_p={a_p},p={p}: disc={a_p*a_p-4*p}>0 eigs={np.round(ev,3)} PD={np.all(ev>0)} (INDEFINITE=off-circle)")
    print("\nRH-for-E = 'the discriminant a_p^2-4p stays negative', guaranteed by Hasse (geometry).")
