#!/usr/bin/env python3
"""Cycle 10 (2026-08-26): characterize Delta_prime(a), the prime-form eigenvalue LIFT, in the screw frame.
From cycle 9: lambda_a = mu1 - log a + Delta_prime(a),  RH <=> Delta_prime(a) >= log a - mu1  forall a.
Build_Q pieces: Q_arch = pole+core+tail (arch-only), Qprime = von Mangoldt form.
(A) lambda_arch = eig0(Q_arch)/a -- confirm = mu1 - log a (cross-check cycle 9's separate build).
(B) FIRST-ORDER perturbation lift: pair = <c_arch|Qprime|c_arch>/<c_arch|c_arch>/a  (explicit prime sum
    against the KNOWN scale-invariant arch eigenfunction). RH <=> pair(a) [+ higher order] >= log a - mu1.
(C) actual Delta_prime = lambda_full - lambda_arch. Compare pair (1st order) vs actual vs (log a - mu1)."""
import numpy as np, suzuki_exact_Q as EQ
mu1=0.23676; K=80
print(f"{'a':>5} {'lam_arch':>9} {'mu1-log a':>10} {'pair(1st)':>10} {'Dprime(act)':>11} {'log a-mu1':>10}")
for a in [0.7,1.0,1.5,2.0,2.5,3.0,3.5,4.0]:
    Q,pc=EQ.build_Q(a,K,return_pieces=True)
    def sym(M): return (M+M.T)/2
    Qarch=sym(pc["pole"]+pc["core"]+pc["tail"]); Qprime=sym(pc["prime"]); Qfull=sym(Q)
    ea,Va=np.linalg.eigh(Qarch); c=Va[:,0]; lam_arch=ea[0]/a
    pair=float(c@Qprime@c)/float(c@c)/a
    lam_full=np.linalg.eigvalsh(Qfull)[0]/a
    Dact=lam_full-lam_arch
    print(f"{a:5.1f} {lam_arch:+9.4f} {mu1-np.log(a):+10.4f} {pair:+10.4f} {Dact:+11.4f} {np.log(a)-mu1:+10.4f}")
print("\nlam_arch should match mu1-log a (col2). If pair(1st order) tracks log a-mu1 (col6) => the prime")
print("lift IS the arch deficit at 1st order => explicit Guinand-Weil pairing; RH <=> pair>=log a-mu1.")
print("Gap (Dprime_actual - pair) = higher-order/level-repulsion correction.")
