#!/usr/bin/env python3
"""On-site moments of the a-field dL = Lambda(n)-1, deterministic sieve (no estimator).
Banked per Gate's ruling 2026-07-28: E[dL^2] (Gram diagonal scale, Hilbert's ask) and
m4/m2^2 (the mislabeled 'rho', kept for the record). build: rt6-afield-moments-v1"""
import json, numpy as np
NMAX = 4_000_000
lam = np.zeros(NMAX + 1)
comp = np.zeros(NMAX + 1, dtype=bool)
for p in range(2, NMAX + 1):
    if not comp[p]:
        comp[2*p::p] = True
        lp = np.log(p); pk = p
        while pk <= NMAX:
            lam[pk] = lp; pk *= p
rows = []
for N in (10**4, 10**5, 10**6, 2*10**6, 4*10**6):
    dL = lam[2:N+1] - 1.0
    m2 = float((dL**2).mean()); m4 = float((dL**4).mean())
    rows.append(dict(N=N, m2=m2, m4=m4, m4_over_m2sq=m4/m2**2, logN=float(np.log(N))))
    print(rows[-1])
json.dump(dict(build="rt6-afield-moments-v1-2026-07-28", field="dL = Lambda(n)-1",
               note="m2 = on-site variance E[dL^2] ~ logN (Gram diagonal scale); m4/m2^2 was mislabeled rho",
               rows=rows), open("results/rt6_afield_moments.json","w"), indent=1)
print("banked -> results/rt6_afield_moments.json")
