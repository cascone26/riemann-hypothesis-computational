#!/usr/bin/env python3
"""§361 run: the k=3 NEXT-ORDER probe (criteria pre-registered in LEDGER §361, unchanged).
S3(H) to HMAX=8000, ~20 checkpoints; fit a H^2 + b H logH + c H; residual vs sqrt(H)."""
import numpy as np

HMAX = 8000
print(f"sieving spf to {HMAX} ...", flush=True)
spf = np.zeros(HMAX + 1, dtype=np.int64)
for p in range(2, int(HMAX ** 0.5) + 1):
    if spf[p] == 0:
        sl = spf[p * p::p]
        sl[sl == 0] = p
        spf[p * p::p] = sl

def odd_prime_factors(m):
    fs = []
    while m % 2 == 0:
        m //= 2
    while m > 1:
        p = spf[m] if spf[m] != 0 else m
        fs.append(int(p))
        while m % p == 0:
            m //= p
    return fs

sieve = np.ones(1_000_001, dtype=bool)
sieve[:2] = False
for p in range(2, 1001):
    if sieve[p]:
        sieve[p * p::p] = False
primes = np.nonzero(sieve)[0]
op5 = primes[primes >= 5].astype(np.float64)
C3 = float(np.prod((1.0 - 3.0 / op5) / (1.0 - 1.0 / op5) ** 3))
print(f"  C3 (p>=5) = {C3:.10f}", flush=True)

def S3(d1, d2):
    r1, r2 = d1 % 3, d2 % 3
    nu3 = len({0, r1, r2})
    if nu3 == 3:
        return 0.0
    val = 4.0 * ((1.0 - nu3 / 3.0) / (1.0 - 1.0 / 3.0) ** 3) * C3
    ps = set(odd_prime_factors(d1)) | set(odd_prime_factors(d2)) | set(odd_prime_factors(abs(d1 - d2)))
    ps.discard(3)
    for p in ps:
        nu = len({0, d1 % p, d2 % p})
        val *= ((1.0 - nu / p) / (1.0 - 1.0 / p) ** 3) / ((1.0 - 3.0 / p) / (1.0 - 1.0 / p) ** 3)
    return val

print("accumulating shells ...", flush=True)
shell = np.zeros(HMAX + 1)
evens = list(range(2, HMAX + 1, 2))
for d1 in evens:
    row = 0.0
    for d2 in evens:
        if d2 > d1:
            break
        if d2 == d1:
            continue
        v = S3(d1, d2)
        shell[d1] += 2.0 * v          # ordered pairs: (d1,d2) and (d2,d1), max = d1
    if d1 % 800 == 0:
        print(f"  d1={d1}", flush=True)
CUM = np.cumsum(shell)

Hs = [250, 300, 350, 450, 500, 600, 700, 850, 1000, 1200, 1400, 1700, 2000, 2400, 2800, 3400, 4000, 4800, 5600, 6800, 8000]
S3v = np.array([CUM[H] for H in Hs], dtype=float)
H = np.array(Hs, dtype=float)
A = np.column_stack([H * H, H * np.log(H), H])
coef, *_ = np.linalg.lstsq(A, S3v, rcond=None)
fitv = A @ coef
resid = S3v - fitv
print(f"\n3-param fit: a={coef[0]:.6f}  b={coef[1]:+.4f}  c={coef[2]:+.3f}")
print(f"{'H':>6} {'resid':>12} {'resid/sqrtH':>12}")
for i, Hi in enumerate(Hs):
    print(f"{Hi:6d} {resid[i]:+12.2f} {resid[i]/np.sqrt(Hi):+12.4f}")
jitter = np.diff(resid)
noise = float(np.sqrt(np.mean(jitter * jitter) / 2))
beta, *_ = np.linalg.lstsq(np.sqrt(H)[:, None], resid[:, None], rcond=None)
beta = float(beta[0][0])
pred = beta * np.sqrt(H)
upper = pred[len(Hs) // 2:]
sign_stable = np.all(np.sign(resid[len(Hs) // 2:]) == np.sign(beta)) and abs(beta * np.sqrt(H[-1])) > 2 * noise
print(f"\nnoise band (successive-jitter rms) = {noise:.2f}")
print(f"beta (sqrtH fit) = {beta:+.4f}; |beta*sqrt(Hmax)| = {abs(beta)*np.sqrt(H[-1]):.2f} vs 2*noise = {2*noise:.2f}")
if np.max(np.abs(resid)) < 2 * noise:
    print("VERDICT R1: no visible k=3 remainder at H<=8000 (bound, not disproof)")
elif sign_stable:
    print(f"VERDICT R2: M-S remainder measured, beta = {beta:+.4f}")
else:
    print("VERDICT R3: systematic but not sqrtH-shaped — flag, no story")
