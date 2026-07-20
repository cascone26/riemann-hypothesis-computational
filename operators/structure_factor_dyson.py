#!/usr/bin/env python3
"""
structure_factor_dyson.py -- Dyson quasicrystal test on the unfolded zeta zeros.

Dyson (2009): the zeros, suitably unfolded to unit mean spacing, should behave like
a 1D QUASICRYSTAL -> pure-point ("Bragg") diffraction. The structure factor
    S(k) = (1/N) | sum_n e^{i k x_n} |^2       (x_n = unfolded zeros, unit density)
should show SHARP PEAKS at special k, over a DIFFUSE background that -> 0 with N
(hyperuniformity, Zhang-Martelli-Torquato 2018). The explicit formula pins the peak
locations to k = log(p^m) (prime powers), in the UNUNFOLDED variable.

This script runs BOTH views on real Odlyzko zeros:
 (A) Raw-variable form factor F(u)=|sum e^{i u gamma_n}|^2/N -> peaks at u=log(p^m).
     (Weil explicit formula made visible; guaranteed, but a correctness check.)
 (B) Unfolded structure factor S(k): quantify Bragg-peak height vs diffuse floor,
     and how the diffuse floor scales with N (the genuinely-quantitative Dyson test
     that the literature left underspecified).
"""
import numpy as np
import json, time, sys

NZ = int(sys.argv[1]) if len(sys.argv) > 1 else 500000
ZFILE = "/Users/scones/projects/riemann/data/zeros_odlyzko_2M"

print(f"[sf] loading first {NZ} zeros ...", flush=True)
g = np.loadtxt(ZFILE, max_rows=NZ)
print(f"[sf] {len(g)} zeros, range {g[0]:.3f}..{g[-1]:.1f}", flush=True)

# ---------- (A) raw form factor: peaks at log(prime powers) ----------
def primes_up_to(n):
    sieve = np.ones(n+1, bool); sieve[:2] = False
    for i in range(2, int(n**0.5)+1):
        if sieve[i]: sieve[i*i::i] = False
    return np.where(sieve)[0]

pp = []
for p in primes_up_to(200):
    q = p
    while q <= 200:
        pp.append((np.log(q), np.log(p)/np.sqrt(q), p, q)); q *= p
pp.sort()

print("\n" + "="*70)
print("(A) RAW FORM FACTOR  F(u)=|sum e^{i u gamma}|/N  at u=log(p^m)")
print("    explicit-formula prediction: sharp peaks at these u, weight ~ log p / sqrt(p^m)")
print("="*70)
formA = []
for (u, wt, p, q) in pp[:15]:
    F = np.abs(np.sum(np.exp(1j * u * g))) / len(g)
    # compare to off-peak baseline at nearby non-log-prime u
    u_off = u + 0.05
    Foff = np.abs(np.sum(np.exp(1j * u_off * g))) / len(g)
    ratio = F / max(Foff, 1e-12)
    tag = f"p={p}" if q == p else f"{p}^{int(round(np.log(q)/np.log(p)))}={q}"
    print(f"  u=log({tag})={u:.4f}  |F|={F:.4e}  off-peak={Foff:.4e}  ratio={ratio:6.1f}")
    formA.append({"u": float(u), "prime_power": int(q), "F": float(F), "F_off": float(Foff), "ratio": float(ratio)})

# ---------- (B) unfolded structure factor + diffuse floor scaling ----------
# unfold: N(T) ~ (T/2pi)(log(T/2pi)-1); use running Riemann-von Mangoldt smooth count
def unfold(gam):
    T = gam / (2*np.pi)
    return T * (np.log(T) - 1) + 7.0/8.0   # smooth zero-counting function (unit density)

def diffuse_floor_and_peaks(gam, label):
    x = unfold(gam); N = len(x)
    # sample S(k) on a fine grid away from k=0; unit density -> first Bragg at k=2pi
    ks = np.linspace(0.05, 4*np.pi, 4000)
    # S(k) = |sum e^{i k x}|^2 / N   (compute in chunks to control memory)
    S = np.empty_like(ks)
    for i, k in enumerate(ks):
        S[i] = np.abs(np.sum(np.exp(1j * k * x)))**2 / N
    floor = np.median(S)                    # diffuse background level
    peak = np.max(S)                        # strongest Bragg-like peak
    kpeak = ks[np.argmax(S)]
    # fraction of grid above 10x floor = "how crystalline"
    frac_bragg = float(np.mean(S > 10*floor))
    print(f"  [{label}] N={N}: diffuse floor(median S)={floor:.3f}  max peak={peak:.1f} at k={kpeak:.3f}"
          f"  (k/2pi={kpeak/(2*np.pi):.3f})  frac>10xfloor={frac_bragg:.4f}")
    return {"N": N, "floor_median_S": float(floor), "max_peak": float(peak),
            "k_peak": float(kpeak), "k_peak_over_2pi": float(kpeak/(2*np.pi)),
            "frac_above_10x_floor": frac_bragg}

print("\n" + "="*70)
print("(B) UNFOLDED STRUCTURE FACTOR S(k): diffuse floor scaling with N")
print("    hyperuniform/quasicrystal => floor stays ~O(1) or shrinks, sharp Bragg peaks")
print("="*70)
formB = []
for frac in [0.1, 0.25, 0.5, 1.0]:
    m = int(len(g) * frac)
    formB.append(diffuse_floor_and_peaks(g[:m], f"first {m}"))

out = {"n_zeros": len(g), "raw_form_factor_A": formA, "unfolded_structure_B": formB}
p = "/Users/scones/projects/riemann/results/structure_factor_dyson.json"
json.dump(out, open(p, "w"), indent=2)
print(f"\nsaved -> {p}")
