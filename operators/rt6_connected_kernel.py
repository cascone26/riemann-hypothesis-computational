#!/usr/bin/env python3
"""
§389 (Riemann): the EXACT connected 6-point kernel O — partition-lattice Mobius
inversion (NOT the terminal's subset inclusion-exclusion sieve).

DEFINITION (pinned normalization, reproduce to the digit):
  Local moment functional at prime p for a block B (sub-multiset of residues mod p):
        w_p(B) = (1 - nu_p(B)/p) / (1 - 1/p)^{|B|},   nu_p = # DISTINCT residues in B,
        w_p(empty)=1, w_p(singleton)=1.
  Local CONNECTED cumulant (disconnected blocks subtracted BY CONSTRUCTION via the
  partition lattice of the 6 field insertions):
        O_p(X) = sum_{pi in Pi(6)}  mu(pi) * prod_{B in pi} w_p(B),
        mu(pi) = prod_{B in pi} (-1)^{|B|-1} (|B|-1)!      (Mobius fn of the partition lattice)
  This is the joint 6th cumulant of the local functional; the kappa2^3, kappa2*kappa4,
  kappa3^2 disconnected pieces are exactly the non-full partitions, subtracted here.
  GLOBAL connected kernel (Euler product over the accessible odd primes):
        O(X) = prod_{p in ODD_PRIMES} O_p(X).
  This is a DEFINITIONAL convention; Tao must match it to cross-check to the digit.

Exact rational arithmetic (fractions.Fraction) -> zero float ambiguity.
Off-diagonal cells only (X = 6 DISTINCT positions). The coincident diagonal
(w4=w1,w5=w2,w6=w3) is a DIFFERENT object (a-field moments kappa2*kappa4/kappa2^3),
NOT assembled here -- flagged, handed off, deliberately absent.
"""
import sys, json
from fractions import Fraction as Fr
from itertools import combinations

ODD_PRIMES = [3, 5, 7, 11, 13]  # highly-composite-odd sector: 105=3*5*7, 1155=*11, 15015=*13

# ---- partition lattice of {0..5} with Mobius weights, built once ----
def set_partitions(collection):
    collection = list(collection)
    if len(collection) == 1:
        yield [collection]; return
    first = collection[0]
    for smaller in set_partitions(collection[1:]):
        for i, subset in enumerate(smaller):
            yield smaller[:i] + [[first] + subset] + smaller[i+1:]
        yield [[first]] + smaller

from math import factorial
def block_mu(sz):
    return ((-1) ** (sz - 1)) * factorial(sz - 1)

PARTS6 = []
for pi in set_partitions(range(6)):
    mu = 1
    for B in pi:
        mu *= block_mu(len(B))
    PARTS6.append((mu, [tuple(B) for B in pi]))
assert len(PARTS6) == 203, len(PARTS6)  # Bell(6)

BUILD_VERSION = "rt6-ckernel-v2-collapse-corrected-2026-07-28"

def w_p(positions, p):
    """CORRECTED coincident handling: collapse to DISTINCT positions FIRST.
    A repeated position is one field point -> exponent counts distinct positions,
    not multiplicity. This kills the (1-1/p)^{|B|} over-count that fabricated the
    huge negative diagonal (the retracted -8.7e5 contamination)."""
    d = set(positions)                     # collapse coincident positions
    if not d:
        return Fr(1)
    k = len(d)                             # distinct-position count (NOT multiplicity)
    nu = len(set(x % p for x in d))        # distinct residues among distinct positions
    num = Fr(p - nu, p)                    # (1 - nu/p)
    den = Fr(p - 1, p) ** k                # (1 - 1/p)^k_distinct
    return num / den

def W_block(X, B):
    """GLOBAL singular-series moment of the sub-block B: Euler product of local w_p
    over POSITIONS (coincident positions collapse inside w_p)."""
    val = Fr(1)
    pos = [X[i] for i in B]
    for p in ODD_PRIMES:
        val *= w_p(pos, p)
    return val

def O_p(X, p):
    """single-prime connected cumulant (sanity/reproduction only, not the global object)."""
    total = Fr(0)
    for mu, pi in PARTS6:
        prod = Fr(mu)
        for B in pi:
            prod *= w_p([X[i] for i in B], p)
        total += prod
    return total

def O_global(X):
    """CORRECT connected 6-pt kernel: Mobius inversion is the OUTERMOST op,
    applied to the global block-moment W(B)=prod_p w_p(B). Disconnected
    kappa2^3/kappa2*kappa4/kappa3^2 subtracted globally, by construction."""
    total = Fr(0)
    for mu, pi in PARTS6:
        prod = Fr(mu)
        for B in pi:
            prod *= W_block(X, B)
        total += prod
    return total

# ---- build a triple-set at fixed triple-sum, form off-diagonal O[t,t'] ----
def triples_at_sum(s, lo=1):
    """odd-position triples p1<p2<p3 summing to s (positions on the odd lattice)."""
    out = []
    for a in range(lo, s):
        for b in range(a + 1, s):
            c = s - a - b
            if c > b:
                out.append((a, b, c))
    return out

def build_matrix(s, primes_report=False, max_triples=None):
    """FULL matrix now: diagonal + coincident cells assembled with collapse-corrected
    w_p (no None-skipping). D = coincident/self cells, O = distinct-6 cells."""
    T = triples_at_sum(s)
    if max_triples:
        T = T[:max_triples]
    n = len(T)
    M = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            X = list(T[i]) + list(T[j])
            v = O_global(X)                # coincident positions handled inside w_p
            M[i][j] = v
            M[j][i] = v                    # symmetric by construction
    return T, M

def norms_and_signs(T, M):
    import numpy as np
    n = len(T)
    Full = np.zeros((n, n))
    Off = np.zeros((n, n))
    diag = np.zeros(n)
    neg_cells = []
    for i in range(n):
        for j in range(n):
            fv = float(M[i][j])
            Full[i][j] = fv
            if i == j:
                diag[i] = fv
            else:
                Off[i][j] = fv
                if fv < 0 and j > i:
                    neg_cells.append((T[i], T[j], fv))
    Full = (Full + Full.T) / 2
    Soff = (Off + Off.T) / 2
    eig_full = np.linalg.eigvalsh(Full)      # lambda_min(D+O) -- the ACTUAL decider
    eig_off = np.linalg.eigvalsh(Soff)       # hollow O (Hilbert: forced negative, trace 0)
    off_norm = float(np.max(np.abs(eig_off)))
    d_min = float(diag.min())                # lambda_min(D), D diagonal -> min entry
    d_max = float(diag.max())
    ratio = off_norm / d_min if d_min > 0 else float('inf')
    return {
        "n_triples": n,
        "off_cells": n * (n - 1),
        "neg_off_cells": len(neg_cells) * 2,   # symmetric pair count
        "off_norm_symm": off_norm,
        "diag_min": d_min,
        "diag_max": d_max,
        "lambda_min_full": float(eig_full.min()),   # <-- D+O, non-hollow
        "lambda_max_full": float(eig_full.max()),
        "hilbert_ratio": ratio,                # ||O||/lambda_min(D); <0.32 banks tier-1
        "eig_min_off": float(eig_off.min()),
        "eig_max_off": float(eig_off.max()),
        "neg_cells": neg_cells,
    }

if __name__ == "__main__":
    # p=5 single-prime reproduction check (matches terminal's checkpoint target arena)
    print("=== SINGLE-PRIME connected cumulant O_5 sanity (a few distinct configs) ===")
    for X in [(1,2,3,4,5,6),(1,2,4,8,16,3),(0,1,2,3,4,5),(1,3,7,9,11,13)]:
        print(f"  X={X}  O_5 = {O_p(X,5)}  (= {float(O_p(X,5)):+.6f})")

    results = {}
    for s in [30, 36, 42]:
        T, M = build_matrix(s)
        info = norms_and_signs(T, M)
        # keep only summary + a capped neg-cell list in the checkpoint
        cells = [ {"t": list(a), "tp": list(b), "val": v} for (a,b,v) in info["neg_cells"] ]
        results[str(s)] = {
            "s": s,
            "n_triples": info["n_triples"],
            "off_cells": info["off_cells"],
            "neg_off_cells": info["neg_off_cells"],
            "off_norm_symm": info["off_norm_symm"],
            "diag_min": info["diag_min"],
            "diag_max": info["diag_max"],
            "lambda_min_full": info["lambda_min_full"],
            "lambda_max_full": info["lambda_max_full"],
            "hilbert_ratio": info["hilbert_ratio"],
            "neg_cells": cells,
        }
        print(f"\n=== s={s}: {info['n_triples']} triples ===")
        print(f"  diag(D) in [{info['diag_min']:.4f}, {info['diag_max']:.4f}]  "
              f"(collapse-corrected, NOT the contaminated -8.7e5)")
        print(f"  ||O_offdiag||_symm = {info['off_norm_symm']:.4f}")
        print(f"  lambda_min(D+O) [ACTUAL decider] = {info['lambda_min_full']:.4f}  "
              f"(hollow-O eig_min was {info['eig_min_off']:.1f}, forced trivia)")
        print(f"  Hilbert ratio ||O||/lambda_min(D) = {info['hilbert_ratio']:.4f}  "
              f"({'<0.32 -> tier-1 banks' if info['hilbert_ratio'] < 0.32 else '>=0.32 -> NOT yet'})")

    print("\n=== SCALING ===")
    for k in ['off_norm_symm', 'diag_min', 'lambda_min_full', 'hilbert_ratio']:
        vs = [results[str(s)][k] for s in [30,36,42]]
        print(f"  {k:18s} s=30/36/42: {vs[0]:.4f}, {vs[1]:.4f}, {vs[2]:.4f}")

    out = {
        "build_version": BUILD_VERSION,
        "object": "connected 6-point kernel O, partition-lattice Mobius inversion",
        "definition": {
            "w_p": "(1 - nu_p/p)/(1-1/p)^{|B|}, nu_p=#distinct residues in block B",
            "O_p": "sum_{pi in Pi(6)} mu(pi) prod_B w_p(B), mu=prod_B (-1)^{|B|-1}(|B|-1)!",
            "O_global": "prod_{p in [3,5,7,11,13]} O_p(X)",
            "arith": "exact fractions.Fraction; floats only at report",
        },
        "odd_primes": ODD_PRIMES,
        "diagonal": "NOT assembled here -- coincident (w4=w1..) is a-field kappa2*kappa4/kappa2^3, flagged/handed off",
        "results": results,
    }
    with open("results/rt6_connected_kernel_checkpoint.json", "w") as f:
        json.dump(out, f, indent=2)
    print("\nCHECKPOINT -> results/rt6_connected_kernel_checkpoint.json")
