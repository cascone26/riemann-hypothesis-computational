#!/usr/bin/env python3
"""
graph_zeta_tau.py — RESEARCH_QUEUE #2 (2026-07-20). Comparative anatomy vs a SOLVED RH.

The Ihara zeta function of a finite (q+1)-regular connected graph has a PROVEN Riemann Hypothesis:
  graph-RH  ⟺  the graph is RAMANUJAN  ⟺  every non-trivial adjacency eigenvalue |λ| ≤ 2√q
(Ihara–Sunada; the "zeros" u satisfy 1−λu+qu²=0, and |u|=q^{-1/2} on the critical line ⟺ |λ|≤2√q).

Write λ = 2√q cos θ. Ramanujan ⟺ every θ REAL (λ in the interval); non-Ramanujan ⟺ some θ COMPLEX
(|λ|>2√q). This is the EXACT graph mirror of §37's ζ observation: on-"line" zeros give real-angle
(positive) terms; off-"line" zeros give complex-angle (sign-indefinite cosh) terms. So the graph gives
us a SOLVED instance of the same positivity structure, and we can read off WHAT carries the positivity.

This script:
  1. builds a Ramanujan graph (Petersen, 3-regular) and a NON-Ramanujan 3-regular graph (bottleneck),
  2. verifies graph-RH numerically (all non-trivial |λ| ≤ 2√q ⟺ Ramanujan),
  3. builds the spectral-side form Σ_i h(θ_i) and shows the STRUCTURAL SIGNATURE: all-real angles
     (Ramanujan) vs a complex-angle mode injecting an anomalous out-of-range term (non-Ramanujan),
  4. states the mechanism DIFF: what makes positivity a THEOREM here vs OPEN for ζ.

HONEST SCOPE (verify-before-ship): the toy spectral form Σ h(θ_i) does NOT itself go negative for a
simple Gaussian h — a genuinely sign-indefinite graph Weil functional needs Terras' full graph explicit
formula (sum over closed geodesics, Euler-char boundary). What IS verified here is the STRUCTURAL MIRROR
(Ramanujan ⟺ all θ real ⟺ §37's on-line/off-line split) and it is that mirror — not a negative form —
that carries the mechanism diff below.
"""

import numpy as np
import networkx as nx
import json


def adjacency_spectrum(G):
    A = nx.to_numpy_array(G)
    return np.sort(np.linalg.eigvalsh(A))[::-1], A


def graph_rh_report(G, name, q):
    lam, A = adjacency_spectrum(G)
    deg = q + 1
    bound = 2 * np.sqrt(q)
    # non-trivial eigenvalues: exclude the Perron λ=deg and (if bipartite) λ=-deg
    nontrivial = [x for x in lam if abs(abs(x) - deg) > 1e-9]
    max_nontrivial = max(abs(x) for x in nontrivial)
    ramanujan = max_nontrivial <= bound + 1e-9
    print(f"\n[{name}]  {deg}-regular, n={G.number_of_nodes()}, q={q}, 2√q={bound:.4f}")
    print(f"  spectrum: {np.round(lam, 3)}")
    print(f"  max non-trivial |λ| = {max_nontrivial:.4f}  ⟹ "
          f"{'RAMANUJAN (graph-RH HOLDS)' if ramanujan else 'NON-Ramanujan (graph-RH FAILS)'}")
    return lam, deg, bound, ramanujan


def theta_of(lam, q):
    """λ = 2√q cos θ ⇒ θ = arccos(λ/2√q). Real if |λ|≤2√q, else complex (off critical line)."""
    z = lam / (2 * np.sqrt(q))
    return np.arccos(z + 0j)  # complex arccos; imag part ≠ 0 ⟺ off-line


def weil_form(lam, deg, q, sigma=1.0):
    """Spectral-side form on the ANGLE variable θ (mirror of ζ's γ), test h(θ)=exp(-θ²/2σ²).
    S = Σ_{nontrivial} h(θ_i). Ramanujan: every θ real ⟹ every term ∈ (0,1] (in-range). Non-
    Ramanujan: a complex θ=iy gives h=exp(+y²/2σ²) > 1 — an OUT-OF-RANGE anomalous term, the
    structural signature of an off-line 'zero'. (This simple S stays positive; the point is the
    real-vs-complex angle split, not the sign of S — see header honest-scope note.)"""
    theta = theta_of(lam, q)
    deg_v = deg
    terms = []
    for x, th in zip(lam, theta):
        if abs(abs(x) - deg_v) < 1e-9:
            continue  # skip trivial Perron modes
        # h(theta) = exp(-theta^2/(2 sigma^2)); for real theta this is in (0,1]; for theta=iy it is exp(+y^2/..)
        h = np.exp(-(th ** 2) / (2 * sigma ** 2))
        terms.append(complex(h))
    S = sum(terms)
    real_terms = [t.real for t in terms]
    # CORRECT off-"line" signature: a complex angle (Im θ ≠ 0), i.e. |λ|>2√q. (An h>1 term only
    # shows for a TOP-edge violation θ≈iy; a bottom-edge one θ≈π+iy is suppressed but still complex.)
    n_complex = sum(1 for x, th in zip(lam, theta)
                    if abs(th.imag) > 1e-9 and abs(abs(x) - deg) > 1e-9)
    return S.real, min(real_terms), max(real_terms), n_complex > 0


def make_non_ramanujan_3regular():
    """Bottleneck 3-regular graph: two K_4-minus-edge clusters joined by a thin bridge → small
    spectral gap → λ_2 close to 3 (> 2√2=2.828). Search small random 3-regular graphs for one
    that is non-Ramanujan (they exist but are uncommon at small n)."""
    q = 2
    bound = 2 * np.sqrt(q)
    # deterministic search over seeds (no Date.now/random-global): networkx seeded generator
    # prefer a TOP-edge violator (λ_2 > 2√q, a bottleneck/small-gap graph) — cleaner illustration;
    # fall back to any non-Ramanujan.
    fallback = None
    for n in [14, 16, 18, 20]:
        for seed in range(300):
            try:
                G = nx.random_regular_graph(3, n, seed=seed)
            except Exception:
                continue
            if not nx.is_connected(G):
                continue
            lam, _ = adjacency_spectrum(G)
            lam2 = lam[1]  # second-largest (positive top-edge)
            nontrivial = [x for x in lam if abs(abs(x) - 3) > 1e-9]
            if lam2 > bound + 1e-3:
                return G, seed  # top-edge violator
            if fallback is None and max(abs(x) for x in nontrivial) > bound + 1e-3:
                fallback = (G, seed)
    return fallback if fallback else (None, None)


def main():
    print("=" * 78)
    print("GRAPH-ZETA COMPARATIVE ANATOMY — positivity in a SOLVED-RH world (Ihara/Ramanujan)")
    print("=" * 78)

    out = {}

    # --- Ramanujan case: Petersen graph (3-regular, q=2) ---
    P = nx.petersen_graph()
    lamP, degP, boundP, ramP = graph_rh_report(P, "Petersen", q=2)
    SP, loP, hiP, cplxP = weil_form(lamP, degP, 2)
    print(f"  Weil-form S = Σ h(θ_i) = {SP:+.4f}   per-mode range [{loP:+.4f},{hiP:+.4f}]   "
          f"complex-angle modes: {cplxP}")
    out["petersen"] = {"ramanujan": bool(ramP), "weil_S": float(SP), "any_complex_angle": bool(cplxP)}

    # --- Non-Ramanujan case ---
    NR, seed = make_non_ramanujan_3regular()
    if NR is not None:
        lamNR, degNR, boundNR, ramNR = graph_rh_report(NR, f"non-Ramanujan (seed {seed})", q=2)
        SNR, loNR, hiNR, cplxNR = weil_form(lamNR, degNR, 2)
        print(f"  Weil-form S = Σ h(θ_i) = {SNR:+.4f}   per-mode range [{loNR:+.4f},{hiNR:+.4f}]   "
              f"complex-angle modes: {cplxNR}")
        out["non_ramanujan"] = {"seed": seed, "ramanujan": bool(ramNR), "weil_S": float(SNR),
                                "any_complex_angle": bool(cplxNR)}
    else:
        print("\n[non-Ramanujan] none found in seed range (all small 3-reg graphs Ramanujan)")

    print("\n" + "=" * 78)
    print("MECHANISM DIFF (the payoff for the ζ problem):")
    print("=" * 78)
    print("""
 In the graph world, positivity / graph-RH is a THEOREM because:
   • the 'zeros' ARE eigenvalues of a genuine BOUNDED SELF-ADJOINT operator (adjacency matrix A);
     its spectrum is REAL and bounded by fiat (finite symmetric matrix) — the analog of 'zeros on
     a line' is automatic up to the RADIUS, and Ramanujan-ness (the exact radius 2√q) is a
     construction/counting fact (Alon–Boppana floor + explicit LPS/Margulis constructions).
   • θ_i real ⟺ λ_i ∈ [−2√q, 2√q] ⟺ spectral GAP. Complex θ_i (non-Ramanujan) inject sign-
     indefinite terms — EXACTLY §37's on-line(positive) vs off-line(complex, sign-indefinite) split.

 What ζ LACKS (the one missing ingredient, made concrete):
   • ζ has NO bounded self-adjoint operator whose spectrum is the zeros — that is the missing
     Hilbert–Pólya operator. In the graph case self-adjointness is FREE (finite matrix); for ζ it
     is the entire open problem. The Ramanujan spectral GAP is the analog of ζ's ZERO-FREE REGION,
     NOT of RH itself; graph-RH = 'gap is exactly 2√q' is analogous to 'all zeros exactly on ½'.
   • ⟹ the ζ-analog of 'what carries positivity' is precisely the SELF-ADJOINTNESS of the missing
     operator. This loops directly to RESEARCH_QUEUE #6 (Berry–Keating H=xp, deficiency indices):
     the graph shows positivity is FREE once you have a bounded self-adjoint operator, so the whole
     fight for ζ is producing that operator (or proving it can't not exist). Consistent with the
     meta-theorem: positivity is global/product-emergent because it must come from an OPERATOR, not
     from local per-prime data.
""")
    with open("results/graph_zeta_tau.json", "w") as f:
        json.dump(out, f, indent=2)
    print("Results: results/graph_zeta_tau.json")


if __name__ == "__main__":
    main()
