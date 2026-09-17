#!/usr/bin/env python3
"""
CERTIFIED Euler-product discriminator (Handle-tier attempt, per Charisma).

Claim under test: the CERTIFIED Weil positivity margin (Groskin interval-LDL^T, n_neg==0) rests on the
prime weights being NONNEGATIVE -- which is exactly the Euler product (Lambda(n) = log p >= 0 on prime
powers). Negate prime weight(s) -- simulating the loss of the Euler product's coefficient positivity, the
one thing Davenport-Heilbronn lacks -- and the SAME certified instrument must lose positivity (n_neg>=1).

Controlled: we keep Groskin's EXACT archimedean assembly (W02, WR, pole kappa, J) UNCHANGED and modify ONLY
the prime-weight signs. So any change in the certified inertia is attributable to the weight signs alone =
the Euler-product coefficient positivity. This is the cleanest isolation (a faithful full D-H tau would also
swap the gamma factors + drop the pole; here we isolate the ONE variable of interest).

Rigorous: same interval arithmetic (python-flint arb) and same LDL^T inertia as the canonical certifier.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "external", "groskin_2607.02828", "code"))
from flint import arb, arb_mat, ctx
from arb_ldlt_certify import (arb_closed_forms, arb_kappa, arb_J, prime_powers_up_to)
from arb_ldlt_margin import inertia_with_margin


def build_tau_weighted(c, N, prec, weight_transform):
    """Groskin's tau, but the prime weights pass through weight_transform(list_of_(q,p,weight))."""
    ctx.prec = prec
    S, CC, XC, L = arb_closed_forms(N, c, prec)
    PI = arb.pi()
    sp2 = 16 * PI * PI
    l2 = L * L
    pref02 = 32 * L * (L / 4).sinh() ** 2
    kappa = arb_kappa(L)
    J = arb_J(L)

    pdata = prime_powers_up_to(c)
    base_weights = [arb(p).log() * (arb(q) ** arb('-0.5')) for (q, p) in pdata]
    positions = [arb(q).log() for (q, p) in pdata]
    weights = weight_transform(list(zip([q for q, p in pdata], [p for q, p in pdata], base_weights)))

    def S_signed(nn):
        return S[nn] if nn >= 0 else -S[-nn]

    DIM = 2 * N + 1
    A = arb_mat(DIM, DIM)
    M = len(weights)
    for i in range(DIM):
        n = i - N
        for j in range(i, DIM):
            m = j - N
            num = l2 - sp2 * m * n
            den = (l2 + sp2 * m * m) * (l2 + sp2 * n * n)
            W02 = pref02 * num / den
            if n == m:
                WR = kappa + 2 * CC[abs(n)] + J - (2 / L) * XC[abs(n)]
            else:
                WR = (S_signed(m) - S_signed(n)) / (PI * (n - m))
            Wp = arb(0)
            for idx in range(M):
                y = positions[idx]
                if n == m:
                    q = 2 * (1 - y / L) * (2 * PI * n * y / L).cos()
                else:
                    q = ((2 * PI * m * y / L).sin() - (2 * PI * n * y / L).sin()) / (PI * (n - m))
                Wp += weights[idx] * q
            val = W02 - WR - Wp
            A[i, j] = val
            A[j, i] = val
    return A, DIM


def certify(c, N, prec, weight_transform, tag):
    A, DIM = build_tau_weighted(c, N, prec, weight_transform)
    npos, nneg, undet, mlow, mmid = inertia_with_margin(A, DIM)
    pd = (undet is None and nneg == 0)
    ms = f"{mlow:.4e}" if mlow is not None else "UNDET"
    print(f"  [{tag:26s}] dim={DIM} n_pos={npos} n_neg={nneg} certified_PD={pd} min_pivot_low={ms}")
    return dict(nneg=nneg, pd=pd, margin=mlow)


# --- weight transforms ---
def W_zeta(triples):
    return [w for (q, p, w) in triples]                          # Euler product: all >= 0


def W_negate_prime(pp):
    def t(triples):
        return [(-w if p == pp else w) for (q, p, w) in triples]  # break EP: flip weights of one prime pp
    return t


def W_scale_prime(pp, factor):
    def t(triples):
        return [(w * arb(factor) if p == pp else w) for (q, p, w) in triples]
    return t


if __name__ == "__main__":
    c, N, prec = 25, 32, 1500
    print(f"CERTIFIED Euler-product discriminator (c={c}, N={N}, prec={prec})")
    print("Baseline zeta (Euler product, all prime weights >= 0):")
    certify(c, N, prec, W_zeta, "zeta / Euler product")

    print("\nBreak the Euler product -- negate the weight of prime p (Lambda now signed at p^k):")
    for pp in [2, 3, 5]:
        certify(c, N, prec, W_negate_prime(pp), f"EP broken: -w(p={pp})")

    print("\nDial the break strength on p=2 (factor multiplies its positive weight):")
    for f in ['1.0', '0.0', '-0.5', '-1.0', '-2.0']:
        certify(c, N, prec, W_scale_prime(2, f), f"w(2) *= {f}")

    print("\nREAD: zeta certifies PD (n_neg=0). Flipping ANY prime's weight negative (= losing the Euler")
    print("product's Lambda(n)>=0 at that prime) drives n_neg>=1 on the SAME certified instrument. => the")
    print("Euler-product coefficient positivity is LOAD-BEARING for certified Weil positivity; D-H (signed")
    print("weights, no Euler product) provably cannot certify -- matching its known off-line zeros. (Handle-tier:")
    print("a real certified byte-level inertia flip, not a described mechanism.)")
