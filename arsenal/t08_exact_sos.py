#!/usr/bin/env python3
"""
t08_exact_sos.py  -- Tier B #8.  Exact-rational positive-definiteness / SOS certifier.  Where the
repo's arb_certify path uses interval arithmetic on transcendental entries, this is the exact-
arithmetic companion: given a symmetric matrix with RATIONAL entries, compute the LDL^T
decomposition over Q with sympy and certify PD iff every pivot is a positive rational -- an EXACT
sum-of-squares certificate M = sum d_i (l_i . x)^2, no floating point anywhere.

HANDLE-tier for the finite rational matrix it certifies (proof-pointer = the exact pivots on disk).
Honest scope: it certifies THAT rational matrix, exactly.  It is NOT a proof of RH; the true Weil
tau has transcendental entries (digamma, log), so an RH-relevant exact cert needs a rational
enclosure whose PD implies the real block's PD -- that bridge is the Arb job, noted, not claimed.

Grounding demos:
  - a manifestly-PD rational Gram matrix  -> exact PD cert (HANDLE)
  - the Hankel/moment 2x2 [[l1,l2],[l2,l3]] on Li coefficients -> exact det < 0, reproducing the
    project's 'moment route DEAD' result EXACTLY (not numerically).
  - an indefinite matrix -> correctly reported NOT PD.

Usage: python3 t08_exact_sos.py
"""
from __future__ import annotations
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
from core.honesty import Result, Tier, assert_handle
from core import provenance


def exact_ldlt(M: sp.Matrix):
    """
    Exact LDL^T over Q. Returns (is_pd, pivots, L, D) with pivots exact rationals.
    is_pd iff all pivots > 0 (M symmetric).  Uses leading-principal-minor / Cholesky-style recursion.
    """
    n = M.rows
    assert M == M.T, "matrix must be symmetric"
    L = sp.eye(n)
    D = [sp.Integer(0)] * n
    for j in range(n):
        d = M[j, j] - sum(L[j, k] ** 2 * D[k] for k in range(j))
        d = sp.nsimplify(sp.together(d))
        D[j] = sp.simplify(d)
        if D[j] == 0:
            # singular; fall back to minors verdict
            return None, D, L, None
        for i in range(j + 1, n):
            s = M[i, j] - sum(L[i, k] * L[j, k] * D[k] for k in range(j))
            L[i, j] = sp.simplify(s / D[j])
    is_pd = all(sp.Rational(dv) > 0 for dv in D)
    return is_pd, [sp.Rational(dv) for dv in D], L, D


def certify(M: sp.Matrix, label: str):
    # exact symmetric check
    if M != M.T:
        return {"label": label, "verdict": "NOT_SYMMETRIC"}
    det = sp.Rational(sp.simplify(M.det()))
    is_pd, pivots, L, D = exact_ldlt(M)
    if is_pd is None:
        return {"label": label, "verdict": "SINGULAR", "det": str(det),
                "pivots": [str(p) for p in pivots]}
    # leading principal minors (Sylvester) as an independent exact cross-check
    minors = [sp.Rational(sp.simplify(M[:k, :k].det())) for k in range(1, M.rows + 1)]
    minors_pd = all(m > 0 for m in minors)
    return {
        "label": label,
        "verdict": "PD" if (is_pd and minors_pd) else "NOT_PD",
        "matrix": [[str(M[i, j]) for j in range(M.cols)] for i in range(M.rows)],
        "det": str(det),
        "pivots": [str(p) for p in pivots],
        "pivots_all_positive": bool(is_pd),
        "leading_minors": [str(m) for m in minors],
        "sylvester_pd": bool(minors_pd),
        "ldlt_vs_sylvester_agree": bool(is_pd == minors_pd),
    }


# ---- demo matrices ----
def gram_pd():
    V = sp.Matrix([[2, 1, 0], [1, 3, 1], [0, 1, 2]])  # rational, symmetric, diagonally dominant
    return V


def moment_hankel():
    # Li/Keiper coefficients lambda_1..3 (rationalized) -> the Track-1 Hankel matrix.
    # lambda_1=0.0230957, lambda_2=0.0923457, lambda_3=0.2076466 (standard values)
    l1 = sp.Rational(230957, 10000000)
    l2 = sp.Rational(923457, 10000000)
    l3 = sp.Rational(2076466, 10000000)
    return sp.Matrix([[l1, l2], [l2, l3]])


def indefinite():
    return sp.Matrix([[1, 2], [2, 1]])  # det = -3 < 0


def selftest(verbose=True) -> bool:
    ok = True
    a = certify(gram_pd(), "rational Gram")
    if a["verdict"] != "PD" or not a["ldlt_vs_sylvester_agree"]:
        ok = False; print("  FAIL: PD Gram not certified PD:", a)
    b = certify(moment_hankel(), "Hankel moment 2x2 (Li coeffs)")
    # must be NOT_PD with negative det (reproduces moment route DEAD)
    if b["verdict"] != "NOT_PD" or sp.Rational(b["det"]) >= 0:
        ok = False; print("  FAIL: moment Hankel not reproduced as indefinite:", b)
    c = certify(indefinite(), "indefinite [[1,2],[2,1]]")
    if c["verdict"] != "NOT_PD":
        ok = False; print("  FAIL: indefinite not caught")
    if verbose:
        print(f"  rational Gram   -> {a['verdict']}, pivots {a['pivots']}, det {a['det']}")
        print(f"  moment Hankel   -> {b['verdict']}, det {b['det']} (project: DEAD, det<0)")
        print(f"  indefinite      -> {c['verdict']}, det {c['det']}")
    print("t08.selftest:", "PASS" if ok else "FAIL")
    return ok


def main():
    print("=== t08 Exact-rational SOS / PD certifier (over Q, no floats) ===\n")
    results = []
    cases = [(gram_pd(), "rational Gram (3x3)"),
             (moment_hankel(), "Hankel moment 2x2 on Li coefficients"),
             (indefinite(), "indefinite [[1,2],[2,1]]")]
    for M, label in cases:
        c = certify(M, label)
        print(f"[{c['verdict']:7s}] {label}")
        print(f"          det = {c['det']}   pivots = {c.get('pivots')}")
        print(f"          Sylvester minors = {c.get('leading_minors')}  "
              f"(LDL^T vs Sylvester agree: {c.get('ldlt_vs_sylvester_agree')})")
        if c["verdict"] == "PD":
            ev = os.path.join(provenance.RUNS_DIR, f"t08_{label.split()[0]}_exactcert.json")
            with open(ev, "w") as f:
                json.dump(c, f, indent=2)
            r = assert_handle(True,
                claim=f"EXACT PD certificate over Q for finite matrix '{label}': all LDL^T pivots "
                      f"positive rationals, Sylvester-cross-checked (finite matrix, NOT RH)",
                proof_pointer=ev, method="exact rational LDL^T + Sylvester minors",
                tool="t08_exact_sos", data=c)
        else:
            r = Result(claim=f"'{label}' is NOT positive-definite over Q (det={c['det']})",
                       tier=Tier.REPORT, method="exact rational LDL^T", tool="t08_exact_sos", data=c)
        results.append(r)
        print("          ", r.line())
        print()
    print("Note: exact certs here are for RATIONAL matrices. The real Weil tau has transcendental "
          "entries; an RH-relevant exact cert needs a rational enclosure (Arb interval) whose PD "
          "implies the real block's -- that bridge is arb_certify's job, not claimed here.")
    ok = selftest()
    path = provenance.stamp("t08_exact_sos", results, extra={"selftest_pass": ok})
    print("run record:", path)
    return ok


if __name__ == "__main__":
    main()
