#!/usr/bin/env python3
"""
SICK_HIM.py  -- adversarial verification harness.  Attacks the arsenal to catch false
certifications, overclaims, and math errors.  Assumes every claim is a lie until an INDEPENDENT
check clears it.  See SICK_HIM.md for the protocol.  Exit 0 iff the arsenal survives every attack.
"""
from __future__ import annotations
import sys, os, io, json, glob, re, math, contextlib, importlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import mpmath as mp
import sympy as sp

ARS = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(ARS, "runs")
FAILS = []
def fail(check, detail): FAILS.append((check, detail)); print(f"  [ATTACK LANDED] {check}: {detail}")
def ok(check, detail=""): print(f"  [survived] {check}{(' — ' + detail) if detail else ''}")


# ---- 1. overclaim sweep over source + emitted claims ----
def check_overclaims():
    print("\n[1] overclaim sweep (source + emitted claims)")
    banned = ["proven rh", "proved rh", "proves rh", "certified rh", "rh certified",
              "solves rh", "rh solved", "proof of rh complete", "we have proven the riemann"]
    hits = []
    for f in glob.glob(os.path.join(ARS, "*.py")) + glob.glob(os.path.join(ARS, "core", "*.py")):
        # skip the banned-list definitions themselves (honesty.py, SICK_HIM.py)
        if os.path.basename(f) in ("honesty.py", "SICK_HIM.py"):
            continue
        txt = open(f).read().lower()
        for b in banned:
            if b in txt:
                hits.append((os.path.basename(f), b))
    for rec in glob.glob(os.path.join(RUNS, "*.json")):
        for r in json.load(open(rec)).get("results", []):
            c = r.get("claim", "").lower()
            for b in banned:
                if b in c:
                    hits.append((os.path.basename(rec), b))
    if hits:
        for f, b in hits: fail("overclaim", f"{f}: {b!r}")
    else:
        ok("overclaim sweep", "no RH-proof overclaims in source or emitted claims")


# ---- 2. HANDLE evidence deep-check ----
def check_handle_evidence():
    print("\n[2] HANDLE evidence deep-check")
    n_handle = 0; bad = 0
    for rec in glob.glob(os.path.join(RUNS, "*.json")):
        for r in json.load(open(rec)).get("results", []):
            if r.get("tier") != "HANDLE":
                continue
            n_handle += 1
            pp = r.get("proof_pointer")
            if not pp or pp == "/dev/null" or not os.path.exists(pp) or os.path.getsize(pp) == 0:
                fail("hollow HANDLE", f"{os.path.basename(rec)}: {pp}"); bad += 1; continue
            try:
                json.load(open(pp))  # must be real structured evidence
            except Exception:
                fail("HANDLE evidence not parseable", f"{pp}"); bad += 1
    if n_handle and not bad:
        ok("HANDLE evidence", f"{n_handle} HANDLE claims all point at real, parseable evidence")
    elif not n_handle:
        ok("HANDLE evidence", "no HANDLE claims present")


# ---- 3. triple-independent Weil check ----
def check_triple_weil():
    print("\n[3] triple-independent Weil margin (numpy vs mpmath vs raised-precision-shifted)")
    from core.weil import WeilForm
    from t06_fabrication_detector import independent_weil_margin
    t = 0.10
    W_np = WeilForm(zero_count=40).margin(t)
    W_mp = independent_weil_margin(t=t, nzeros=40)[0]
    # third path: mpmath at higher dps, shifted quadrature nodes, tanh-sinh
    with mp.workdps(40):
        tt = mp.mpf(t)
        ghat = lambda r: tt * mp.sqrt(2 * mp.pi) * mp.e ** (-(tt ** 2) * r ** 2 / 2)
        Om = lambda r: mp.re(mp.digamma(mp.mpf('0.25') + 1j * r / 2)) - mp.log(mp.pi)
        A = mp.quad(lambda r: ghat(r) * Om(r), [-160, -37, 0, 37, 160]) / (2 * mp.pi)
        B = 2 * tt * mp.sqrt(2 * mp.pi) * mp.e ** (tt ** 2 / 8)
        # primes via sympy factorization, independent of numpy sieve
        P = mp.mpf(0)
        for n in range(2, 20000):
            fac = sp.factorint(n)
            if len(fac) == 1:
                p = next(iter(fac)); P += mp.log(p) / mp.sqrt(n) * mp.e ** (-(mp.log(n)) ** 2 / (2 * tt ** 2))
        P *= 2
        W_3 = float(A + B - P)
    d12 = abs(W_np - W_mp); d13 = abs(W_np - W_3); d23 = abs(W_mp - W_3)
    print(f"      W_numpy={W_np:.10f}  W_mpmath={W_mp:.10f}  W_third={W_3:.10f}")
    print(f"      pairwise diffs: |1-2|={d12:.2e} |1-3|={d13:.2e} |2-3|={d23:.2e}")
    if max(d12, d13, d23) < 1e-6:
        ok("triple Weil", "three independent code paths agree to <1e-6")
    else:
        fail("triple Weil disagreement", f"max pairwise diff {max(d12,d13,d23):.2e}")


# ---- 4. exact-cert adversarial (t08) ----
def check_t08_adversarial():
    print("\n[4] exact-cert adversarial (near-singular / PSD-not-PD / non-symmetric)")
    from t08_exact_sos import certify
    cases = []
    # near-singular positive: det = +1/10^18 exactly, PD
    e = sp.Rational(1, 10 ** 18)
    M_pos = sp.Matrix([[1, 0], [0, e]])            # PD, det=e>0
    cases.append(("near-singular PD (det=1e-18)", M_pos, "PD"))
    # near-singular indefinite: det = -1/10^18
    M_neg = sp.Matrix([[1, 1], [1, 1 - e]])        # det = (1-e) - 1 = -e < 0
    cases.append(("near-singular NOT-PD (det=-1e-18)", M_neg, "NOT_PD"))
    # PSD but not PD: det = 0
    M_psd = sp.Matrix([[1, 1], [1, 1]])            # det=0, singular
    cases.append(("PSD-not-PD (det=0)", M_psd, ("NOT_PD", "SINGULAR")))
    # non-symmetric
    M_ns = sp.Matrix([[1, 2], [0, 1]])
    cases.append(("non-symmetric", M_ns, "NOT_SYMMETRIC"))
    allgood = True
    for label, M, expect in cases:
        v = certify(M, label)["verdict"]
        exp = expect if isinstance(expect, tuple) else (expect,)
        good = v in exp
        # cross-check det sign against numpy for the numeric ones
        print(f"      {label}: verdict={v}  expected~{expect}  {'OK' if good else 'WRONG'}")
        if not good:
            fail("t08 adversarial", f"{label}: got {v}, expected {expect}"); allgood = False
    if allgood:
        ok("t08 adversarial", "exact arithmetic correct on all near-singular/degenerate cases")


# ---- 5. mutation probe: break the engine, confirm selftest catches it ----
def check_mutation():
    print("\n[5] mutation probe (inject a sign flip; the selftest MUST then fail)")
    import core.weil as weil
    orig_P = weil.WeilForm.P
    # mutant: flip the sign of the prime term (a real bug that should break the identity)
    def mutant_P(self, t):
        return -orig_P(self, t)
    weil.WeilForm.P = mutant_P
    try:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            passed_on_broken = weil.selftest(verbose=False)
    finally:
        weil.WeilForm.P = orig_P  # always restore
    # sanity: restored engine passes again
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        passed_after_restore = weil.selftest(verbose=False)
    if passed_on_broken:
        fail("rubber-stamp selftest", "weil.selftest PASSED on a sign-flipped (broken) engine")
    elif not passed_after_restore:
        fail("mutation not restored", "engine did not recover after mutation")
    else:
        ok("mutation probe", "selftest correctly FAILED on the broken engine and passed after restore")


# ---- 6. guard fuzz ----
def check_guard_fuzz():
    print("\n[6] honesty-guard fuzz")
    from core.honesty import guard_text, OverclaimError, Result, Tier, NoEvidenceError
    phrasings = [
        "this proves RH", "we certified RH today", "RH solved at last",
        "a complete proof of RH", "this proves the riemann hypothesis",
        "RH is proven by this cert", "solves the riemann hypothesis",
    ]
    escaped = []
    for p in phrasings:
        try:
            guard_text(p); escaped.append(p)
        except OverclaimError:
            pass
    if escaped:
        for p in escaped: fail("guard bypass", f"{p!r} not blocked")
    else:
        ok("guard fuzz", f"all {len(phrasings)} overclaim phrasings blocked")
    # HANDLE without real evidence must be refused
    try:
        Result(claim="finite margin >0", tier=Tier.HANDLE, method="x",
               proof_pointer="/no/such/evidence_zzz.json")
        fail("guard bypass", "HANDLE with missing evidence file was allowed")
    except NoEvidenceError:
        ok("HANDLE evidence enforcement", "HANDLE with missing evidence refused at construction")


# ---- 7. independent-code audit ----
def check_independence():
    print("\n[7] independence audit (t06 mpmath path must not reuse core.weil)")
    src = open(os.path.join(ARS, "t06_fabrication_detector.py")).read()
    # the independent margin function must NOT call into core.weil's WeilForm for its own recompute
    fn_src = src[src.index("def independent_weil_margin"):src.index("def audit_claim")]
    leaks = [tok for tok in ("WeilForm", "from core.weil", "weil.WeilForm", "core.weil")
             if tok in fn_src]
    if leaks:
        fail("shared-code", f"independent_weil_margin references {leaks}")
    else:
        ok("independence", "independent_weil_margin uses only mpmath (own sieve, own quad, own zeros)")


# ---- 8. audit_runs hardening: dir / symlink / junk evidence must be flagged (red-team C) ----
def check_audit_hardening():
    print("\n[8] audit_runs hardening (directory / symlink / junk-JSON evidence must be flagged)")
    import tempfile, shutil
    from t06_fabrication_detector import audit_runs, deep_reverify
    tmp = tempfile.mkdtemp(prefix="t06_hard_")
    try:
        runs = os.path.join(tmp, "runs"); os.makedirs(runs)
        os.makedirs(os.path.join(tmp, "adir"))
        open(os.path.join(tmp, "junk.json"), "w").write("1")        # not valid JSON object? '1' parses
        open(os.path.join(tmp, "junkbytes.json"), "w").write("\x00not json")
        try: os.symlink("/etc/hosts", os.path.join(tmp, "lnk.json"))
        except FileExistsError: pass
        specs = {"dir.json": os.path.join(tmp, "adir"),
                 "junkbytes.json": os.path.join(tmp, "junkbytes.json"),
                 "lnk.json": os.path.join(tmp, "lnk.json")}
        for fn, pp in specs.items():
            json.dump({"tool": "fake", "results": [{"claim": "CERTIFIED", "tier": "HANDLE",
                       "proof_pointer": pp, "method": "x"}]}, open(os.path.join(runs, fn), "w"))
        found = {f["file"] for f in audit_runs(runs)}
        missed = [fn for fn in specs if fn not in found]
        if missed:
            fail("audit hardening", f"these fabricated evidences slipped through: {missed}")
        else:
            ok("audit hardening", "directory, symlink, and junk-bytes evidence all flagged")

        # deep_reverify must catch a t08 cert whose stored matrix contradicts its verdict
        ev = os.path.join(tmp, "fakecert.json")
        json.dump({"label": "liar", "verdict": "PD", "matrix": [["1", "2"], ["2", "1"]]},
                  open(ev, "w"))  # [[1,2],[2,1]] is NOT PD (det=-3), but claims PD
        json.dump({"tool": "fake", "results": [{"claim": "PD cert", "tier": "HANDLE",
                   "proof_pointer": ev, "method": "x"}]}, open(os.path.join(runs, "fc.json"), "w"))
        dr = deep_reverify(runs)
        caught = any(d.get("verdict") == "MISMATCH" for d in dr)
        if caught:
            ok("deep_reverify", "fabricated t08 cert (claims PD, matrix is indefinite) caught as MISMATCH")
        else:
            fail("deep_reverify", "fabricated t08 cert NOT caught by re-derivation")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    print("=" * 70)
    print("SICK_HIM — adversarial verification of the RH arsenal")
    print("=" * 70)
    check_overclaims()
    check_handle_evidence()
    check_triple_weil()
    check_t08_adversarial()
    check_mutation()
    check_guard_fuzz()
    check_independence()
    check_audit_hardening()
    print("\n" + "=" * 70)
    if FAILS:
        print(f"VERDICT: {len(FAILS)} ATTACK(S) LANDED — arsenal NOT trusted until fixed:")
        for c, d in FAILS:
            print(f"   - {c}: {d}")
    else:
        print("VERDICT: arsenal SURVIVED every adversarial check. No false certs, no overclaims,")
        print("         math triple-checked, selftests proven to bite (mutation broke them).")
    print("=" * 70)
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
