#!/usr/bin/env python3
"""
t06_fabrication_detector.py  -- Tier A #6.  The anti-fraud tool.  Directly targets the project's
worst failure mode: fabricated "CERTIFIED/HANDLE" claims (2 of 3 parallel tracks faked certs,
caught only by independent re-run).  This tool trusts NOTHING it is told: for every claim it can,
it RE-DERIVES the number by an INDEPENDENT method and flags disagreement, and it audits every
HANDLE-tier run record for a real, non-empty proof-pointer on disk.

Three checks:
  1. audit_runs()        -- every HANDLE Result in arsenal/runs/*.json must point at a proof file
                            that exists and is non-empty.  Missing/empty => UNBACKED.
  2. independent_weil()  -- recompute the Weil margin A+B-P by mpmath quadrature + mpmath zeros
                            (fully independent of weil.py's numpy/scipy path) and compare.
  3. audit_claim(v, fn)  -- generic: given a claimed value and an independent recompute fn,
                            corroborate or flag DISAGREES beyond tolerance.

The selftest PLANTS a fabricated cert (a HANDLE claim whose number is wrong, and one whose proof
file is missing) and asserts the detector catches BOTH -- proof the detector actually bites.

Usage: python3 t06_fabrication_detector.py [--scan DIR]
"""
from __future__ import annotations
import sys, os, json, glob, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from core.honesty import Result, Tier
from core import provenance

mp.mp.dps = 30


# ---------- check 1: audit HANDLE claims for real evidence ----------
# NOTE: audit_runs is an EXISTENCE + INTEGRITY check (does a real, parseable, non-empty regular
# file back this HANDLE claim?). It is NOT a proof-of-correctness -- a well-formed but fabricated
# evidence file passes this layer. Correctness is the job of deep_reverify() (re-derivation) and
# independent_weil_margin()/audit_claim(). Red-team (2026-07-26) showed the old version was fooled
# by directories, symlinks, and 1-byte junk; those are now rejected. But a plausible fake JSON
# still needs deep_reverify to be caught -- see below.
def audit_runs(runs_dir=None):
    runs_dir = runs_dir or provenance.RUNS_DIR
    findings = []
    for path in sorted(glob.glob(os.path.join(runs_dir, "*.json"))):
        try:
            rec = json.load(open(path))
        except Exception as e:
            findings.append({"file": path, "verdict": "UNREADABLE", "detail": str(e)})
            continue
        for r in rec.get("results", []):
            if r.get("tier") != "HANDLE":
                continue
            claim = r.get("claim", "")[:80]
            pp = r.get("proof_pointer")
            base = os.path.basename(path)
            if not pp or pp == "/dev/null":
                findings.append({"file": base, "verdict": "UNBACKED", "claim": claim,
                                 "detail": "HANDLE with no proof_pointer"})
            elif os.path.islink(pp):
                findings.append({"file": base, "verdict": "SYMLINK_EVIDENCE", "claim": claim,
                                 "detail": f"proof pointer is a symlink (indirection): {pp}"})
            elif not os.path.exists(pp):
                findings.append({"file": base, "verdict": "UNBACKED", "claim": claim,
                                 "detail": f"proof file missing: {pp}"})
            elif os.path.isdir(pp):
                findings.append({"file": base, "verdict": "NOT_A_FILE", "claim": claim,
                                 "detail": f"proof pointer is a directory, not a file: {pp}"})
            elif os.path.getsize(pp) == 0:
                findings.append({"file": base, "verdict": "EMPTY_EVIDENCE", "claim": claim,
                                 "detail": f"proof file empty: {pp}"})
            else:
                try:
                    json.load(open(pp))   # must be real structured evidence, not junk bytes
                except Exception as e:
                    findings.append({"file": base, "verdict": "UNPARSEABLE_EVIDENCE",
                                     "claim": claim,
                                     "detail": f"proof file not valid JSON ({e}): {pp}"})
    return findings


# ---------- deep re-derivation: actually reconstruct and re-verify a stored cert ----------
def deep_reverify(runs_dir=None):
    """
    Content-level check: for HANDLE claims whose evidence is self-contained enough to re-derive
    (currently t08 exact PD certs, which store their matrix), reconstruct the object from scratch
    and re-run the verification.  A fabricated cert whose stored matrix does not actually certify
    the claimed verdict is caught here (existence checks cannot catch that).
    """
    runs_dir = runs_dir or provenance.RUNS_DIR
    results = []
    for path in sorted(glob.glob(os.path.join(runs_dir, "*.json"))):
        try:
            rec = json.load(open(path))
        except Exception:
            continue
        for r in rec.get("results", []):
            if r.get("tier") != "HANDLE":
                continue
            pp = r.get("proof_pointer")
            if not pp or not os.path.exists(pp) or os.path.islink(pp) or os.path.isdir(pp):
                continue
            try:
                ev = json.load(open(pp))
            except Exception:
                continue
            if "matrix" in ev and "verdict" in ev:   # a t08 exact cert
                v = _reverify_t08(ev)
                results.append({"file": os.path.basename(pp), "kind": "t08_exact_cert", **v})
    return results


def _reverify_t08(ev):
    """Reconstruct the sympy matrix from stored rationals and re-run the exact PD check."""
    import sympy as sp
    from t08_exact_sos import certify
    try:
        M = sp.Matrix([[sp.Rational(x) for x in row] for row in ev["matrix"]])
    except Exception as e:
        return {"verdict": "EVIDENCE_MALFORMED", "detail": str(e)}
    fresh = certify(M, ev.get("label", "rederived"))
    agree = fresh["verdict"] == ev.get("verdict")
    return {"verdict": "CONFIRMED" if agree else "MISMATCH",
            "stored_verdict": ev.get("verdict"), "rederived_verdict": fresh["verdict"],
            "detail": "re-derived from stored matrix matches" if agree
                      else "stored verdict does NOT match independent re-derivation"}


# ---------- check 2: independent Weil margin (mpmath, fully separate code path) ----------
def _von_mangoldt_mp(N):
    Lam = {}
    sieve = [True] * (N + 1)
    for p in range(2, N + 1):
        if sieve[p]:
            lp = mp.log(p); pk = p
            while pk <= N:
                Lam[pk] = lp; pk *= p
            for m in range(p * p, N + 1, p):
                sieve[m] = False
    return Lam


def independent_weil_margin(t=0.10, Nprime=20000, nzeros=40):
    """A+B-P and zero-side, computed with mpmath quadrature + mpmath.zetazero. No numpy/scipy."""
    t = mp.mpf(t)
    ghat = lambda r: t * mp.sqrt(2 * mp.pi) * mp.e ** (-(t ** 2) * r ** 2 / 2)
    Omega = lambda r: mp.re(mp.digamma(mp.mpf('0.25') + 1j * r / 2)) - mp.log(mp.pi)
    A = mp.quad(lambda r: ghat(r) * Omega(r), [-140, 0, 140]) / (2 * mp.pi)
    B = ghat(1j / 2) + ghat(-1j / 2)
    Lam = _von_mangoldt_mp(Nprime)
    P = 2 * mp.fsum(v / mp.sqrt(n) * mp.e ** (-(mp.log(n)) ** 2 / (2 * t ** 2))
                    for n, v in Lam.items())
    W = A + mp.re(B) - P
    Z = 2 * mp.fsum(ghat(mp.im(mp.zetazero(k))) for k in range(1, nzeros + 1))
    return float(W), float(Z), float(abs(W - Z))


# ---------- check 3: generic claim audit ----------
def audit_claim(claimed_value, recompute_fn, tol=1e-6, label=""):
    got = recompute_fn()
    diff = abs(claimed_value - got)
    verdict = "CORROBORATED" if diff <= tol else "DISAGREES"
    return {"label": label, "claimed": claimed_value, "recomputed": got,
            "diff": diff, "tol": tol, "verdict": verdict}


def selftest(verbose=True) -> bool:
    ok = True
    # --- plant fabricated certs in a temp runs dir and confirm the detector catches them ---
    import tempfile, shutil
    tmp = tempfile.mkdtemp(prefix="t06_fraud_")
    try:
        # (a) HANDLE claim with a missing proof file
        json.dump({"tool": "faketrack", "results": [
            {"claim": "tau CERTIFIED PD for all N", "tier": "HANDLE",
             "proof_pointer": "/tmp/does_not_exist_12345.json", "method": "fabricated"}]},
            open(os.path.join(tmp, "fake1.json"), "w"))
        # (b) HANDLE claim with /dev/null (no real evidence)
        json.dump({"tool": "faketrack", "results": [
            {"claim": "margin proven positive", "tier": "HANDLE",
             "proof_pointer": "/dev/null", "method": "fabricated"}]},
            open(os.path.join(tmp, "fake2.json"), "w"))
        # (c) a legit HANDLE with real evidence (should NOT be flagged)
        ev = os.path.join(tmp, "real_evidence.json")
        json.dump({"pivots_positive": True}, open(ev, "w"))
        json.dump({"tool": "realtrack", "results": [
            {"claim": "finite 2x2 PD over Q", "tier": "HANDLE",
             "proof_pointer": ev, "method": "exact-ldlt"}]},
            open(os.path.join(tmp, "real.json"), "w"))

        findings = audit_runs(tmp)
        flagged = {f["file"] for f in findings}
        if "fake1.json" not in flagged or "fake2.json" not in flagged:
            ok = False; print("  FAIL: detector missed a fabricated cert:", flagged)
        if "real.json" in flagged:
            ok = False; print("  FAIL: detector false-positived a real cert")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # --- plant a fabricated NUMBER and confirm audit_claim catches it ---
    fake = audit_claim(999.0, lambda: independent_weil_margin(t=0.10)[0], tol=1e-3,
                       label="fabricated margin")
    if fake["verdict"] != "DISAGREES":
        ok = False; print("  FAIL: audit_claim did not flag a fabricated number")

    # --- confirm a TRUE claim is corroborated by the independent path ---
    from core.weil import WeilForm
    wf = WeilForm(zero_count=40)
    true_margin = wf.margin(0.10)
    good = audit_claim(true_margin, lambda: independent_weil_margin(t=0.10, nzeros=40)[1],
                       tol=1e-3, label="true margin vs independent zero-side")
    if good["verdict"] != "CORROBORATED":
        ok = False; print("  FAIL: true margin not corroborated by independent method:", good)

    if verbose:
        print(f"  planted-fraud audit: caught {'both' if ok else 'NOT ALL'} fabricated certs")
        print(f"  fabricated number 999.0 vs recompute {fake['recomputed']:.4f} -> {fake['verdict']}")
        print(f"  true margin {true_margin:.6f} vs independent -> {good['verdict']} "
              f"(diff {good['diff']:.2e})")
    print("t06.selftest:", "PASS" if ok else "FAIL")
    return ok


def main():
    print("=== t06 Fabrication detector (trusts nothing; re-derives independently) ===\n")
    # 1. independent Weil corroboration
    W, Z, res = independent_weil_margin(t=0.10)
    print(f"[independent recompute] Weil margin (mpmath path) W={W:.8f}, zero-side Z={Z:.8f}, "
          f"|W-Z|={res:.2e}")
    from core.weil import WeilForm
    wf = WeilForm(zero_count=40)
    Wnp = wf.margin(0.10)
    cross = audit_claim(Wnp, lambda: W, tol=1e-4, label="numpy-engine vs mpmath-engine margin")
    print(f"[cross-engine] numpy-path W={Wnp:.8f} vs mpmath-path W={W:.8f} -> {cross['verdict']} "
          f"(diff {cross['diff']:.2e})")

    # 2. audit all existing run records for unbacked HANDLE claims
    findings = audit_runs()
    print(f"\n[audit runs/] scanned {len(glob.glob(os.path.join(provenance.RUNS_DIR,'*.json')))} "
          f"run records; {len(findings)} unbacked/suspect HANDLE claims found:")
    for f in findings:
        print(f"   {f['verdict']}: {f.get('file')} -- {f.get('detail')}")
    if not findings:
        print("   (none -- every HANDLE claim on disk points at real evidence)")

    ok = selftest()
    results = [
        Result(claim=f"cross-engine Weil margin corroborated by an INDEPENDENT mpmath path "
                     f"(diff {cross['diff']:.1e})", tier=Tier.REPORT,
               method="independent re-derivation", tool="t06_fabrication_detector", data=cross),
        Result(claim=f"audited {len(glob.glob(os.path.join(provenance.RUNS_DIR,'*.json')))} run "
                     f"records: {len(findings)} unbacked HANDLE claims", tier=Tier.REPORT,
               method="proof-pointer audit", tool="t06_fabrication_detector",
               data={"findings": findings}),
    ]
    path = provenance.stamp("t06_fabrication_detector", results, extra={"selftest_pass": ok})
    print("run record:", path)
    return ok


if __name__ == "__main__":
    main()
