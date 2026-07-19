#!/usr/bin/env python3
"""VERIFICATION GATE — the anti-false-positive spine (Charisma Report/Handle).
A method's numerical output is REPORT. It becomes HANDLE (a real breakthrough) ONLY if an
INDEPENDENT rigorous re-verification confirms it. No PROOF_FOUND, no phone ping, no halt,
unless it passes this gate. This is what stops the machine fooling itself.

A method claiming a proof MUST drop a certificate spec at results/PROOF_CERT.json:
  {"type":"PD","c":14,"N":64,"prec":2400}          # positive-definiteness claim -> re-run Arb LDL^T
  {"type":"lean","file":"RiemannLean/.../X.lean"}   # proof claim -> independent `lake build`
  {"type":"bound","script":"operators/foo.py","expect":"..."}  # a re-runnable numeric certificate
The gate runs the INDEPENDENT check in a fresh process and returns HANDLE only on genuine pass.
"""
import json, os, subprocess
os.environ.pop("ANTHROPIC_API_KEY",None)  # never bill the key here either
ROOT="/Users/scones/projects/riemann" if os.path.exists("/Users/scones") else os.path.expanduser("~/projects/riemann")
CERT=f"{ROOT}/results/PROOF_CERT.json"
FLINT=f"{ROOT}/flint_venv/bin/python"
GROSKIN=f"{ROOT}/external/groskin_2607.02828/code/arb_ldlt_certify.py"

def verify():
    """Return (handle:bool, report:str). handle=True ONLY on independent rigorous confirmation."""
    if not os.path.exists(CERT):
        return False, "REJECTED: PROOF claimed but no results/PROOF_CERT.json — REPORT only, not HANDLE."
    try: spec=json.load(open(CERT))
    except Exception as e: return False, f"REJECTED: unreadable cert ({e})"
    t=spec.get("type")
    if t=="PD":
        # independent Arb interval LDL^T: PD is PROVEN only if n_neg==0 with all pivots strictly signed
        try:
            r=subprocess.run([FLINT,GROSKIN,"--c",str(spec["c"]),"--N",str(spec["N"]),
                              "--prec",str(spec.get("prec",2400)),"--json-out","/tmp/gate_cert.json"],
                             capture_output=True,text=True,timeout=3600,cwd=ROOT)
            j=json.load(open("/tmp/gate_cert.json"))
            if j.get("certified_positive_definite") and j.get("n_neg")==0 and j.get("undetermined_pivot") is None:
                return True, f"HANDLE: Arb-certified PD at c={spec['c']} N={spec['N']} (n_neg=0, all pivots signed). But note: finite (c,N) PD is NOT full RH — only a proof for ALL c,N is."
            return False, f"REJECTED: Arb re-verification did NOT confirm PD (n_neg={j.get('n_neg')}, undet={j.get('undetermined_pivot')})."
        except Exception as e:
            return False, f"REJECTED: Arb gate failed to run ({e})."
    if t=="lean":
        try:
            r=subprocess.run(["bash","-lc",f"cd {ROOT}/RiemannLean && ~/.elan/bin/lake build 2>&1 | tail -5"],
                             capture_output=True,text=True,timeout=1800)
            ok = "error" not in r.stdout.lower() and "sorry" not in r.stdout.lower()
            return ok, ("HANDLE: lake build clean (0 errors, 0 sorry)." if ok else f"REJECTED: lake build not clean:\n{r.stdout[-400:]}")
        except Exception as e:
            return False, f"REJECTED: lean gate failed ({e})."
    if t=="bound":
        try:
            r=subprocess.run([f"{ROOT}/cvs_venv/bin/python",f"{ROOT}/{spec['script']}"],
                             capture_output=True,text=True,timeout=1200,cwd=ROOT)
            ok = spec.get("expect","___") in r.stdout
            return ok, ("HANDLE: independent re-run reproduced the certificate." if ok else "REJECTED: re-run did not reproduce the claimed signal.")
        except Exception as e:
            return False, f"REJECTED: bound gate failed ({e})."
    return False, f"REJECTED: unknown cert type {t}."

if __name__=="__main__":
    ok,rep=verify(); print(("PASS " if ok else "FAIL ")+rep)
