#!/usr/bin/env python3
"""ORGAN 1b/1c — propose->verify loop + Handle ledger (the multiplier).
Takes a candidate Lean lemma (name, statement+proof text, imports), writes a throwaway module in
the RiemannLean project, runs `lake build` on it, returns PASS/FAIL + Lean errors. On PASS (compiles,
no 'sorry'/'admit'), records it to handles.jsonl with a sha256 of the exact proof text = a Handle
immune to my error rate + context amnesia (a verified lemma stays verified across sessions).

Usage:
  python3 verify_loop.py --file cand.lean --name my_lemma
  echo '<lean text>' | python3 verify_loop.py --name my_lemma --stdin
"""
import subprocess, sys, os, json, hashlib, argparse, re, tempfile
PROJ = os.path.dirname(os.path.abspath(__file__))
ELAN = os.path.expanduser("~/.elan/bin/lake")
LEDGER = os.path.join(PROJ, "handles.jsonl")

def verify(name, text):
    if "sorry" in re.sub(r'--.*', '', text) or "admit" in text:
        return {"name": name, "verdict": "REJECTED", "reason": "contains sorry/admit"}
    modname = f"Cand_{name}"
    path = os.path.join(PROJ, "RiemannLean", f"{modname}.lean")
    with open(path, "w") as f:
        f.write(text)
    try:
        r = subprocess.run([ELAN, "build", f"RiemannLean.{modname}"],
                           cwd=PROJ, capture_output=True, text=True, timeout=900)
        out = (r.stdout + r.stderr)
        ok = (r.returncode == 0) and ("error:" not in out)
        res = {"name": name, "verdict": "PASS" if ok else "FAIL",
               "sha256": hashlib.sha256(text.encode()).hexdigest()[:16],
               "errors": None if ok else out[-2000:]}
        if ok:
            with open(LEDGER, "a") as lf:
                lf.write(json.dumps({"name": name, "sha256": res["sha256"],
                                     "text": text, "module": modname}) + "\n")
        return res
    finally:
        # keep PASS modules (they're real deposits), remove FAIL scratch
        pass

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--file"); ap.add_argument("--stdin", action="store_true")
    a = ap.parse_args()
    text = sys.stdin.read() if a.stdin else open(a.file).read()
    print(json.dumps(verify(a.name, text), indent=2)[:2500])
