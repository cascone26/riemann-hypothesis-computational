#!/usr/bin/env python3
"""ORGAN 1b/1c — propose->verify loop + Handle ledger (the multiplier). HARDENED per superprism
skeptic audit (2026-08-05): the gate now checks the AXIOM DEPENDENCY SET, not just a sorry regex.

Gate (a PASS must clear ALL):
  1. text has no `sorry`/`admit`/`axiom`/`native_decide`/`@[implemented_by]` after stripping BOTH
     line (`--`) and block (`/- -/`) comments (cheap pre-filter);
  2. `lake env lean <file>` compiles with no `error:`;
  3. for EVERY top-level theorem/lemma, `#print axioms <name>` reports an axiom set that is a
     SUBSET of the trusted base {propext, Classical.choice, Quot.sound}. Any `sorryAx`, custom
     axiom, or `Lean.ofReduceBool`/`Lean.trustCompiler` (native_decide) => FAIL.
This is the real Lean trust boundary: a proof is a Handle iff it reduces to the trusted kernel axioms.

Usage: python3 verify_loop.py --file cand.lean --name my_lemma
Provenance recorded: name, sha256, timestamp, lean_version, machine, axioms, build_time."""
import subprocess, sys, os, json, hashlib, argparse, re, time, datetime
PROJ = os.path.dirname(os.path.abspath(__file__))
LAKE = os.path.expanduser("~/.elan/bin/lake")
LEDGER = os.path.join(PROJ, "handles.jsonl")
TRUSTED = {"propext", "Classical.choice", "Quot.sound"}
BANNED_TOK = ["sorry", "admit", "axiom", "native_decide", "implemented_by", "trustCompiler"]

def strip_comments(text):
    t = re.sub(r'--.*', '', text)
    t = re.sub(r'/-.*?-/', '', t, flags=re.DOTALL)
    return t

def lean_version():
    try:
        return subprocess.run([LAKE, "env", "lean", "--version"], cwd=PROJ,
                              capture_output=True, text=True, timeout=60).stdout.strip()
    except Exception: return "unknown"

def verify(name, text):
    clean = strip_comments(text)
    for tok in BANNED_TOK:
        if tok in clean:
            return {"name": name, "verdict": "REJECTED", "reason": f"contains '{tok}'"}
    names = re.findall(r'^\s*(?:theorem|lemma)\s+([A-Za-z_][A-Za-z0-9_\'\.]*)', text, re.MULTILINE)
    if not names:
        return {"name": name, "verdict": "REJECTED", "reason": "no theorem/lemma found"}
    modname = f"Cand_{name}"
    path = os.path.join(PROJ, "RiemannLean", f"{modname}.lean")
    body = text + "\n" + "\n".join(f"#print axioms {nm}" for nm in names) + "\n"
    with open(path, "w") as f: f.write(body)
    t0 = time.time()
    r = subprocess.run([LAKE, "env", "lean", path], cwd=PROJ, capture_output=True, text=True, timeout=900)
    dt = round(time.time() - t0, 1)
    out = r.stdout + r.stderr
    if r.returncode != 0 or "error:" in out:
        return {"name": name, "verdict": "FAIL", "build_time": dt, "errors": out[-2000:]}
    # parse #print axioms results. Lean prints "'name' depends on axioms: [a, b, c]"
    # or "'name' does not depend on any axioms"
    used = set()
    for m in re.finditer(r"depends on axioms:\s*\[([^\]]*)\]", out):
        used |= {a.strip() for a in m.group(1).split(",") if a.strip()}
    extra = used - TRUSTED
    if extra:
        return {"name": name, "verdict": "FAIL", "build_time": dt,
                "reason": f"untrusted axioms: {sorted(extra)}", "axioms": sorted(used)}
    rec = {"name": name, "verdict": "PASS", "sha256": hashlib.sha256(text.encode()).hexdigest()[:16],
           "timestamp": datetime.datetime.now().isoformat(timespec="seconds"),
           "lean_version": lean_version(), "machine": "Sam", "axioms": sorted(used),
           "build_time": dt, "theorems": names, "module": modname}
    with open(LEDGER, "a") as lf: lf.write(json.dumps({**rec, "text": text}) + "\n")
    return rec

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True); ap.add_argument("--file"); ap.add_argument("--stdin", action="store_true")
    a = ap.parse_args()
    text = sys.stdin.read() if a.stdin else open(a.file).read()
    print(json.dumps(verify(a.name, text), indent=2)[:2500])
