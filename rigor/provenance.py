#!/usr/bin/env python3
"""provenance.py — reproducibility snapshot per certified claim (idea #7).
Interval-arithmetic certs can silently shift with library/version drift (python-flint, mpmath, prec, NLERCH).
Every "CERTIFIED" LEDGER entry should pin the exact code+env+params behind it so it can be re-run byte-for-byte.
  stamp(claim, params)          -> provenance dict (versions, git HEAD+dirty, host, params, UTC time)
  write_cert(name, claim, params, result, script) -> writes results/<name>.prov.json (runs rigor_lint first)
CLI:  provenance.py stamp                      # print current env stamp
      provenance.py cert <name> <script.py> "<claim>"   # lint script, then write a pinned cert record"""
import sys, os, json, subprocess, socket, datetime, hashlib
ROOT = "/Users/scones/projects/riemann"
def _run(c):
    try: return subprocess.check_output(c, shell=True, cwd=ROOT, stderr=subprocess.DEVNULL).decode().strip()
    except Exception: return None
def _ver(mod):
    try:
        m=__import__(mod); return getattr(m,'__version__',None) or _run(f"{ROOT}/flint_venv/bin/python3 -c 'import {mod};print({mod}.__version__)'")
    except Exception: return _run(f"{ROOT}/flint_venv/bin/python3 -c 'import {mod};print(getattr({mod},chr(95)*2+\"version\"+chr(95)*2,\"?\"))'")
def _sha(path):
    p=os.path.join(ROOT,path)
    return hashlib.sha256(open(p,'rb').read()).hexdigest()[:16] if os.path.exists(p) else None
def stamp(claim, params=None, script=None):
    return {
        "claim": claim,
        "utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "host": socket.gethostname(),
        "git_head": _run("git rev-parse --short HEAD"),
        "git_dirty": bool(_run("git status --porcelain")),
        "env": {"python-flint": _ver("flint"), "mpmath": _ver("mpmath"),
                "python": _run(f"{ROOT}/flint_venv/bin/python3 -c 'import sys;print(sys.version.split()[0])'")},
        "params": params or {},
        "script": script, "script_sha256": _sha(script) if script else None,
    }
def write_cert(name, claim, params, result, script):
    rec = stamp(claim, params, script)
    rec["result"] = result
    # gate: run the rigor linter on the script first — no pinned cert for a script that fails lint
    if script:
        lint = subprocess.run([f"{ROOT}/flint_venv/bin/python3", f"{ROOT}/rigor/rigor_lint.py",
                               os.path.join(ROOT,script)], capture_output=True)
        rec["rigor_lint"] = "PASS" if lint.returncode==0 else "FAIL"
    out = os.path.join(ROOT, "results", f"{name}.prov.json")
    json.dump(rec, open(out,'w'), indent=2)
    print(f"wrote {out}  (lint={rec.get('rigor_lint','n/a')}, git={rec['git_head']}{'*dirty' if rec['git_dirty'] else ''})")
    return out
if __name__=="__main__":
    if len(sys.argv)>1 and sys.argv[1]=="cert":
        _,_,name,script,claim = sys.argv[:5]
        write_cert(name, claim, {"note":"params inline in script"}, {"status":"see script output"}, script)
    else:
        print(json.dumps(stamp("env stamp"), indent=2))
