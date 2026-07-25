#!/usr/bin/env python3
"""Compute scatter dispatcher (room 2026-07-25, Purple/Yellow's ask).
Correction on record: dbn_dynamics jobs are mpmath, arbitrary-precision, CPU-bound —
GPUs buy nothing. The win is CORE COUNT. Olympus has 16 cores sitting mostly idle;
this scatters independent (c,N,dps,T) jobs across them over SSH and gathers results
back, falling back to local multiprocessing if Olympus is unreachable.

Usage:
  python3 tools/dispatch.py --jobs jobs.json [--local] [--sync]
  python3 tools/dispatch.py --script _moment_ladder.py --args "7 40 120" --args "15 40 120"

jobs.json: [{"script": "_moment_ladder.py", "args": ["7","40","120"], "label": "c7"}, ...]

Each job's stdout/stderr land in results/dispatch/<label or index>.log plus a
manifest.json with exit codes + wall time, so a claim's proof_pointer can cite
"results/dispatch/<run_id>/manifest.json" instead of a scrollback screenshot.
"""
import argparse, json, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.path.expanduser("~/projects/riemann")
DBN = "operators/dbn_dynamics"
REMOTE_HOST = "olympus"
REMOTE_ROOT = "~/projects/riemann"
# Only what a dbn_dynamics job actually needs remotely — not the 288M cvs_venv.
SYNC_PATHS = [
    f"{DBN}/",
    "cvs_venv/lib/python3.14/site-packages/connes_cvs",
]


def remote_alive(timeout=6):
    try:
        r = subprocess.run(["ssh", "-o", f"ConnectTimeout={timeout}", REMOTE_HOST, "nproc"],
                            capture_output=True, text=True, timeout=timeout + 4)
        return r.returncode == 0, int(r.stdout.strip() or 0)
    except Exception:
        return False, 0


def sync_to_remote():
    subprocess.run(["ssh", REMOTE_HOST, f"mkdir -p {REMOTE_ROOT}/{DBN} "
                     f"{REMOTE_ROOT}/cvs_venv/lib/python3.14/site-packages"], check=True)
    subprocess.run(["rsync", "-az", "--exclude=__pycache__", "--exclude=_cvs_cache",
                     f"{ROOT}/{DBN}/", f"{REMOTE_HOST}:{REMOTE_ROOT}/{DBN}/"], check=True)
    subprocess.run(["rsync", "-az",
                     f"{ROOT}/cvs_venv/lib/python3.14/site-packages/connes_cvs",
                     f"{REMOTE_HOST}:{REMOTE_ROOT}/cvs_venv/lib/python3.14/site-packages/"], check=True)


def run_job(job, idx, out_dir, remote):
    label = job.get("label", f"job{idx}")
    script = job["script"]
    args = job.get("args", [])
    log_path = os.path.join(out_dir, f"{label}.log")
    t0 = time.time()
    if remote:
        remote_cmd = (f"cd {REMOTE_ROOT} && PYTHONPATH={REMOTE_ROOT}/cvs_venv/lib/python3.14/site-packages "
                       f"python3 {DBN}/{script} {' '.join(args)}")
        cmd = ["ssh", REMOTE_HOST, remote_cmd]
    else:
        env = dict(os.environ)
        cmd = [sys.executable, os.path.join(ROOT, DBN, script)] + list(args)
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=job.get("timeout", 3600), cwd=ROOT)
        ok = r.returncode == 0
        out = r.stdout + ("\n--STDERR--\n" + r.stderr if r.stderr else "")
    except subprocess.TimeoutExpired as e:
        ok = False
        out = f"TIMEOUT after {job.get('timeout', 3600)}s\n{e.stdout or ''}"
    with open(log_path, "w") as f:
        f.write(out)
    return {"label": label, "script": script, "args": args, "ok": ok,
            "wall_s": round(time.time() - t0, 1), "log": log_path, "remote": remote}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", help="path to jobs.json")
    ap.add_argument("--script", help="single script name under operators/dbn_dynamics/")
    ap.add_argument("--args", action="append", default=[], help="space-separated arg string, repeatable")
    ap.add_argument("--local", action="store_true", help="force local execution, skip Olympus")
    ap.add_argument("--sync", action="store_true", help="rsync dbn_dynamics + connes_cvs to Olympus first")
    ap.add_argument("--run-id", default=None)
    a = ap.parse_args()

    if a.jobs:
        jobs = json.load(open(a.jobs))
    elif a.script:
        jobs = [{"script": a.script, "args": s.split(), "label": f"{a.script.rstrip('.py')}_{i}"}
                for i, s in enumerate(a.args)]
    else:
        ap.error("need --jobs or --script/--args")

    remote = False
    if not a.local:
        alive, cores = remote_alive()
        if alive:
            remote = True
            print(f"[dispatch] Olympus alive, {cores} cores.")
            if a.sync:
                print("[dispatch] syncing dbn_dynamics + connes_cvs to Olympus...")
                sync_to_remote()
        else:
            print("[dispatch] Olympus unreachable, falling back to local multiprocessing.")

    run_id = a.run_id or f"run_{int(time.time())}"
    out_dir = os.path.join(ROOT, "results", "dispatch", run_id)
    os.makedirs(out_dir, exist_ok=True)

    workers = 16 if remote else max(1, os.cpu_count() - 1)
    print(f"[dispatch] {len(jobs)} jobs, {workers} workers, remote={remote}")

    results = []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(run_job, job, i, out_dir, remote): job for i, job in enumerate(jobs)}
        for fut in as_completed(futs):
            r = fut.result()
            status = "OK" if r["ok"] else "FAIL"
            print(f"  [{status}] {r['label']} ({r['wall_s']}s)")
            results.append(r)

    manifest = {"run_id": run_id, "remote": remote, "n_jobs": len(jobs), "results": results}
    manifest_path = os.path.join(out_dir, "manifest.json")
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"[dispatch] manifest: results/dispatch/{run_id}/manifest.json")


if __name__ == "__main__":
    main()
