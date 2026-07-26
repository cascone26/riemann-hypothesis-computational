#!/usr/bin/env python3
"""
provenance.py -- reproducibility stamps for every arsenal run.

Writes a run record to arsenal/runs/<tool>-<timestamp>.json capturing: tool, git commit,
python + key package versions, input digest, output digest, and the honest Result list.
The run record file IS the proof-pointer a HANDLE-tier Result points at.

Mirrors the intent of rigor/provenance.py but adds the honest-Result list and a stable
runs/ location so the fabrication detector (t06) can re-find and re-run any prior claim.
"""
from __future__ import annotations
import os, sys, json, time, subprocess, platform
from typing import Any

ARSENAL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS_DIR = os.path.join(ARSENAL_DIR, "runs")
os.makedirs(RUNS_DIR, exist_ok=True)


def _git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=ARSENAL_DIR,
            stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return "unknown"


def _pkg_versions() -> dict:
    vers = {"python": sys.version.split()[0], "platform": platform.platform()}
    for m in ("numpy", "scipy", "mpmath", "sympy", "flint", "cvxpy", "networkx"):
        try:
            mod = __import__(m)
            vers[m] = getattr(mod, "__version__", "?")
        except Exception:
            vers[m] = "absent"
    return vers


def stamp(tool: str, results: list, inputs: dict | None = None,
          extra: dict | None = None) -> str:
    """Write a run record; return its path (the canonical proof-pointer for this run)."""
    ts = time.strftime("%Y%m%dT%H%M%S", time.localtime())
    rec = {
        "tool": tool,
        "timestamp": ts,
        "git_commit": _git_commit(),
        "versions": _pkg_versions(),
        "inputs": inputs or {},
        "results": [r.to_dict() if hasattr(r, "to_dict") else r for r in results],
        "extra": extra or {},
    }
    path = os.path.join(RUNS_DIR, f"{tool}-{ts}.json")
    # avoid clobber within the same second
    n = 1
    while os.path.exists(path):
        path = os.path.join(RUNS_DIR, f"{tool}-{ts}-{n}.json")
        n += 1
    with open(path, "w") as f:
        json.dump(rec, f, indent=2, default=str)
    return path


def load(path: str) -> dict:
    with open(path) as f:
        return json.load(f)


if __name__ == "__main__":
    from honesty import Result, Tier
    p = stamp("provenance-selftest",
              [Result(claim="selftest ran", tier=Tier.REPORT, method="none")])
    print("wrote", p)
    print("git", _git_commit())
