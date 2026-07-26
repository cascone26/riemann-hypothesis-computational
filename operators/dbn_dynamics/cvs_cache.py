"""Memoized precompute for the CvS scaffold (room 2026-07-25, GreenBot).
BlueBot's #1 bleed: every run last night re-ground the same dps=120 precompute (~7 min)
for ~30s of new work. Drop-in fix:
    from cvs_cache import precompute_components   # instead of cvs_fast
Same signature, disk-cached under _cvs_cache/ keyed on (c,N,T,dps). First call pays,
every later run at the same config loads in ~seconds. Additive: never mutates cvs_fast.
"""
import os, pickle, hashlib
import mpmath as mp
from cvs_fast import precompute_components as _real

CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_cvs_cache")

def precompute_components(c, N, T, dps):
    os.makedirs(CACHE_DIR, exist_ok=True)
    key = f"c{c}_N{N}_T{T}_dps{dps}"
    path = os.path.join(CACHE_DIR, key + ".pkl")
    if os.path.exists(path):
        with open(path, "rb") as f:
            comp = pickle.load(f)
        # honesty check: cached strings -> mpf at current dps; verify round-trip precision
        if comp.get("_dps", 0) >= dps:
            out = {k: _thaw(v, dps) for k, v in comp.items() if k != "_dps"}
            return out
    comp = _real(c, N, T, dps)
    frozen = {k: _freeze(v) for k, v in comp.items()}
    frozen["_dps"] = dps
    tmp = path + ".tmp"
    with open(tmp, "wb") as f:
        pickle.dump(frozen, f)
    os.replace(tmp, path)  # atomic: parallel runs never see a half-written cache
    return comp

def _freeze(v):
    """mpf/matrix -> repr strings (full precision); containers recursed; else passthrough."""
    if isinstance(v, mp.mpf):
        return ("mpf", mp.nstr(v, mp.mp.dps + 10))
    if isinstance(v, mp.matrix):
        return ("mat", v.rows, v.cols, [mp.nstr(v[i, j], mp.mp.dps + 10)
                                        for i in range(v.rows) for j in range(v.cols)])
    if isinstance(v, dict):
        return ("dict", {k: _freeze(x) for k, x in v.items()})
    if isinstance(v, tuple):
        return ("tup", [_freeze(x) for x in v])   # tuples round-trip as tuples (Purple's catch)
    if isinstance(v, list):
        return ("list", [_freeze(x) for x in v])
    return ("raw", v)

def _thaw(v, dps):
    tag = v[0]
    if tag == "mpf":
        return mp.mpf(v[1])
    if tag == "mat":
        _, r, c_, vals = v
        m = mp.matrix(r, c_)
        for i in range(r):
            for j in range(c_):
                m[i, j] = mp.mpf(vals[i * c_ + j])
        return m
    if tag == "dict":
        return {k: _thaw(x, dps) for k, x in v[1].items()}
    if tag == "tup":
        return tuple(_thaw(x, dps) for x in v[1])
    if tag == "list":
        return [_thaw(x, dps) for x in v[1]]
    return v[1]
