#!/usr/bin/env python3
"""rh_ledger_dedup.py -- similarity fingerprint for "is this RH idea already on record."

Adapts Base Layer Eyes' six-mark fusion idea (addr+skel+parts+reach/spine/craft, summed,
three-way same/drifted/different verdict -- see [[baselayer-console-project]]) down to
what's actually checkable from LEDGER.md prose: no MinHash lib dependency, plain stdlib.

Two marks, fused:
  - skel  = Jaccard similarity of word 3-shingles (overall shape/wording overlap)
  - parts = containment ratio (does the SMALLER text's shingle set sit almost entirely
            inside the LARGER one -- catches "this is a subset/rewrite of that", which
            plain Jaccard undercounts when the two texts are very different lengths)
Fused by average (this project's analog of BLE's "SUM never product, so one weak mark
can't zero out the verdict"). Verdict bands: >=0.55 SAME, 0.25-0.55 RELATED, <0.25 DIFFERENT.

This is NOT a semantic checker -- it catches near-duplicate wording/paraphrase, not "the
same mathematical idea expressed completely differently." Use alongside rh_ledger_graph.py
search (title keyword search) and t12_roundtable_gate.py (the real anti-rediscovery gate,
which requires a logged literature check) -- this tool is the fast local pre-check before
those, not a replacement for either.

Usage:
    python3 tools/rh_ledger_dedup.py check "<new idea text, a paragraph or so>"
    python3 tools/rh_ledger_dedup.py against S866 "<candidate text>"   # one specific section
    python3 tools/rh_ledger_dedup.py --selftest
"""
from __future__ import annotations
import os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rh_ledger_graph as G

WORD_RE = re.compile(r"[a-z0-9]+")


def shingles(text: str, n: int = 3) -> set:
    words = WORD_RE.findall(text.lower())
    if len(words) < n:
        return {" ".join(words)} if words else set()
    return {" ".join(words[i : i + n]) for i in range(len(words) - n + 1)}


def fuse_score(a: set, b: set) -> tuple:
    if not a or not b:
        return 0.0, 0.0, 0.0
    inter = len(a & b)
    union = len(a | b)
    skel = inter / union if union else 0.0
    smaller = min(len(a), len(b))
    parts = inter / smaller if smaller else 0.0
    fused = (skel + parts) / 2
    return fused, skel, parts


def verdict(fused: float) -> str:
    if fused >= 0.55:
        return "SAME"
    if fused >= 0.25:
        return "RELATED"
    return "DIFFERENT"


def section_texts():
    """Re-parse with bodies kept (rh_ledger_graph's stored graph drops raw body text)."""
    with open(G.LEDGER, encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()
    headers = [(i, G.HEADER_RE.match(line.rstrip("\n"))) for i, line in enumerate(lines)]
    headers = [(i, m) for i, m in headers if m]
    out = []
    for idx, (line_no, m) in enumerate(headers):
        start = line_no
        end = headers[idx + 1][0] if idx + 1 < len(headers) else len(lines)
        body = "".join(lines[start + 1 : end])
        title = m.group("title").strip()
        sid = G.section_id(m.group("lo"), m.group("hi"), m.group("letter"))
        out.append({"id": sid, "line": start + 1, "title": title, "text": title + " " + body})
    return out


def check(candidate_text: str, top_n: int = 10):
    cand_sh = shingles(candidate_text)
    results = []
    for s in section_texts():
        fused, skel, parts = fuse_score(cand_sh, shingles(s["text"]))
        results.append((fused, skel, parts, s))
    results.sort(key=lambda r: -r[0])
    return results[:top_n]


def cmd_check(text):
    hits = check(text)
    print(f"top matches for candidate text ({len(text)} chars):\n")
    for fused, skel, parts, s in hits:
        v = verdict(fused)
        print(f"  [{v:9}] fused={fused:.3f} (skel={skel:.3f} parts={parts:.3f})  {s['id']:>10} L{s['line']:<6} {s['title'][:80]}")
    top = hits[0][0] if hits else 0.0
    if top >= 0.55:
        print("\n>> HIGH similarity to an existing section -- read it before pursuing this as new.")
    elif top >= 0.25:
        print("\n>> Some overlap with existing work -- worth a skim, not necessarily a duplicate.")
    else:
        print("\n>> No strong lexical match found in the ledger (still run t12's literature gate for real prior-art -- this only checks THIS project's own text).")
    return 0


def cmd_against(sid, text):
    target = next((s for s in section_texts() if s["id"] == sid), None)
    if not target:
        print(f"no section '{sid}' found")
        return 1
    fused, skel, parts = fuse_score(shingles(text), shingles(target["text"]))
    print(f"{sid}: {target['title'][:100]}")
    print(f"fused={fused:.3f} (skel={skel:.3f} parts={parts:.3f}) -> {verdict(fused)}")
    return 0


def selftest():
    real = section_texts()
    # pick a real, meaty section as ground truth
    sample = next(s for s in real if s["id"] == "S866")
    exact = sample["text"]
    # near-duplicate: same text with ~30% of characters chopped off the tail (paraphrase proxy)
    near_dup = exact[: int(len(exact) * 0.7)]
    # unrelated: a distant section by id, different topic
    unrelated_id = "S27" if any(s["id"] == "S27" for s in real) else real[0]["id"]
    unrelated = next(s for s in real if s["id"] == unrelated_id)["text"]

    fused_dup, *_ = fuse_score(shingles(exact), shingles(near_dup))
    fused_unrel, *_ = fuse_score(shingles(exact), shingles(unrelated))
    fused_self, *_ = fuse_score(shingles(exact), shingles(exact))

    ok = True
    checks = [
        ("exact self-match scores 1.0", abs(fused_self - 1.0) < 1e-9, fused_self),
        ("truncated near-dup verdict is SAME or RELATED, not DIFFERENT", verdict(fused_dup) != "DIFFERENT", fused_dup),
        ("unrelated section scores lower than the near-dup", fused_unrel < fused_dup, (fused_unrel, fused_dup)),
        ("unrelated section verdict is DIFFERENT or RELATED, not SAME", verdict(fused_unrel) != "SAME", fused_unrel),
    ]
    for name, cond, val in checks:
        status = "PASS" if cond else "FAIL"
        if not cond:
            ok = False
        print(f"[{status}] {name}  (value={val})")
    return 0 if ok else 1


def main():
    if "--selftest" in sys.argv:
        return selftest()
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    if sys.argv[1] == "check":
        return cmd_check(" ".join(sys.argv[2:]))
    if sys.argv[1] == "against":
        return cmd_against(sys.argv[2], " ".join(sys.argv[3:]))
    print("unknown command")
    return 1


if __name__ == "__main__":
    sys.exit(main())
