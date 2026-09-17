#!/usr/bin/env python3
"""rh_report_handle_check.py -- flags prose closure-claims in LEDGER.md that have no
backing evidence marker, i.e. Report-tier language dressed as Handle-tier fact.

Why this specific gap: the ledger ALREADY has a strong Handle-tier convention for Lean
work (EXIT0 + axiom_audit.lean/axiom_gate.sh, see rh-forge-axiom-gate memory) -- that
catches hidden `sorryAx` in compiled .lean files. It does NOT catch a pure-PROSE claim
of "PROVEN"/"unconditional"/"CLOSED" for a numerical/analytic result that was never
Lean-checked and never had its own certified/numeric backing. That is exactly the shape
of the real historical bug on this project: the retracted "F_n>0 for n=86..3732" claim
(LEDGER.md line ~914-963, session 4f, 2026-07-11) -- asserted as an analytic Lemma,
turned out false on direct computation, and was caught only by a LATER re-derivation
session, not by any automated gate. This tool is that gate, run going forward.

Backing evidence (any ONE is enough to clear a section):
  - a `.lean` file reference (Lean-checked, covered by the separate axiom gate)
  - "EXIT0" (this project's own zero-sorry marker)
  - "arb-certified" / "certified" (interval-arithmetic backed)
  - a script/data path reference (operators/xxx.py, results/xxx.md, a plotted N=... run)
  - explicit "RETRACTED"/"FALSE"/strikethrough (~~...~~) in the SAME section -- already
    self-corrected, not a live risk, so it's cleared rather than flagged.

Anything using closure language with NONE of the above is flagged CLAIMED-NOT-HANDLED.
This does not mean the claim is wrong -- it means the ledger text alone doesn't prove it,
same spirit as [[feedback_charisma_report_handle_protocol]]: Report is real work, but it
is not Handle until something outside the prose confirms it.

Usage:
    python3 tools/rh_report_handle_check.py             # full scan, print flagged sections
    python3 tools/rh_report_handle_check.py --selftest   # run against known synthetic cases
"""
from __future__ import annotations
import os, re, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rh_ledger_graph import parse_ledger, LEDGER  # reuse the same section parser

CLOSURE_RE = re.compile(
    r"\b(PROVEN|CLOSED\b|QED|confirmed rigorously|"
    r"rigorously (?:proven|established)|fully (?:proven|closed|resolved))\b",
    re.IGNORECASE,
)
# NOTE: "unconditional(ly)" was dropped 2026-09-14 after the first real scan (267 hits,
# ~17% of all sections) showed it's ordinary analytic-number-theory vocabulary in this
# ledger ("unconditional result" = "doesn't assume RH/GRH"), not a closure claim on its
# own -- keeping it produced overwhelming false-positive noise. Left out rather than
# hand-tuned into a longer regex, since no lexical rule reliably tells "unconditional
# method" from "unconditionally proven"; that distinction needs a human/LLM reader.
# Deliberately still excludes bare "DONE"/"VERIFIED"/"closes" for the same over-common
# reason.

BACKING_RE = re.compile(
    r"(EXIT0|arb-certified|\bcertified\b|`[\w/]+\.py`|operators/\w+\.py|results/[\w.-]+\.md|"
    r"axiom_audit|axiom_gate|\bLean\s*4\b|\bmathlib\b|\blean4\b)",
    re.IGNORECASE,
)
# NOTE: deliberately does NOT include a bare "N=<number>" pattern -- this ledger uses
# lowercase "n=<number>" constantly as ordinary math-index notation (e.g. "n=86..3732",
# the exact retracted claim below), which under re.IGNORECASE collided with an intended
# "N=<computed value>" evidence signal and produced a false negative on the real
# historical bad case. Caught by --selftest 2026-09-14; removed rather than patched
# narrower, since no reliable regex distinguishes "index n" from "parameter N" in prose.
SELF_CORRECTED_RE = re.compile(r"(RETRACTED|~~.*~~|IS FALSE\b|⚠️\s*CORRECTION)", re.IGNORECASE | re.DOTALL)


def check_section(title: str, body: str, lean_files: list) -> str | None:
    """Returns a flag reason string, or None if the section is clean."""
    text = title + "\n" + body
    if not CLOSURE_RE.search(text):
        return None
    if lean_files:
        return None
    if BACKING_RE.search(text):
        return None
    if SELF_CORRECTED_RE.search(text):
        return None
    m = CLOSURE_RE.search(text)
    return f"closure language '{m.group(0)}' with no Lean file / EXIT0 / certified / script-path backing in-section"


def scan(path=LEDGER):
    sections = parse_ledger(path)
    flagged = []
    for s in sections:
        # need the raw body back for the check; re-derive from file since parse_ledger
        # doesn't keep it (graph tool only stores derived fields) -- cheap re-read.
        reason = check_section(s["title"], "", s["lean_files"])
        if reason:
            flagged.append((s, reason))
    return flagged


def scan_with_bodies(path=LEDGER):
    """Like scan() but re-parses with full body text kept, for accurate in-body backing checks."""
    import rh_ledger_graph as G

    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()
    headers = [(i, G.HEADER_RE.match(line.rstrip("\n"))) for i, line in enumerate(lines)]
    headers = [(i, m) for i, m in headers if m]

    flagged = []
    for idx, (line_no, m) in enumerate(headers):
        start = line_no
        end = headers[idx + 1][0] if idx + 1 < len(headers) else len(lines)
        body = "".join(lines[start + 1 : end])
        title = m.group("title").strip()
        lean_files = G.LEAN_RE.findall(title + " " + body)
        reason = check_section(title, body, lean_files)
        if reason:
            sid = G.section_id(m.group("lo"), m.group("hi"), m.group("letter"))
            flagged.append({"id": sid, "line": start + 1, "title": title[:120], "reason": reason})
    return flagged


def selftest():
    import rh_ledger_graph as G

    cases = [
        ("no closure language", "some prose about a partial result", ""),
        ("closure + lean file", "the bound is PROVEN", "see `Foo.lean` for the proof"),
        ("closure + EXIT0", "PROVEN unconditionally", "compile check: EXIT0, no sorries"),
        ("closure + certified", "the region is CLOSED", "arb-certified over [0,1]"),
        ("closure + already retracted", "F_n > 0 is PROVEN for n=86..3732", "~~RETRACTED, IS FALSE~~"),
        ("BAD: closure, no backing, real historical shape",
         "F_n > 0 for n=86..3732 (Lemma C, analytic) is PROVEN unconditionally",
         "no reference given, just asserted"),
        ("BAD: ordinary math index must NOT count as backing",
         "the bound holds for n=86..3732 and is fully proven", "n ranges over integers, no citation"),
    ]
    ok = True
    for name, title, body in cases:
        lean_files = G.LEAN_RE.findall(title + " " + body)  # mirror real scan_with_bodies extraction
        r = check_section(title, body, lean_files)
        expect_flag = name.startswith("BAD")
        got_flag = r is not None
        status = "PASS" if got_flag == expect_flag else "FAIL"
        if status == "FAIL":
            ok = False
        print(f"[{status}] {name}: {'flagged' if got_flag else 'clean'}" + (f" ({r})" if r else ""))
    return 0 if ok else 1


def main():
    if "--selftest" in sys.argv:
        return selftest()
    flagged = scan_with_bodies()
    print(f"{len(flagged)} section(s) flagged CLAIMED-NOT-HANDLED (closure language, no in-section backing):\n")
    for f in flagged:
        print(f"  {f['id']:>10}  L{f['line']:<6}  {f['title']}")
        print(f"             -> {f['reason']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
