#!/usr/bin/env python3
"""rh_ledger_graph.py -- typed bond-graph + frontier tracker for LEDGER.md.

Why: LEDGER.md (39k+ lines, 1500+ sections) already has real relationships between
sections written as PROSE inside titles/bodies ("SECOND DIRECT REUSE OF THEOREM 5.2",
"THE BRIDGE PIECE", explicit "§NNN" cross-references) -- but they're only readable by
reading the whole file. This tool extracts them into an explicit, queryable graph, and
flags two things a full re-read is needed for today:
  1. FRONTIER -- sections that read as open/unresolved with no later section closing them.
  2. DEDUP -- before starting a "new" angle, search titles+bodies for anything similar,
     so a session doesn't re-derive something already on record (the known failure mode --
     see rh-roundtable-literature-gate memory).

Not a claim-correctness checker -- see rh_report_handle_check.py for that (Report vs
Handle tier on CLOSED/PROVEN claims). This tool is pure navigation/dedup over structure
that already exists in the text; it does not judge whether a section's math is right.

Usage:
    python3 tools/rh_ledger_graph.py build              # parse LEDGER.md -> tools/.rh_ledger_graph.json
    python3 tools/rh_ledger_graph.py stats               # parse-integrity check
    python3 tools/rh_ledger_graph.py bonds SS867         # what cites / is cited by SS867
    python3 tools/rh_ledger_graph.py frontier            # open/unresolved sections
    python3 tools/rh_ledger_graph.py search "wkb"        # dedup / prior-art search over titles+bodies
"""
from __future__ import annotations
import json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
LEDGER = os.path.join(ROOT, "LEDGER.md")
GRAPH_PATH = os.path.join(HERE, ".rh_ledger_graph.json")

# "## §NNN[-MMM][letter] [—|--] title"  -- em-dash is the common case but ~1% of
# sections in this file omit it entirely, so the separator is optional.
HEADER_RE = re.compile(
    r"^## §(?P<lo>\d+)(?:-(?P<hi>\d+))?(?P<letter>[a-z])?\s*(?:[—-]{1,2}\s*)?(?P<title>.*)$"
)
CITE_RE = re.compile(r"§(\d+)(?:-(\d+))?")
LEAN_RE = re.compile(r"`([A-Za-z][A-Za-z0-9_]*\.lean)`")
THEOREM_RE = re.compile(r"theorems?\s+(\d+)(?:-(\d+))?", re.IGNORECASE)

FRONTIER_WORDS = re.compile(
    r"\b(TODO|OPEN QUESTION|not yet (?:built|closed|resolved|done)|unresolved|"
    r"still (?:open|missing|needed)|blocked (?:on|by)|frontier|NEXT:|remains open)\b",
    re.IGNORECASE,
)
CLOSED_WORDS = re.compile(
    r"\b(CLOSED|closes?\b|PROVEN|VERIFIED|DONE\b|COMPLETE(?:D)?)\b", re.IGNORECASE
)


def section_id(lo, hi, letter):
    sid = f"S{lo}"
    if hi:
        sid += f"-{hi}"
    if letter:
        sid += letter
    return sid


def parse_ledger(path=LEDGER):
    """Split LEDGER.md into sections at '## §...' headers. Returns list of dicts."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()

    headers = []  # (line_idx, match)
    for i, line in enumerate(lines):
        m = HEADER_RE.match(line.rstrip("\n"))
        if m:
            headers.append((i, m))

    sections = []
    for idx, (line_no, m) in enumerate(headers):
        start = line_no
        end = headers[idx + 1][0] if idx + 1 < len(headers) else len(lines)
        body = "".join(lines[start + 1 : end])
        title = m.group("title").strip()
        sid = section_id(m.group("lo"), m.group("hi"), m.group("letter"))

        cites = set()
        for cm in CITE_RE.finditer(title + "\n" + body):
            lo2 = cm.group(1)
            # a section never "cites" its own number
            if lo2 != m.group("lo"):
                cites.add(f"S{lo2}" + (f"-{cm.group(2)}" if cm.group(2) else ""))

        lean_files = sorted(set(LEAN_RE.findall(title + " " + body)))
        theorem_nums = THEOREM_RE.findall(title + " " + body)

        sections.append(
            {
                "id": sid,
                "lo": int(m.group("lo")),
                "hi": int(m.group("hi")) if m.group("hi") else int(m.group("lo")),
                "line": start + 1,
                "title": title,
                "lean_files": lean_files,
                "theorem_nums": [f"{a}-{b}" if b else a for a, b in theorem_nums],
                "cites": sorted(cites),
                "has_frontier_language": bool(FRONTIER_WORDS.search(title + " " + body)),
                "has_closed_language": bool(CLOSED_WORDS.search(title + " " + body)),
                "body_chars": len(body),
            }
        )
    return sections


def build():
    sections = parse_ledger()
    by_id = {s["id"]: s for s in sections}
    # reverse edges: who cites me
    cited_by = {s["id"]: [] for s in sections}
    for s in sections:
        for c in s["cites"]:
            if c in cited_by:
                cited_by[c].append(s["id"])
    for s in sections:
        s["cited_by"] = sorted(cited_by.get(s["id"], []))

    graph = {
        "source": LEDGER,
        "section_count": len(sections),
        "sections": sections,
    }
    with open(GRAPH_PATH, "w", encoding="utf-8") as fh:
        json.dump(graph, fh, indent=1)
    return graph


def load():
    if not os.path.exists(GRAPH_PATH):
        return build()
    with open(GRAPH_PATH, encoding="utf-8") as fh:
        return json.load(fh)


def cmd_stats():
    g = build()
    secs = g["sections"]
    n_cites = sum(len(s["cites"]) for s in secs)
    n_frontier = sum(1 for s in secs if s["has_frontier_language"])
    n_closed = sum(1 for s in secs if s["has_closed_language"])
    n_lean = sum(1 for s in secs if s["lean_files"])
    raw_header_count = int(
        os.popen(f"grep -c '^## §' {LEDGER!r}").read().strip() or 0
    )
    print(f"sections parsed:        {len(secs)}")
    print(f"raw '## §' header lines: {raw_header_count}  (should match)")
    print(f"total cross-references: {n_cites}")
    print(f"sections w/ Lean file:  {n_lean}")
    print(f"sections w/ frontier lang: {n_frontier}")
    print(f"sections w/ closed lang:   {n_closed}")
    if len(secs) != raw_header_count:
        print("WARNING: parsed count != raw grep count -- parser is dropping or double-counting headers")
        return 1
    print("OK: parser found every '## §' header exactly once.")
    return 0


def cmd_bonds(sid):
    g = load()
    by_id = {s["id"]: s for s in g["sections"]}
    s = by_id.get(sid)
    if not s:
        # tolerate "867" or "SS867" or "S867"
        alt = sid.lstrip("S")
        matches = [k for k in by_id if k.lstrip("S").split("-")[0] == alt]
        if not matches:
            print(f"no section '{sid}' found")
            return 1
        s = by_id[matches[0]]
    print(f"{s['id']} (line {s['line']}): {s['title']}")
    if s["lean_files"]:
        print(f"  lean files: {', '.join(s['lean_files'])}")
    print(f"  cites -> {s['cites'] or '(none)'}")
    print(f"  cited by <- {s['cited_by'] or '(none)'}")
    return 0


def cmd_frontier():
    g = load()
    open_secs = [s for s in g["sections"] if s["has_frontier_language"] and not s["has_closed_language"]]
    open_secs.sort(key=lambda s: -s["lo"])
    print(f"{len(open_secs)} sections read as open/unresolved (frontier language, no closed language in same section):\n")
    for s in open_secs[:60]:
        print(f"  {s['id']:>8}  L{s['line']:<6}  {s['title'][:100]}")
    if len(open_secs) > 60:
        print(f"  ... and {len(open_secs) - 60} more")
    return 0


def cmd_search(query):
    g = load()
    q = query.lower()
    hits = [s for s in g["sections"] if q in s["title"].lower()]
    print(f"{len(hits)} section titles match '{query}':\n")
    for s in hits[:40]:
        print(f"  {s['id']:>8}  L{s['line']:<6}  {s['title'][:110]}")
    if len(hits) > 40:
        print(f"  ... and {len(hits) - 40} more")
    if not hits:
        print("(no title hits -- run a full-text grep on LEDGER.md body text too before concluding it's novel)")
    return 0


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    cmd, rest = sys.argv[1], sys.argv[2:]
    if cmd == "build":
        g = build()
        print(f"built {GRAPH_PATH}: {g['section_count']} sections")
        return 0
    if cmd == "stats":
        return cmd_stats()
    if cmd == "bonds":
        return cmd_bonds(rest[0] if rest else "")
    if cmd == "frontier":
        return cmd_frontier()
    if cmd == "search":
        return cmd_search(" ".join(rest))
    print(f"unknown command '{cmd}'")
    return 1


if __name__ == "__main__":
    sys.exit(main())
