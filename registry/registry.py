#!/usr/bin/env python3
"""Structured claims registry (room 2026-07-25). Sits beside LEDGER.md, doesn't replace it.

LEDGER.md is the narrative/provenance record (append-only prose, per feedback_process_provenance).
This registry is the QUERYABLE substrate the room kept reaching for by hand tonight and getting
wrong: BlueBot misremembered §254's status, YellowBot's prediction got misattributed to itself,
GreenBot repeated the misattribution one message later. "What's the current status of claim X" and
"who predicted X" should be a lookup, not scrollback archaeology.

Append-only JSONL, same discipline as LEDGER.md: retracting or updating a claim APPENDS a new
record for the same id, never edits an old line. Current status = the most recent record for that id.

Usage:
  registry.py add --id C255-K --claim "..." --author BlueBot --status live \
      --config '{"c":[9,11,13],"N":40,"dps":120}' --falsifier "..." \
      [--script operators/dbn_dynamics/foo.py] [--proof-pointer "LEDGER.md#..."] [--depends-on C255-J]
  registry.py update --id C255-K --status retracted --falsifier "..."   # append a status change
  registry.py list [--status live] [--author BlueBot]
  registry.py show <id>                  # full history for one id
  registry.py scoreboard                 # per-author record: live/confirmed/retracted/wounded counts
"""
import argparse, json, os, sys
from collections import defaultdict, OrderedDict

REG_DIR = os.path.dirname(os.path.abspath(__file__))
CLAIMS = os.path.join(REG_DIR, "claims.jsonl")

FIELDS = ["id", "ts", "claim", "author", "status", "config", "falsifier",
          "script", "proof_pointer", "depends_on", "notes"]

STATUSES = {"live", "confirmed", "retracted", "wounded", "superseded"}


def load():
    if not os.path.exists(CLAIMS):
        return []
    with open(CLAIMS) as f:
        return [json.loads(line) for line in f if line.strip()]


def append(rec):
    with open(CLAIMS, "a") as f:
        f.write(json.dumps(rec, sort_keys=False) + "\n")


def current_state(records):
    """Most recent record per id, in first-seen order."""
    order = []
    latest = {}
    for r in records:
        if r["id"] not in latest:
            order.append(r["id"])
        latest[r["id"]] = r
    return [latest[i] for i in order]


def cmd_add(a):
    if a.status not in STATUSES:
        sys.exit(f"status must be one of {sorted(STATUSES)}")
    rec = OrderedDict()
    rec["id"] = a.id
    rec["ts"] = a.ts
    rec["claim"] = a.claim
    rec["author"] = a.author
    rec["status"] = a.status
    rec["config"] = json.loads(a.config) if a.config else {}
    rec["falsifier"] = a.falsifier or ""
    rec["script"] = a.script or ""
    rec["proof_pointer"] = a.proof_pointer or ""
    rec["depends_on"] = a.depends_on or []
    rec["notes"] = a.notes or ""
    append(rec)
    print(f"added {a.id} status={a.status}")


def cmd_update(a):
    records = load()
    ids = {r["id"] for r in records}
    if a.id not in ids:
        sys.exit(f"unknown id {a.id} — use add for a new claim")
    prev = [r for r in records if r["id"] == a.id][-1]
    rec = OrderedDict(prev)
    rec["ts"] = a.ts
    if a.status:
        if a.status not in STATUSES:
            sys.exit(f"status must be one of {sorted(STATUSES)}")
        rec["status"] = a.status
    if a.falsifier:
        rec["falsifier"] = a.falsifier
    if a.proof_pointer:
        rec["proof_pointer"] = a.proof_pointer
    if a.notes:
        rec["notes"] = (prev.get("notes", "") + " | " + a.notes).strip(" |")
    append(rec)
    print(f"updated {a.id} -> status={rec['status']}")


def cmd_list(a):
    state = current_state(load())
    if a.status:
        state = [r for r in state if r["status"] == a.status]
    if a.author:
        state = [r for r in state if r["author"] == a.author]
    if not state:
        print("(no matching claims)")
        return
    for r in state:
        cfg = ",".join(f"{k}={v}" for k, v in r.get("config", {}).items())
        print(f"[{r['status']:>10}] {r['id']:<10} ({r['author']:<10}) {r['claim'][:70]}"
              f"{'...' if len(r['claim']) > 70 else ''}  {{{cfg}}}")


def cmd_show(a):
    records = [r for r in load() if r["id"] == a.claim_id]
    if not records:
        sys.exit(f"no such id {a.claim_id}")
    for r in records:
        print(json.dumps(r, indent=2))
        print("-" * 60)


def cmd_scoreboard(a):
    state = current_state(load())
    by_author = defaultdict(lambda: defaultdict(int))
    for r in state:
        by_author[r["author"]][r["status"]] += 1
    print(f"{'author':<12} {'live':>5} {'confirmed':>10} {'retracted':>10} {'wounded':>8} {'superseded':>11}")
    for author, counts in sorted(by_author.items()):
        print(f"{author:<12} {counts['live']:>5} {counts['confirmed']:>10} {counts['retracted']:>10} "
              f"{counts['wounded']:>8} {counts['superseded']:>11}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("add")
    p.add_argument("--id", required=True)
    p.add_argument("--ts", required=True, help="e.g. 2026-07-25T02:07")
    p.add_argument("--claim", required=True)
    p.add_argument("--author", required=True)
    p.add_argument("--status", required=True)
    p.add_argument("--config", default=None, help="JSON string, e.g. '{\"c\":[9,11,13]}'")
    p.add_argument("--falsifier", default=None)
    p.add_argument("--script", default=None)
    p.add_argument("--proof-pointer", default=None)
    p.add_argument("--depends-on", action="append", default=[])
    p.add_argument("--notes", default=None)
    p.set_defaults(func=cmd_add)

    p = sub.add_parser("update")
    p.add_argument("--id", required=True)
    p.add_argument("--ts", required=True)
    p.add_argument("--status", default=None)
    p.add_argument("--falsifier", default=None)
    p.add_argument("--proof-pointer", default=None)
    p.add_argument("--notes", default=None)
    p.set_defaults(func=cmd_update)

    p = sub.add_parser("list")
    p.add_argument("--status", default=None)
    p.add_argument("--author", default=None)
    p.set_defaults(func=cmd_list)

    p = sub.add_parser("show")
    p.add_argument("claim_id")
    p.set_defaults(func=cmd_show)

    p = sub.add_parser("scoreboard")
    p.set_defaults(func=cmd_scoreboard)

    a = ap.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
