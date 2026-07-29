#!/usr/bin/env python3
"""
lit_recheck.py -- scheduled live literature re-pull for the RH lit-gate.

t01_route_prefilter's PRIOR_ART dict is a snapshot taken at build time (2026-07-26). An engine
meant to run indefinitely can't rely on a snapshot -- new papers publish continuously, and a
static dict is exactly the kind of "standing constraint" that goes silently stale over a long
run (the same governance-decay shape as a rule dropped by context compaction, just at the
literature layer instead of the prompt layer). This re-queries arXiv per tracked topic since its
last-checked timestamp and stages new hits for human/session review -- it does NOT auto-merge
into PRIOR_ART, the same way orient-before-act says two plausible candidates is a stop-and-ask,
not a guess.

Run:  python3 engine/lit_recheck.py            # check only topics due for a refresh
      python3 engine/lit_recheck.py --force     # ignore due-dates, check everything now
State:  engine/lit_state.json   (per-topic last-checked timestamp + seen arXiv ids)
Output: engine/LIT_WATCH.md     (append-only staging log of new hits needing human review)
"""
from __future__ import annotations
import json, os, sys, datetime, urllib.request, urllib.parse
import xml.etree.ElementTree as ET

ROOT = "/Users/scones/projects/riemann" if os.path.exists("/Users/scones") else os.path.expanduser("~/projects/riemann")
ENG = os.path.join(ROOT, "engine")
STATE = os.path.join(ENG, "lit_state.json")
WATCH = os.path.join(ENG, "LIT_WATCH.md")

# Topic -> arXiv search query. Seeded from t01_route_prefilter.PRIOR_ART keys + this project's
# live keystones (LEDGER 329-332). Add a topic here whenever a new route/keystone gets adopted.
TOPICS = {
    "jensen_laguerre_polya": 'all:"Jensen polynomial" AND all:"Riemann hypothesis"',
    "de_bruijn_newman": 'all:"de Bruijn-Newman"',
    "beurling_selberg": 'all:"Beurling-Selberg" AND all:"Riemann"',
    "crystalline_quasicrystal": 'all:"Fourier quasicrystal" AND all:"Riemann zeta"',
    "connes_consani": 'all:"Connes" AND all:"Weil positivity"',
    "arithmetic_hodge_index": 'all:"arithmetic Hodge index" OR all:"Faltings-Hriljac"',
    "de_branges_spaces": 'all:"de Branges space" AND all:"Riemann hypothesis"',
    "mollifier_levinson": 'all:"mollifier" AND all:"zeros of the Riemann zeta"',
}
REFRESH_DAYS = 14  # cheap, bounded re-check cadence -- never spams arXiv, never goes stale for long


def now_iso():
    return datetime.datetime.now().isoformat(timespec="seconds")


def load_state():
    try:
        return json.load(open(STATE))
    except Exception:
        return {}


def save_state(s):
    json.dump(s, open(STATE, "w"), indent=2)


def arxiv_query(query, max_results=15):
    url = "http://export.arxiv.org/api/query?" + urllib.parse.urlencode({
        "search_query": query, "sortBy": "submittedDate", "sortOrder": "descending",
        "max_results": max_results})
    with urllib.request.urlopen(url, timeout=30) as r:
        raw = r.read()
    root = ET.fromstring(raw)
    ns = {"a": "http://www.w3.org/2005/Atom"}
    out = []
    for entry in root.findall("a:entry", ns):
        id_ = entry.find("a:id", ns).text.strip()
        title = " ".join(entry.find("a:title", ns).text.split())
        published = entry.find("a:published", ns).text.strip()
        out.append({"id": id_, "title": title, "published": published})
    return out


def is_due(topic, state, force):
    if force:
        return True
    last = state.get(topic, {}).get("checked_at")
    if not last:
        return True
    days = (datetime.datetime.now() - datetime.datetime.fromisoformat(last)).days
    return days >= REFRESH_DAYS


def main():
    force = "--force" in sys.argv
    state = load_state()
    due_topics = [t for t in TOPICS if is_due(t, state, force)]
    new_hits = []
    for topic in due_topics:
        q = TOPICS[topic]
        try:
            entries = arxiv_query(q)
        except Exception as e:
            print(f"[{now_iso()}] {topic}: query failed ({e})")
            continue
        seen = set(state.get(topic, {}).get("seen_ids", []))
        fresh = [e for e in entries if e["id"] not in seen]
        state[topic] = {"checked_at": now_iso(),
                         "seen_ids": list(seen | {e["id"] for e in entries})[-200:]}
        if fresh:
            new_hits.append((topic, fresh))
        print(f"[{now_iso()}] {topic}: {len(fresh)} new / {len(entries)} fetched")
    save_state(state)
    if new_hits:
        with open(WATCH, "a") as f:
            f.write(f"\n## Lit re-check {now_iso()}\n\n")
            for topic, fresh in new_hits:
                f.write(f"### {topic}\n")
                for e in fresh:
                    f.write(f"- [{e['published'][:10]}] {e['title']} -- {e['id']}\n")
        print(f"[{now_iso()}] {len(new_hits)} topic(s) with new hits -- staged in {WATCH}, "
              f"NOT auto-merged into PRIOR_ART (needs human/session review before folding in)")
    else:
        print(f"[{now_iso()}] no new hits ({len(due_topics)} topic(s) checked, "
              f"{len(TOPICS) - len(due_topics)} not yet due)")


if __name__ == "__main__":
    main()
