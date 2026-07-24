#!/usr/bin/env python3
"""
arxiv_watch.py — daily check for new RH/Weil-positivity-adjacent arXiv papers.

Queries a fixed set of targeted searches against the arXiv API, diffs against
a seen-IDs store, and reports/ntfy's anything new. Zero Claude tokens to run.

Run manually:  python3 arxiv_watch.py
Scheduled by:  ~/Library/LaunchAgents/com.jacob.riemann-arxiv-watch.plist
"""
import json
import os
import time
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET

ROOT = os.path.dirname(os.path.abspath(__file__))
STATE_FILE = os.path.join(ROOT, "arxiv_watch_state.json")
LOG_FILE = os.path.join(ROOT, "arxiv_watch_log.md")
NTFY_TOPIC = "jacob-cobo-80c6d2b9e9c4"
NTFY_URL = f"https://ntfy.sh/{NTFY_TOPIC}"

QUERIES = [
    'abs:"Weil positivity"',
    'abs:"Weil quadratic form"',
    'abs:"Guinand-Weil"',
    'abs:"de Bruijn-Newman"',
    'abs:"Riemann hypothesis" AND abs:"explicit formula"',
    'abs:"Hilbert-Polya"',
    "abs:\"Li's criterion\"",
    'abs:"screw function" AND abs:"zeta"',
    'abs:"trace formula" AND abs:"Riemann hypothesis"',
    'abs:"positivity" AND abs:"critical line"',
]

NS = {"a": "http://www.w3.org/2005/Atom"}
LOW_SIGNAL_CATS = {"math.GM"}  # still reported, just flagged


def fetch(query, max_results=10):
    q = urllib.parse.quote(query)
    url = (
        f"https://export.arxiv.org/api/query?search_query={q}"
        f"&sortBy=submittedDate&sortOrder=descending&max_results={max_results}"
    )
    with urllib.request.urlopen(url, timeout=20) as resp:
        data = resp.read()
    root = ET.fromstring(data)
    out = []
    for e in root.findall("a:entry", NS):
        arxid = e.find("a:id", NS).text.rsplit("/", 1)[-1]
        title = " ".join(e.find("a:title", NS).text.split())
        summary = " ".join(e.find("a:summary", NS).text.split())
        authors = ", ".join(a.find("a:name", NS).text for a in e.findall("a:author", NS))
        cats = [c.get("term") for c in e.findall("a:category", NS)]
        published = e.find("a:published", NS).text[:10]
        out.append(
            {
                "id": arxid,
                "title": title,
                "authors": authors,
                "cats": cats,
                "published": published,
                "summary": summary,
            }
        )
    return out


def load_state():
    if os.path.exists(STATE_FILE):
        with open(STATE_FILE) as f:
            return json.load(f)
    return {"seen_ids": []}


def save_state(state):
    with open(STATE_FILE, "w") as f:
        json.dump(state, f, indent=2)


def notify(title, body):
    payload = body.encode("utf-8")
    req = urllib.request.Request(NTFY_URL, data=payload, method="POST")
    req.add_header("Title", title)
    req.add_header("Priority", "default")
    try:
        urllib.request.urlopen(req, timeout=10)
    except Exception as ex:
        print(f"ntfy failed: {ex}")


def main():
    state = load_state()
    seen = set(state["seen_ids"])
    new_papers = {}

    for q in QUERIES:
        try:
            entries = fetch(q)
        except Exception as ex:
            print(f"query failed [{q}]: {ex}")
            continue
        for e in entries:
            if e["id"] not in seen:
                new_papers[e["id"]] = e
        time.sleep(1)  # be polite to arXiv API

    if not new_papers:
        print("No new papers.")
        return

    lines = [f"\n## {time.strftime('%Y-%m-%d')} — {len(new_papers)} new paper(s)\n"]
    ntfy_lines = []
    for arxid, e in sorted(new_papers.items()):
        flag = " [GM/low-signal category]" if set(e["cats"]) & LOW_SIGNAL_CATS else ""
        lines.append(
            f"- **{e['id']}** ({e['published']}){flag} — {e['title']}\n"
            f"  Authors: {e['authors']}\n"
            f"  Cats: {', '.join(e['cats'])}\n"
            f"  https://arxiv.org/abs/{arxid}\n"
            f"  > {e['summary'][:300]}\n"
        )
        ntfy_lines.append(f"{e['id']}{flag}: {e['title']}")

    with open(LOG_FILE, "a") as f:
        f.write("\n".join(lines))

    notify(
        f"RH literature watch: {len(new_papers)} new paper(s)",
        "\n".join(ntfy_lines)[:3800],  # ntfy body size guard
    )

    state["seen_ids"] = list(seen | set(new_papers.keys()))
    save_state(state)
    print(f"Found {len(new_papers)} new paper(s), logged + notified.")


if __name__ == "__main__":
    main()
