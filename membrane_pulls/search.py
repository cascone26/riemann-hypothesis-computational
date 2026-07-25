#!/usr/bin/env python3
"""Full-text search over the pulled prolate/explicit-formula literature (room 2026-07-25).
GreenBot's grep for the edge-asymptotics law in 2310.18423 came up empty at 02-ish AM because
the paper doesn't contain it -- the real source (Fuchs '64 / Widom / des Cloizeaux-Mehta) isn't
pulled yet, and there was no fast way to tell "not in our corpus" from "not searched right."
This makes "who proved X" a query against what's ACTUALLY here (Charisma: Report-tier, since
it's text search, not verified proof-reading) instead of a dead grep against one file.

Build/refresh the index (run after adding new pulls):
  python3 search.py index
Query:
  python3 search.py query "prolate spheroidal asymptotics"
  python3 search.py query "edge asymptotics" --context 200
List what's actually indexed (so "not found" can be trusted):
  python3 search.py sources
"""
import argparse, os, re, sqlite3, sys, tarfile

ROOT = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(ROOT, "corpus.db")
EXTS = (".tex", ".txt")


def iter_source_files():
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [d for d in dirnames if d not in ("__pycache__",)]
        for fn in filenames:
            if fn.endswith(EXTS):
                yield os.path.join(dirpath, fn)


def iter_tarballs():
    for fn in os.listdir(ROOT):
        if fn.endswith((".tar.gz", ".tgz")):
            yield os.path.join(ROOT, fn)


def build_index():
    con = sqlite3.connect(DB)
    con.execute("DROP TABLE IF EXISTS corpus")
    con.execute("""CREATE VIRTUAL TABLE corpus USING fts5(
        source, path, para_id UNINDEXED, body, tokenize='porter'
    )""")

    n_docs = 0
    n_paras = 0
    seen_basenames = set()

    def index_text(source, path, text):
        nonlocal n_paras
        # split on blank lines -> paragraph-ish chunks; keeps FTS snippet context tight
        paras = [p.strip() for p in re.split(r"\n\s*\n", text) if p.strip()]
        for i, p in enumerate(paras):
            if len(p) < 20:
                continue
            con.execute("INSERT INTO corpus (source, path, para_id, body) VALUES (?,?,?,?)",
                        (source, path, i, p[:8000]))
            n_paras += 1

    for path in iter_source_files():
        rel = os.path.relpath(path, ROOT)
        try:
            text = open(path, encoding="utf-8", errors="ignore").read()
        except Exception as e:
            print(f"  skip {rel}: {e}")
            continue
        index_text(rel, rel, text)
        seen_basenames.add(os.path.basename(path))
        n_docs += 1
        print(f"  indexed {rel} ({len(text)} chars)")

    for tb in iter_tarballs():
        rel = os.path.relpath(tb, ROOT)
        try:
            with tarfile.open(tb) as tf:
                for member in tf.getmembers():
                    if member.isfile() and member.name.endswith(EXTS):
                        basename = os.path.basename(member.name)
                        if basename in seen_basenames:
                            continue  # already covered by an extracted copy
                        try:
                            text = tf.extractfile(member).read().decode("utf-8", errors="ignore")
                        except Exception:
                            continue
                        inner = f"{rel}::{member.name}"
                        index_text(inner, inner, text)
                        seen_basenames.add(basename)
                        n_docs += 1
        except Exception as e:
            print(f"  skip tarball {rel}: {e}")

    con.commit()
    con.close()
    print(f"[index] {n_docs} source docs, {n_paras} paragraphs indexed -> {DB}")


def query(q, context, limit):
    if not os.path.exists(DB):
        sys.exit("no index yet -- run: python3 search.py index")
    con = sqlite3.connect(DB)
    rows = con.execute(
        """SELECT source, para_id, snippet(corpus, 3, '>>>', '<<<', ' ... ', ?) , rank
           FROM corpus WHERE corpus MATCH ? ORDER BY rank LIMIT ?""",
        (max(10, context // 8), q, limit)
    ).fetchall()
    if not rows:
        print(f"NO MATCHES for {q!r} in the current corpus.")
        print("Run 'python3 search.py sources' to see what's actually indexed --")
        print("a miss here means 'not in our pulled papers,' not 'doesn't exist.'")
        return
    for source, para_id, snippet, rank in rows:
        print(f"\n[{source}] para {para_id} (rank {rank:.2f})")
        print("  " + snippet.replace("\n", " "))


def sources():
    if not os.path.exists(DB):
        sys.exit("no index yet -- run: python3 search.py index")
    con = sqlite3.connect(DB)
    rows = con.execute("SELECT source, COUNT(*) FROM corpus GROUP BY source ORDER BY source").fetchall()
    for source, n in rows:
        print(f"  {n:>5} paras  {source}")
    print(f"\n{len(rows)} sources indexed.")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("index")
    p = sub.add_parser("query")
    p.add_argument("q")
    p.add_argument("--context", type=int, default=200)
    p.add_argument("--limit", type=int, default=15)
    sub.add_parser("sources")
    a = ap.parse_args()
    if a.cmd == "index":
        build_index()
    elif a.cmd == "query":
        query(a.q, a.context, a.limit)
    elif a.cmd == "sources":
        sources()


if __name__ == "__main__":
    main()
