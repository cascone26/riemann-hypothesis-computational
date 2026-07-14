#!/usr/bin/env python3
"""
scan-db.py — G1: SQLite events database for Riemann dense_sys scans.
Deploy to Olympus: /tmp/scan-db.py

Usage:
  python3 scan-db.py ingest /tmp/dense_sys_AS_report.json     # add JSON report to DB
  python3 scan-db.py ingest-log /tmp/dense_sys_AS.log         # parse log directly into DB
  python3 scan-db.py records                                   # show all-time records
  python3 scan-db.py bc                                        # list all BC events
  python3 scan-db.py spacing                                   # BC inter-arrival gaps (C1)
  python3 scan-db.py alc [--min 99]                           # ND events above ALC threshold
  python3 scan-db.py blobs                                     # blob summary table
  python3 scan-db.py zones                                     # zone summary table
  python3 scan-db.py predict                                   # next blob/zone ETA (G4)
  python3 scan-db.py drift                                     # blob peak trend (C5)
  python3 scan-db.py next-bc                                   # AR(5) next-BC estimate (A1)

DB lives at: /tmp/riemann_events.db (persistent across scans)
"""

import sys, json, sqlite3, re, math
from pathlib import Path
from datetime import datetime

DB_PATH       = '/tmp/riemann_events.db'
RECORDS_PATH  = '/tmp/riemann_records.json'
DEEP_THRESH   = 1e-10
RC_THRESH     = 1.5e-12
ND_THRESH     = 99.0

def get_db():
    db = sqlite3.connect(DB_PATH)
    db.row_factory = sqlite3.Row
    db.executescript("""
        CREATE TABLE IF NOT EXISTS events (
            lsq INTEGER PRIMARY KEY,
            l1 REAL, l2 REAL, alc REAL,
            tag TEXT, region TEXT,
            ingested_from TEXT
        );
        CREATE TABLE IF NOT EXISTS bc_events (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            lsq INTEGER UNIQUE, l1 REAL, l2 REAL, alc REAL, nd INTEGER,
            scan_label TEXT
        );
        CREATE TABLE IF NOT EXISTS nd_events (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            lsq INTEGER UNIQUE, l1 REAL, l2 REAL, alc REAL,
            scan_label TEXT
        );
        CREATE TABLE IF NOT EXISTS blobs (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            start_lsq INTEGER UNIQUE, end_lsq INTEGER, span INTEGER,
            peak_l1 REAL, peak_lsq INTEGER,
            entry_mc REAL, exit_ratio REAL,
            hump_count INTEGER, open INTEGER DEFAULT 0,
            scan_label TEXT
        );
        CREATE TABLE IF NOT EXISTS zones (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            start_lsq INTEGER UNIQUE, end_lsq INTEGER, span INTEGER,
            bc_count INTEGER, nd_count INTEGER, max_alc REAL,
            open INTEGER DEFAULT 0,
            scan_label TEXT
        );
        CREATE TABLE IF NOT EXISTS ingested_scans (
            label TEXT PRIMARY KEY,
            ingested_at TEXT,
            lsq_start INTEGER, lsq_end INTEGER
        );
    """)
    db.commit()
    return db

def alc(l1, l2):
    mn, mx = min(abs(l1), abs(l2)), max(abs(l1), abs(l2))
    return (mn / mx * 100.0) if mx > 0 else 0.0

def parse_log_rows(path):
    rows = []
    with open(path) as f:
        for line in f:
            m = re.match(r'\s+(\d+)\s+(-[\d.e+\-]+)\s+(-[\d.e+\-]+)\s+\d+s\s*(.*)', line.rstrip())
            if not m: continue
            lsq = int(m.group(1))
            l1, l2 = float(m.group(2)), float(m.group(3))
            tag = m.group(4).strip()
            rows.append({'lsq': lsq, 'l1': l1, 'l2': l2, 'tag': tag,
                         'alc': alc(l1, l2)})
    return rows

def ingest_json(db, json_path, label=None):
    data = json.loads(Path(json_path).read_text())
    if label is None:
        label = Path(json_path).stem.replace('_report', '').replace('dense_sys_', '')

    if db.execute("SELECT 1 FROM ingested_scans WHERE label=?", (label,)).fetchone():
        print(f"  {label} already ingested — skipping")
        return

    for e in data.get('bc_events', []):
        db.execute("""INSERT OR IGNORE INTO bc_events (lsq,l1,l2,alc,nd,scan_label)
                      VALUES (?,?,?,?,?,?)""",
                   (e['lsq'], e['l1'], e['l2'], e['alc'], int(e.get('nd', False)), label))

    for e in data.get('nd_events', []):
        db.execute("""INSERT OR IGNORE INTO nd_events (lsq,l1,l2,alc,scan_label)
                      VALUES (?,?,?,?,?)""",
                   (e['lsq'], e['l1'], e['l2'], e['alc'], label))

    for b in data.get('blobs', []):
        db.execute("""INSERT OR IGNORE INTO blobs
                      (start_lsq,end_lsq,span,peak_l1,peak_lsq,entry_mc,exit_ratio,hump_count,open,scan_label)
                      VALUES (?,?,?,?,?,?,?,?,?,?)""",
                   (b['start'], b['end'], b['span'], b['peak_l1'], b['peak_lsq'],
                    b.get('entry_mc'), b.get('exit_drop_ratio'), b.get('hump_count', 1),
                    int(b.get('open', False)), label))

    for z in data.get('zones', []):
        db.execute("""INSERT OR IGNORE INTO zones
                      (start_lsq,end_lsq,span,bc_count,nd_count,max_alc,open,scan_label)
                      VALUES (?,?,?,?,?,?,?,?)""",
                   (z['start'], z['end'], z['span'], z['bc_count'], z['nd_count'],
                    z['max_alc'], int(z.get('open', False)), label))

    t = data.get('totals', {})
    db.execute("""INSERT OR REPLACE INTO ingested_scans (label,ingested_at,lsq_start,lsq_end)
                  VALUES (?,?,?,?)""",
               (label, datetime.now().isoformat(), t.get('lsq_start'), t.get('lsq_end')))
    db.commit()
    print(f"  Ingested {label}: {t.get('bc_count',0)} BCs, {t.get('nd_count',0)} NDs, "
          f"{t.get('blob_count',0)} blobs, {t.get('zone_count',0)} zones")

def ingest_log_direct(db, log_path):
    label = Path(log_path).stem.replace('dense_sys_', '')
    if db.execute("SELECT 1 FROM ingested_scans WHERE label=?", (label,)).fetchone():
        print(f"  {label} already ingested — skipping")
        return

    # Run scan-analyzer to get JSON first
    import subprocess, tempfile
    tmp = tempfile.NamedTemporaryFile(suffix='.json', delete=False)
    tmp.close()
    result = subprocess.run(['python3', '/tmp/scan-analyzer.py', log_path, '--json', tmp.name],
                            capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  Analyzer failed: {result.stderr[:200]}")
        return
    ingest_json(db, tmp.name, label)
    Path(tmp.name).unlink()

def cmd_records(db):
    print("\n=== ALL-TIME RECORDS ===")
    r = {}

    row = db.execute("SELECT lsq, l1, l2 FROM bc_events ORDER BY ABS(l1) ASC LIMIT 1").fetchone()
    if row: r['min_l1'] = f"λ₁={row['l1']:.3e} at lsq {row['lsq']}"

    row = db.execute("SELECT MAX(id) as streak FROM bc_events").fetchone()
    if row: r['bc_streak'] = f"{row['streak']} total BCs"

    row = db.execute("SELECT lsq, alc FROM nd_events ORDER BY alc DESC LIMIT 1").fetchone()
    if row: r['max_alc'] = f"{row['alc']:.4f}% at lsq {row['lsq']}"

    row = db.execute("SELECT COUNT(*) as cnt FROM nd_events WHERE alc >= 99").fetchone()
    if row: r['nd_bc_count'] = db.execute("SELECT COUNT(*) as cnt FROM bc_events WHERE nd=1").fetchone()['cnt']

    row = db.execute("SELECT start_lsq, span, peak_l1 FROM blobs WHERE open=0 ORDER BY ABS(peak_l1) DESC LIMIT 1").fetchone()
    if row: r['max_peak'] = f"λ₁={row['peak_l1']:.3e} (blob start {row['start_lsq']})"

    row = db.execute("SELECT start_lsq, span FROM blobs WHERE open=0 ORDER BY span DESC LIMIT 1").fetchone()
    if row: r['longest_blob'] = f"{row['span']} lsq (start {row['start_lsq']})"

    row = db.execute("SELECT start_lsq, exit_ratio FROM blobs WHERE exit_ratio IS NOT NULL ORDER BY exit_ratio DESC LIMIT 1").fetchone()
    if row: r['max_exit'] = f"÷{row['exit_ratio']:.0f} (blob start {row['start_lsq']})"

    row = db.execute("SELECT start_lsq, bc_count, span FROM zones ORDER BY bc_count DESC LIMIT 1").fetchone()
    if row: r['max_bc_zone'] = f"{row['bc_count']} BCs in {row['span']} lsq (start {row['start_lsq']})"

    row = db.execute("SELECT start_lsq, nd_count, span FROM zones WHERE span > 0 ORDER BY CAST(nd_count AS REAL)/span DESC LIMIT 1").fetchone()
    if row: r['max_nd_density'] = f"{row['nd_count']}/{row['span']} = {row['nd_count']/row['span']*100:.1f}% (start {row['start_lsq']})"

    for k, v in r.items():
        print(f"  {k:20s}: {v}")
    update_records_json(db)

def update_records_json(db):
    records = {}
    row = db.execute("SELECT lsq, l1, l2 FROM bc_events ORDER BY ABS(l1) ASC LIMIT 1").fetchone()
    if row: records['min_l1'] = {'value': row['l1'], 'lsq': row['lsq']}

    row = db.execute("SELECT COUNT(*) as cnt FROM bc_events").fetchone()
    if row: records['bc_total'] = row['cnt']

    row = db.execute("SELECT COUNT(*) as cnt FROM bc_events WHERE nd=1").fetchone()
    if row: records['nd_bc_count'] = row['cnt']

    row = db.execute("SELECT lsq, alc FROM nd_events ORDER BY alc DESC LIMIT 1").fetchone()
    if row: records['max_alc'] = {'value': row['alc'], 'lsq': row['lsq']}

    rows = db.execute("SELECT start_lsq, span, peak_l1, peak_lsq FROM blobs WHERE open=0 ORDER BY ABS(peak_l1) DESC LIMIT 10").fetchall()
    records['top10_peaks'] = [{'start': r['start_lsq'], 'peak_l1': r['peak_l1'], 'peak_lsq': r['peak_lsq']} for r in rows]

    row = db.execute("SELECT start_lsq, span FROM blobs ORDER BY span DESC LIMIT 1").fetchone()
    if row: records['longest_blob'] = {'start': row['start_lsq'], 'span': row['span']}

    row = db.execute("SELECT start_lsq, exit_ratio FROM blobs WHERE exit_ratio IS NOT NULL ORDER BY exit_ratio DESC LIMIT 1").fetchone()
    if row: records['max_exit_ratio'] = {'start': row['start_lsq'], 'ratio': row['exit_ratio']}

    row = db.execute("SELECT start_lsq, bc_count, span FROM zones ORDER BY bc_count DESC LIMIT 1").fetchone()
    if row: records['max_bc_zone'] = {'start': row['start_lsq'], 'bc_count': row['bc_count'], 'span': row['span']}

    records['updated'] = datetime.now().isoformat()
    Path(RECORDS_PATH).write_text(json.dumps(records, indent=2))
    print(f"  records.json updated → {RECORDS_PATH}")

def cmd_bc(db):
    rows = db.execute("SELECT * FROM bc_events ORDER BY lsq").fetchall()
    print(f"\n{'BC#':>4} {'lsq':>6} {'λ₁':>14} {'λ₂':>14} {'ALC%':>8} {'ND':>4}")
    print("-" * 58)
    for i, r in enumerate(rows, 1):
        nd = '✓' if r['nd'] else ''
        print(f"{i:>4} {r['lsq']:>6} {r['l1']:>14.3e} {r['l2']:>14.3e} {r['alc']:>8.3f}% {nd:>4}")
    print(f"\nTotal: {len(rows)} BC events")

def cmd_spacing(db):
    rows = db.execute("SELECT lsq FROM bc_events ORDER BY lsq").fetchall()
    lsqs = [r['lsq'] for r in rows]
    if len(lsqs) < 2:
        print("Need ≥2 BC events"); return
    gaps = [lsqs[i+1]-lsqs[i] for i in range(len(lsqs)-1)]
    avg = sum(gaps)/len(gaps)
    variance = sum((g-avg)**2 for g in gaps)/len(gaps)
    sd = math.sqrt(variance)
    print(f"\n=== BC INTER-ARRIVAL SPACING ===")
    print(f"  Events:  {len(lsqs)}")
    print(f"  Mean gap: {avg:.1f} lsq")
    print(f"  Std dev:  {sd:.1f} lsq")
    print(f"  Min gap:  {min(gaps)} lsq")
    print(f"  Max gap:  {max(gaps)} lsq")
    # Histogram (text)
    buckets = [0]*10
    bw = max(1, int((max(gaps)-min(gaps))/10))
    for g in gaps:
        b = min(9, (g - min(gaps)) // bw)
        buckets[b] += 1
    print(f"\n  Gap histogram (bucket width ≈{bw} lsq):")
    for i, cnt in enumerate(buckets):
        lo = min(gaps) + i*bw
        bar = '█' * (cnt * 30 // max(1, max(buckets)))
        print(f"  {lo:4d}-{lo+bw-1:4d}: {bar} {cnt}")

def cmd_alc(db, min_alc=99.0):
    rows = db.execute("SELECT * FROM nd_events WHERE alc >= ? ORDER BY alc DESC LIMIT 30", (min_alc,)).fetchall()
    print(f"\n=== ND EVENTS (ALC ≥ {min_alc}%) — top 30 ===")
    print(f"{'lsq':>6} {'λ₁':>14} {'λ₂':>14} {'ALC%':>10}")
    print("-" * 50)
    for r in rows:
        print(f"{r['lsq']:>6} {r['l1']:>14.3e} {r['l2']:>14.3e} {r['alc']:>10.4f}%")

def cmd_blobs(db):
    rows = db.execute("SELECT * FROM blobs ORDER BY start_lsq").fetchall()
    print(f"\n{'#':>3} {'start':>7} {'end':>7} {'span':>5} {'peak λ₁':>14} {'peak lsq':>9} {'exit÷':>8}")
    print("-" * 60)
    for i, r in enumerate(rows, 1):
        ex = f"÷{r['exit_ratio']:.0f}" if r['exit_ratio'] else '(open)'
        print(f"{i:>3} {r['start_lsq']:>7} {r['end_lsq'] or '?':>7} {r['span']:>5} "
              f"{r['peak_l1']:>14.3e} {r['peak_lsq']:>9} {ex:>8}")

def cmd_zones(db):
    rows = db.execute("SELECT * FROM zones ORDER BY start_lsq").fetchall()
    print(f"\n{'#':>3} {'start':>7} {'end':>7} {'span':>5} {'BCs':>4} {'NDs':>4} {'maxALC%':>8}")
    print("-" * 50)
    for i, r in enumerate(rows, 1):
        print(f"{i:>3} {r['start_lsq']:>7} {r['end_lsq'] or '?':>7} {r['span']:>5} "
              f"{r['bc_count']:>4} {r['nd_count']:>4} {r['max_alc']:>8.3f}%")

def cmd_predict(db):
    """G4: Predict next blob/zone entry based on historical span statistics."""
    blob_spans = [r['span'] for r in db.execute("SELECT span FROM blobs WHERE open=0").fetchall()]
    zone_spans = [r['span'] for r in db.execute("SELECT span FROM zones WHERE open=0").fetchall()]

    def stats(arr):
        if not arr: return None, None
        avg = sum(arr)/len(arr)
        sd = math.sqrt(sum((x-avg)**2 for x in arr)/len(arr))
        return avg, sd

    b_avg, b_sd = stats(blob_spans)
    z_avg, z_sd = stats(zone_spans)

    last_lsq = db.execute("SELECT MAX(lsq) FROM bc_events").fetchone()[0] or 0

    print(f"\n=== ZONE/BLOB PREDICTOR (G4) ===")
    if b_avg: print(f"  Blob span: mean={b_avg:.1f} lsq  σ={b_sd:.1f}  n={len(blob_spans)}")
    if z_avg: print(f"  Zone span: mean={z_avg:.1f} lsq  σ={z_sd:.1f}  n={len(zone_spans)}")

    open_blob = db.execute("SELECT start_lsq FROM blobs WHERE open=1 ORDER BY start_lsq DESC LIMIT 1").fetchone()
    open_zone = db.execute("SELECT start_lsq FROM zones WHERE open=1 ORDER BY start_lsq DESC LIMIT 1").fetchone()

    if open_blob and b_avg:
        expected_end = open_blob['start_lsq'] + int(b_avg)
        print(f"\n  Current open blob started at {open_blob['start_lsq']}")
        print(f"  Expected blob exit: lsq ~{expected_end} (±{b_sd:.0f})")
        if z_avg:
            print(f"  Expected zone entry after that: lsq ~{expected_end + int(z_avg)}")
    elif open_zone and z_avg:
        expected_end = open_zone['start_lsq'] + int(z_avg)
        print(f"\n  Current open zone started at {open_zone['start_lsq']}")
        print(f"  Expected zone exit (blob entry): lsq ~{expected_end} (±{z_sd:.0f})")

def cmd_drift(db):
    """C5: Linear regression on blob peak λ₁ vs lsq — is the spectrum drifting?"""
    rows = db.execute("SELECT peak_lsq, peak_l1 FROM blobs WHERE open=0 ORDER BY peak_lsq").fetchall()
    if len(rows) < 3:
        print("Need ≥3 complete blobs for drift analysis")
        return

    xs = [r['peak_lsq'] for r in rows]
    ys = [abs(r['peak_l1']) for r in rows]
    log_ys = [math.log10(y) for y in ys]  # log scale

    n = len(xs)
    sx = sum(xs); sy = sum(log_ys)
    sxx = sum(x*x for x in xs); sxy = sum(x*y for x, y in zip(xs, log_ys))
    slope = (n*sxy - sx*sy) / (n*sxx - sx*sx)
    intercept = (sy - slope*sx) / n

    print(f"\n=== BLOB PEAK DRIFT ANALYSIS (C5) ===")
    print(f"  Blobs: {n} complete")
    print(f"  log₁₀(|λ₁ peak|) = {slope:.6f} × lsq + {intercept:.3f}")
    print(f"  Slope: {slope:+.2e} per lsq")
    if slope < 0:
        print(f"  Trend: blob peaks are DECREASING — {abs(slope)*1000:.4f} decades per 1000 lsq")
    elif slope > 0:
        print(f"  Trend: blob peaks are INCREASING — {slope*1000:.4f} decades per 1000 lsq")
    else:
        print(f"  Trend: FLAT — no significant drift")

    # Show top and bottom peaks
    sorted_by_peak = sorted(rows, key=lambda r: abs(r['peak_l1']), reverse=True)
    print(f"\n  Top 5 blob peaks:")
    for r in sorted_by_peak[:5]:
        print(f"    lsq {r['peak_lsq']}: λ₁={r['peak_l1']:.3e}")

def cmd_next_bc(db):
    """A1: AR(5) model — predict approximate lsq of next BC event."""
    rows = db.execute("SELECT lsq FROM bc_events ORDER BY lsq").fetchall()
    lsqs = [r['lsq'] for r in rows]
    if len(lsqs) < 6:
        print("Need ≥6 BC events for AR(5) model")
        return

    gaps = [lsqs[i+1]-lsqs[i] for i in range(len(lsqs)-1)]
    # Simple AR(5): predict next gap as weighted average of last 5 gaps
    last5 = gaps[-5:]
    weights = [0.1, 0.15, 0.2, 0.25, 0.3]  # more recent = higher weight
    predicted_gap = sum(w*g for w, g in zip(weights, last5))
    next_lsq = lsqs[-1] + predicted_gap

    avg_gap = sum(gaps)/len(gaps)
    print(f"\n=== NEXT-BC PREDICTOR (AR-5) ===")
    print(f"  Last BC: lsq {lsqs[-1]}")
    print(f"  Last 5 gaps: {last5}")
    print(f"  Predicted next gap: {predicted_gap:.1f} lsq")
    print(f"  Predicted next BC: ~lsq {next_lsq:.0f}")
    print(f"  (Long-run avg gap: {avg_gap:.1f} lsq)")

def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(0)

    cmd = sys.argv[1]
    db = get_db()

    if cmd == 'ingest':
        for path in sys.argv[2:]:
            print(f"Ingesting {path}...")
            ingest_json(db, path)
        cmd_records(db)

    elif cmd == 'ingest-log':
        for path in sys.argv[2:]:
            print(f"Ingesting log {path}...")
            ingest_log_direct(db, path)
        cmd_records(db)

    elif cmd == 'ingest-all':
        # Ingest all existing JSON reports on Olympus
        import glob
        reports = sorted(glob.glob('/tmp/dense_sys_*_report.json'))
        print(f"Found {len(reports)} reports to ingest...")
        for r in reports:
            ingest_json(db, r)
        cmd_records(db)

    elif cmd == 'records':   cmd_records(db)
    elif cmd == 'bc':        cmd_bc(db)
    elif cmd == 'spacing':   cmd_spacing(db)
    elif cmd == 'blobs':     cmd_blobs(db)
    elif cmd == 'zones':     cmd_zones(db)
    elif cmd == 'predict':   cmd_predict(db)
    elif cmd == 'drift':     cmd_drift(db)
    elif cmd == 'next-bc':   cmd_next_bc(db)
    elif cmd == 'alc':
        min_alc = float(sys.argv[2]) if len(sys.argv) > 2 else ND_THRESH
        cmd_alc(db, min_alc)
    else:
        print(f"Unknown command: {cmd}")
        sys.exit(1)

if __name__ == '__main__':
    main()
