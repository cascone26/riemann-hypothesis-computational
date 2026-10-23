#!/usr/bin/env python3
"""
scan-analyzer.py — Autonomous Riemann scan log parser
Usage: python3 scan-analyzer.py /tmp/dense_sys_AS.log [--ntfy topic] [--json out.json] [--md out.md]

Detects: BC events, ND events, blob entries/exits, zone transitions, all-time records.
Generates: structured JSON + markdown summary. Zero Claude tokens needed.
Deploy to Olympus: scp this to olympus:/tmp/scan-analyzer.py
"""

import sys, re, json, math, argparse, subprocess, os
from pathlib import Path

DEEP_THRESH = 1e-10
RC_THRESH   = 1.5e-12
ND_THRESH   = 99.0   # ALC% for near-degenerate

def alc(l1, l2):
    mn, mx = min(abs(l1), abs(l2)), max(abs(l1), abs(l2))
    return (mn / mx * 100.0) if mx > 0 else 0.0

def parse_log(path):
    rows = []
    done = False
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('DONE'):
                done = True
                continue
            # Match: "  22689  -1.282913e-12  -9.847857e-13  0s  <--RECORD-CLASS"
            m = re.match(r'\s+(\d+)\s+(-[\d.e+\-]+)\s+(-[\d.e+\-]+)\s+\d+s\s*(.*)', line)
            if not m:
                continue
            lsq  = int(m.group(1))
            l1   = float(m.group(2))
            l2   = float(m.group(3))
            tag  = m.group(4).strip()
            rows.append({'lsq': lsq, 'l1': l1, 'l2': l2, 'tag': tag})
    return rows, done

def analyze(rows):
    bc_events   = []
    nd_events   = []
    blobs       = []
    zones       = []

    # State machine
    in_blob = False
    blob_start = None
    blob_rows = []
    zone_start = None
    zone_rows = []

    bc_count = 0
    nd_count = 0

    prev = None

    for r in rows:
        lsq = r['lsq']
        l1, l2, tag = r['l1'], r['l2'], r['tag']
        a = alc(l1, l2)

        is_deep = (abs(l1) < DEEP_THRESH and abs(l2) < DEEP_THRESH)
        is_rc   = (abs(l1) < RC_THRESH   and abs(l2) < RC_THRESH)
        is_nd   = (a >= ND_THRESH)

        # BC event
        if is_rc:
            bc_count += 1
            bc_events.append({
                'bc_num': bc_count,
                'lsq': lsq,
                'l1': l1,
                'l2': l2,
                'alc': round(a, 3),
                'nd': is_nd
            })

        # ND event
        if is_nd and is_deep:
            nd_count += 1
            nd_events.append({'lsq': lsq, 'l1': l1, 'l2': l2, 'alc': round(a, 3)})

        # Blob/zone state machine
        if in_blob:
            blob_rows.append(r)
            if is_deep:
                # Blob exit
                exit_lsq = lsq
                blobs.append({
                    'start': blob_start,
                    'end': exit_lsq - 1,
                    'span': exit_lsq - blob_start,
                    'peak_l1': max(blob_rows, key=lambda x: abs(x['l1']))['l1'],
                    'peak_lsq': max(blob_rows, key=lambda x: abs(x['l1']))['lsq'],
                    'exit_drop_ratio': abs(prev['l1']) / abs(l1) if l1 != 0 else None,
                    'exit_lsq': exit_lsq,
                })
                in_blob = False
                blob_rows = []
                zone_start = lsq
                zone_rows = [r]
        else:
            # In zone
            if zone_start is not None:
                zone_rows.append(r)

            if not is_deep:
                # Blob entry
                if zone_start is not None:
                    zones.append({
                        'start': zone_start,
                        'end': lsq - 1,
                        'span': lsq - zone_start,
                        'bc_count': sum(1 for e in bc_events if zone_start <= e['lsq'] < lsq),
                        'nd_count': sum(1 for e in nd_events if zone_start <= e['lsq'] < lsq),
                        'max_alc': max((alc(x['l1'], x['l2']) for x in zone_rows), default=0),
                    })
                blob_start = lsq
                blob_rows  = [r]
                in_blob    = True
                zone_rows  = []
            else:
                if zone_start is None:
                    zone_start = lsq
                    zone_rows = [r]

        prev = r

    # Handle open blob/zone at end of log
    if in_blob and blob_rows:
        blobs.append({
            'start': blob_start,
            'end': blob_rows[-1]['lsq'],
            'span': blob_rows[-1]['lsq'] - blob_start + 1,
            'peak_l1': max(blob_rows, key=lambda x: abs(x['l1']))['l1'],
            'peak_lsq': max(blob_rows, key=lambda x: abs(x['l1']))['lsq'],
            'exit_drop_ratio': None,
            'exit_lsq': None,
            'open': True,
        })
    elif zone_start and zone_rows:
        zones.append({
            'start': zone_start,
            'end': zone_rows[-1]['lsq'],
            'span': zone_rows[-1]['lsq'] - zone_start + 1,
            'bc_count': sum(1 for e in bc_events if zone_start <= e['lsq']),
            'nd_count': sum(1 for e in nd_events if zone_start <= e['lsq']),
            'max_alc': max((alc(x['l1'], x['l2']) for x in zone_rows), default=0),
            'open': True,
        })

    return {
        'bc_events': bc_events,
        'nd_events': nd_events,
        'blobs': blobs,
        'zones': zones,
        'totals': {
            'lsq_start': rows[0]['lsq'] if rows else None,
            'lsq_end':   rows[-1]['lsq'] if rows else None,
            'bc_count':  len(bc_events),
            'nd_count':  len(nd_events),
            'blob_count': len(blobs),
            'zone_count': len(zones),
        }
    }

def fmt_e(v):
    return f'{v:.3e}'

def render_md(data, log_path):
    lines = [f'# Scan Analysis: {Path(log_path).name}', '']
    t = data['totals']
    lines += [
        f'**lsq range:** {t["lsq_start"]} – {t["lsq_end"]}',
        f'**BC events:** {t["bc_count"]}  |  **ND events:** {t["nd_count"]}',
        f'**Blobs:** {t["blob_count"]}  |  **Zones:** {t["zone_count"]}',
        '',
    ]

    if data['bc_events']:
        lines.append('## BC Events (RECORD-CLASS)')
        lines.append('| BC# | lsq | λ₁ | λ₂ | ALC% | ND? |')
        lines.append('|-----|-----|----|----|------|-----|')
        for e in data['bc_events']:
            nd = '✓' if e['nd'] else ''
            lines.append(f"| {e['bc_num']} | {e['lsq']} | {fmt_e(e['l1'])} | {fmt_e(e['l2'])} | {e['alc']:.2f}% | {nd} |")
        lines.append('')

    if data['blobs']:
        lines.append('## Blobs')
        lines.append('| # | start | end | span | peak λ₁ | peak lsq | exit ÷ |')
        lines.append('|---|-------|-----|------|---------|----------|--------|')
        for i, b in enumerate(data['blobs']):
            ex = f"÷{b['exit_drop_ratio']:.0f}" if b.get('exit_drop_ratio') else '(open)'
            lines.append(f"| {i+1} | {b['start']} | {b['end']} | {b['span']} | {fmt_e(b['peak_l1'])} | {b['peak_lsq']} | {ex} |")
        lines.append('')

    if data['zones']:
        lines.append('## Zones')
        lines.append('| # | start | end | span | BCs | NDs | max ALC% |')
        lines.append('|---|-------|-----|------|-----|-----|---------|')
        for i, z in enumerate(data['zones']):
            lines.append(f"| {i+1} | {z['start']} | {z['end']} | {z['span']} | {z['bc_count']} | {z['nd_count']} | {z['max_alc']:.2f}% |")
        lines.append('')

    return '\n'.join(lines)

def ntfy_send(topic, msg, title='Riemann Scan'):
    try:
        subprocess.run(
            ['node', os.path.expanduser('~/tools/notify/notify.js'), msg,
             '--title', title, '--priority', 'high'],
            capture_output=True, timeout=15
        )
    except Exception as e:
        print(f'ntfy failed: {e}', file=sys.stderr)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('log', help='Path to scan log file')
    ap.add_argument('--ntfy', help='ntfy topic for notifications')
    ap.add_argument('--json', help='JSON output path')
    ap.add_argument('--md',   help='Markdown output path')
    args = ap.parse_args()

    print(f'Parsing {args.log}...')
    rows, done = parse_log(args.log)
    print(f'  Rows: {len(rows)}, Complete: {done}')

    data = analyze(rows)
    t = data['totals']
    print(f'  BC: {t["bc_count"]}, ND: {t["nd_count"]}, Blobs: {t["blob_count"]}, Zones: {t["zone_count"]}')

    if args.json:
        Path(args.json).write_text(json.dumps(data, indent=2))
        print(f'  JSON → {args.json}')

    md = render_md(data, args.log)
    if args.md:
        Path(args.md).write_text(md)
        print(f'  MD   → {args.md}')
    else:
        print('\n' + md)

    if args.ntfy:
        summary = (
            f"Scan done: {Path(args.log).name}\n"
            f"lsq {t['lsq_start']}-{t['lsq_end']}\n"
            f"BC: {t['bc_count']} | ND: {t['nd_count']} | Blobs: {t['blob_count']}"
        )
        ntfy_send(args.ntfy, summary, f'Scan Complete: {Path(args.log).stem}')
        print(f'  ntfy sent → {args.ntfy}')

if __name__ == '__main__':
    main()
