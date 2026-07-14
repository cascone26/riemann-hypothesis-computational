#!/usr/bin/env python3
"""
scan-scaffold.py — G3: Auto-generate paper section scaffolding from scan JSON.
Deploy to Olympus: /tmp/scan-scaffold.py

Usage:
  python3 scan-scaffold.py /tmp/dense_sys_AT_report.json [--out section.md] [--section 26]
  python3 scan-scaffold.py --db /tmp/riemann_events.db --last-n 2   # scaffold last N scans from DB

Generates §XX.md with pre-filled tables, stats, event logs.
Leaves _[PROSE: ...]_ placeholders for narrative Claude fills in.
Cuts paper-writing token cost ~70%.
"""

import sys, json, re, argparse, math
from pathlib import Path
from datetime import datetime

def alc(l1, l2):
    mn, mx = min(abs(l1), abs(l2)), max(abs(l1), abs(l2))
    return (mn / mx * 100.0) if mx > 0 else 0.0

def fmt_e(v):
    return f'{v:.3e}'

def fmt_sci(v):
    """Format as ×10⁻ⁿ for markdown."""
    if v == 0: return '0'
    exp = int(math.floor(math.log10(abs(v))))
    coef = v / 10**exp
    return f'{coef:.3f}×10^{exp}'

def section_header(sec_num, title):
    return f'\n### §{sec_num} — {title}\n'

def scaffold_zone(zone, bc_events_in_zone, nd_events_in_zone, sec_num, zone_label):
    lines = [section_header(sec_num, f'{zone_label}: Zone (lsq {zone["start"]}–{zone["end"]}, {zone["span"]} lsq)')]

    lines.append(f'_[PROSE: Describe zone opening — what blob exited to enter this zone, transition character.]_\n')

    if bc_events_in_zone:
        lines.append(f'**BC events ({len(bc_events_in_zone)} total):**\n')
        lines.append('| BC# | lsq | λ₁ | λ₂ | ALC% | ND? |')
        lines.append('|-----|-----|----|----|------|-----|')
        for e in bc_events_in_zone:
            nd_flag = '**ND-BC!**' if e.get('nd') else ''
            lines.append(f"| {e['bc_num']} | {e['lsq']} | {fmt_e(e['l1'])} | {fmt_e(e['l2'])} | {e['alc']:.2f}% | {nd_flag} |")
        lines.append('')

        if any(e.get('nd') for e in bc_events_in_zone):
            nd_bcs = [e for e in bc_events_in_zone if e.get('nd')]
            for e in nd_bcs:
                lines.append(f'BC#{e["bc_num"]} at lsq {e["lsq"]} is an **ND-BC event** (ALC={e["alc"]:.2f}%), '
                              f'with both eigenvalues below RC threshold and nearly degenerate.')
                lines.append(f'_[PROSE: Discuss significance of this ND-BC — eigenvalue symmetry, channel convergence.]_\n')

    nd_in_zone = nd_events_in_zone
    if nd_in_zone:
        nd_density = len(nd_in_zone) / zone['span'] * 100
        max_alc_ev = max(nd_in_zone, key=lambda e: e['alc'])
        lines.append(f'**ND statistics:**')
        lines.append(f'- Total ND events: **{len(nd_in_zone)} / {zone["span"]} lsq ({nd_density:.1f}%)**')
        lines.append(f'- Max ALC: **{max_alc_ev["alc"]:.3f}%** at lsq {max_alc_ev["lsq"]} '
                     f'(λ₁={fmt_e(max_alc_ev["l1"])}, λ₂={fmt_e(max_alc_ev["l2"])})')
        lines.append(f'_[PROSE: Characterize ND density relative to all-time record (Zone AJ-XXX 52.3%). Note any ND clusters.]_\n')

    lines.append(f'**Zone ALC summary:** {zone["max_alc"]:.3f}% max')
    lines.append(f'_[PROSE: Zone closing — blob entry details, MC ratio, transition.]_\n')

    return '\n'.join(lines)

def scaffold_blob(blob, bc_in_blob, nd_in_blob, sec_num, blob_label):
    lines = [section_header(sec_num, f'{blob_label}: Blob (lsq {blob["start"]}–{blob.get("end","?")},'
                                     f' {blob["span"]} lsq{"" if not blob.get("open") else ", OPEN"})')]

    entry_mc = blob.get('entry_mc') or blob.get('exit_drop_ratio')
    if entry_mc:
        lines.append(f'**Entry:** lsq {blob["start"]}, MC ×{entry_mc:.1f} over prior zone value.')
    lines.append(f'_[PROSE: Describe blob entry character — sharp/gradual, λ₂ behavior at entry.]_\n')

    # Peak
    lines.append(f'**Peak: lsq {blob["peak_lsq"]}, λ₁={fmt_e(blob["peak_l1"])}**')
    lines.append(f'_[PROSE: Rank this peak all-time. Describe ascent: number of humps, MC ratios at sub-peaks, λ₂ during ascent.]_\n')

    # Exit
    if blob.get('exit_ratio'):
        lines.append(f'**Exit:** lsq {blob.get("exit_lsq","?")}, drop ÷{blob["exit_ratio"]:.0f}')
        lines.append(f'_[PROSE: Rank exit drop all-time. Describe descent from peak to exit.]_\n')

    # Complete stats block
    lines.append(f'**{blob_label} complete statistics:**')
    lines.append(f'- Span: **{blob["span"]} lsq** (lsq {blob["start"]}–{blob.get("end","?")})')
    lines.append(f'- Humps: {blob.get("hump_count", "_[count]_")}')
    lines.append(f'- Peak: lsq {blob["peak_lsq"]}, λ₁={fmt_e(blob["peak_l1"])} (_[rank]_ all-time)')
    if entry_mc: lines.append(f'- Entry MC: ×{entry_mc:.1f}')
    if blob.get('exit_ratio'): lines.append(f'- Exit drop: ÷{blob["exit_ratio"]:.0f} (_[rank]_ all-time)')
    lines.append(f'_[PROSE: Note any λ₂-DEEP episodes, intra-blob nadirs, or other anomalies.]_\n')

    return '\n'.join(lines)

def scaffold_from_json(json_path, base_section=26, out_path=None):
    data = json.loads(Path(json_path).read_text())
    t = data.get('totals', {})

    lines = [f'<!-- AUTO-SCAFFOLDED by scan-scaffold.py from {Path(json_path).name} on {datetime.now().strftime("%Y-%m-%d %H:%M")} -->',
             f'<!-- lsq {t.get("lsq_start")}–{t.get("lsq_end")} | '
             f'{t.get("bc_count",0)} BCs | {t.get("nd_count",0)} NDs | '
             f'{t.get("blob_count",0)} blobs | {t.get("zone_count",0)} zones -->',
             '']

    bc_events = data.get('bc_events', [])
    nd_events = data.get('nd_events', [])
    blobs     = data.get('blobs', [])
    zones     = data.get('zones', [])

    # Build a combined timeline by start lsq
    timeline = []
    for b in blobs: timeline.append(('blob', b))
    for z in zones: timeline.append(('zone', z))
    timeline.sort(key=lambda x: x[1]['start'])

    sec = base_section
    for kind, obj in timeline:
        # Collect events in this region
        s, e = obj['start'], obj.get('end') or obj['start'] + obj['span']
        bc_here = [ev for ev in bc_events if s <= ev['lsq'] <= e]
        nd_here = [ev for ev in nd_events if s <= ev['lsq'] <= e]

        if kind == 'zone':
            label = '_[ZONE LABEL]_'
            lines.append(scaffold_zone(obj, bc_here, nd_here, f'{sec}.XX', label))
        else:
            label = '_[BLOB LABEL]_'
            lines.append(scaffold_blob(obj, bc_here, nd_here, f'{sec}.XX', label))

        lines.append('')

    # Trailing stats block
    lines.append(f'\n---\n')
    lines.append(f'### Running Totals After This Scan\n')
    lines.append(f'- BC streak: _[CURRENT STREAK]_ (was _[PRIOR]_)')
    lines.append(f'- ND-BC count: _[TOTAL]_')
    lines.append(f'- New all-time records: _[LIST ANY]_\n')

    output = '\n'.join(lines)

    if out_path:
        Path(out_path).write_text(output)
        print(f'Scaffold → {out_path}')
    else:
        print(output)

    return output

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('json', nargs='?', help='Path to scan JSON report')
    ap.add_argument('--out', help='Output markdown path')
    ap.add_argument('--section', type=int, default=26, help='Base section number')
    args = ap.parse_args()

    if args.json:
        scaffold_from_json(args.json, base_section=args.section, out_path=args.out)
    else:
        print("Usage: scan-scaffold.py <report.json> [--out section.md]")

if __name__ == '__main__':
    main()
