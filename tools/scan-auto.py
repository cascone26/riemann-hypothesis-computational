#!/usr/bin/env python3
"""
scan-auto.py — Autonomous scan sequencer for Riemann dense_sys scans.
Monitors a running scan log; on completion: analyzes, notifies, launches next scan.
Deploy to Olympus: scp this to olympus:/tmp/scan-auto.py
Usage: nohup python3 /tmp/scan-auto.py AS AT AU AV > /tmp/scan-auto.log 2>&1 &

Sequence is: scan each label in order. For each, it:
  1. Tails /tmp/dense_sys_XX.log for "DONE XX"
  2. Runs scan-analyzer.py to extract structured data
  3. Sends ntfy notification with summary
  4. Creates /tmp/dense_sys_YY.py from prior scan script
  5. Launches nohup python3 /tmp/dense_sys_YY.py > /tmp/dense_sys_YY.log 2>&1 &
  6. Moves to next scan label

BC events during active scan are detected in real-time and notified immediately.
"""

import sys, re, time, subprocess, os
from pathlib import Path
from datetime import datetime

NTFY_TOPIC   = 'jacob-cobo-80c6d2b9e9c4'
ANALYZER     = '/tmp/scan-analyzer.py'
DB_TOOL      = '/tmp/scan-db.py'
SCAFFOLD     = '/tmp/scan-scaffold.py'
HURST        = '/tmp/scan-hurst.py'
GUE          = '/tmp/scan-gue.py'
LPPL         = '/tmp/scan-lppl.py'
SCRIPT_TMPL  = '/tmp/dense_sys_{}.py'
LOG_TMPL     = '/tmp/dense_sys_{}.log'
REPORT_TMPL  = '/tmp/dense_sys_{}_report.json'
SCAFFOLD_TMPL= '/tmp/dense_sys_{}_scaffold.md'
HURST_TMPL   = '/tmp/dense_sys_{}_hurst.json'
GUE_TMPL     = '/tmp/dense_sys_{}_gue.json'
LPPL_TMPL    = '/tmp/dense_sys_{}_lppl.json'
RC_THRESH   = 1.5e-12
DEEP_THRESH = 1e-10
POLL_SEC    = 3

def ts():
    return datetime.now().strftime('%H:%M:%S')

def log(msg):
    print(f'[{ts()}] {msg}', flush=True)

def ntfy(msg, title='Riemann'):
    try:
        subprocess.run(
            ['curl', '-s', '-d', msg, f'https://ntfy.sh/{NTFY_TOPIC}',
             '-H', f'Title: {title}', '-H', 'Priority: high'],
            capture_output=True, timeout=10
        )
    except:
        pass

def alc(l1, l2):
    mn, mx = min(abs(l1), abs(l2)), max(abs(l1), abs(l2))
    return (mn / mx * 100.0) if mx > 0 else 0.0

def make_next_script(prev_label, next_label):
    """Create next scan script from previous by substituting labels and ranges."""
    prev_script = SCRIPT_TMPL.format(prev_label)
    next_script = SCRIPT_TMPL.format(next_label)
    if not Path(prev_script).exists():
        log(f'ERROR: {prev_script} not found — cannot create {next_script}')
        return False
    content = Path(prev_script).read_text()
    # Update label references
    content = content.replace(f'{prev_label} scan', f'{next_label} scan')
    content = content.replace(f'DONE {prev_label}', f'DONE {next_label}')
    # Update range — find "range(START, END)" and shift by 500
    m = re.search(r'range\((\d+),\s*(\d+)\)', content)
    if m:
        old_start, old_end = int(m.group(1)), int(m.group(2))
        new_start = old_start + 500
        new_end   = old_end   + 500
        content = content.replace(m.group(0), f'range({new_start}, {new_end})')
        log(f'  Range: {old_start}-{old_end-1} → {new_start}-{new_end-1}')
    Path(next_script).write_text(content)
    log(f'  Created {next_script}')
    return True

def launch_scan(label):
    log_path = LOG_TMPL.format(label)
    script   = SCRIPT_TMPL.format(label)
    cmd = f'nohup python3 {script} > {log_path} 2>&1 &'
    subprocess.run(['bash', '-c', cmd])
    time.sleep(2)
    # Verify it started
    result = subprocess.run(['pgrep', '-f', f'dense_sys_{label}.py'], capture_output=True, text=True)
    if result.stdout.strip():
        log(f'  Launched {label} scan (PID {result.stdout.strip()})')
        return True
    else:
        log(f'  WARNING: Could not verify {label} scan started')
        return False

def tail_scan(label, start_line=0):
    """
    Tail a scan log. Yields (line_no, line) tuples.
    Returns True when "DONE XX" is detected, False if interrupted.
    """
    log_path = LOG_TMPL.format(label)
    bc_this_scan = 0
    line_no = 0

    while True:
        if not Path(log_path).exists():
            time.sleep(POLL_SEC)
            continue

        lines = Path(log_path).read_text().splitlines()

        for i in range(start_line, len(lines)):
            line = lines[i]
            line_no = i + 1

            if f'DONE {label}' in line:
                log(f'  DONE {label} detected at line {line_no}')
                return True, bc_this_scan

            # Real-time BC detection
            m = re.match(r'\s+(\d+)\s+(-[\d.e+\-]+)\s+(-[\d.e+\-]+)\s+\d+s\s*(.*)', line)
            if m and 'RECORD-CLASS' in m.group(4):
                lsq = m.group(1)
                l1  = float(m.group(2))
                l2  = float(m.group(3))
                bc_this_scan += 1
                a = alc(l1, l2)
                nd_flag = ' [ND!]' if a >= 99.0 else ''
                log(f'  BC at lsq {lsq}: λ₁={l1:.3e} λ₂={l2:.3e} ALC={a:.1f}%{nd_flag}')
                # ntfy only if λ₁ goes positive (disproves RH) or hits an extraordinary threshold
                if l1 > 0 or l2 > 0:
                    ntfy(f'*** POSITIVE EIGENVALUE ***\nlsq {lsq}: λ₁={l1:.3e} λ₂={l2:.3e}\nRH MAY BE DISPROVED', 'CRITICAL — Positive λ')

        start_line = len(lines)
        time.sleep(POLL_SEC)

def analyze_scan(label):
    log_path      = LOG_TMPL.format(label)
    report_path   = REPORT_TMPL.format(label)
    scaffold_path = SCAFFOLD_TMPL.format(label)

    # Step 1: JSON report
    if not Path(ANALYZER).exists():
        log(f'  WARNING: {ANALYZER} not found — skipping analysis')
        return
    result = subprocess.run(
        ['python3', ANALYZER, log_path, '--json', report_path],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        log(f'  Analyzer error: {result.stderr[:200]}')
        return
    log(f'  Analysis complete → {report_path}')

    # Step 2: Ingest into SQLite DB + update records.json (G1+G2)
    if Path(DB_TOOL).exists():
        r2 = subprocess.run(['python3', DB_TOOL, 'ingest', report_path],
                            capture_output=True, text=True)
        if r2.returncode == 0:
            log(f'  DB ingested → /tmp/riemann_events.db + records.json updated')
        else:
            log(f'  DB ingest error: {r2.stderr[:100]}')

    # Step 3: Generate paper scaffolding (G3)
    if Path(SCAFFOLD).exists():
        r3 = subprocess.run(['python3', SCAFFOLD, report_path, '--out', scaffold_path],
                            capture_output=True, text=True)
        if r3.returncode == 0:
            log(f'  Scaffold → {scaffold_path}')
        else:
            log(f'  Scaffold error: {r3.stderr[:100]}')

    # Step 4: Hurst exponent on λ₁ sequence
    hurst_path = HURST_TMPL.format(label)
    if Path(HURST).exists():
        r4 = subprocess.run(['python3', HURST, log_path, '--json', hurst_path],
                            capture_output=True, text=True)
        if r4.returncode == 0:
            log(f'  Hurst → {hurst_path}')
            # Extract H value from output
            for line in r4.stdout.splitlines():
                if 'Hurst (DFA)' in line or 'Hurst (R/S)' in line:
                    log(f'  {line.strip()}')
        else:
            log(f'  Hurst error: {r4.stderr[:100]}')

    # Step 5: GUE spacing test on BC events
    gue_path = GUE_TMPL.format(label)
    db_path = '/tmp/riemann_events.db'
    if Path(GUE).exists() and Path(db_path).exists():
        r5 = subprocess.run(['python3', GUE, db_path, '--json', gue_path],
                            capture_output=True, text=True)
        if r5.returncode == 0:
            log(f'  GUE → {gue_path}')
            for line in r5.stdout.splitlines():
                if 'Best fit' in line or 'variance' in line:
                    log(f'  {line.strip()}')
        else:
            log(f'  GUE error: {r5.stderr[:100]}')

    # Step 6: LPPL log-periodic test on blob humps
    lppl_path = LPPL_TMPL.format(label)
    if Path(LPPL).exists():
        r6 = subprocess.run(['python3', LPPL, log_path, '--json', lppl_path],
                            capture_output=True, text=True)
        if r6.returncode == 0:
            log(f'  LPPL → {lppl_path}')
            for line in r6.stdout.splitlines():
                if 'LOG-PERIODIC' in line or 'blobs show' in line:
                    log(f'  {line.strip()}')
        else:
            log(f'  LPPL error: {r6.stderr[:100]}')

    log(f'  Post-processing complete for {label}')

def main():
    labels = sys.argv[1:] if len(sys.argv) > 1 else []
    if not labels:
        print('Usage: scan-auto.py LABEL1 [LABEL2 ...] [--start-running]')
        print('Example: scan-auto.py AS AT AU AV')
        sys.exit(1)

    log(f'scan-auto.py starting. Sequence: {" → ".join(labels)}')

    for i, label in enumerate(labels):
        log_path = LOG_TMPL.format(label)
        log(f'\n=== SCAN {label} ===')

        # If log doesn't exist, need to launch (unless it's the first and presumably already running)
        if not Path(log_path).exists():
            if i == 0:
                log(f'  {log_path} not found — waiting for it to appear...')
                while not Path(log_path).exists():
                    time.sleep(5)
            else:
                # Create from previous and launch
                prev_label = labels[i-1]
                log(f'  Creating {label} scan from {prev_label}...')
                if not make_next_script(prev_label, label):
                    log(f'  ERROR: Cannot create {label} scan. Stopping.')
                    sys.exit(1)
                launch_scan(label)

        # Tail until done
        log(f'  Tailing {log_path}...')
        done, bc_count = tail_scan(label)

        if done:
            log(f'  {label} scan complete. {bc_count} BCs detected this scan.')
            analyze_scan(label)
        else:
            log(f'  {label} scan monitoring interrupted')

    log('All scans in sequence complete.')

if __name__ == '__main__':
    main()
