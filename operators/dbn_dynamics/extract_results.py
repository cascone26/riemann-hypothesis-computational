#!/usr/bin/env python3
"""Extract certification results from output."""
import sys
import re

def extract_results(output_text):
    """Parse certification output and extract results."""
    lines = output_text.split('\n')

    results = {
        'certified': False,
        'runtime': None,
        'pivots': [],
        'all_positive': True,
    }

    for line in lines:
        # Extract runtime
        if 'runtime=' in line:
            m = re.search(r'runtime=([0-9.]+)s', line)
            if m:
                results['runtime'] = float(m.group(1))

        # Extract pivot lines
        if 'pivot' in line and ':' in line:
            # Format: "  pivot i: m.mmmme+ee +/- r.rre-ee  lower=l.lllle+ee  >0 OK"
            m = re.search(r'pivot (\d+):[^l]*lower=([0-9.e+-]+)\s+([>]0 OK|FAIL)', line)
            if m:
                idx = int(m.group(1))
                lower_bound = float(m.group(2))
                status = m.group(3)
                results['pivots'].append({
                    'index': idx,
                    'lower_bound': lower_bound,
                    'status': status,
                    'positive': status == '>0 OK'
                })
                if status != '>0 OK':
                    results['all_positive'] = False

        # Check for certification
        if 'CERTIFIED PD' in line:
            results['certified'] = True
        elif 'not certified' in line:
            results['certified'] = False

    return results

if __name__ == '__main__':
    if len(sys.argv) > 1:
        with open(sys.argv[1]) as f:
            output = f.read()
    else:
        output = sys.stdin.read()

    results = extract_results(output)
    print(f"Certified: {results['certified']}")
    print(f"Runtime: {results['runtime']}s")
    print(f"All pivots positive: {results['all_positive']}")
    print(f"Pivots:")
    for p in sorted(results['pivots'], key=lambda x: x['index']):
        print(f"  [{p['index']}] {p['lower_bound']:.6e} ({p['status']})")
