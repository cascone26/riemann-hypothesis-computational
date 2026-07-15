# C=100 Analysis — Pre-Run Checklist

Run through this checklist before executing `analyze_c100.py`.

## Environment Setup (5 min)

- [ ] Python 3.8+ installed
  ```bash
  python3 --version
  ```

- [ ] mpmath installed
  ```bash
  python3 -c "import mpmath; print(mpmath.__version__)"
  ```
  *If missing: `pip install mpmath`*

- [ ] SSH configured and working
  ```bash
  ssh olympus echo "OK"
  ```
  *If fails: Check Tailscale (`ping olympus`) or SSH config*

- [ ] Tailscale connected
  ```bash
  ping olympus
  ```
  *If fails: Start Tailscale or check network*

- [ ] Script files in place
  ```bash
  ls -l /Users/scones/projects/riemann/operators/analyze_c100.py
  ```

## Remote Status (2 min)

- [ ] Result file directory exists
  ```bash
  ssh olympus ls -ld /home/cobo/riemann/results/
  ```

- [ ] Check computation status
  ```bash
  ssh olympus tail -20 /home/cobo/riemann/results/chain_log.txt
  ```

- [ ] Result file exists (optional, can run even if not ready)
  ```bash
  ssh olympus ls -lh /home/cobo/riemann/results/cvs_highprec_c100.json
  ```

## Local Output Directory (1 min)

- [ ] Findings output file exists
  ```bash
  ls -l /Users/scones/projects/riemann/results/SESSION_FINDINGS_2026_07_11.md
  ```
  *If missing: Create empty file or let script create it*

- [ ] Write permission on results directory
  ```bash
  test -w /Users/scones/projects/riemann/results/ && echo "OK" || echo "FAIL"
  ```

## Execution (1 min)

**Choose one method:**

### Method A: Shell Launcher (Recommended)
```bash
cd /Users/scones/projects/riemann/
./operators/RUN_C100_ANALYSIS.sh
```

### Method B: Direct Python
```bash
cd /Users/scones/projects/riemann/
python3 operators/analyze_c100.py
```

## Verification (2 min)

After script completes:

- [ ] Check console output
  ```
  ====================================================================
  ANALYSIS COMPLETE
  ====================================================================
  ```

- [ ] Verify findings were written
  ```bash
  tail -50 /Users/scones/projects/riemann/results/SESSION_FINDINGS_2026_07_11.md
  ```
  *Should see "§49. CvS High-Precision Analysis"*

- [ ] Check for section §49
  ```bash
  grep -n "§49" /Users/scones/projects/riemann/results/SESSION_FINDINGS_2026_07_11.md
  ```

- [ ] Review key metrics
  ```bash
  grep -A 20 "lam_min" /Users/scones/projects/riemann/results/SESSION_FINDINGS_2026_07_11.md | head -30
  ```

## Troubleshooting Decision Tree

```
Does script run?
├─ YES → Go to "Verify findings" above
└─ NO → Check error message:
    ├─ "SSH timeout" → Check Tailscale: ping olympus
    ├─ "No such file or directory" → Result not ready yet; wait or check chain_log
    ├─ "ModuleNotFoundError: mpmath" → pip install mpmath
    ├─ "Could not extract lam_min" → JSON structure unknown; debug with:
    │   ssh olympus cat /home/cobo/riemann/results/cvs_highprec_c100.json | jq . | head -50
    └─ Other error → Check ANALYZE_C100_README.md "Troubleshooting" section
```

## Timing Expectations

| Task | Time |
|------|------|
| Pre-flight checks | 5 min |
| Script execution | 5–10 sec |
| Findings review | 2 min |
| **Total** | **~8 min** |

## Quick Command Summary

```bash
# All at once:
cd /Users/scones/projects/riemann/ && \
python3 -c "import mpmath; print('✓ mpmath')" && \
ssh olympus echo "✓ SSH" && \
python3 operators/analyze_c100.py && \
echo "✓ Done. Check results/SESSION_FINDINGS_2026_07_11.md"
```

---

**Status**: Ready to execute. Script is production-ready and idempotent.
