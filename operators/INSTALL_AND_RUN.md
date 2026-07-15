# Installation & Execution Guide — analyze_c100.py

## Quick Start (30 seconds)

```bash
# 1. Ensure you're in the riemann project
cd /Users/scones/projects/riemann/

# 2. Ensure mpmath is installed
pip install mpmath

# 3. Run analysis (when remote results are ready)
python3 operators/analyze_c100.py
```

## What This Does

After the Olympus c=100 CvS computation completes:

1. **Fetches** results via SSH: `/home/cobo/riemann/results/cvs_highprec_c100.json`
2. **Extracts** ground-state eigenvalue (lam_min) and Riemann zeros
3. **Compares** each zero against 10 known reference values
4. **Computes** digit precision for each zero (high-precision arithmetic)
5. **Detects** outliers (zeros with precision >2σ below mean)
6. **Analyzes** k=7 zero: |sin(γ₇ × log(100) / 2)|
7. **Verifies** chain log for T=1200 and extraction job launches
8. **Appends** findings to `results/SESSION_FINDINGS_2026_07_11.md` (section §49)

## Prerequisites

### Required
- **Python 3.8+** (included on Mac)
- **mpmath** library
- **SSH access** to Olympus
- **Tailscale** connection (or network access to Olympus)

### Optional (for convenience)
- **bash** (included on Mac, for RUN_C100_ANALYSIS.sh)

## Installation

### Step 1: Install mpmath
```bash
# Using pip (recommended)
pip install mpmath

# Or using conda
conda install mpmath

# Or using uv (Jacob's preferred)
uv pip install mpmath
```

**Verify**:
```bash
python3 -c "import mpmath; print(f'✓ mpmath {mpmath.__version__}')"
```

### Step 2: Verify SSH to Olympus
```bash
# Test SSH connection
ssh olympus echo "OK"

# If that fails, check Tailscale
ping olympus

# If that fails, ensure Tailscale is running
sudo launchctl start com.tailscale.ipn.macsys
```

### Step 3: Copy Script (already done, but verify)
```bash
ls -l /Users/scones/projects/riemann/operators/analyze_c100.py
# Should show: -rwxr-xr-x ... 13K ... analyze_c100.py
```

## Usage

### Option A: Using Shell Launcher (Recommended)
```bash
cd /Users/scones/projects/riemann/
./operators/RUN_C100_ANALYSIS.sh
```

This script:
- ✓ Checks SSH access
- ✓ Checks if result file exists
- ✓ Installs mpmath if missing
- ✓ Runs the analysis
- ✓ Reports where findings were written

### Option B: Direct Python
```bash
cd /Users/scones/projects/riemann/
python3 operators/analyze_c100.py
```

This runs the analysis directly (assumes SSH and mpmath are ready).

## Execution Flow

```
[Step 1] SSH fetch remote result file
    ↓
[Step 2] Inspect JSON structure (print top-level keys)
    ↓
[Step 3] Extract lam_min and zeros
    ↓
[Step 4] Compute precision metrics (mpmath dps=200)
    ↓
[Step 5] Detect outliers (>2σ below mean)
    ↓
[Step 6] Analyze k=7 zero special case
    ↓
[Step 7] Fetch chain log for verification
    ↓
[Step 8] Format findings as markdown section §49
    ↓
[Output] Append to SESSION_FINDINGS_2026_07_11.md
    ↓
[Report] Summary metrics to console
```

**Total time**: ~5–10 seconds (network-dependent)

## Expected Output

### Console
```
[2026-07-12T01:54:30] Starting c=100 analysis...
[1/5] Fetching results from olympus:/home/cobo/riemann/results/cvs_highprec_c100.json
[2/5] Inspecting JSON structure...
JSON structure:
  lam_min: float = -1.234e-10
  zeros: list (len=10)
[3/5] Computing precision metrics...
  ✓ Analyzed 10 zeros
[4/5] Detecting outliers...
  ✓ Found 2 outliers
[5/5] Fetching chain log for verification...
[6/6] Appending to results/SESSION_FINDINGS_2026_07_11.md
  ✓ Findings appended

====================================================================
ANALYSIS COMPLETE
====================================================================
lam_min (sign): -1.234e-10 (-)
Zeros extracted: 10
Mean digit precision: 12.34
Outliers: 2
Findings written to: /Users/scones/projects/riemann/results/SESSION_FINDINGS_2026_07_11.md
====================================================================
```

### File Output
New section appended to `results/SESSION_FINDINGS_2026_07_11.md`:

```markdown
## §49. CvS High-Precision Analysis (c=100, T=800, dps=150)

**Date**: 2026-07-12T...

### lam_min
- **Value**: -1.234e-10
- **Sign**: -

### Extracted Riemann Zeros
| k | γ_extracted | γ_ref | error | digit_precision | sign |
...
```

## Monitoring the Remote Computation

While waiting for the c=100 job to finish on Olympus:

### Check computation progress
```bash
ssh olympus tail -f /home/cobo/riemann/results/chain_log.txt
```

### Check if result file exists yet
```bash
ssh olympus ls -lh /home/cobo/riemann/results/cvs_highprec_c100.json
```

### Peek at partial result (if in progress)
```bash
ssh olympus head -c 500 /home/cobo/riemann/results/cvs_highprec_c100.json
```

## Error Handling

| Error | Diagnosis | Fix |
|-------|-----------|-----|
| `SSH timeout` | Tailscale issue | `ping olympus` and `ssh olympus echo OK` |
| `No such file or directory` | Result file not ready | Wait for computation; check chain log |
| `Invalid JSON` | Corrupted result file | Rerun computation on Olympus |
| `ModuleNotFoundError: mpmath` | mpmath not installed | `pip install mpmath` |
| `Could not extract lam_min` | JSON structure unexpected | Print JSON: `ssh olympus cat ... \| jq . \| head -50` |

## Files

| File | Size | Purpose |
|------|------|---------|
| `analyze_c100.py` | 13 KB | Main analysis engine (executable) |
| `RUN_C100_ANALYSIS.sh` | 2.3 KB | Shell launcher with pre-flight checks |
| `ANALYZE_C100_README.md` | ~5 KB | Comprehensive documentation |
| `C100_ANALYSIS_SUITE.md` | ~8 KB | Overview and integration guide |
| `INSTALL_AND_RUN.md` | This file | Installation & execution guide |

**All files located at**: `/Users/scones/projects/riemann/operators/`

## Next Steps

1. **Wait for remote computation** to complete (monitor via `ssh olympus tail -f /home/cobo/riemann/results/chain_log.txt`)
2. **Run analysis** when ready: `./operators/RUN_C100_ANALYSIS.sh`
3. **Review findings** in `results/SESSION_FINDINGS_2026_07_11.md` (section §49)
4. **Archive** findings if needed (markdown is self-contained, easy to copy/paste)

## Support

For issues:
1. Check **Troubleshooting** section above
2. Read **ANALYZE_C100_README.md** for detailed docs
3. Review **C100_ANALYSIS_SUITE.md** for architecture overview

---

**Ready to run!** The script is production-ready and idempotent.
You can safely re-run it multiple times (findings will append, not duplicate).
