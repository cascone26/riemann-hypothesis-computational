"""
Fetch non-trivial zeros of the Riemann zeta function from LMFDB.
Stores them locally as high-precision decimal strings and numpy-ready floats.
"""

import urllib.request
import time
import os
import json

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
os.makedirs(DATA_DIR, exist_ok=True)

LMFDB_URL = "https://www.lmfdb.org/zeros/zeta/list"
BATCH_SIZE = 1000  # max per request
TARGET = 100_000


def fetch_batch(start_n, limit=BATCH_SIZE):
    """Fetch a batch of zeros starting at index start_n."""
    url = f"{LMFDB_URL}?N={start_n}&limit={limit}"
    req = urllib.request.Request(url, headers={"User-Agent": "RH-Research-Pipeline/1.0"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        text = resp.read().decode("utf-8")

    zeros = []
    for line in text.strip().split("\n"):
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) == 2:
            idx, value = parts
            try:
                int(idx)
                zeros.append(value)
            except ValueError:
                continue
    return zeros


def fetch_all(target=TARGET):
    """Fetch all zeros up to target count, with resume support."""
    outfile = os.path.join(DATA_DIR, "zeros_100k.txt")

    # Check for existing progress
    existing = 0
    if os.path.exists(outfile):
        with open(outfile, "r") as f:
            existing = sum(1 for line in f if line.strip())
        print(f"Found {existing} existing zeros, resuming...")

    if existing >= target:
        print(f"Already have {existing} zeros. Done.")
        return outfile

    with open(outfile, "a") as f:
        current = existing + 1
        while current <= target:
            batch_size = min(BATCH_SIZE, target - current + 1)
            print(f"Fetching zeros {current}-{current + batch_size - 1}...")

            try:
                zeros = fetch_batch(current, batch_size)
                for z in zeros:
                    f.write(z + "\n")
                f.flush()
                current += len(zeros)

                if len(zeros) < batch_size:
                    print(f"Got fewer zeros than expected ({len(zeros)} vs {batch_size}), retrying...")
                    time.sleep(2)
                    continue

                # Be nice to LMFDB servers
                time.sleep(0.5)

            except Exception as e:
                print(f"Error at N={current}: {e}")
                print("Waiting 5s and retrying...")
                time.sleep(5)

    print(f"Done. {target} zeros saved to {outfile}")
    return outfile


if __name__ == "__main__":
    fetch_all()
