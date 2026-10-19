#!/usr/bin/env python3
# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# Funded by the European Space Agency (ESA).

"""Verify fresh B4 results against archived original CSVs."""

import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXP = ROOT / "experiments/B4_dynamic_vs_fixed"

REF = EXP / "reference/B4_dt150_original"
NEW = EXP / "results/B4_dt150_fresh_v2"

RTOL = 2e-6
ATOL = 1e-8

def read(path):
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        return reader.fieldnames, list(reader)

if not REF.is_dir() or not NEW.is_dir():
    raise SystemExit("STOP: missing reference or fresh B4 directory")

reference_files = {p.name: p for p in REF.glob("B4_*.csv")}
fresh_files = {p.name: p for p in NEW.glob("B4_*.csv")}

if not reference_files or reference_files.keys() != fresh_files.keys():
    raise SystemExit(
        "STOP: reference and fresh CSV file sets differ.\n"
        f"Reference: {sorted(reference_files)}\n"
        f"Fresh: {sorted(fresh_files)}"
    )

report = {}
failed = False

for name in sorted(reference_files):
    cols_a, rows_a = read(reference_files[name])
    cols_b, rows_b = read(fresh_files[name])

    if cols_a != cols_b or len(rows_a) != len(rows_b):
        print(f"FAIL {name}: columns or row counts differ")
        failed = True
        continue

    mismatches = 0
    max_abs = 0.0
    examples = []

    for index, (a, b) in enumerate(zip(rows_a, rows_b), 1):
        for col in cols_a:
            va, vb = a[col], b[col]

            try:
                x, y = float(va), float(vb)
            except (ValueError, TypeError):
                equal = va == vb
            else:
                if math.isnan(x) and math.isnan(y):
                    equal = True
                else:
                    equal = math.isclose(
                        x, y, rel_tol=RTOL, abs_tol=ATOL
                    )
                if math.isfinite(x) and math.isfinite(y):
                    max_abs = max(max_abs, abs(x-y))

            if not equal:
                mismatches += 1
                if len(examples) < 3:
                    examples.append(f"row {index}, {col}: {va} vs {vb}")

    passed = mismatches == 0
    failed |= not passed

    report[name] = {
        "rows": len(rows_a),
        "max_absolute_difference": max_abs,
        "mismatches": mismatches,
        "passed": passed,
        "examples": examples,
    }

    print(
        f"{'PASS' if passed else 'FAIL'}: {name} | "
        f"rows={len(rows_a)} | max_abs={max_abs:.6g} | "
        f"mismatches={mismatches}"
    )

out = ROOT / "results/B4_dt150_fresh_v2_verification.json"
out.write_text(json.dumps(report, indent=2) + "\n")
print("\nVerification saved:", out)

if failed:
    raise SystemExit("FAIL: B4 numerical reproduction needs investigation")

print("\nSUCCESS: All archived B4 CSVs reproduced within tolerance.")
