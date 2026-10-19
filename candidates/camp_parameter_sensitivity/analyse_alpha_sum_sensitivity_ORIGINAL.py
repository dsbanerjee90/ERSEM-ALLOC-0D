#!/usr/bin/env python3
# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# ERSEM-ALLOC / ESA-CAMP: exploratory P1 parameter sensitivity.

from pathlib import Path
from netCDF4 import Dataset
import numpy as np
import csv

ROOT = Path("/work/dba/ERSEM-ALLOC-0D")
EXP = ROOT / "experiments"
OUT = ROOT / "paper/analysis/alpha_sum_sensitivity_pilot.csv"

CASES = {
    "No external nitrate supply": {
        "reference": "14_ALLOC_LOW_INITIAL_N3",
        "sum": (
            "18_SENS_NOSUPPLY_SUM_M10",
            "19_SENS_NOSUPPLY_SUM_P10"
        ),
        "alpha": (
            "22_SENS_NOSUPPLY_ALPHA_M10",
            "23_SENS_NOSUPPLY_ALPHA_P10"
        ),
    },
    "External nitrate supply": {
        "reference": "16_ALLOC_N3_SUPPLY",
        "sum": (
            "20_SENS_SUPPLY_SUM_M10",
            "21_SENS_SUPPLY_SUM_P10"
        ),
        "alpha": (
            "24_SENS_SUPPLY_ALPHA_M10",
            "25_SENS_SUPPLY_ALPHA_P10"
        ),
    },
}


def diagnostics(case):
    folder = EXP / case

    with Dataset(folder / "output.nc") as ds:

        def read(name):
            x = np.asarray(
                np.ma.filled(ds[name][:], np.nan),
                dtype=float
            ).ravel()
            assert np.all(np.isfinite(x)), (case, name)
            return x

        t = read("time")
        c = read("P1_c")
        gross = read("P1_fO3PIc")
        nitrate = read("P1_fN3PIn")

        assert np.all(np.diff(t) > 0), case
        assert np.all(c > 0), case
        assert np.isclose(t[-1] - t[0], 7 * 86400), case

        gross_specific = np.trapezoid(gross / c, t) / np.ptp(t)
        nitrate_specific = np.trapezoid(nitrate / c, t) / np.ptp(t)

    with open(folder / "allocation_trajectory.csv", newline="") as f:
        rows = list(csv.DictReader(f))

    tf = np.array([float(r["time_seconds"]) for r in rows])
    fc = np.array([float(r["P1_fC"]) for r in rows])

    assert np.all(np.diff(tf) > 0), case
    assert np.all(np.isfinite(fc)), case
    assert np.isclose(tf[-1] - tf[0], 7 * 86400), case

    mean_fc = np.trapezoid(fc, tf) / np.ptp(tf)

    return {
        "Mean fC": mean_fc,
        "Gross fixation / C": gross_specific,
        "Nitrate uptake / C": nitrate_specific,
    }


# Read each case only once.
all_cases = set()

for conditions in CASES.values():
    all_cases.add(conditions["reference"])
    for parameter in ("alpha", "sum"):
        all_cases.update(conditions[parameter])

values = {case: diagnostics(case) for case in sorted(all_cases)}

results = []

for environment, conditions in CASES.items():

    print("\n" + "=" * 78)
    print(environment.upper())
    print("=" * 78)

    reference = values[conditions["reference"]]

    for parameter in ("alpha", "sum"):

        minus_case, plus_case = conditions[parameter]
        minus = values[minus_case]
        plus = values[plus_case]

        print(f"\nPARAMETER: {parameter}")
        print(
            f"{'Diagnostic':23s} "
            f"{'Reference':>12s} "
            f"{'-10%':>12s} "
            f"{'+10%':>12s} "
            f"{'Sensitivity':>13s} "
            f"{'Asymmetry %':>13s}"
        )

        for metric in reference:

            y0 = reference[metric]
            ym = minus[metric]
            yp = plus[metric]

            assert abs(y0) > 1.e-12, (environment, metric)

            sensitivity = (yp - ym) / (0.2 * y0)

            # Deviation from a locally symmetric response, relative to baseline.
            asymmetry = 100 * (yp + ym - 2*y0) / abs(y0)

            print(
                f"{metric:23s} "
                f"{y0:12.6f} "
                f"{ym:12.6f} "
                f"{yp:12.6f} "
                f"{sensitivity:+13.4f} "
                f"{asymmetry:+13.3f}"
            )

            results.append({
                "environment": environment,
                "parameter": parameter,
                "diagnostic": metric,
                "reference": y0,
                "minus_10": ym,
                "plus_10": yp,
                "normalised_sensitivity": sensitivity,
                "asymmetry_percent": asymmetry,
            })

with OUT.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=results[0].keys())
    writer.writeheader()
    writer.writerows(results)

print("\nSAVED:", OUT)
print("All ten simulations analysed successfully.")
