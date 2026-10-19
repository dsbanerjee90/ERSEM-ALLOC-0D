#!/usr/bin/env python3
# Author: Deep S. Banerjee, Plymouth Marine Laboratory (PML).
# ERSEM-ALLOC / ESA-CAMP: photosynthetic parameter sensitivity.

from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

ROOT = Path("/work/dba/ERSEM-ALLOC-0D/paper")
INPUT = ROOT / "analysis/alpha_sum_sensitivity_pilot.csv"
OUTPUT = ROOT / "figures"
OUTPUT.mkdir(parents=True, exist_ok=True)

environments = [
    ("No external nitrate supply", "No external supply"),
    ("External nitrate supply", "With external supply"),
]

parameters = ["alpha", "sum"]

diagnostics = [
    ("Mean fC", r"Mean allocation, $f_C$"),
    ("Gross fixation / C", "Gross fixation / C"),
    ("Nitrate uptake / C", r"NO$_3$ uptake / C"),
]

data = {}

with INPUT.open(newline="") as f:
    for row in csv.DictReader(f):
        key = (
            row["environment"],
            row["parameter"],
            row["diagnostic"]
        )
        assert key not in data, f"Duplicate: {key}"
        data[key] = float(row["normalised_sensitivity"])

matrices = []

for environment, _ in environments:
    matrix = np.array([
        [
            data[(environment, parameter, diagnostic)]
            for parameter in parameters
        ]
        for diagnostic, _ in diagnostics
    ])
    assert np.all(np.isfinite(matrix))
    matrices.append(matrix)

limit = max(np.max(np.abs(x)) for x in matrices)
limit = max(0.1, np.ceil(limit * 10) / 10)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 7.5,
    "axes.titlesize": 8,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

fig = plt.figure(figsize=(180/25.4, 82/25.4))

gs = fig.add_gridspec(
    1, 3,
    width_ratios=[1, 1, 0.055],
    left=0.25,
    right=0.91,
    bottom=0.18,
    top=0.80,
    wspace=0.22
)

axes = [fig.add_subplot(gs[0, i]) for i in range(2)]
cax = fig.add_subplot(gs[0, 2])

norm = TwoSlopeNorm(
    vmin=-limit,
    vcenter=0.0,
    vmax=limit
)

for index, (ax, matrix) in enumerate(zip(axes, matrices)):

    im = ax.imshow(
        matrix,
        cmap="RdBu_r",
        norm=norm,
        interpolation="nearest",
        aspect="auto"
    )

    ax.set_xticks(
        [0, 1],
        [r"$\alpha$", r"$\mathrm{sum}$"]
    )

    ax.set_yticks(range(3))

    if index == 0:
        ax.set_yticklabels([d[1] for d in diagnostics])
    else:
        ax.set_yticklabels([])

    ax.set_title(
        environments[index][1],
        pad=12
    )

    ax.tick_params(length=0)

    # Separate cells without obscuring the data.
    ax.set_xticks(np.arange(-0.5, 2, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, 3, 1), minor=True)
    ax.grid(
        which="minor",
        color="white",
        linewidth=1.6
    )
    ax.tick_params(which="minor", length=0)

    # Display signed sensitivity values.
    for i in range(3):
        for j in range(2):
            value = matrix[i, j]
            ax.text(
                j, i,
                f"{value:+.3f}",
                ha="center",
                va="center",
                fontsize=8,
                fontweight="medium",
                color="black"
            )

    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.text(
        -0.12, 1.09,
        "ab"[index],
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold"
    )

cb = fig.colorbar(im, cax=cax)
cb.set_label("Normalised sensitivity", fontsize=7.5)
cb.ax.tick_params(labelsize=7, width=0.5, length=2)
cb.outline.set_linewidth(0.5)

stem = OUTPUT / "fig_alpha_sum_sensitivity"

fig.savefig(stem.with_suffix(".pdf"))
fig.savefig(stem.with_suffix(".png"), dpi=450)

plt.close(fig)

print("SAVED:", stem.with_suffix(".pdf"))
print("SAVED:", stem.with_suffix(".png"))
print("Sensitivity range:", -limit, "to", limit)
