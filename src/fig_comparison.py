#!/usr/bin/env python3
"""Figure 6: arm 1 against arm 2, per ligand-receptor cell.

The claim this figure has to carry is a negative one -- that changing the
protonation protocol did not change the outcome -- and a negative is easy to
assert and hard to show. Plotting each cell's count under one arm against the
same cell under the other puts every point on the identity line when nothing
changed, so agreement is visible rather than asserted, and any exception is
immediately off the line and labelled.

Left panel: top-1 counts. Right panel: any-pose counts. Same axes, same scale,
so the two panels are directly comparable.

Usage: fig_comparison.py <arm1_dir> <arm2_dir> <out_dir>
"""

import csv
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ORDER = ["8UTB", "8V80", "8V82", "8V88"]
NAMES = {"8UTB": "NS-1738", "8V80": "(-)-TQS",
         "8V82": "PNU-120596", "8V88": "GAT107"}

A1, A2, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 300, "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "savefig.bbox": "tight",
})


def counts(path):
    top = defaultdict(int)
    any_ = defaultdict(int)
    n = 0
    with open(os.path.join(path, "matrix_summary.csv")) as fh:
        for row in csv.DictReader(fh):
            key = (row["ligand"], row["receptor"])
            top[key] += row["top1_success"] == "True"
            any_[key] += row["any_pose_success"] == "True"
            n += 1
    assert n == 80, "expected 80 runs in %s, found %d" % (path, n)
    return top, any_


t1, y1 = counts(A1)
t2, y2 = counts(A2)

fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.7))

for ax, d1, d2, title in (
        (axes[0], t1, t2, "Top-1 recovery"),
        (axes[1], y1, y2, "Any-pose recovery")):

    moved = []
    for lig in ORDER:
        for rec in ORDER:
            k = (lig, rec)
            x, y = d1[k], d2[k]
            cognate = lig == rec
            if x != y:
                moved.append((lig, rec, x, y))
            ax.scatter(x, y,
                       s=74 if cognate else 46,
                       c="#d1495b" if x != y else ("#1f6fb4" if cognate
                                                   else "#8fb8de"),
                       marker="D" if cognate else "o",
                       edgecolor="white", linewidth=0.7,
                       zorder=4 if x != y else 3)

    ax.plot([-0.3, 5.3], [-0.3, 5.3], color="#999999",
            linewidth=1, linestyle="--", zorder=1)

    for lig, rec, x, y in moved:
        ax.annotate("%s$\\rightarrow$%s" % (lig, rec), (x, y),
                    textcoords="offset points", xytext=(7, -3),
                    fontsize=7.5, color="#d1495b")

    ax.set_xlim(-0.3, 5.6)
    ax.set_ylim(-0.3, 5.6)
    ax.set_xticks(range(6))
    ax.set_yticks(range(6))
    ax.set_xlabel("Arm 1 (PDBFixer) seeds recovering")
    ax.set_ylabel("Arm 2 (PDB2PQR/PROPKA) seeds recovering")
    ax.set_title("%s\n%d of 16 cells off the identity line"
                 % (title, len(moved)))
    ax.set_aspect("equal")

# one legend for both panels, identity carried by shape not colour alone
handles = [
    plt.Line2D([], [], marker="D", color="none", markerfacecolor="#1f6fb4",
               markeredgecolor="white", markersize=8, label="cognate cell"),
    plt.Line2D([], [], marker="o", color="none", markerfacecolor="#8fb8de",
               markeredgecolor="white", markersize=7, label="cross cell"),
    plt.Line2D([], [], marker="o", color="none", markerfacecolor="#d1495b",
               markeredgecolor="white", markersize=7, label="changed between arms"),
]
fig.legend(handles=handles, frameon=False, ncol=3,
           loc="lower center", bbox_to_anchor=(0.5, -0.07))

fig.savefig(os.path.join(OUT, "fig6_arm1_vs_arm2.png"))
plt.close(fig)

print("top-1  cells off identity:", sum(1 for lig in ORDER for rec in ORDER
                                        if t1[(lig, rec)] != t2[(lig, rec)]))
print("any    cells off identity:", sum(1 for lig in ORDER for rec in ORDER
                                        if y1[(lig, rec)] != y2[(lig, rec)]))
for lig in ORDER:
    for rec in ORDER:
        k = (lig, rec)
        if y1[k] != y2[k]:
            print("  any-pose %s -> %s : %d -> %d"
                  % (lig, rec, y1[k], y2[k]))
print("wrote", os.path.join(OUT, "fig6_arm1_vs_arm2.png"))
