#!/usr/bin/env python3
"""Extend the PoseBusters stage from a 2-pose-per-run subset to every pose.

The original stage checked each run's rank-1 pose and its best-sampled pose,
which is 130 distinct poses once the runs where those coincide are counted
once. This patch checks all 9 poses of all 80 runs (720), so the
physical-validity claim covers the whole ensemble rather than a subset.

Writes to analysis_posebusters_all.csv so the original 130-pose file is kept
rather than overwritten.

Changes exactly three things and nothing else:
  1. the pose loop iterates over every pose instead of a 2-element dict
  2. the role label falls back to "other" for poses that are neither
  3. the output filename
"""

import pathlib
import sys

src = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "/content/analysis.py")
dst = pathlib.Path(sys.argv[2] if len(sys.argv) > 2 else "/content/analysis_all.py")

text = src.read_text()
lines = text.split("\n")
out = []
n_loop = n_name = n_log = 0

for line in lines:
    stripped = line.strip()

    # 1 + 2: iterate over every pose, label the rest "other"
    if stripped == "for rank, role in sorted(wanted.items()):":
        indent = line[: len(line) - len(line.lstrip())]
        out.append(indent + "for rank in range(1, len(poses) + 1):")
        out.append(indent + '    role = wanted.get(rank, "other")')
        n_loop += 1
        continue

    # 3: do not overwrite the original subset file
    if 'DRIVE_OUT / "analysis_posebusters.csv"' in line:
        out.append(line.replace('"analysis_posebusters.csv"',
                                '"analysis_posebusters_all.csv"'))
        n_name += 1
        continue

    # cosmetic: the progress line would otherwise fire 36 times
    if "posebusters:" in line and "len(rows) % 20" in line:
        out.append(line.replace("% 20", "% 100"))
        n_log += 1
        continue
    if stripped == "if len(rows) % 20 == 0:":
        out.append(line.replace("% 20", "% 100"))
        n_log += 1
        continue

    out.append(line)

assert n_loop == 1, "expected exactly 1 pose loop, patched %d" % n_loop
assert n_name == 1, "expected exactly 1 output filename, patched %d" % n_name

dst.write_text("\n".join(out))

# the patched file must still be valid Python
import ast
ast.parse(dst.read_text())

print("patched pose loop:      %d" % n_loop)
print("patched output name:    %d" % n_name)
print("patched log cadence:    %d" % n_log)
print("wrote %s (%d bytes), syntax OK" % (dst, dst.stat().st_size))
