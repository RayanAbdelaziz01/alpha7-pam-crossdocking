#!/usr/bin/env python3
"""Determine which preparation arm's poses are currently in results/.

Both runners default to the same RUN_DIR, so results/ holds one arm's poses at
a time. The 720-pose validity sweep reads arm-1 checkpoints, so it is only
correct if results/ holds arm-1 poses. File timestamps cannot settle this
(an archive extraction rewrote them all at once), so this compares the rank-1
RMSD computed from the pose file on disk against the value each arm recorded in
its own checkpoint. Whichever arm matches is the arm on disk.

Usage: pose_provenance.py <arm1_drive_dir> <arm2_drive_dir> [work_dir]
"""

import json
import os
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolAlign

RDLogger.DisableLog("rdApp.*")

PAMS = {"8UTB": "XG3", "8V80": "YLR", "8V82": "I34", "8V88": "YLI"}

D1 = sys.argv[1]
A2 = sys.argv[2]
W = sys.argv[3] if len(sys.argv) > 3 else "/content/alpha7_matrix"

# Cells chosen because the two arms recorded different outcomes there, so a
# match is diagnostic rather than coincidental.
PROBES = [
    ("8V82", "8UTB", 2027),
    ("8V82", "8UTB", 2028),
    ("8UTB", "8UTB", 2029),
    ("8V80", "8V80", 2026),
    ("8V88", "8V88", 2030),
]


def checkpoint_rmsds(base, tag):
    path = os.path.join(base, "checkpoints", tag + ".json")
    if not os.path.exists(path):
        return None
    rec = json.load(open(path))
    return [round(r["rmsd_A"], 3) for r in rec["rows"]]


match1 = match2 = ambiguous = missing = 0

for lig, rec, seed in PROBES:
    tag = "%s_%s__into__%s_seed%d" % (lig, PAMS[lig], rec, seed)
    sdf = os.path.join(W, "results", tag + ".sdf")
    if not os.path.exists(sdf):
        print(tag, "| NO_SDF")
        missing += 1
        continue

    poses = [Chem.RemoveHs(m) for m in Chem.SDMolSupplier(sdf, removeHs=False)
             if m is not None]
    native = Chem.SDMolSupplier(
        os.path.join(W, "work", "%s_%s_native_heavy.sdf" % (lig, PAMS[lig])),
        removeHs=True)[0]
    disk = [round(rdMolAlign.CalcRMS(p, native, maxMatches=10000), 3)
            for p in poses]

    a1 = checkpoint_rmsds(D1, tag)
    a2 = checkpoint_rmsds(A2, tag)

    # compare the full ranked RMSD list, not just rank 1
    eq1 = a1 is not None and len(a1) == len(disk) and all(
        abs(x - y) < 0.01 for x, y in zip(disk, a1))
    eq2 = a2 is not None and len(a2) == len(disk) and all(
        abs(x - y) < 0.01 for x, y in zip(disk, a2))

    if eq1 and not eq2:
        verdict = "ARM1"
        match1 += 1
    elif eq2 and not eq1:
        verdict = "ARM2"
        match2 += 1
    elif eq1 and eq2:
        verdict = "BOTH (arms identical here, not diagnostic)"
        ambiguous += 1
    else:
        verdict = "NEITHER"
        ambiguous += 1

    print("%s\n   disk r1=%s  arm1 r1=%s  arm2 r1=%s  -> %s"
          % (tag, disk[0] if disk else None,
             a1[0] if a1 else None, a2[0] if a2 else None, verdict))

print("\nSUMMARY: arm1_matches=%d arm2_matches=%d ambiguous=%d missing=%d"
      % (match1, match2, ambiguous, missing))
if match1 and not match2:
    print("VERDICT: results/ holds ARM 1 poses. The arm-1 sweep is correct.")
elif match2 and not match1:
    print("VERDICT: results/ holds ARM 2 poses. The arm-1 sweep is MISMATCHED.")
else:
    print("VERDICT: INCONCLUSIVE - do not rely on the sweep until resolved.")
