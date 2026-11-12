#!/usr/bin/env python3
"""Assemble alpha7_artifacts.zip: the data file the notebook ships with.

The notebook resolves its inputs from this archive, so the archive has to
contain everything the analysis reads and nothing that depends on a particular
Google Drive. Layout inside the zip mirrors what the notebook expects:

    matrix_summary.csv              arm 1, top level
    matrix_all_poses.csv
    analysis_*.csv
    structural_correlates.json
    matrix_manifest.json
    regression_check.json
    checkpoints/*.json
    figures/*.png
    arm2_pdb2pqr/...                arm 2, same shape

A manifest with a SHA-256 for every entry is written into the archive so a
reviewer can verify the contents are the ones the manuscript reports.

Usage: build_archive.py <arm1_dir> <arm2_dir> <out_zip>
"""

import hashlib
import json
import os
import sys
import zipfile

A1, A2, OUT = sys.argv[1], sys.argv[2], sys.argv[3]

# Files whose absence means the archive is not usable, so a silent partial
# build is impossible.
REQUIRED_ARM1 = [
    "matrix_summary.csv", "matrix_all_poses.csv", "matrix_manifest.json",
    "regression_check.json", "analysis_contacts.csv", "analysis_clusters.csv",
    "analysis_posebusters_all.csv", "structural_correlates.json",
]
REQUIRED_ARM2 = [
    "matrix_summary.csv", "matrix_all_poses.csv", "matrix_manifest.json",
    "regression_check.json", "analysis_contacts.csv", "analysis_clusters.csv",
]

SKIP_SUFFIX = (".log", ".zip")


def collect(root, prefix):
    """(arcname, abspath) for everything worth shipping under root."""
    items = []
    for dirpath, dirnames, filenames in os.walk(root):
        # arm 2 lives inside arm 1's folder; it is collected separately
        dirnames[:] = [d for d in dirnames
                       if d not in ("arm2_pdb2pqr", "pb_tmp")]
        for fn in sorted(filenames):
            if fn.endswith(SKIP_SUFFIX):
                continue
            abs_ = os.path.join(dirpath, fn)
            rel = os.path.relpath(abs_, root)
            items.append((os.path.join(prefix, rel) if prefix else rel, abs_))
    return items


missing = [f for f in REQUIRED_ARM1 if not os.path.exists(os.path.join(A1, f))]
missing += ["arm2/" + f for f in REQUIRED_ARM2
            if not os.path.exists(os.path.join(A2, f))]
if missing:
    sys.exit("REFUSING to build a partial archive. Missing:\n  "
             + "\n  ".join(missing))

entries = collect(A1, "") + collect(A2, "arm2_pdb2pqr")

manifest = {"generated_utc": __import__("datetime").datetime.utcnow()
            .strftime("%Y-%m-%dT%H:%M:%SZ"),
            "n_files": len(entries), "files": {}}

total = 0
with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED) as zf:
    for arcname, abs_ in entries:
        data = open(abs_, "rb").read()
        total += len(data)
        manifest["files"][arcname] = {
            "bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest(),
        }
        zf.writestr(arcname, data)
    zf.writestr("ARCHIVE_MANIFEST.json", json.dumps(manifest, indent=1))

print("files archived     :", len(entries))
print("uncompressed bytes :", total)
print("archive size       :", os.path.getsize(OUT))
print("wrote              :", OUT)

# report what went in, by category, so a short count is checkable at a glance
cats = {}
for arcname, _ in entries:
    top = arcname.split(os.sep)[0]
    if top.endswith(".csv") or top.endswith(".json"):
        top = "(top-level files)"
    cats[top] = cats.get(top, 0) + 1
for k in sorted(cats):
    print("  %-28s %d" % (k, cats[k]))
