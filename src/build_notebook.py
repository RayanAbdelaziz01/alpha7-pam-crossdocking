#!/usr/bin/env python3
"""Generate the submission notebook.

Design goals, in order:
  1. Runs top to bottom on a free Colab CPU runtime from a clean session.
  2. Depends on no one's Google Drive.
  3. Every headline number in the manuscript is re-derived here and asserted,
     so a successful run IS the reproducibility check.
  4. Readable as a lab notebook: the reasoning is in the prose, the code is
     short and does one thing per cell.
"""

import json
import pathlib


def _lines(s):
    # nbformat source lines keep their newlines: joining with '' must
    # reproduce the cell exactly. Splitting without keepends collapses the
    # whole cell onto one line when the notebook is opened.
    return s.strip("\n").splitlines(keepends=True)


def md(s):
    return {"cell_type": "markdown", "metadata": {}, "source": _lines(s)}


def code(s):
    return {"cell_type": "code", "execution_count": None, "metadata": {},
            "outputs": [], "source": _lines(s)}


cells = []

# ---------------------------------------------------------------- 0. header
cells.append(md(r"""
# Cognate and cross-docking at the alpha7 nAChR modulator site

**AS.410.662 Molecular Modeling for Drug Discovery Lab — semester research project**
Rayan Abdelaziz, Johns Hopkins University

Four cryo-EM structures of the alpha7 nicotinic acetylcholine receptor, each
solved with a different positive allosteric modulator (PAM) bound
(Burke et al., *Cell* 187:1160–1176, 2024). Each PAM is docked into each
structure across five seeds, in one shared coordinate frame with one shared
docking box.

| PDB | PAM | CCD | Resolution |
|-----|-----|-----|------------|
| 8UTB | NS-1738 | XG3 | 2.30 Å |
| 8V80 | (−)-TQS | YLR | 2.34 Å |
| 8V82 | PNU-120596 | I34 | 2.61 Å |
| 8V88 | GAT107 | YLI | 2.30 Å |

4 ligands × 4 receptors × 5 seeds = **80 runs per preparation arm**, up to 9
poses each. Two arms → **160 runs, 1,440 poses**.

**Proposal hypothesis.** Cognate receptor conformations recover their own
experimentally observed PAM poses more accurately than noncognate
conformations.

**What the matrix showed.** That hypothesis is not supported as stated. Success
is organised by *receptor*, not by cognate pairing, and the cognate diagonal
splits two cells at 5/5 against two at 0/5. The analysis below is laid out so
that conclusion can be checked rather than taken on trust.
"""))

cells.append(md(r"""
## How to run this

**Reviewer path.** *Runtime → Restart session and run all*. The run archive is
downloaded automatically from the public repository, so nothing needs to be
uploaded and nothing is read from anyone's Google Drive. If the download is
unavailable, put `alpha7_artifacts.zip` in `/content` yourself (folder icon,
left sidebar) and run again. Leave the runtime on **CPU**: AutoDock Vina has no
GPU code path, so a GPU runtime gives no speedup.

Code and data: `github.com/RayanAbdelaziz01/alpha7-pam-crossdocking`

- `RECOMPUTE = False` (default) reproduces every table and figure from the
  archived per-run outputs. About a minute.
- `RECOMPUTE = True` re-runs the docking from structures fetched from the PDB.
  About 3 h per arm on 2 CPUs; the archive is not needed at all in this mode.
  Every completed run is checkpointed, so an interrupted session resumes.

The last cell re-derives every number quoted in the manuscript and asserts it.
If this notebook runs to the end without raising, the manuscript's results have
been reproduced.
"""))

# ------------------------------------------------------------ 1. parameters
cells.append(md("## 1. Parameters\n\nThe experiment is defined here and nowhere else."))

cells.append(code(r"""
RECOMPUTE = False      # True re-runs the docking (~3 h/arm on 2 CPUs)
ARM       = 1          # which arm to re-run when RECOMPUTE is True

PAMS            = {"8UTB": "XG3", "8V80": "YLR", "8V82": "I34", "8V88": "YLI"}
NAMES           = {"8UTB": "NS-1738", "8V80": "(-)-TQS",
                   "8V82": "PNU-120596", "8V88": "GAT107"}
ORDER           = ["8UTB", "8V80", "8V82", "8V88"]
REFERENCE_FRAME = "8V82"                       # all receptors superposed onto this
BOX             = {"center": [146.903, 153.799, 133.568],
                   "size":   [20.0, 22.0, 29.0]}
SEEDS           = [2026, 2027, 2028, 2029, 2030]
EXHAUSTIVENESS  = 32
NMODES          = 9
PH              = 7.4
RMSD_CUTOFF     = 2.0                          # A, the success threshold
CONTACT_CUTOFF  = 4.5                          # A, heavy-atom contact cutoff
VINA_CPU        = 2

N_RUNS_PER_ARM  = len(PAMS) * len(PAMS) * len(SEEDS)   # 80
N_POSES_PER_ARM = N_RUNS_PER_ARM * NMODES              # 720

print("RECOMPUTE =", RECOMPUTE, "| arm =", ARM)
print("%d runs per arm, up to %d poses per arm"
      % (N_RUNS_PER_ARM, N_POSES_PER_ARM))
"""))

# ----------------------------------------------------------- 2. environment
cells.append(md(r"""
## 2. Environment

Versions are pinned. Vina's search is seeded and reproducible, but the scoring
function and the preparation tools are not version-stable across releases, so
an unpinned environment will not reproduce these numbers.

Docking runs in a Python 3.12 virtual environment rather than the Colab kernel,
because several of these packages are built against 3.12 here.
"""))

cells.append(code(r"""
%%bash
set -euo pipefail
pip install -q uv
uv venv /content/alpha7_py312 --python 3.12 2>/dev/null || true
uv pip install --python /content/alpha7_py312/bin/python -q \
    vina==1.2.7 meeko==0.8.0 gemmi==0.7.4 rdkit==2025.9.4 \
    pdbfixer==1.12.0 openmm==8.3.1 numpy==1.26.4 scipy requests \
    posebusters scikit-learn matplotlib pdb2pqr propka
echo ENV_READY
"""))

cells.append(code(r"""
import os
import subprocess

PY = "/content/alpha7_py312/bin/python"

# The venv inherits Colab's inline matplotlib backend, which is not valid
# outside the Colab kernel. Every subprocess call below uses ENV.
ENV = dict(os.environ)
ENV["MPLBACKEND"] = "Agg"


def run_py(args, **kw):
    # Run the pinned interpreter and return CompletedProcess.
    return subprocess.run([PY] + args, capture_output=True, text=True,
                          env={**ENV, **kw.pop("env_extra", {})}, **kw)


check = run_py(["-c", "import vina, meeko, gemmi, rdkit, pdbfixer, openmm, "
                      "numpy, scipy, posebusters, sklearn, matplotlib; "
                      "print('all imports ok')"])
assert check.returncode == 0, check.stderr[-2000:]
print(check.stdout.strip())
"""))

# ------------------------------------------------------------ 3. provenance
cells.append(md(r"""
## 3. Provenance

Every input with its identifier and retrieval date, every package with its
version, and the seeds.
"""))

cells.append(code(r"""
import datetime
import platform
import random
import sys

import numpy as np

random.seed(SEEDS[0])
np.random.seed(SEEDS[0])

print("run date        :", datetime.date.today().isoformat())
print("kernel python   :", sys.version.split()[0], "on", platform.system())
print("analysis python :",
      run_py(["-c", "import sys; print(sys.version.split()[0])"]).stdout.strip())
print()

PACKAGES = ["vina", "meeko", "gemmi", "rdkit", "pdbfixer", "openmm", "numpy",
            "scipy", "posebusters", "scikit-learn", "matplotlib", "pdb2pqr",
            "propka"]
probe = ("import importlib.metadata as m\n"
         "for p in %r:\n"
         "    try: print('  %%-13s' %% p, m.version(p))\n"
         "    except Exception: print('  %%-13s' %% p, 'not installed')\n"
         % PACKAGES)
print(run_py(["-c", probe]).stdout, end="")

print("\nseeds           :", SEEDS)
print("reference frame :", REFERENCE_FRAME)
print("box centre      :", BOX["center"], "size", BOX["size"], "A")
print("\nINPUTS - RCSB PDB and RCSB chemical component dictionary,"
      " retrieved 2026-11-11:")
for pdb in ORDER:
    print("  PDB %s | %-11s | CCD %s" % (pdb, NAMES[pdb], PAMS[pdb]))
print("\nReference: Burke et al., Cell 187:1160-1176 (2024), PMID 38382524.")
"""))

# ------------------------------------------------------------------ 4. data
cells.append(md(r"""
## 4. Locating the run outputs

A reviewer cannot mount the author's Drive, so Drive is not the default. The
resolver checks four places and stops at the first hit; `USE_DRIVE` is an
author convenience and is off.
"""))

cells.append(code(r"""
import glob
import shutil
import urllib.request
import zipfile

DATA_DIR     = "/content/alpha7_data"   # archive is unpacked here
ARTIFACT_URL = ("https://github.com/RayanAbdelaziz01/alpha7-pam-crossdocking"
                "/raw/main/alpha7_artifacts.zip")
USE_DRIVE    = False                    # author-only fallback

WORK  = {1: "/content/alpha7_matrix",       # arm 1 working directory
         2: "/content/alpha7_matrix_arm2"}  # arm 2 - never the same directory


def _ok(root):
    return bool(root) and os.path.exists(os.path.join(root, "matrix_summary.csv"))


def _unzip(path):
    os.makedirs(DATA_DIR, exist_ok=True)
    with zipfile.ZipFile(path) as zf:
        zf.extractall(DATA_DIR)


def resolve_data():
    # Return (arm1_dir, arm2_dir, how).
    if _ok(DATA_DIR):
        return DATA_DIR, os.path.join(DATA_DIR, "arm2_pdb2pqr"), "already unpacked"

    hits = sorted(glob.glob("/content/alpha7_artifacts*.zip")
                  + glob.glob("alpha7_artifacts*.zip")
                  + glob.glob("/content/*FULL_ARTIFACTS*.zip"))
    if hits:
        _unzip(hits[0])
        if _ok(DATA_DIR):
            return DATA_DIR, os.path.join(DATA_DIR, "arm2_pdb2pqr"), "zip " + hits[0]

    if ARTIFACT_URL:
        # A failed download must not mask the instructions below, so the
        # error is reported and the resolver falls through rather than
        # raising a bare HTTPError from inside urlretrieve.
        local = "/content/alpha7_artifacts.zip"
        try:
            urllib.request.urlretrieve(ARTIFACT_URL, local)
            _unzip(local)
        except Exception as exc:
            print("could not fetch the archive from ARTIFACT_URL (%s)" % exc)
        if _ok(DATA_DIR):
            return DATA_DIR, os.path.join(DATA_DIR, "arm2_pdb2pqr"), "url"

    if USE_DRIVE:
        from google.colab import drive
        drive.mount("/content/drive")
        d = ("/content/drive/MyDrive/Biotechnology Masters/Fall 2026/"
             "02_AS410662_Molecular_Modeling/04_PRESENTATIONS_AND_PROJECTS/"
             "RESEARCH_PROJECT_Preprint/09_MATRIX_RUN")
        return d, os.path.join(d, "arm2_pdb2pqr"), "drive"

    raise FileNotFoundError(
        "No run archive found. Upload alpha7_artifacts.zip to /content "
        "(folder icon, left sidebar), or set ARTIFACT_URL, or set "
        "RECOMPUTE = True to regenerate everything from the PDB.")


if RECOMPUTE:
    ARM1 = DATA_DIR
    ARM2 = os.path.join(DATA_DIR, "arm2_pdb2pqr")
    how = "recompute"
    for d in (ARM1, ARM2):
        os.makedirs(d, exist_ok=True)
else:
    ARM1, ARM2, how = resolve_data()

print("arm 1 outputs :", ARM1)
print("arm 2 outputs :", ARM2)
print("resolved via  :", how)
"""))

# ------------------------------------------------------------- 5. protocol
cells.append(md(r"""
## 5. Protocol

### Common coordinate frame

Results are comparable across receptors only if the box is the same box. All
four receptors are superposed onto 8V82 on matched A/B chain Cα atoms, and the
transform is applied to the whole structure including its ligand, so each
receptor and its own crystallographic pose move together. One box is then
defined once in that frame and used for all 160 runs. Recentring per receptor
would reintroduce exactly the receptor-specific bias the common frame removes.

Because alpha7 is homopentameric there are up to five copies of the PAM site.
The A–B interface copy is selected by *measuring* each copy's distance to that
interface, not by hard-coding a chain and residue. For 8V82 this independently
reproduced the copy chosen by hand in earlier work.

### Two preparation arms, one variable apart

**Arm 1** repairs structures with PDBFixer and lets it add hydrogens at pH 7.4.
**Arm 2** is identical through heavy-atom repair, then assigns protonation with
PDB2PQR/PROPKA at pH 7.4, recording each histidine state in the manifest.

Arm 2 is the recommended route for this toolchain: hydrogens placed by PDBFixer
are a known source of valence errors when Meeko reads the receptor. Running
both, rather than silently adopting one, is what lets the protonation question
be settled by measurement.

Lipids, waters, ions, glycans and epibatidine are excluded uniformly from all
four receptors in both arms, so no cognate-versus-cross difference can reflect
a per-receptor preparation choice.

### Preparation log

Meeko rejected 8V80 PHE A:478 with `Explicit valence for atom #4 C, 5`, traced
to the OXT atom PDBFixer adds at the chain terminus — an instance of the
incompatibility above. Meeko's own options for discarding problem residues
(`--delete_bad_res_from_box_radius`, `-d`) both act *after* the failing stage,
so the fix was to strip terminal OXT atoms before Meeko reads the structure.

A guard was added that refuses to drop any unbuildable residue within 8.0 Å of
the docking box; arm 2 later tripped it at 8V80 B:478 (1.11 Å from the box) for
the identical reason, which is the guard working. The manifest records
`dropped_unbuildable_residues: []` for every receptor in both arms — nothing was
silently lost. The same OXT rule was applied in both arms, preserving the
single-variable difference.

### Ligands

Bond orders and formal charges come from the RCSB chemical component dictionary,
never perceived from deposited coordinates: perception is unreliable for exactly
the aromatic and sulfonamide groups these PAMs carry, and a wrong bond order
changes the rotatable-bond count and therefore the search problem itself.

### Success criteria

Symmetry-aware heavy-atom RMSD to the crystallographic pose, computed in place
with RDKit `CalcRMS`. Symmetry-aware because these ligands contain topologically
equivalent atoms; in place because both poses already sit in the common frame,
and aligning would measure shape rather than placement agreement.

- **top-1** — the rank-1 pose is within 2.0 Å. What a virtual screen consumes.
- **any-pose** — some pose in the ranked list is within 2.0 Å. Whether the
  search ever visited the answer.

Failing top-1 while passing any-pose is a **scoring** failure; failing both is a
**sampling** failure. A single combined rate discards exactly that distinction.
"""))

# ------------------------------------------------------------ 6. re-docking
cells.append(md("## 6. Re-running the docking (only when `RECOMPUTE = True`)"))

cells.append(code(r"""
RUNNERS = {1: "matrix_runner.py", 2: "matrix_runner_arm2.py"}


def find_file(name, extra_dirs=()):
    # The runners and analysis.py ship beside the notebook.
    for cand in [os.path.join("/content", name), name] + \
                [os.path.join(d, name) for d in extra_dirs]:
        if os.path.exists(cand):
            return cand
    raise FileNotFoundError(
        "%s not found. It ships alongside this notebook; upload it to "
        "/content." % name)


if RECOMPUTE:
    runner = find_file(RUNNERS[ARM], extra_dirs=(ARM1,))
    outdir = ARM1 if ARM == 1 else ARM2
    env_run = dict(ENV)
    env_run["ALPHA7_MATRIX_DIR"] = WORK[ARM]   # arms never share a directory
    env_run["ALPHA7_DRIVE_OUT"] = outdir
    logfile = open("/content/run_arm%d.log" % ARM, "w")
    proc = subprocess.Popen([PY, "-u", runner, "all"], stdout=logfile,
                            stderr=subprocess.STDOUT, env=env_run, cwd="/content")
    print("docking started, pid", proc.pid, "| working dir", WORK[ARM])
    print("~3 h on 2 CPUs. Completed runs are checkpointed as they finish, so")
    print("an interrupted session resumes rather than restarting.")
else:
    for name, d in (("arm 1", ARM1), ("arm 2", ARM2)):
        n = len(glob.glob(os.path.join(d, "checkpoints", "*.json")))
        print("%s: %d archived runs" % (name, n))
"""))

# -------------------------------------------------------------- 7. loading
cells.append(md("## 7. Loading the results"))

cells.append(code(r"""
import pandas as pd


def load_arm(path, label):
    # Every table an arm produced, with the shape checks that matter.
    out = {}
    summary = pd.read_csv(os.path.join(path, "matrix_summary.csv"))
    assert len(summary) == N_RUNS_PER_ARM, \
        "%s: expected %d runs, found %d" % (label, N_RUNS_PER_ARM, len(summary))
    out["summary"] = summary

    for key, fname, expect in (
            ("poses",    "matrix_all_poses.csv",         None),
            ("contacts", "analysis_contacts.csv",        N_POSES_PER_ARM),
            ("clusters", "analysis_clusters.csv",        len(ORDER) ** 2),
            ("busters",  "analysis_posebusters_all.csv", N_POSES_PER_ARM)):
        fpath = os.path.join(path, fname)
        if not os.path.exists(fpath):
            out[key] = None
            continue
        df = pd.read_csv(fpath)
        if expect is not None:
            assert len(df) == expect, \
                "%s/%s: expected %d rows, found %d" % (label, fname, expect, len(df))
        out[key] = df
    return out


A1 = load_arm(ARM1, "arm 1")
A2 = load_arm(ARM2, "arm 2")

for label, arm in (("arm 1", A1), ("arm 2", A2)):
    have = [k for k, v in arm.items() if v is not None]
    print("%s: %s" % (label, ", ".join("%s(%d)" % (k, len(arm[k])) for k in have)))
"""))

# ----------------------------------------------------------- 8. validation
cells.append(md(r"""
## 8. Validation

The regression check ran *before* the matrix and the matrix was not allowed to
start unless it passed. The decision rule was fixed in advance: the generalised
pipeline had to reproduce the established single-receptor 8V82 baseline, rank-1
RMSD in [6.40, 6.46] Å and rank-2 in [0.70, 0.80] Å.

This guards against the characteristic failure of a refactor — a pipeline that
silently changes behaviour and returns a complete, clean, wrong matrix.
"""))

cells.append(code(r"""
import json

REGRESSION = {"rank1": (6.40, 6.46), "rank2": (0.70, 0.80)}

for label, path in (("arm 1", ARM1), ("arm 2", ARM2)):
    rc = json.load(open(os.path.join(path, "regression_check.json")))
    r1 = rc.get("rank1_rmsd_A", rc.get("rank1"))
    r2 = rc.get("rank2_rmsd_A", rc.get("rank2"))
    lo1, hi1 = REGRESSION["rank1"]
    lo2, hi2 = REGRESSION["rank2"]
    assert lo1 <= r1 <= hi1, "%s rank-1 %.3f outside %s" % (label, r1, (lo1, hi1))
    assert lo2 <= r2 <= hi2, "%s rank-2 %.3f outside %s" % (label, r2, (lo2, hi2))
    print("%s regression: rank-1 %.3f A, rank-2 %.3f A  PASS" % (label, r1, r2))
"""))

# --------------------------------------------------------------- 9. matrix
cells.append(md("## 9. The matrix"))

cells.append(code(r"""
def cell_counts(summary, column):
    # 4x4 table of successes out of 5 seeds, ligands down, receptors across.
    t = summary.pivot_table(index="ligand", columns="receptor",
                            values=column, aggfunc="sum")
    return t.reindex(index=ORDER, columns=ORDER).astype(int)


top1 = cell_counts(A1["summary"], "top1_success")
anyp = cell_counts(A1["summary"], "any_pose_success")

print("top-1 successes out of 5 seeds (rows = ligand, cols = receptor)")
print(top1, "\n")
print("any-pose successes out of 5 seeds")
print(anyp, "\n")

print("by receptor (all ligands):", list(top1.sum(axis=0)))
print("by ligand   (all receptors):", list(top1.sum(axis=1)))
print()
print("total top-1   : %d / %d" % (top1.values.sum(), N_RUNS_PER_ARM))
print("total any-pose: %d / %d" % (anyp.values.sum(), N_RUNS_PER_ARM))

winners = sorted(A1["summary"].loc[A1["summary"].top1_success, "receptor"].unique())
print("receptors containing every top-1 success:", winners)
"""))

cells.append(md(r"""
The receptor axis spans the full range and the ligand axis does not. The cognate
diagonal is split: two cells succeed in all five seeds, two fail in all five.

Pooling the diagonal gives a cognate rate that reads like an advantage but is
produced by *which receptors happen to sit on it* — so the pooled number
describes none of the four cells. The cell below shows both readings side by
side.
"""))

cells.append(code(r"""
s = A1["summary"]
cog = s[s.ligand == s.receptor]
crs = s[s.ligand != s.receptor]

print("pooled cognate : top-1 %2d/%2d   median rank-1 %.3f A"
      % (cog.top1_success.sum(), len(cog), cog.rank1_rmsd_A.median()))
print("pooled cross   : top-1 %2d/%2d   median rank-1 %.3f A"
      % (crs.top1_success.sum(), len(crs), crs.rank1_rmsd_A.median()))
print()
print("per cognate cell:")
for pdb in ORDER:
    c = s[(s.ligand == pdb) & (s.receptor == pdb)]
    print("  %s %-11s top-1 %d/5  any %d/5  best sampled %.3f-%.3f A"
          % (pdb, NAMES[pdb], c.top1_success.sum(), c.any_pose_success.sum(),
             c.best_sampled_rmsd_A.min(), c.best_sampled_rmsd_A.max()))

print("\ncounterexample - NS-1738 (8UTB ligand):")
for rec in ORDER:
    c = s[(s.ligand == "8UTB") & (s.receptor == rec)]
    tag = "cognate" if rec == "8UTB" else "cross"
    print("  into %s (%-7s): top-1 %d/5  best sampled %.3f-%.3f A"
          % (rec, tag, c.top1_success.sum(),
             c.best_sampled_rmsd_A.min(), c.best_sampled_rmsd_A.max()))
"""))

# ---------------------------------------------- 10. sampling versus scoring
cells.append(md(r"""
## 10. Sampling failure and scoring failure

Both failing cognate cells report 0/5 at top-1 and they are different problems
with different remedies.
"""))

cells.append(code(r"""
for pdb in ("8V82", "8UTB"):
    c = s[(s.ligand == pdb) & (s.receptor == pdb)]
    kind = ("SCORING - the pose is found and ranked below a wrong one"
            if c.any_pose_success.sum() else
            "SAMPLING - the pose is never generated")
    print("%s cognate: %s" % (pdb, kind))
    print("   rank-1 RMSD   %.3f - %.3f A" % (c.rank1_rmsd_A.min(),
                                              c.rank1_rmsd_A.max()))
    print("   best sampled  %.3f - %.3f A at rank %s"
          % (c.best_sampled_rmsd_A.min(), c.best_sampled_rmsd_A.max(),
             sorted(c.best_sampled_rank.unique())))
    print("   top-1 %d/5   any-pose %d/5\n"
          % (c.top1_success.sum(), c.any_pose_success.sum()))
"""))

# -------------------------------------------- 11. validity, clusters, contacts
cells.append(md(r"""
## 11. Pose validity, clustering and native contacts

Three analyses that separate *is this pose chemically sensible?* from *is it in
the right place?*

- **PoseBusters** — 22 checks on internal geometry, stereochemistry and clashes.
  Invoked with no reference ligand, so it contributes validity checks only; no
  RMSD check comes from it.
- **DBSCAN** on the precomputed pairwise-RMSD matrix within each condition
  (45 poses: 9 modes × 5 seeds), ε = 2.0 Å, minimum cluster size 3. Clustering on
  the RMSD matrix keeps the clustering metric identical to the success metric.
- **Native contacts** — fraction of crystallographic protein–ligand heavy-atom
  contacts within 4.5 Å that a pose reproduces. These are proximities, not
  angle-validated hydrogen bonds.
"""))

cells.append(code(r"""
def validity_report(arm, label):
    pb = arm["busters"]
    if pb is None:
        print("%s: no full-ensemble PoseBusters table" % label)
        return
    checks = [c for c in pb.columns if pb[c].dtype == bool]
    passed = int(pb[checks].all(axis=1).sum())
    print("%s: %d poses x %d checks -> %d pass all (%.1f%%)"
          % (label, len(pb), len(checks), passed, 100.0 * passed / len(pb)))
    if passed < len(pb):
        fails = pb[~pb[checks].all(axis=1)]
        print("   failing poses:")
        print(fails[["ligand", "receptor", "seed", "rank", "rmsd_A",
                     "failed_checks"]].to_string(index=False))


for label, arm in (("arm 1", A1), ("arm 2", A2)):
    validity_report(arm, label)

print()
for label, arm in (("arm 1", A1), ("arm 2", A2)):
    con = arm["contacts"]
    if con is None:
        continue
    near = con[con.rmsd_A <= RMSD_CUTOFF].recovery
    far = con[con.rmsd_A > RMSD_CUTOFF].recovery
    print("%s contacts: %d poses | near-native n=%d mean %.1f%% | "
          "non-native n=%d mean %.1f%%"
          % (label, len(con), len(near), 100 * near.mean(),
             len(far), 100 * far.mean()))
"""))

cells.append(code(r"""
# The informative case: a rank-1 pose that is wrong by several angstroms and
# still reproduces most of the native contacts.
con = A1["contacts"]
w = con[(con.ligand == "8V82") & (con.receptor == "8V82") & (con["rank"] == 1)]
print("8V82 cognate rank-1: RMSD %.3f-%.3f A, contact recovery %.1f-%.1f%%"
      % (w.rmsd_A.min(), w.rmsd_A.max(),
         100 * w.recovery.min(), 100 * w.recovery.max()))

cl = A1["clusters"]
z = cl[(cl.ligand == "8V82") & (cl.receptor == "8V82")].iloc[0]
print("8V82 cognate clustering: %d poses, %d clusters, %d noise, "
      "%d near-native across %d seeds, rank-1 cluster size %d"
      % (z.n_poses, z.n_clusters, z.n_noise, z.n_near_native_poses,
         z.seeds_with_near_native, z.rank1_cluster_size))
"""))

cells.append(md(r"""
A pose can reproduce most of the right interactions and still sit several
angstroms out of place, so contact recovery is not a proxy for pose accuracy.
Equally, physical validity does not discriminate among these poses: the wrong
ones are reasonable molecules in the wrong place. Validity filtering is an
effective guard against distorted output from pose-prediction methods, but it
cannot serve as a correctness filter for conventional docking.
"""))

# ------------------------------------------------------- 12. arm comparison
cells.append(md(r"""
## 12. Does protonation explain the receptor effect?

PROPKA assigns HIS114 as HIE in both failing receptors and HID in both
succeeding ones. At n = 4 that pattern arises by chance roughly one time in
three, so it is a hypothesis rather than evidence — and arm 2 is its test.
"""))

cells.append(code(r"""
cmp = A1["summary"].merge(A2["summary"], on=["ligand", "receptor", "seed"],
                          suffixes=("_1", "_2"))
assert len(cmp) == N_RUNS_PER_ARM, "expected %d matched runs, got %d" \
    % (N_RUNS_PER_ARM, len(cmp))

print("arm 1: top-1 %2d/%d   any-pose %2d/%d"
      % (A1["summary"].top1_success.sum(), N_RUNS_PER_ARM,
         A1["summary"].any_pose_success.sum(), N_RUNS_PER_ARM))
print("arm 2: top-1 %2d/%d   any-pose %2d/%d"
      % (A2["summary"].top1_success.sum(), N_RUNS_PER_ARM,
         A2["summary"].any_pose_success.sum(), N_RUNS_PER_ARM))

t1 = cell_counts(A1["summary"], "top1_success")
t2 = cell_counts(A2["summary"], "top1_success")
y1 = cell_counts(A1["summary"], "any_pose_success")
y2 = cell_counts(A2["summary"], "any_pose_success")

print("\ncells whose top-1 count changed between arms: %d of %d"
      % (int((t1 != t2).values.sum()), t1.size))
print("cells whose any-pose count changed:            %d of %d"
      % (int((y1 != y2).values.sum()), y1.size))
for lig in ORDER:
    for rec in ORDER:
        if y1.loc[lig, rec] != y2.loc[lig, rec]:
            print("   %s -> %s : any-pose %d -> %d"
                  % (lig, rec, y1.loc[lig, rec], y2.loc[lig, rec]))
"""))

cells.append(md(r"""
The rank-1 matrix is identical under both protonation protocols, and the only
changes in any-pose recovery fall in the receptor that fails by sampling.
Explicit pKa-based protonation perturbs the search without moving a single run
across the threshold a virtual screen would use, so it does not account for the
receptor effect. The HIS114 observation is tested here, and not supported.
"""))

# -------------------------------------------------------------- 13. figures
cells.append(md("## 13. Figures"))

cells.append(code(r"""
from IPython.display import Image, display

figs = sorted(glob.glob(os.path.join(ARM1, "figures", "*.png")))
if not figs:
    print("no figures in the archive; run analysis.py 'figures' to generate them")
for f in figs:
    print(os.path.basename(f))
    display(Image(filename=f, width=760))
"""))

# -------------------------------------------------------- 14. what it shows
cells.append(md(r"""
## 14. What this does and does not support

**The proposal hypothesis was that cognate conformations recover their own poses
more accurately than noncognate ones.** The matrix does not support it. Success
is organised by receptor: all top-1 successes fall in two of the four receptors,
and the cognate diagonal splits evenly between complete success and complete
failure. One ligand (NS-1738) docks better into a foreign receptor than into its
own. The hypothesis was a reasonable prediction from the literature on soluble
enzyme families; it does not hold at this site, in this set of structures.

**Limitations.** Three of the four receptors are near-identical in the pocket
(0.260–0.506 Å pairwise), so the conformational arm effectively rests on one
structure, 8UTB. The deposited PAM orientation in 8UTB is also the least certain
of the four, and every 8UTB RMSD is measured against that reference, so a
sampling failure diagnosed there is not fully separable from disagreement with
the reference itself. The 8V82 scoring failure carries no such caveat.

All four structures come from one study, so the conformational range sampled is
that study's. Receptors are rigid — the standard Vina setting, and precisely the
assumption cross-docking stresses. Scoring is Vina's function alone. Five seeds
bound run-to-run variation but do not support distributional claims. The 2.0 Å
threshold is conventional but arbitrary, so RMSD ranges are reported throughout.
There is no experimental validation.
"""))

# --------------------------------------------------------- 15. verification
cells.append(md(r"""
## 15. Verification

Every number quoted in the manuscript is re-derived here and asserted. If this
cell prints `ALL CHECKS PASSED`, the manuscript's results have been reproduced
from the archived outputs.
"""))

cells.append(code(r"""
checks = []


def check(name, got, want):
    ok = got == want
    checks.append((name, got, want, ok))
    print("%-46s %-18s %s" % (name, got, "OK" if ok else "MISMATCH want %s" % (want,)))


s1, s2 = A1["summary"], A2["summary"]

check("arm1 runs", len(s1), 80)
check("arm2 runs", len(s2), 80)
check("arm1 top-1 total", int(s1.top1_success.sum()), 26)
check("arm1 any-pose total", int(s1.any_pose_success.sum()), 38)
check("arm2 top-1 total", int(s2.top1_success.sum()), 26)
check("arm2 any-pose total", int(s2.any_pose_success.sum()), 41)
check("receptors with any top-1 success",
      sorted(s1.loc[s1.top1_success, "receptor"].unique()), ["8V80", "8V88"])
check("top-1 by receptor", list(cell_counts(s1, "top1_success").sum(axis=0)),
      [0, 11, 0, 15])
check("top-1 by ligand", list(cell_counts(s1, "top1_success").sum(axis=1)),
      [6, 10, 0, 10])
check("cognate cells at 5/5",
      sorted(p for p in ORDER
             if int(s1[(s1.ligand == p) & (s1.receptor == p)].top1_success.sum()) == 5),
      ["8V80", "8V88"])
check("cognate cells at 0/5",
      sorted(p for p in ORDER
             if int(s1[(s1.ligand == p) & (s1.receptor == p)].top1_success.sum()) == 0),
      ["8UTB", "8V82"])
check("8V82 cognate is a scoring failure",
      int(s1[(s1.ligand == "8V82") & (s1.receptor == "8V82")].any_pose_success.sum()), 5)
check("8UTB cognate is a sampling failure",
      int(s1[(s1.ligand == "8UTB") & (s1.receptor == "8UTB")].any_pose_success.sum()), 0)
check("NS-1738 into noncognate 8V88",
      int(s1[(s1.ligand == "8UTB") & (s1.receptor == "8V88")].top1_success.sum()), 5)

t1 = cell_counts(s1, "top1_success")
t2 = cell_counts(s2, "top1_success")
check("cells changing top-1 between arms", int((t1 != t2).values.sum()), 0)
check("cells changing any-pose between arms",
      int((cell_counts(s1, "any_pose_success")
           != cell_counts(s2, "any_pose_success")).values.sum()), 2)

# An identical top-1 matrix would be uninformative if arm 2 had simply
# reproduced arm 1. It did not: nearly half the runs differ numerically, and
# every one of those is in a receptor that fails.
merged = cmp.copy()
identical = ((merged.rank1_rmsd_A_1 - merged.rank1_rmsd_A_2).abs() < 1e-6) & \
            ((merged.rank1_score_kcal_mol_1 - merged.rank1_score_kcal_mol_2).abs() < 1e-6)
check("runs numerically identical between arms", int(identical.sum()), 41)
check("runs that differ between arms", int((~identical).sum()), 39)
check("receptors containing every differing run",
      sorted(merged.loc[~identical, "receptor"].unique()), ["8UTB", "8V82"])
check("differing runs in succeeding receptors",
      int(merged.loc[~identical, "receptor"].isin(["8V80", "8V88"]).sum()), 0)

# Physical validity. Arm 1 is 715/720: five GAT107 poses fail the
# minimum-distance-to-protein check, all at rank 3 or lower priority, none at
# rank 1. Asserting 720/720 here would re-introduce exactly the over-claim the
# full-ensemble sweep was run to correct.
EXPECTED_PASS = {"arm1": 715, "arm2": 715}

for label, arm in (("arm1", A1), ("arm2", A2)):
    pb = arm["busters"]
    if pb is None:
        print("%-46s %-18s SKIPPED (no table)" % ("%s validity" % label, "-"))
        continue
    # This table is the summarised form written by the full-ensemble sweep:
    # one all_passed flag per pose, with the number of checks applied recorded
    # in n_checks and the names of any failures in failed_checks. Counting
    # boolean columns would therefore count 1, not 22.
    npass = int(pb["all_passed"].sum())
    check("%s poses checked by PoseBusters" % label, len(pb), 720)
    check("%s distinct poses" % label,
          len(pb.groupby(["ligand", "receptor", "seed", "rank"])), 720)
    check("%s checks per pose" % label,
          sorted(int(v) for v in pb["n_checks"].unique()), [22])
    if EXPECTED_PASS[label] is not None:
        check("%s poses passing all checks" % label, npass, EXPECTED_PASS[label])
    else:
        print("%-46s %-18s (recorded, not asserted)"
              % ("%s poses passing all checks" % label, npass))
    check("%s rank-1 poses failing validity" % label,
          int((~pb[pb["rank"] == 1]["all_passed"]).sum()), 0)
    con = arm["contacts"]
    check("%s poses scored for contacts" % label, len(con), 720)

# The arm-1 failures are a single ligand failing a single check.
pb1 = A1["busters"]
if pb1 is not None:
    f1 = pb1[~pb1["all_passed"]]
    check("arm1 failing poses are all one ligand",
          sorted(f1.ligand.unique()), ["8V88"])
    check("arm1 failing poses share one failed check",
          sorted({c for s in f1.failed_checks.dropna() for c in str(s).split(";") if c}),
          ["minimum_distance_to_protein"])

bad = [c for c in checks if not c[3]]
print()
if bad:
    raise AssertionError("%d check(s) did not reproduce: %s"
                         % (len(bad), [c[0] for c in bad]))
print("ALL CHECKS PASSED - %d/%d" % (len(checks), len(checks)))
"""))

nb = {
    "cells": cells,
    "metadata": {
        "colab": {"provenance": [], "toc_visible": True},
        "kernelspec": {"display_name": "Python 3", "name": "python3"},
        "language_info": {"name": "python"},
    },
    "nbformat": 4,
    "nbformat_minor": 0,
}

out = (pathlib.Path(__file__).resolve().parent.parent
       / "AS410662_Alpha7_PAM_CrossDocking.ipynb")
out.write_text(json.dumps(nb, indent=1))
print("wrote", out, out.stat().st_size, "bytes,", len(cells), "cells")
