# α7 nAChR PAM cross-docking — submission package

AS.410.662 Molecular Modeling for Drug Discovery Lab
Rayan Abdelaziz, Johns Hopkins University

This is everything needed to reproduce the manuscript's numbers. It does not
depend on the author's Google Drive.

---

## What is here

| File | What it is |
|------|------------|
| `AS410662_Alpha7_PAM_CrossDocking.ipynb` | The notebook. Runs top to bottom on a free Colab CPU runtime. |
| `alpha7_artifacts.zip` | All run outputs: per-run checkpoints, pose tables, analysis CSVs, figures, manifests. Contains `ARCHIVE_MANIFEST.json` with a SHA-256 for every file. |
| `matrix_runner.py` | Arm 1 docking runner (PDBFixer hydrogens). |
| `matrix_runner_arm2.py` | Arm 2 docking runner (PDB2PQR/PROPKA protonation). Differs in one variable. |
| `analysis.py` | Post-docking stages: posebusters, cluster, contacts, figures. |
| `main.tex`, `references.bib` | The manuscript, on the course Report Template. |
| `talk_script.md` | Speaking notes for the slide deck, section by section, with anticipated questions. |

The slide deck is at `https://claude.ai/artifact/GpyeZ7TyTSXG2TVxjfq6XP`.
Code and data are public at
`https://github.com/RayanAbdelaziz01/alpha7-pam-crossdocking`.

## Running it

1. Open the notebook in Google Colab. **Leave the runtime on CPU** — AutoDock
   Vina has no GPU code path, so a GPU runtime gives no speedup.
2. Runtime → Restart session and run all.

The run archive is downloaded automatically from the public repository, so
nothing needs to be uploaded. If that download is unavailable, put
`alpha7_artifacts.zip` in `/content` yourself (folder icon in the left sidebar)
and run again.

Code and data: `https://github.com/RayanAbdelaziz01/alpha7-pam-crossdocking`

With `RECOMPUTE = False` (the default) this reproduces every table and figure
from the archived outputs in about a minute. The docking itself is not re-run.

With `RECOMPUTE = True` the notebook fetches the four structures from the PDB
and re-runs the docking from scratch — roughly 3 h per arm on 2 CPUs. In that
mode the archive is not needed at all. Every completed run is checkpointed, so
an interrupted session resumes rather than restarting.

Nothing reads from any Google Drive. `USE_DRIVE` exists in the resolver as an
author convenience and defaults to `False`.

## What the experiment is

Four cryo-EM structures of α7, each solved with a different positive allosteric
modulator bound (Burke et al., *Cell* 187:1160–1176, 2024). Each modulator is
docked into each structure across five seeds, in one shared coordinate frame
with one shared docking box.

- 4 ligands × 4 receptors × 5 seeds = **80 runs per arm**, up to 9 poses each
- Two preparation arms = **160 runs, 1,440 poses**

Success is symmetry-aware heavy-atom RMSD to the deposited pose, computed in
place with no prior alignment, at a 2.0 Å threshold. Two measures are reported
separately throughout: **top-1** (the rank-1 pose is correct — what a virtual
screen consumes) and **any-pose** (some pose in the ranked list is correct —
whether the search ever found it). Failing top-1 while passing any-pose is a
scoring failure; failing both is a sampling failure.

## Headline result

All 26 top-1 successes across the 80 arm-1 runs occur in receptor 8V80 or 8V88.
Receptors 8UTB and 8V82 return none, for any ligand, including their own.
Success varies more across receptors (0, 11, 0, 15 of 20) than across ligands
(0, 6, 10, 10 of 20), and the cognate diagonal splits two cells at 5/5 against
two at 0/5 — so cognate docking was not an upper bound here.

Tested at the level of the sixteen cells rather than the eighty correlated runs, by exact enumeration of all 48,048 distinct arrangements: the receptor axis departs from the permutation null (p = 0.029) and the ligand axis does not (p = 0.50).

The second arm, differing only in protonation protocol, reproduced the rank-1
matrix cell for cell.

## Two things to know about the layout

**The arms use separate working directories.** Arm 1 runs in
`alpha7_matrix/`, arm 2 in `alpha7_matrix_arm2/`, each with its own `results/`,
`work/` and checkpoints. They must not share a directory — a shared one lets
the second arm overwrite the first arm's pose files, which is silent and would
mismatch poses to checkpoints. The notebook sets `ALPHA7_MATRIX_DIR` per arm
explicitly.

**Which arm a pose file belongs to is checkable.** `pose_provenance.py`
recomputes a pose's rank-1 RMSD from the file and matches the full ranked list
against each arm's checkpoints, reporting which arm it came from and flagging
cells where the arms agree and so cannot discriminate.

## Verifying the archive

```python
import json, hashlib, zipfile
with zipfile.ZipFile('alpha7_artifacts.zip') as zf:
    man = json.loads(zf.read('ARCHIVE_MANIFEST.json'))
    bad = [name for name, rec in man['files'].items()
           if hashlib.sha256(zf.read(name)).hexdigest() != rec['sha256']]
print('files:', man['n_files'], '| mismatched:', bad)
```

## Environment

AutoDock Vina 1.2.7, Meeko 0.8.0, gemmi 0.7.4, RDKit 2025.9.4, PDBFixer 1.12.0,
OpenMM 8.3.1, NumPy 1.26.4, PoseBusters 0.6.5, PDB2PQR 3.8.0 and PROPKA 3.5.1,
under Python 3.12. Seeds 2026–2030. The notebook carries a provenance cell that
prints all of this at run time, together with every input accession and its
retrieval date.

The patch version is deliberately not pinned. A clean-room run on a fresh Colab
VM resolved Python 3.12.3 rather than the 3.12.15 of the authoring machine,
while every package version above reproduced exactly; the environment is
reproducible at the package level, not at the interpreter patch level.

The docking runs in a pinned Python 3.12 virtual environment rather than the
Colab kernel, because several of these packages are built against 3.12 there.
The notebook creates that environment in its setup cell.
