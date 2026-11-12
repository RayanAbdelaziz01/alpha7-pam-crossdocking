# Cognate and cross-docking at the α7 nAChR modulator site

Code and data for *Receptor Identity, Not Cognate Pairing, Governs Docking
Success at the α7 Nicotinic Receptor Modulator Site*.

Four cryo-EM structures of the α7 nicotinic acetylcholine receptor, each solved
with a different positive allosteric modulator (PAM) bound (Burke et al., *Cell*
187:1160–1176, 2024), are used to run the full 4 × 4 cognate/cross-docking
matrix. Each PAM is docked into each receptor across five seeds, in one shared
coordinate frame and one shared docking box: 80 runs and up to 720 poses per
preparation arm, two arms, 160 runs and 1,440 poses in total.

## Result in one line

Docking success is organised by **receptor**, not by cognate pairing. All 26
top-1 successes in 80 runs fall in receptors 8V80 or 8V88; 8UTB and 8V82 return
none, including for their own ligands. The proposal's hypothesis — that cognate
conformations recover their own poses best — is not supported.

## Reproducing the analysis

Open `AS410662_Alpha7_PAM_CrossDocking.ipynb` in Google Colab on a **CPU**
runtime (AutoDock Vina has no GPU code path, so a GPU runtime gives no speedup).

**Default path (`RECOMPUTE = False`), about a minute.** Every table, figure and
number in the manuscript is re-derived from the archived per-run outputs. The
notebook resolves the archive in this order and stops at the first hit:

1. `/content/alpha7_data` already unpacked
2. `alpha7_artifacts.zip` present in `/content`
3. the `ARTIFACT_URL` set in the notebook
4. Google Drive — **off by default** (`USE_DRIVE = False`); nothing is read from
   anyone's personal cloud storage

If none of these is available the notebook raises `FileNotFoundError` with
instructions rather than failing silently.

**Full re-run (`RECOMPUTE = True`), about 3 h per arm on 2 CPUs.** Structures are
fetched from the PDB and the docking is re-run from scratch; the archive is not
needed at all in this mode. Every completed run is checkpointed as it finishes,
so an interrupted session resumes rather than restarting. Set `ARM` to 1 or 2.

The last cell re-derives every number quoted in the manuscript and asserts it. If
the notebook runs to the end without raising, the results have been reproduced.

## Environment

The notebook builds an isolated Python 3.12 virtual environment and pins:

| Package | Version | | Package | Version |
|---|---|---|---|---|
| AutoDock Vina | 1.2.7 | | PoseBusters | 0.6.5 |
| Meeko | 0.8.0 | | scikit-learn | 1.9.1 |
| gemmi | 0.7.4 | | SciPy | 1.17.1 |
| RDKit | 2025.9.4 | | NumPy | 1.26.4 |
| PDBFixer | 1.12.0 | | PDB2PQR | 3.8.0 |
| OpenMM | 8.3.1 | | PROPKA | 3.5.1 |

Seeds are 2026–2030. A provenance cell prints the resolved versions at runtime.

## The two preparation arms

The arms differ in exactly one variable, protonation:

- **Arm 1** — PDBFixer repairs heavy atoms *and* adds hydrogens at pH 7.4, then
  Meeko converts to PDBQT.
- **Arm 2** — identical through heavy-atom repair, after which PDB2PQR (AMBER)
  with PROPKA assigns protonation and titration states at pH 7.4. This is the
  recommended route for this toolchain; hydrogens placed by PDBFixer are a known
  source of Meeko valence errors.

Arm 2 reproduces arm 1's rank-1 matrix cell for cell, which is what excludes
protonation as the explanation for the receptor effect.

## Layout

```
AS410662_Alpha7_PAM_CrossDocking.ipynb   analysis notebook (start here)
src/matrix_runner.py                     arm 1 docking runner
src/matrix_runner_arm2.py                arm 2 docking runner (PDB2PQR/PROPKA)
src/analysis.py                          tables, clustering, contacts, figures
src/structural_correlates.py             Cα RMSDs, lipid and ion counts
src/patch_pb_all.py                      PoseBusters sweep over all 720 poses
src/fig_comparison.py                    arm 1 vs arm 2 per-run comparison
src/build_archive.py                     assembles alpha7_artifacts.zip
src/build_notebook.py                    generates the notebook
src/lint_notebook.py                     static execution-order checker
src/pose_provenance.py                   identifies which arm a pose set came from
main.tex, references.bib                 manuscript source
```

`alpha7_artifacts.zip` is attached to the repository's release. It carries the
per-run checkpoints, all 1,440 poses, the analysis tables, the figures and a
manifest with a SHA-256 for every entry.

## Data sources

Structures from the Protein Data Bank (8UTB, 8V80, 8V82, 8V88) and ligand
definitions from the RCSB chemical component dictionary (XG3, YLR, I34, YLI),
both retrieved 2026-11-11. Ligand bond orders and formal charges are taken from
the chemical component dictionary and never perceived from deposited geometry.

## Citation

Burke, S. M. et al. Structural mechanisms of α7 nicotinic receptor allosteric
modulation and activation. *Cell* **187**, 1160–1176 (2024).
doi:10.1016/j.cell.2024.01.032

Course project for AS.410.662 Molecular Modeling for Drug Discovery Lab,
Johns Hopkins University.
