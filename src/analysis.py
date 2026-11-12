#!/usr/bin/env python
"""
AS.410.662 - alpha7 nAChR PAM validation study
Post-docking analysis of the 4 x 4 matrix.

Stages
------
posebusters : physical plausibility of the rank-1 and best-sampled pose of
              every run, checked against the prepared receptor.
cluster     : DBSCAN pose clustering per condition, on a precomputed matrix of
              symmetry-aware in-place heavy-atom RMSDs.
contacts    : native-contact recovery - which receptor residues contacting the
              deposited pose are reproduced by each docked pose.
figures     : the manuscript figure set.

Every stage writes a CSV to Drive and is independent of the others, so a stage
can be re-run without repeating the rest.
"""

import argparse
import csv
import json
import os
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolAlign

RDLogger.DisableLog("rdApp.*")

RUN_DIR = Path(os.environ.get("ALPHA7_MATRIX_DIR", "/content/alpha7_matrix"))
DRIVE_OUT = Path(os.environ["ALPHA7_DRIVE_OUT"])
WORK = RUN_DIR / "work"
RESULTS = RUN_DIR / "results"
CKPT = DRIVE_OUT / "checkpoints"
FIGS = DRIVE_OUT / "figures"

PAMS = {"8UTB": "XG3", "8V80": "YLR", "8V82": "I34", "8V88": "YLI"}
NAMES = {"8UTB": "NS-1738", "8V80": "(-)-TQS",
         "8V82": "PNU-120596", "8V88": "GAT107"}
ORDER = ["8UTB", "8V80", "8V82", "8V88"]
SEEDS = [2026, 2027, 2028, 2029, 2030]

RMSD_CUTOFF = 2.0          # pose-recovery criterion, Module 7 convention
CONTACT_CUTOFF = 4.5       # heavy-atom contact cutoff, Module 8 convention
DBSCAN_EPS = 2.0           # same 2.0 A as the pose-recovery criterion
DBSCAN_MIN_SAMPLES = 3     # a cluster must recur beyond a single seed's noise


def log(*parts):
    print("[%s] %s" % (time.strftime("%H:%M:%S"),
                       " ".join(str(p) for p in parts)), flush=True)


def tag_for(lig, rec, seed):
    return "%s_%s__into__%s_seed%d" % (lig, PAMS[lig], rec, seed)


def load_records():
    records = [json.loads(p.read_text()) for p in sorted(CKPT.glob("*.json"))]
    assert records, "no checkpoints found"
    return records


def load_poses(lig, rec, seed):
    sdf = RESULTS / (tag_for(lig, rec, seed) + ".sdf")
    if not sdf.exists():
        return []
    return [Chem.RemoveHs(m) for m in
            Chem.SDMolSupplier(str(sdf), removeHs=False) if m is not None]


def native_mol(lig):
    path = WORK / ("%s_%s_native_heavy.sdf" % (lig, PAMS[lig]))
    return Chem.SDMolSupplier(str(path), removeHs=True)[0]


def write_csv(path, rows):
    if not rows:
        log("no rows for", path.name)
        return
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    log("wrote", path.name, "(%d rows)" % len(rows))


# --------------------------------------------------------------------------
# stage: posebusters
# --------------------------------------------------------------------------

def stage_posebusters():
    """Physical plausibility of the rank-1 and best-sampled pose of each run.

    A pose that recovers the crystallographic geometry but fails a physical
    check is not evidence of a working protocol, and a top-ranked pose that is
    physically impossible is a different kind of scoring failure from one that
    is merely in the wrong place.
    """
    from posebusters import PoseBusters

    buster = PoseBusters(config="dock")
    tmp = RUN_DIR / "pb_tmp"
    tmp.mkdir(exist_ok=True)

    rows = []
    for record in load_records():
        lig, rec, seed = record["ligand"], record["receptor"], record["seed"]
        poses = load_poses(lig, rec, seed)
        if not poses:
            continue
        protein = WORK / ("%s_AB_H.pdb" % rec)
        best_rank = record["best_sampled_rank"]
        wanted = {1: "rank1", best_rank: "best_sampled"}

        for rank, role in sorted(wanted.items()):
            mol = poses[rank - 1]
            single = tmp / ("%s_r%d.sdf" % (tag_for(lig, rec, seed), rank))
            writer = Chem.SDWriter(str(single))
            writer.write(mol)
            writer.close()
            try:
                frame = buster.bust([str(single)], None, str(protein))
            except Exception as exc:                      # noqa: BLE001
                log("  posebusters failed on", single.name, repr(exc))
                continue
            checks = frame.iloc[0].to_dict()
            passed = {k: bool(v) for k, v in checks.items()
                      if isinstance(v, (bool, np.bool_))}
            row = {"ligand": lig, "receptor": rec,
                   "condition": record["condition"], "seed": seed,
                   "rank": rank, "role": role,
                   "rmsd_A": record["rows"][rank - 1]["rmsd_A"],
                   "score_kcal_mol": record["rows"][rank - 1]["score_kcal_mol"],
                   "n_checks": len(passed),
                   "n_failed": sum(1 for v in passed.values() if not v),
                   "all_passed": all(passed.values()),
                   "failed_checks": ";".join(sorted(k for k, v in passed.items()
                                                    if not v))}
            rows.append(row)
        if len(rows) % 20 == 0:
            log("  posebusters:", len(rows), "poses checked")

    write_csv(DRIVE_OUT / "analysis_posebusters.csv", rows)
    return rows


# --------------------------------------------------------------------------
# stage: clustering
# --------------------------------------------------------------------------

def stage_cluster():
    """DBSCAN clustering of the 45 poses (5 seeds x 9 modes) per condition.

    Distance is the same symmetry-aware in-place heavy-atom RMSD used for pose
    recovery, supplied to DBSCAN as a precomputed matrix, so clusters are in the
    same units as the 2.0 A criterion. Label -1 marks noise: a pose that no
    other pose reproduces.
    """
    from sklearn.cluster import DBSCAN

    rows = []
    for lig in ORDER:
        for rec in ORDER:
            mols, meta = [], []
            for seed in SEEDS:
                poses = load_poses(lig, rec, seed)
                for rank, mol in enumerate(poses, start=1):
                    mols.append(mol)
                    meta.append((seed, rank))
            if not mols:
                continue

            n = len(mols)
            dist = np.zeros((n, n), dtype=float)
            for i in range(n):
                for j in range(i):
                    value = rdMolAlign.CalcRMS(mols[i], mols[j],
                                               maxMatches=10000)
                    dist[i, j] = dist[j, i] = value

            labels = DBSCAN(eps=DBSCAN_EPS, min_samples=DBSCAN_MIN_SAMPLES,
                            metric="precomputed").fit(dist).labels_

            native = native_mol(lig)
            to_native = np.array([rdMolAlign.CalcRMS(m, native, maxMatches=10000)
                                  for m in mols])

            sizes = defaultdict(int)
            for label in labels:
                sizes[int(label)] += 1
            rank1_idx = [k for k, (s, r) in enumerate(meta) if r == 1]
            rank1_labels = [int(labels[k]) for k in rank1_idx]

            near = to_native < RMSD_CUTOFF
            near_labels = sorted({int(labels[k]) for k in range(n) if near[k]})

            rows.append({
                "ligand": lig, "receptor": rec,
                "condition": "cognate" if lig == rec else "cross",
                "n_poses": n,
                "n_clusters": len({int(l) for l in labels if l != -1}),
                "n_noise": sizes.get(-1, 0),
                "largest_cluster": max((v for k, v in sizes.items() if k != -1),
                                       default=0),
                "rank1_cluster_labels": ";".join(str(x) for x in rank1_labels),
                "rank1_cluster_size": sizes.get(rank1_labels[0], 0)
                                      if rank1_labels else 0,
                "n_near_native_poses": int(near.sum()),
                "near_native_cluster_labels": ";".join(str(x) for x in near_labels),
                "min_rmsd_A": round(float(to_native.min()), 3),
                "seeds_with_near_native": len({meta[k][0] for k in range(n)
                                               if near[k]}),
            })
            log("  cluster %s->%s: %d clusters, %d noise, %d near-native"
                % (lig, rec, rows[-1]["n_clusters"], rows[-1]["n_noise"],
                   rows[-1]["n_near_native_poses"]))

    write_csv(DRIVE_OUT / "analysis_clusters.csv", rows)
    return rows


# --------------------------------------------------------------------------
# stage: native contacts
# --------------------------------------------------------------------------

def receptor_heavy_atoms(rec):
    """Heavy atoms of the prepared receptor, as (chain, resnum, resname, xyz)."""
    path = WORK / ("%s_AB_H.pdb" % rec)
    keys, coords = [], []
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name = line[12:16].strip()
        if name.startswith("H") or line[76:78].strip() == "H":
            continue
        keys.append((line[21], line[22:26].strip(), line[17:20].strip()))
        coords.append([float(line[30:38]), float(line[38:46]),
                       float(line[46:54])])
    return keys, np.asarray(coords)


def contact_set(keys, coords, mol):
    conf = mol.GetConformer()
    lig = np.array([list(conf.GetAtomPosition(i))
                    for i in range(mol.GetNumAtoms())])
    diff = coords[:, None, :] - lig[None, :, :]
    close = (np.linalg.norm(diff, axis=-1) < CONTACT_CUTOFF).any(axis=1)
    return {keys[i] for i in np.nonzero(close)[0]}


def stage_contacts():
    """Fraction of the deposited pose's receptor contacts that each pose reproduces.

    Contacts are heavy-atom proximities at 4.5 A, the Module 8 convention. They
    are proximities, not validated hydrogen bonds: no donor-acceptor angle is
    tested here.
    """
    cache = {rec: receptor_heavy_atoms(rec) for rec in ORDER}
    rows = []
    for record in load_records():
        lig, rec, seed = record["ligand"], record["receptor"], record["seed"]
        poses = load_poses(lig, rec, seed)
        if not poses:
            continue
        keys, coords = cache[rec]
        native_contacts = contact_set(keys, coords, native_mol(lig))
        if not native_contacts:
            continue
        for rank, mol in enumerate(poses, start=1):
            found = contact_set(keys, coords, mol)
            shared = native_contacts & found
            union = native_contacts | found
            rows.append({
                "ligand": lig, "receptor": rec,
                "condition": record["condition"], "seed": seed, "rank": rank,
                "rmsd_A": record["rows"][rank - 1]["rmsd_A"],
                "score_kcal_mol": record["rows"][rank - 1]["score_kcal_mol"],
                "n_native_contacts": len(native_contacts),
                "n_pose_contacts": len(found),
                "n_shared": len(shared),
                "recovery": round(len(shared) / len(native_contacts), 4),
                "jaccard": round(len(shared) / len(union), 4) if union else 0.0,
            })
        log("  contacts %s->%s seed %d" % (lig, rec, seed))

    write_csv(DRIVE_OUT / "analysis_contacts.csv", rows)
    return rows


# --------------------------------------------------------------------------
# stage: figures
# --------------------------------------------------------------------------

def stage_figures():
    """Manuscript figure set. Axes carry units; everything saved at 300 dpi."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIGS.mkdir(parents=True, exist_ok=True)
    records = load_records()

    plt.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 300, "font.size": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
        "savefig.bbox": "tight",
    })

    def cell(lig, rec, key, reducer):
        vals = [r[key] for r in records
                if r["ligand"] == lig and r["receptor"] == rec]
        return reducer(vals) if vals else np.nan

    # -- Figure 1: top-1 success rate ------------------------------------
    success = np.array([[cell(l, r, "top1_success", lambda v: 100.0 * sum(v) / len(v))
                         for r in ORDER] for l in ORDER])
    fig, ax = plt.subplots(figsize=(4.4, 3.9))
    im = ax.imshow(success, cmap="Blues", vmin=0, vmax=100)
    for i in range(4):
        for j in range(4):
            ax.text(j, i, "%.0f%%" % success[i, j], ha="center", va="center",
                    fontsize=9,
                    color="white" if success[i, j] > 55 else "#222222")
            if i == j:
                ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fill=False,
                                           edgecolor="#d1495b", linewidth=2.2))
    ax.set_xticks(range(4), ORDER)
    ax.set_yticks(range(4), ["%s\n%s" % (NAMES[l], l) for l in ORDER])
    ax.set_xlabel("Receptor structure")
    ax.set_ylabel("Ligand (deposited PAM)")
    ax.set_title("Top-1 pose recovery below 2.0 $\\AA$\n(5 seeds per condition)")
    ax.grid(False)
    fig.colorbar(im, ax=ax, label="Seeds recovering the pose at rank 1 (%)",
                 fraction=0.046)
    fig.savefig(FIGS / "fig1_top1_success.png")
    plt.close(fig)

    # -- Figure 2: best-sampled RMSD -------------------------------------
    best = np.array([[cell(l, r, "best_sampled_rmsd_A", np.median)
                      for r in ORDER] for l in ORDER])
    fig, ax = plt.subplots(figsize=(4.4, 3.9))
    im = ax.imshow(best, cmap="Blues_r", vmin=0, vmax=max(8.0, np.nanmax(best)))
    for i in range(4):
        for j in range(4):
            ax.text(j, i, "%.2f" % best[i, j], ha="center", va="center",
                    fontsize=9,
                    color="white" if best[i, j] < 3 else "#222222")
            if i == j:
                ax.add_patch(plt.Rectangle((j - .5, i - .5), 1, 1, fill=False,
                                           edgecolor="#d1495b", linewidth=2.2))
    ax.set_xticks(range(4), ORDER)
    ax.set_yticks(range(4), ["%s\n%s" % (NAMES[l], l) for l in ORDER])
    ax.set_xlabel("Receptor structure")
    ax.set_ylabel("Ligand (deposited PAM)")
    ax.set_title("Median best-sampled pose RMSD\n(closest of 9 poses, 5 seeds)")
    ax.grid(False)
    fig.colorbar(im, ax=ax, label="RMSD to deposited pose ($\\AA$)",
                 fraction=0.046)
    fig.savefig(FIGS / "fig2_best_sampled_rmsd.png")
    plt.close(fig)

    # -- Figure 3: sampling versus scoring -------------------------------
    fig, ax = plt.subplots(figsize=(5.0, 4.2))
    for cond, colour, marker in (("cognate", "#1f6fb4", "o"),
                                 ("cross", "#e08a1e", "^")):
        xs = [r["best_sampled_rmsd_A"] for r in records if r["condition"] == cond]
        ys = [r["rank1_rmsd_A"] for r in records if r["condition"] == cond]
        ax.scatter(xs, ys, s=34, c=colour, marker=marker, alpha=0.8,
                   edgecolor="white", linewidth=0.6, label=cond, zorder=3)
    lim = 11.0
    ax.plot([0, lim], [0, lim], color="#999999", linewidth=1,
            linestyle="--", zorder=1)
    ax.axhline(RMSD_CUTOFF, color="#d1495b", linewidth=1, zorder=2)
    ax.axvline(RMSD_CUTOFF, color="#d1495b", linewidth=1, zorder=2)
    ax.text(lim * 0.97, RMSD_CUTOFF + 0.2, "2.0 $\\AA$ criterion", ha="right",
            fontsize=8, color="#d1495b")
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("Best-sampled pose RMSD ($\\AA$)  -  search")
    ax.set_ylabel("Rank-1 pose RMSD ($\\AA$)  -  scoring")
    ax.set_title("Sampling succeeds where scoring fails\n(one point per run, n=%d)"
                 % len(records))
    ax.legend(frameon=False, loc="lower right")
    fig.savefig(FIGS / "fig3_sampling_vs_scoring.png")
    plt.close(fig)

    # -- Figure 4: score against RMSD, all poses -------------------------
    fig, ax = plt.subplots(figsize=(5.0, 4.2))
    for cond, colour, marker in (("cognate", "#1f6fb4", "o"),
                                 ("cross", "#e08a1e", "^")):
        xs = [row["rmsd_A"] for r in records if r["condition"] == cond
              for row in r["rows"]]
        ys = [row["score_kcal_mol"] for r in records if r["condition"] == cond
              for row in r["rows"]]
        ax.scatter(xs, ys, s=11, c=colour, marker=marker, alpha=0.45,
                   linewidth=0, label=cond, zorder=3)
    ax.axvline(RMSD_CUTOFF, color="#d1495b", linewidth=1, zorder=2)
    ax.set_xlabel("Pose RMSD to deposited ligand ($\\AA$)")
    ax.set_ylabel("Vina score (kcal/mol)")
    ax.set_title("Score does not track pose accuracy\n(all %d poses)"
                 % sum(len(r["rows"]) for r in records))
    ax.legend(frameon=False, loc="lower right")
    fig.savefig(FIGS / "fig4_score_vs_rmsd.png")
    plt.close(fig)

    # -- Figure 5: contact recovery, if that stage has run ---------------
    contacts_csv = DRIVE_OUT / "analysis_contacts.csv"
    if contacts_csv.exists():
        with contacts_csv.open() as fh:
            crows = list(csv.DictReader(fh))
        if crows:
            fig, ax = plt.subplots(figsize=(5.0, 4.2))
            xs = [float(r["rmsd_A"]) for r in crows]
            ys = [100.0 * float(r["recovery"]) for r in crows]
            ax.scatter(xs, ys, s=9, c="#1f6fb4", alpha=0.35, linewidth=0)
            ax.axvline(RMSD_CUTOFF, color="#d1495b", linewidth=1)
            ax.set_xlabel("Pose RMSD to deposited ligand ($\\AA$)")
            ax.set_ylabel("Native receptor contacts reproduced (%)")
            ax.set_title("Contact recovery against pose accuracy\n"
                         "(heavy-atom proximity, %.1f $\\AA$ cutoff)"
                         % CONTACT_CUTOFF)
            fig.savefig(FIGS / "fig5_contact_recovery.png")
            plt.close(fig)

    log("figures written to", FIGS)


# --------------------------------------------------------------------------

STAGES = {"posebusters": stage_posebusters, "cluster": stage_cluster,
          "contacts": stage_contacts, "figures": stage_figures}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=list(STAGES) + ["all"])
    args = parser.parse_args()
    DRIVE_OUT.mkdir(parents=True, exist_ok=True)
    names = list(STAGES) if args.stage == "all" else [args.stage]
    for name in names:
        log("=== stage:", name, "===")
        STAGES[name]()
    log("analysis complete")


if __name__ == "__main__":
    main()
