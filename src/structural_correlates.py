#!/usr/bin/env python3
"""Recompute the structural correlates cited in the Discussion.

Two claims in the manuscript were asserted without an archived artifact behind
them: the global and pocket Calpha RMSDs, and the absence of bound Ca2+ in 8UTB.
This script recomputes both from the structures on disk and writes a JSON
artifact so the numbers are citable rather than recalled.

Calpha RMSD is computed directly, with no further superposition: the *_AB_raw.pdb
files are already in the common frame (all receptors superposed onto 8V82 on
matched A/B Calpha atoms), so fitting again here would measure fold similarity
rather than the frame displacement the Discussion is talking about.

Ions and lipids are counted from the original downloaded structures, because
receptor preparation strips heteroatoms.
"""

import glob
import json
import os
import sys

import gemmi
import numpy as np

W = os.environ.get("ALPHA7_MATRIX_DIR", "/content/alpha7_matrix")
ORDER = ["8UTB", "8V80", "8V82", "8V88"]
REF = "8V82"
BOX_CENTER = np.array([146.903, 153.799, 133.568])
POCKET_RADIUS = 12.0   # A around the box centre; the pocket subset

out = {"reference_frame": REF, "pocket_radius_A": POCKET_RADIUS, "sources": {}}


def ca_atoms(path):
    """Matched (chain, resnum) -> Calpha position, chains A and B only."""
    st = gemmi.read_structure(path)
    st.setup_entities()
    found = {}
    for chain in st[0]:
        if chain.name not in ("A", "B"):
            continue
        for res in chain:
            atom = res.find_atom("CA", "*")
            if atom is not None:
                found[(chain.name, res.seqid.num)] = np.array(
                    [atom.pos.x, atom.pos.y, atom.pos.z])
    return found


# -- Calpha RMSD in the common frame ---------------------------------------
raw = {}
for pdb in ORDER:
    path = os.path.join(W, "work", "%s_AB_raw.pdb" % pdb)
    if not os.path.exists(path):
        print("MISSING", path)
        continue
    raw[pdb] = ca_atoms(path)
    out["sources"].setdefault("ca_rmsd", []).append(os.path.basename(path))

ref = raw.get(REF)
rmsd = {}
if ref:
    for pdb in ORDER:
        if pdb == REF or pdb not in raw:
            continue
        shared = sorted(set(ref) & set(raw[pdb]))
        if not shared:
            continue
        a = np.array([ref[k] for k in shared])
        b = np.array([raw[pdb][k] for k in shared])
        d2 = ((a - b) ** 2).sum(axis=1)

        near = np.linalg.norm(a - BOX_CENTER, axis=1) < POCKET_RADIUS
        rmsd[pdb] = {
            "n_matched_ca": len(shared),
            "global_ca_rmsd_A": round(float(np.sqrt(d2.mean())), 3),
            "n_pocket_ca": int(near.sum()),
            "pocket_ca_rmsd_A": (round(float(np.sqrt(d2[near].mean())), 3)
                                 if near.any() else None),
        }
out["ca_rmsd_vs_%s" % REF] = rmsd


# -- ions and lipids from the original structures --------------------------
def original_for(pdb):
    for pattern in ("%s.cif", "%s.pdb", "%s*.cif", "%s*.pdb"):
        hits = sorted(glob.glob(os.path.join(W, "input", pattern % pdb)))
        if hits:
            return hits[0]
    hits = sorted(glob.glob(os.path.join(W, "**", "%s*.cif" % pdb), recursive=True))
    return hits[0] if hits else None


ions = {}
for pdb in ORDER:
    path = original_for(pdb)
    if not path:
        ions[pdb] = {"source": None, "note": "original structure not found"}
        continue
    st = gemmi.read_structure(path)
    counts = {}
    for chain in st[0]:
        for res in chain:
            if res.het_flag == "H" and res.name != "HOH":
                counts[res.name] = counts.get(res.name, 0) + 1
    ions[pdb] = {
        "source": os.path.basename(path),
        "CA_ion_copies": counts.get("CA", 0),
        "NA_ion_copies": counts.get("NA", 0),
        "POV_lipid_copies": counts.get("POV", 0),
        "all_het_residues": dict(sorted(counts.items())),
    }
out["heteroatoms"] = ions

dest = sys.argv[1] if len(sys.argv) > 1 else "/content/structural_correlates.json"
with open(dest, "w") as fh:
    json.dump(out, fh, indent=1)

print(json.dumps(out, indent=1))
print("\nwrote", dest)
