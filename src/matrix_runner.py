#!/usr/bin/env python
"""
AS.410.662 - alpha7 nAChR PAM validation study
4 x 4 cognate / cross-docking matrix runner.

Design notes
------------
* Receptors: A-B transmembrane interface of 8UTB, 8V80, 8V82, 8V88.
* Ligands:   the deposited PAM of each structure (XG3, YLR, I34, YLI),
             built from the RCSB chemical component dictionary so that
             bond orders and formal charges are never guessed from geometry.
* Frame:     every receptor is superposed onto the 8V82 A-B interface using
             matched CA atoms. All receptors, all reference ligands and the
             single shared docking box therefore live in one common frame,
             which is what makes an in-place (unsuperposed) RMSD meaningful
             for the twelve cross-docking conditions.
* Baseline:  rigid receptor, polymer only. Epibatidine, PAMs, lipids, waters,
             glycans and ions are excluded uniformly from all four receptors.
             This is the approved Module 7 baseline and it is the only
             pocket-component rule that is symmetric across the four
             structures (8UTB resolves no lipid within 8.7 A of its PAM,
             while the other three each carry two at 2.3-2.9 A).
* Evidence:  one JSON checkpoint per (ligand, receptor, seed) written to
             Drive before the next seed starts, so an interrupted Colab
             runtime costs one seed rather than the run.
"""

import argparse
import csv
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from importlib import metadata
from pathlib import Path

import gemmi
import numpy as np
import requests
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

# --------------------------------------------------------------------------
# configuration
# --------------------------------------------------------------------------

RUN_DIR = Path(os.environ.get("ALPHA7_MATRIX_DIR", "/content/alpha7_matrix"))
DRIVE_OUT = Path(os.environ["ALPHA7_DRIVE_OUT"])

INPUT = RUN_DIR / "input"
WORK = RUN_DIR / "work"
RESULTS = RUN_DIR / "results"
CKPT = DRIVE_OUT / "checkpoints"

PAMS = {
    "8UTB": "XG3",   # NS-1738
    "8V80": "YLR",   # (-)-TQS
    "8V82": "I34",   # PNU-120596
    "8V88": "YLI",   # GAT107
}

REFERENCE_FRAME = "8V82"

# Approved Module 7 docking box, in the 8V82 frame.
BOX = {"center": [146.903, 153.799, 133.568], "size": [20.0, 22.0, 29.0]}

SEEDS = [2026, 2027, 2028, 2029, 2030]
EXHAUSTIVENESS = 32
NMODES = 9
PH = 7.4
RMSD_CUTOFF = 2.0
VINA_CPU = 2

# Regression target: the rank-1 and rank-2 RMSDs recorded in the
# Module 7 checkpoint for PNU-120596 redocked into cognate 8V82.
REGRESSION = {
    "rank1_range": (6.40, 6.46),
    "rank2_range": (0.70, 0.80),
}


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------

def log(*parts):
    msg = " ".join(str(p) for p in parts)
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fetch(url, path):
    if path.exists() and path.stat().st_size > 3000:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    r = requests.get(url, timeout=180,
                     headers={"User-Agent": "AS410662 academic research/1.0"})
    r.raise_for_status()
    path.write_bytes(r.content)
    assert path.stat().st_size > 3000, f"Unexpectedly small response: {url}"
    log("downloaded", path.name, path.stat().st_size, "bytes")


def run_tool(args, log_name):
    args = [str(a) for a in args]
    args[0] = str(Path(sys.executable).parent / args[0])
    proc = subprocess.run(args, text=True, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, cwd=str(RUN_DIR))
    (WORK / log_name).write_text(proc.stdout)
    if proc.returncode:
        print(proc.stdout[-4000:], flush=True)
        raise RuntimeError(f"command failed: {' '.join(args)}")
    return proc.stdout


def run_tool_soft(args, log_name):
    """Run a tool without raising; return (returncode, combined output)."""
    args = [str(a) for a in args]
    args[0] = str(Path(sys.executable).parent / args[0])
    proc = subprocess.run(args, text=True, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, cwd=str(RUN_DIR))
    (WORK / log_name).write_text(proc.stdout)
    return proc.returncode, proc.stdout


def residue_com_distance_from_box(pdb_path, key):
    """Distance from a residue's centre of mass to the docking box region."""
    if ":" not in key:
        return None
    chain, num = key.split(":", 1)
    num = num.strip()
    rows = [l for l in open(pdb_path)
            if l.startswith(("ATOM", "HETATM"))
            and l[21] == chain and l[22:26].strip() == num]
    if not rows:
        return None
    xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                    for l in rows])
    com = xyz.mean(axis=0)
    outside = np.maximum(
        np.abs(com - np.array(BOX["center"])) - np.array(BOX["size"]) / 2.0, 0.0)
    return float(np.linalg.norm(outside))


def meeko_receptor_with_retries(pdb, prepared, prefix, max_tries=12):
    """Prepare a receptor, dropping only unbuildable residues outside the box."""
    dropped = []
    for _ in range(max_tries):
        args = ["mk_prepare_receptor.py", "--read_pdb", prepared,
                "-o", prefix, "-p"]
        code, out = run_tool_soft(args, "meeko_receptor_%s.log" % pdb)
        if code == 0:
            if dropped:
                log("  %s: dropped %d residue(s) outside the box: %s"
                    % (pdb, len(dropped), ", ".join(dropped)))
            return dropped
        found = re.search(r"corresponding to key (\S+)", out)
        if not found:
            print(out[-3000:], flush=True)
            raise RuntimeError(
                "%s: receptor preparation failed without a residue key" % pdb)
        key = found.group(1).strip().rstrip(".,")
        if key in dropped:
            raise RuntimeError("%s: residue %s still failing after deletion"
                               % (pdb, key))
        distance = residue_com_distance_from_box(prepared, key)
        if distance is None:
            raise RuntimeError("%s: cannot locate unbuildable residue %s"
                               % (pdb, key))
        if distance < 8.0:
            raise RuntimeError(
                "%s: unbuildable residue %s lies %.2f A from the docking box; "
                "refusing to silently drop a binding-site residue"
                % (pdb, key, distance))
        log("  %s: residue %s unbuildable, %.1f A outside the box - deleting"
            % (pdb, key, distance))
        _chain, _num = key.split(":", 1)
        _num = _num.strip()
        _lines = Path(prepared).read_text().splitlines(True)
        _kept = [l for l in _lines
                 if not (l.startswith(("ATOM", "HETATM"))
                         and l[21] == _chain and l[22:26].strip() == _num)]
        if len(_kept) == len(_lines):
            raise RuntimeError("%s: no atom records matched residue %s"
                               % (pdb, key))
        Path(prepared).write_text("".join(_kept))
        dropped.append(key)
    raise RuntimeError("%s: too many unbuildable residues: %s" % (pdb, dropped))


def write_status(payload):
    (DRIVE_OUT / "matrix_status.json").write_text(json.dumps(payload, indent=1))


# --------------------------------------------------------------------------
# ligand construction from the chemical component dictionary
# --------------------------------------------------------------------------

def interface_pam_copy(structure, ccd):
    """Return (chain, seqid, dist_A, dist_B) of the PAM copy at the A-B interface.

    The PAM site is a transmembrane intersubunit pocket, so the correct copy
    is the one with atoms close to both chain A and chain B. Chosen by
    measurement rather than hard-coded, so the same rule applies to all four
    deposits.
    """
    model = structure[0]
    ca = {}
    for chain in model:
        pts = [a.pos for res in chain for a in res
               if a.name == "CA" and res.het_flag == "A"]
        if pts:
            ca[chain.name] = pts
    assert "A" in ca and "B" in ca, "missing chain A or B"

    best, best_score = None, None
    for chain in model:
        for res in chain:
            if res.name != ccd:
                continue
            atoms = [a.pos for a in res if a.element.name not in ("H", "D")]
            dist_a = min(p.dist(q) for p in atoms for q in ca["A"])
            dist_b = min(p.dist(q) for p in atoms for q in ca["B"])
            score = max(dist_a, dist_b)
            if best_score is None or score < best_score:
                best_score, best = score, (chain.name, int(res.seqid.num),
                                           round(dist_a, 2), round(dist_b, 2))
    assert best is not None, f"no {ccd} copy found"
    return best


def build_ligand(structure, ccd, chain_name, seqid, ccd_cif):
    """Build an RDKit molecule for one deposited PAM copy.

    Heavy-atom coordinates come from the deposit. Element, formal charge and
    bond order come from the RCSB chemical component record. Nothing is
    inferred from interatomic distances.
    """
    model = structure[0]
    matches = [res for chain in model if chain.name == chain_name
               for res in chain
               if res.name == ccd and int(res.seqid.num) == seqid]
    assert len(matches) == 1, f"expected one {ccd} {chain_name}/{seqid}"
    res = matches[0]

    atoms = [a for a in res if a.element.name not in ("H", "D")]
    names = [a.name.strip() for a in atoms]
    assert len(set(names)) == len(names), "duplicate atom ids in deposit"
    xyz = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in atoms], dtype=float)

    block = gemmi.cif.read_file(str(ccd_cif)).sole_block()
    atom_rows = {row[0]: row for row in block.find(
        ["_chem_comp_atom.atom_id", "_chem_comp_atom.type_symbol",
         "_chem_comp_atom.charge"])}
    bond_rows = list(block.find(
        ["_chem_comp_bond.atom_id_1", "_chem_comp_bond.atom_id_2",
         "_chem_comp_bond.value_order"]))
    assert len(bond_rows) >= 15, f"{ccd} bond table incomplete"

    rw, index = Chem.RWMol(), {}
    for name in names:
        assert name in atom_rows, f"atom {name} absent from {ccd} record"
        row = atom_rows[name]
        atom = Chem.Atom(row[1].capitalize())
        atom.SetFormalCharge(0 if row[2] in (".", "?", "") else int(row[2]))
        index[name] = rw.AddAtom(atom)

    orders = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE,
              "TRIP": Chem.BondType.TRIPLE, "AROM": Chem.BondType.AROMATIC}
    for n1, n2, order in bond_rows:
        if n1 in index and n2 in index:
            btype = orders[order.upper()]
            rw.AddBond(index[n1], index[n2], btype)
            if btype == Chem.BondType.AROMATIC:
                rw.GetAtomWithIdx(index[n1]).SetIsAromatic(True)
                rw.GetAtomWithIdx(index[n2]).SetIsAromatic(True)

    mol = rw.GetMol()
    conf = Chem.Conformer(len(names))
    for name, coord in zip(names, xyz):
        conf.SetAtomPosition(index[name], [float(v) for v in coord])
    mol.AddConformer(conf)
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistryFrom3D(mol)
    mol.SetProp("_Name", f"{ccd} {chain_name}{seqid}")
    return mol


# --------------------------------------------------------------------------
# common coordinate frame
# --------------------------------------------------------------------------

def ab_backbone(structure):
    """Ordered CA atoms of chains A and B, keyed by (chain, residue number)."""
    model = structure[0]
    out = {}
    for chain in model:
        if chain.name not in ("A", "B"):
            continue
        for res in chain:
            if res.het_flag != "A":
                continue
            for atom in res:
                if atom.name == "CA":
                    out[(chain.name, int(res.seqid.num))] = atom.pos
                    break
    return out


def superpose_to_reference(mobile, reference):
    """Superpose `mobile` onto `reference` on matched A-B CA atoms.

    The transform is applied in place to every atom of `mobile`, so the whole
    deposit - protein and ligands together - moves into the reference frame.
    """
    mob_ca, ref_ca = ab_backbone(mobile), ab_backbone(reference)
    shared = sorted(set(mob_ca) & set(ref_ca))
    assert len(shared) > 300, f"only {len(shared)} matched CA atoms"

    p_mob = [mob_ca[k] for k in shared]
    p_ref = [ref_ca[k] for k in shared]

    sup = gemmi.superpose_positions(p_ref, p_mob)
    for model in mobile:
        for chain in model:
            for res in chain:
                for atom in res:
                    atom.pos = gemmi.Position(sup.transform.apply(atom.pos))
    return sup.rmsd, len(shared)


# --------------------------------------------------------------------------
# receptor and ligand preparation
# --------------------------------------------------------------------------

def prepare_receptor(pdb, structure):
    """Rigid A-B polymer receptor, all non-polymer species removed."""
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    st = structure.clone()
    model = st[0]
    for name in [c.name for c in model]:
        if name not in ("A", "B"):
            model.remove_chain(name)
    model.remove_ligands_and_waters()
    model.remove_alternative_conformations()

    counts = {c.name: sum(r.het_flag == "A" for r in c) for c in model}
    assert all(counts.get(c, 0) > 150 for c in ("A", "B")), \
        f"{pdb}: A/B extraction too small: {counts}"

    raw = WORK / f"{pdb}_AB_raw.pdb"
    st.write_pdb(str(raw))

    fixer = PDBFixer(filename=str(raw))
    fixer.findMissingResidues()
    missing_loops = {str(k): list(v) for k, v in fixer.missingResidues.items()}
    fixer.missingResidues = {}          # do not build unresolved loops
    fixer.findNonstandardResidues()
    if fixer.nonstandardResidues:
        raise RuntimeError(f"{pdb}: nonstandard residues {fixer.nonstandardResidues}")
    fixer.findMissingAtoms()
    rebuilt = sum(len(v) for v in fixer.missingAtoms.values())
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=PH)

    prepared = WORK / f"{pdb}_AB_H.pdb"
    with prepared.open("w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)

    # PDBFixer caps chain termini with OXT, which is absent from Meeko's
    # internal-residue templates and triggers a valence error. No OXT atom is
    # anywhere near the transmembrane PAM site.
    _lines = prepared.read_text().splitlines(True)
    _kept = [l for l in _lines
             if not (l.startswith(("ATOM", "HETATM"))
                     and l[12:16].strip() == "OXT")]
    if len(_kept) != len(_lines):
        log(f"  {pdb}: removed {len(_lines) - len(_kept)} OXT cap atom(s)")
    prepared.write_text("".join(_kept))

    prefix = WORK / f"{pdb}_AB_receptor"
    pdbqt = WORK / f"{pdb}_AB_receptor.pdbqt"
    dropped_residues = meeko_receptor_with_retries(pdb, prepared, prefix)
    if not pdbqt.exists():
        candidates = list(WORK.glob(f"{pdb}_AB_receptor*.pdbqt"))
        assert len(candidates) == 1, f"{pdb}: ambiguous receptor pdbqt {candidates}"
        shutil.copyfile(candidates[0], pdbqt)

    coords = []
    for line in pdbqt.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    coords = np.asarray(coords)
    assert len(coords) > 1000, f"{pdb}: prepared receptor too small"

    centre = np.array(BOX["center"])
    near = int((np.linalg.norm(coords - centre, axis=1) < 14.0).sum())
    assert near >= 50, f"{pdb}: shared box is empty after superposition ({near} atoms)"

    return {"pdbqt": str(pdbqt), "atoms": int(len(coords)),
            "dropped_unbuildable_residues": dropped_residues,
            "atoms_near_box_centre": near,
            "rebuilt_sidechain_atoms": rebuilt,
            "unbuilt_missing_loops": missing_loops,
            "sha256": digest(pdbqt)}


def prepare_ligand_pdbqt(pdb, mol):
    """Protonated ligand in deposited geometry, exported to PDBQT."""
    protonated = Chem.AddHs(mol, addCoords=True)
    props = AllChem.MMFFGetMoleculeProperties(protonated)
    if props is not None:
        ff = AllChem.MMFFGetMoleculeForceField(protonated, props)
        for i in range(mol.GetNumAtoms()):
            ff.AddFixedPoint(i)
        ff.Minimize(maxIts=250)

    sdf = WORK / f"{pdb}_{PAMS[pdb]}_native_H.sdf"
    writer = Chem.SDWriter(str(sdf))
    writer.write(protonated)
    writer.close()

    heavy = WORK / f"{pdb}_{PAMS[pdb]}_native_heavy.sdf"
    writer = Chem.SDWriter(str(heavy))
    writer.write(mol)
    writer.close()

    pdbqt = WORK / f"{pdb}_{PAMS[pdb]}_ligand.pdbqt"
    run_tool(["mk_prepare_ligand.py", "-i", sdf, "-o", pdbqt],
             f"meeko_ligand_{pdb}.log")
    assert pdbqt.stat().st_size > 300
    return {"pdbqt": str(pdbqt), "native_heavy_sdf": str(heavy),
            "heavy_atoms": mol.GetNumHeavyAtoms(),
            "smiles": Chem.MolToSmiles(mol), "sha256": digest(pdbqt)}


# --------------------------------------------------------------------------
# stage 1 - preparation
# --------------------------------------------------------------------------

def stage_prepare():
    for folder in (INPUT, WORK, RESULTS, CKPT):
        folder.mkdir(parents=True, exist_ok=True)

    log("fetching deposits and chemical component records")
    for pdb, ccd in PAMS.items():
        fetch(f"https://files.rcsb.org/download/{pdb}.cif", INPUT / f"{pdb}.cif")
        fetch(f"https://files.rcsb.org/ligands/download/{ccd}.cif", INPUT / f"{ccd}.cif")

    structures = {}
    for pdb in PAMS:
        st = gemmi.read_structure(str(INPUT / f"{pdb}.cif"))
        st.setup_entities()
        structures[pdb] = st

    reference = structures[REFERENCE_FRAME]
    manifest = {"reference_frame": REFERENCE_FRAME, "box": BOX,
                "receptors": {}, "ligands": {}, "superposition": {}}

    log("superposing receptors onto the", REFERENCE_FRAME, "frame")
    for pdb, st in structures.items():
        if pdb == REFERENCE_FRAME:
            manifest["superposition"][pdb] = {"rmsd_A": 0.0, "matched_CA": None,
                                              "note": "reference frame"}
            continue
        rmsd, matched = superpose_to_reference(st, reference)
        manifest["superposition"][pdb] = {"rmsd_A": round(float(rmsd), 4),
                                          "matched_CA": matched}
        log(f"  {pdb}: {matched} matched CA, global RMSD {rmsd:.3f} A")

    log("selecting interface PAM copies and building ligands")
    for pdb, ccd in PAMS.items():
        chain, seqid, d_a, d_b = interface_pam_copy(structures[pdb], ccd)
        log(f"  {pdb} {ccd}: interface copy {chain}/{seqid} "
            f"(min dist to A {d_a} A, to B {d_b} A)")
        mol = build_ligand(structures[pdb], ccd, chain, seqid, INPUT / f"{ccd}.cif")
        info = prepare_ligand_pdbqt(pdb, mol)
        info.update({"ccd": ccd, "interface_copy": f"{chain}/{seqid}",
                     "min_dist_chain_A": d_a, "min_dist_chain_B": d_b})
        manifest["ligands"][pdb] = info

        centre = np.array(BOX["center"])
        half = np.array(BOX["size"]) / 2.0
        conf = mol.GetConformer()
        pts = np.array([list(conf.GetAtomPosition(i))
                        for i in range(mol.GetNumAtoms())])
        clearance = float(min((pts - (centre - half)).min(),
                              ((centre + half) - pts).min()))
        info["box_clearance_A"] = round(clearance, 3)
        log(f"    native pose box clearance {clearance:.2f} A")

    log("preparing receptors")
    for pdb in PAMS:
        manifest["receptors"][pdb] = prepare_receptor(pdb, structures[pdb])
        log(f"  {pdb}: {manifest['receptors'][pdb]['atoms']} receptor atoms, "
            f"{manifest['receptors'][pdb]['atoms_near_box_centre']} near box centre")

    manifest["packages"] = {p: metadata.version(p) for p in
                            ["vina", "meeko", "gemmi", "rdkit", "pdbfixer",
                             "openmm", "numpy", "scipy"]}
    manifest["python"] = sys.version
    manifest["platform"] = platform.platform()
    manifest["receptor_composition"] = (
        "Rigid A-B polymer only. Epibatidine, PAMs, lipids, waters, glycans "
        "and ions excluded uniformly from all four receptors.")
    manifest["rmsd_convention"] = (
        "Symmetry-aware, in-place heavy-atom RMSD (RDKit CalcRMS). Poses are "
        "never superposed onto the reference; all receptors and reference "
        "ligands share the 8V82 frame.")
    manifest["retrieved_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

    (RESULTS / "manifest.json").write_text(json.dumps(manifest, indent=2))
    (DRIVE_OUT / "matrix_manifest.json").write_text(json.dumps(manifest, indent=2))
    log("preparation complete")
    return manifest


# --------------------------------------------------------------------------
# stage 2 - docking
# --------------------------------------------------------------------------

def native_mol(ligand_pdb, manifest):
    path = manifest["ligands"][ligand_pdb]["native_heavy_sdf"]
    mol = Chem.SDMolSupplier(path, removeHs=True)[0]
    assert mol is not None
    return mol


def dock_pair(ligand_pdb, receptor_pdb, seed, manifest):
    """One (ligand, receptor, seed) condition."""
    tag = f"{ligand_pdb}_{PAMS[ligand_pdb]}__into__{receptor_pdb}_seed{seed}"
    ckpt = CKPT / f"{tag}.json"
    if ckpt.exists():
        return json.loads(ckpt.read_text())

    from vina import Vina
    started = time.time()

    v = Vina(sf_name="vina", cpu=VINA_CPU, seed=seed, verbosity=0)
    v.set_receptor(manifest["receptors"][receptor_pdb]["pdbqt"])
    v.set_ligand_from_file(manifest["ligands"][ligand_pdb]["pdbqt"])
    v.compute_vina_maps(center=BOX["center"], box_size=BOX["size"])
    v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=NMODES)

    scores = v.energies(n_poses=NMODES, energy_range=30)
    pose_pdbqt = RESULTS / f"{tag}.pdbqt"
    v.write_poses(str(pose_pdbqt), n_poses=NMODES, energy_range=30, overwrite=True)

    pose_sdf = RESULTS / f"{tag}.sdf"
    run_tool(["mk_export.py", pose_pdbqt, "-s", pose_sdf], f"export_{tag}.log")
    poses = [Chem.RemoveHs(m) for m in
             Chem.SDMolSupplier(str(pose_sdf), removeHs=False) if m is not None]
    assert poses, f"{tag}: no poses exported"
    assert len(poses) == len(scores), f"{tag}: score/pose count mismatch"

    native = native_mol(ligand_pdb, manifest)
    rows = []
    for rank, (pose, score) in enumerate(zip(poses, scores), start=1):
        assert pose.GetNumHeavyAtoms() == native.GetNumHeavyAtoms()
        rmsd = float(rdMolAlign.CalcRMS(pose, native, maxMatches=10000))
        rows.append({"ligand": ligand_pdb, "ccd": PAMS[ligand_pdb],
                     "receptor": receptor_pdb,
                     "condition": "cognate" if ligand_pdb == receptor_pdb else "cross",
                     "seed": seed, "rank": rank,
                     "score_kcal_mol": round(float(score[0]), 4),
                     "rmsd_A": round(rmsd, 4),
                     "passes_2A": bool(rmsd < RMSD_CUTOFF)})

    record = {
        "tag": tag, "ligand": ligand_pdb, "receptor": receptor_pdb, "seed": seed,
        "condition": rows[0]["condition"],
        "exhaustiveness": EXHAUSTIVENESS, "n_modes": NMODES,
        "rank1_rmsd_A": rows[0]["rmsd_A"],
        "rank1_score_kcal_mol": rows[0]["score_kcal_mol"],
        "best_sampled_rmsd_A": min(r["rmsd_A"] for r in rows),
        "best_sampled_rank": min(rows, key=lambda r: r["rmsd_A"])["rank"],
        "top1_success": bool(rows[0]["rmsd_A"] < RMSD_CUTOFF),
        "any_pose_success": any(r["passes_2A"] for r in rows),
        "wall_seconds": round(time.time() - started, 1),
        "pose_sdf_sha256": digest(pose_sdf),
        "rows": rows,
    }
    ckpt.write_text(json.dumps(record, indent=1))
    log(f"  {tag}: rank1 {record['rank1_rmsd_A']:.3f} A, "
        f"best {record['best_sampled_rmsd_A']:.3f} A "
        f"(rank {record['best_sampled_rank']}), "
        f"{record['wall_seconds']:.0f}s")
    return record


def stage_regression(manifest):
    """Reproduce the known cognate 8V82 / PNU-120596 result before the matrix.

    If the generalized pipeline cannot recover the Module 7 numbers, nothing
    downstream of it can be trusted, so this runs first and stops on failure.
    """
    log("regression check: PNU-120596 into cognate 8V82, seed 2026")
    rec = dock_pair("8V82", "8V82", 2026, manifest)
    r1, r2 = rec["rows"][0]["rmsd_A"], rec["rows"][1]["rmsd_A"]
    lo1, hi1 = REGRESSION["rank1_range"]
    lo2, hi2 = REGRESSION["rank2_range"]
    ok = (lo1 <= r1 <= hi1) and (lo2 <= r2 <= hi2)
    log(f"  rank1 {r1:.3f} A (expected {lo1}-{hi1}), "
        f"rank2 {r2:.3f} A (expected {lo2}-{hi2}) -> "
        f"{'MATCH' if ok else 'MISMATCH'}")
    (DRIVE_OUT / "regression_check.json").write_text(json.dumps(
        {"rank1_rmsd_A": r1, "rank2_rmsd_A": r2,
         "expected_rank1": list(REGRESSION["rank1_range"]),
         "expected_rank2": list(REGRESSION["rank2_range"]),
         "reproduces_module7_baseline": ok,
         "wall_seconds": rec["wall_seconds"]}, indent=1))
    return ok, rec


def stage_matrix(manifest):
    order = []
    # cognate conditions first, so the diagonal exists even if the run is cut short
    for pdb in PAMS:
        order.append((pdb, pdb))
    for lig in PAMS:
        for rec in PAMS:
            if lig != rec:
                order.append((lig, rec))

    total = len(order) * len(SEEDS)
    done = 0
    started = time.time()
    for lig, rec in order:
        for seed in SEEDS:
            dock_pair(lig, rec, seed, manifest)
            done += 1
            elapsed = time.time() - started
            rate = elapsed / max(done, 1)
            write_status({
                "completed_runs": done, "total_runs": total,
                "elapsed_seconds": round(elapsed, 1),
                "mean_seconds_per_run": round(rate, 1),
                "projected_remaining_seconds": round(rate * (total - done), 1),
                "current": f"{lig}->{rec} seed {seed}",
                "updated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            })
    log(f"matrix complete: {done}/{total} runs")


def stage_collect():
    records = [json.loads(p.read_text()) for p in sorted(CKPT.glob("*.json"))]
    rows = [r for rec in records for r in rec["rows"]]
    if not rows:
        log("no checkpoints to collect")
        return

    csv_path = DRIVE_OUT / "matrix_all_poses.csv"
    with csv_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    summary = []
    for rec in records:
        summary.append({k: rec[k] for k in
                        ("ligand", "receptor", "condition", "seed",
                         "rank1_rmsd_A", "rank1_score_kcal_mol",
                         "best_sampled_rmsd_A", "best_sampled_rank",
                         "top1_success", "any_pose_success", "wall_seconds")})
    summary_path = DRIVE_OUT / "matrix_summary.csv"
    with summary_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(summary[0].keys()))
        writer.writeheader()
        writer.writerows(summary)

    log(f"collected {len(records)} conditions, {len(rows)} poses")
    log(f"  {csv_path}")
    log(f"  {summary_path}")


# --------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=["prepare", "regression", "matrix",
                                          "collect", "all"])
    parser.add_argument("--force", action="store_true",
                        help="continue the matrix even if regression fails")
    args = parser.parse_args()

    manifest_path = RESULTS / "manifest.json"
    if args.stage == "prepare":
        stage_prepare()
        return

    if not manifest_path.exists():
        manifest = stage_prepare()
    else:
        manifest = json.loads(manifest_path.read_text())

    if args.stage in ("regression", "all"):
        ok, _ = stage_regression(manifest)
        if args.stage == "regression":
            return
        if not ok and not args.force:
            log("STOP: regression check did not reproduce the Module 7 baseline.")
            log("      The matrix was not run. Inspect regression_check.json.")
            return

    if args.stage in ("matrix", "all"):
        stage_matrix(manifest)

    stage_collect()


if __name__ == "__main__":
    main()
