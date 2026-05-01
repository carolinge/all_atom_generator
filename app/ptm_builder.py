# -*- coding: utf-8 -*-
"""
ptm_builder.py
==============
Apply post-translational / post-transcriptional modifications to a structure
by placing new atoms using ideal internal coordinates (NeRF algorithm), then
letting OpenMM fill in hydrogens.

Supported modifications
-----------------------
Protein:
  SEP  -- phosphoserine   (SER -> SEP)
  TPO  -- phosphothreonine (THR -> TPO)
  PTR  -- phosphotyrosine  (TYR -> PTR)

RNA:
  8OG  -- 8-oxoguanosine  (G / RG -> 8OG)

Usage (programmatic)
--------------------
    from ptm_builder import apply_ptm
    new_pdb_str = apply_ptm(pdb_path, chain_id, resnum, ptm_code)
"""

from __future__ import annotations
import math
import numpy as np
from pathlib import Path
from typing import Sequence

# ── NeRF: build one atom from 3 reference atoms + IC ─────────────────────────

def _nerf(a: np.ndarray, b: np.ndarray, c: np.ndarray,
          bond: float, angle_deg: float, dihedral_deg: float) -> np.ndarray:
    """
    Natural Extension Reference Frame.
    Given three known atoms a, b, c, return the position of a new atom d such
    that:
        |c-d|  = bond  (Angstrom)
        angle(b,c,d) = angle_deg
        dihedral(a,b,c,d) = dihedral_deg
    """
    angle   = math.radians(angle_deg)
    dihedral = math.radians(dihedral_deg)

    bc = (c - b) / np.linalg.norm(c - b)
    ab = (b - a) / np.linalg.norm(b - a)

    n = np.cross(ab, bc)
    if np.linalg.norm(n) < 1e-8:          # degenerate: a,b,c collinear
        n = np.array([0.0, 0.0, 1.0])
    n /= np.linalg.norm(n)

    m = np.cross(n, bc)

    d = c + bond * (
        -math.cos(angle) * bc
        + math.sin(angle) * math.cos(dihedral) * m
        + math.sin(angle) * math.sin(dihedral) * n
    )
    return d


# ── Patch definitions ─────────────────────────────────────────────────────────
#
# Each patch is a dict with:
#   "rename"   : new residue name
#   "delete"   : atom names to remove from the residue
#   "add"      : list of dicts, each describing one new atom:
#       "name"     : atom name in the new residue
#       "element"  : element symbol
#       "ref"      : (name1, name2, name3) -- three existing atoms used as
#                    reference for NeRF (name3 is the bonded atom)
#       "bond"     : bond length  name3--new  (Angstrom)
#       "angle"    : bond angle   name2-name3-new  (degrees)
#       "dihedral" : dihedral     name1-name2-name3-new (degrees)
#
# IC values from:
#   Phospho: CHARMM36 PRES PHOSPHO / AMBER99SB SEP/TPO/PTR templates
#   8OG:     Aduri 2007 JCTC + AMBER GAFF ideal geometry

PATCHES: dict[str, dict] = {

    # ── Phosphoserine (SER -> SEP) ────────────────────────────────────────────
    "SEP": {
        "rename": "SEP",
        "from":   "SER",
        "delete": ["HG"],          # remove serine hydroxyl H
        "add": [
            # OG -- P
            {
                "name": "P",  "element": "P",
                "ref":  ("CA", "CB", "OG"),
                "bond": 1.610, "angle": 120.0, "dihedral": 180.0,
            },
            # P -- O1P (non-bridging, anti to CB)
            {
                "name": "O1P", "element": "O",
                "ref":  ("CA", "CB", "OG"),
                "bond": 1.503, "angle": 108.0, "dihedral":  60.0,
            },
            # P -- O2P (non-bridging)
            {
                "name": "O2P", "element": "O",
                "ref":  ("CA", "CB", "OG"),
                "bond": 1.503, "angle": 108.0, "dihedral": -60.0,
            },
            # P -- O3P (non-bridging)
            {
                "name": "O3P", "element": "O",
                "ref":  ("CB", "OG", "P"),
                "bond": 1.503, "angle": 108.0, "dihedral": 180.0,
            },
        ],
    },

    # ── Phosphothreonine (THR -> TPO) ─────────────────────────────────────────
    "TPO": {
        "rename": "TPO",
        "from":   "THR",
        "delete": ["HG1"],
        "add": [
            {
                "name": "P",  "element": "P",
                "ref":  ("CA", "CB", "OG1"),
                "bond": 1.610, "angle": 120.0, "dihedral": 180.0,
            },
            {
                "name": "O1P", "element": "O",
                "ref":  ("CA", "CB", "OG1"),
                "bond": 1.503, "angle": 108.0, "dihedral":  60.0,
            },
            {
                "name": "O2P", "element": "O",
                "ref":  ("CA", "CB", "OG1"),
                "bond": 1.503, "angle": 108.0, "dihedral": -60.0,
            },
            {
                "name": "O3P", "element": "O",
                "ref":  ("CB", "OG1", "P"),
                "bond": 1.503, "angle": 108.0, "dihedral": 180.0,
            },
        ],
    },

    # ── Phosphotyrosine (TYR -> PTR) ──────────────────────────────────────────
    "PTR": {
        "rename": "PTR",
        "from":   "TYR",
        "delete": ["HH"],
        "add": [
            {
                "name": "P",  "element": "P",
                "ref":  ("CE1", "CZ", "OH"),
                "bond": 1.610, "angle": 120.0, "dihedral": 180.0,
            },
            {
                "name": "O1P", "element": "O",
                "ref":  ("CE1", "CZ", "OH"),
                "bond": 1.503, "angle": 108.0, "dihedral":  60.0,
            },
            {
                "name": "O2P", "element": "O",
                "ref":  ("CE1", "CZ", "OH"),
                "bond": 1.503, "angle": 108.0, "dihedral": -60.0,
            },
            {
                "name": "O3P", "element": "O",
                "ref":  ("CZ", "OH", "P"),
                "bond": 1.503, "angle": 108.0, "dihedral": 180.0,
            },
        ],
    },

    # ── 8-Oxoguanosine (G/RG -> 8OG) ─────────────────────────────────────────
    # Structural changes vs G:
    #   C8: aromatic C-H  ->  C8=O8 (sp2 carbonyl)
    #   N7: aromatic N    ->  N7-H  (sp2 N-H)
    # Delete H8, add O8 (at C8) and H7 (at N7).
    # IC from Aduri 2007 + AMBER RNA ideal geometry.
    "8OG": {
        "rename": "8OG",
        "from":   "G",         # also matches RG, G3, G5 -- handled in code
        "delete": ["H8"],
        "add": [
            # C8=O8: place O8 using N7-C8 bond direction
            # dihedral N9-N7-C8-O8 ~ 180 deg (O8 trans to N9 across C8)
            {
                "name": "O8",  "element": "O",
                "ref":  ("N9", "N7", "C8"),
                "bond": 1.230, "angle": 128.8, "dihedral": 180.0,
            },
            # N7-H7: place H7 in the imidazole plane
            # dihedral C8-C5-N7-H7 ~ 180 deg (H7 trans to C8 across N7)
            {
                "name": "H7",  "element": "H",
                "ref":  ("C8", "C5", "N7"),
                "bond": 1.010, "angle": 125.2, "dihedral": 180.0,
            },
        ],
    },
}

# ── Full intra-residue bond lists for PTM products ────────────────────────────
# OpenMM's createStandardBonds() ignores residues whose names aren't in its
# built-in table.  For non-standard PTM residues we must emit CONECT records
# for *every* internal bond, otherwise ForceField template matching fails with
# "atoms match X but bonds are different".

# Standard guanosine (G) intra-residue bonds, plus 8OG patch deltas:
#   delete (C8,H8); add (C8,O8) and (N7,H7)
_BONDS_8OG = [
    # Phosphate
    ("P", "OP1"), ("P", "OP2"), ("P", "O5'"),
    # Sugar backbone
    ("O5'", "C5'"), ("C5'", "C4'"), ("C5'", "H5'"), ("C5'", "H5''"),
    ("C4'", "O4'"), ("C4'", "C3'"), ("C4'", "H4'"),
    ("O4'", "C1'"),
    ("C3'", "O3'"), ("C3'", "C2'"), ("C3'", "H3'"),
    ("C2'", "O2'"), ("C2'", "C1'"), ("C2'", "H2'"),
    ("O2'", "HO2'"),
    ("C1'", "N9"), ("C1'", "H1'"),
    # Purine ring (8OG: no H8; instead C8=O8 and N7-H7)
    ("N9", "C8"), ("N9", "C4"),
    ("C8", "N7"), ("C8", "O8"),       # 8OG: O8 replaces H8
    ("N7", "C5"), ("N7", "H7"),       # 8OG: H7 added
    ("C5", "C4"), ("C5", "C6"),
    ("C4", "N3"),
    ("N3", "C2"),
    ("C2", "N1"), ("C2", "N2"),
    ("N2", "H21"), ("N2", "H22"),
    ("N1", "H1"), ("N1", "C6"),
    ("C6", "O6"),
]

# 5' terminal 8OG (8OG5): no incoming P from 5' side, gets H5T on O5'.
# AMBER 5'-terminal RNA uses HO5' instead of phosphate.
_BONDS_8OG5 = [b for b in _BONDS_8OG if "P" not in b and b != ("P", "OP1") and b != ("P", "OP2") and b != ("P", "O5'")]
_BONDS_8OG5 = [b for b in _BONDS_8OG5 if b != ("O5'", "C5'")] + [("O5'", "C5'"), ("O5'", "HO5'")]

# 3' terminal 8OG (8OG3): O3' has HO3' instead of bond to next residue's P
_BONDS_8OG3 = _BONDS_8OG + [("O3'", "HO3'")]

# Standard serine-based phosphorylation (SEP/TPO/PTR) internal bonds.
# These residues are actually in OpenMM's amber14 standard tables, but listed
# here for completeness in case future PTMs need them.
_BONDS_BY_PTM: dict[str, list[tuple[str, str]]] = {
    "8OG":  _BONDS_8OG,
    "8OG5": _BONDS_8OG5,
    "8OG3": _BONDS_8OG3,
}


# Aliases: G5, G3 (terminal RNA), RG (some PDB naming) all map to 8OG patch
_RESNAME_TO_PATCH: dict[str, str] = {
    "SER": "SEP",
    "THR": "TPO",
    "TYR": "PTR",
    "G":   "8OG",
    "RG":  "8OG",
    "G5":  "8OG",
    "G3":  "8OG",
    "GUA": "8OG",
}

# Target residue names for terminal RNA variants
_8OG_TERMINAL: dict[str, str] = {
    "G5": "8OG5",
    "G3": "8OG3",
}


# ── Standard nucleotide names (OpenMM has templates for these) ────────────────
_STANDARD_NA = {
    "A", "C", "G", "U", "DA", "DC", "DG", "DT",        # internal
    "A3", "C3", "G3", "U3", "DA3", "DC3", "DG3", "DT3", # 3' terminal
    "A5", "C5", "G5", "U5", "DA5", "DC5", "DG5", "DT5", # 5' terminal
    "RA", "RC", "RG", "RU",                               # alternate RNA names
}


def _add_backbone_conect(pdb_str: str) -> str:
    """Add CONECT records for non-standard PTM residues.

    Two kinds of CONECT records are emitted:

    1. **Inter-residue backbone bond** (e.g. O3'(prev) -> P(modified)).
       OpenMM's ``createStandardBonds()`` only creates these when *both*
       residues have standard templates; if one neighbour is non-standard
       (e.g. 8OG), the bond is silently dropped.

    2. **Full intra-residue bonds** for non-standard residues (looked up in
       ``_BONDS_BY_PTM``).  Without these, ForceField template matching
       fails because the atoms match the template but no bonds are seen.

    OpenMM's PDBFile reader honours CONECT records.
    """
    # Collect per-residue info in original PDB order.
    # key = (chain, resseq), value = {atom_name: serial, resname, atoms: {name: serial}}
    residues: dict[tuple, dict] = {}
    order: list[tuple] = []

    for line in pdb_str.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        # PDB fixed-width columns
        serial = int(line[6:11])
        aname  = line[12:16].strip()
        chain  = line[21]
        try:
            resseq = int(line[22:26])
        except ValueError:
            continue
        resname = line[17:20].strip()
        if not resname:
            resname = line[17:21].strip()

        key = (chain, resseq)
        if key not in residues:
            residues[key] = {"resname": resname, "atoms": {}}
            order.append(key)
        residues[key]["atoms"][aname] = serial

    # ── (1) inter-residue backbone bonds across non-standard residues ─────────
    conect_pairs: list[tuple[int, int]] = []
    for i in range(len(order) - 1):
        k1, k2 = order[i], order[i + 1]
        if k1[0] != k2[0]:          # different chains
            continue
        r1, r2 = residues[k1], residues[k2]
        n1 = r1["resname"].upper()
        n2 = r2["resname"].upper()

        if n1 in _STANDARD_NA and n2 in _STANDARD_NA:
            continue

        s_o3 = r1["atoms"].get("O3'")
        s_p  = r2["atoms"].get("P")
        if s_o3 is not None and s_p is not None:
            conect_pairs.append((s_o3, s_p))

    # ── (2) intra-residue bonds for non-standard PTM residues ─────────────────
    for key in order:
        rname = residues[key]["resname"].upper()
        bond_list = _BONDS_BY_PTM.get(rname)
        if not bond_list:
            continue
        atoms = residues[key]["atoms"]
        for a1, a2 in bond_list:
            s1, s2 = atoms.get(a1), atoms.get(a2)
            if s1 is not None and s2 is not None:
                conect_pairs.append((s1, s2))

    if not conect_pairs:
        return pdb_str

    # Build CONECT lines (pair-per-line for clarity & compatibility)
    conect_lines = []
    seen = set()
    for s1, s2 in conect_pairs:
        pair = (min(s1, s2), max(s1, s2))
        if pair in seen:
            continue
        seen.add(pair)
        conect_lines.append(f"CONECT{s1:5d}{s2:5d}")
        conect_lines.append(f"CONECT{s2:5d}{s1:5d}")

    lines = pdb_str.rstrip().splitlines()
    result = []
    for line in lines:
        if line.startswith("END"):
            result.extend(conect_lines)
        result.append(line)
    # If no END line was found, append CONECT at the end
    if not any(l.startswith("END") for l in lines):
        result.extend(conect_lines)

    return "\n".join(result) + "\n"


# ── BioPython-based structure modification ────────────────────────────────────

def apply_ptm(pdb_path: str | Path,
              chain_id: str,
              resnum: int | str,
              ptm_code: str) -> str:
    """
    Apply a PTM to a single residue and return the modified structure as a
    PDB-format string.

    Parameters
    ----------
    pdb_path : path to input PDB file
    chain_id : chain identifier (single character)
    resnum   : residue sequence number (int or str)
    ptm_code : one of SEP, TPO, PTR, 8OG

    Returns
    -------
    PDB string of the modified structure.

    Raises
    ------
    ValueError  if the residue is not found or the patch is incompatible.
    KeyError    if ptm_code is not recognised.
    """
    import Bio.PDB as bpdb
    from io import StringIO

    ptm_code = ptm_code.upper()
    if ptm_code not in PATCHES:
        raise KeyError(f"Unknown PTM code: {ptm_code!r}. "
                       f"Known: {list(PATCHES)}")

    patch = PATCHES[ptm_code]
    resnum = int(resnum)

    # -- parse structure --
    parser = bpdb.PDBParser(QUIET=True)
    structure = parser.get_structure("mol", str(pdb_path))
    model = structure[0]

    # locate target residue
    target_res = None
    for chain in model:
        if chain.id != chain_id:
            continue
        for res in chain:
            if res.id[1] == resnum:
                target_res = res
                break
        if target_res:
            break

    if target_res is None:
        raise ValueError(
            f"Residue {chain_id}{resnum} not found in {pdb_path}")

    orig_name = target_res.resname.strip()

    # -- validate compatible residue --
    expected = patch.get("from", "")
    # For 8OG accept all G variants
    if ptm_code == "8OG":
        if orig_name not in _RESNAME_TO_PATCH or _RESNAME_TO_PATCH[orig_name] != "8OG":
            raise ValueError(
                f"8OG patch requires a guanosine residue (G/RG/G3/G5), "
                f"got {orig_name!r}")
    else:
        if orig_name != expected:
            raise ValueError(
                f"{ptm_code} patch requires {expected!r}, got {orig_name!r}")

    # -- build atom coordinate lookup --
    atom_coords: dict[str, np.ndarray] = {
        a.name: np.array(a.coord, dtype=float)
        for a in target_res.get_atoms()
    }

    # -- delete atoms --
    for del_name in patch["delete"]:
        if del_name in atom_coords:
            del atom_coords[del_name]
            if target_res.has_id(del_name):
                target_res.detach_child(del_name)

    # -- add new atoms via NeRF --
    serial_offset = max(
        (a.serial_number for a in structure.get_atoms()), default=0)

    for atom_def in patch["add"]:
        name = atom_def["name"]
        ref1, ref2, ref3 = atom_def["ref"]

        # Check references exist
        missing = [r for r in (ref1, ref2, ref3) if r not in atom_coords]
        if missing:
            raise ValueError(
                f"Reference atoms {missing} not found in residue "
                f"{chain_id}{resnum} ({orig_name}) for building {name}. "
                f"Available: {sorted(atom_coords)}")

        pos = _nerf(
            atom_coords[ref1], atom_coords[ref2], atom_coords[ref3],
            atom_def["bond"], atom_def["angle"], atom_def["dihedral"],
        )

        serial_offset += 1
        new_atom = bpdb.Atom.Atom(
            name=name,
            coord=pos,
            bfactor=0.0,
            occupancy=1.0,
            altloc=" ",
            fullname=f" {name:<3s}",
            serial_number=serial_offset,
            element=atom_def["element"],
        )
        target_res.add(new_atom)
        atom_coords[name] = pos

    # -- rename residue --
    new_name = patch["rename"]
    # Handle terminal RNA variants
    if ptm_code == "8OG" and orig_name in _8OG_TERMINAL:
        new_name = _8OG_TERMINAL[orig_name]
    target_res.resname = new_name

    # -- serialise back to PDB string --
    io = bpdb.PDBIO()
    io.set_structure(structure)
    buf = StringIO()
    io.save(buf)
    pdb_str = buf.getvalue()

    # -- fix backbone bonds: add CONECT records for O3'->P across modified
    #    residues.  OpenMM's createStandardBonds() skips non-standard residue
    #    names, so the O3'(prev)->P(modified) bond gets lost.  CONECT records
    #    are honoured by OpenMM's PDBFile reader and restore these bonds.
    pdb_str = _add_backbone_conect(pdb_str)
    return pdb_str


# ── Convenience: list what PTMs are applicable to a residue name ──────────────

def applicable_ptms(resname: str) -> list[str]:
    """Return list of PTM codes that can be applied to this residue."""
    resname = resname.strip().upper()
    code = _RESNAME_TO_PATCH.get(resname)
    return [code] if code else []


# ── CLI entry point (for testing) ─────────────────────────────────────────────

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 5:
        print("Usage: python ptm_builder.py <pdb> <chain> <resnum> <ptm>")
        sys.exit(1)
    result = apply_ptm(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    print(result)
