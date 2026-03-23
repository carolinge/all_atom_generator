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
    return buf.getvalue()


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
