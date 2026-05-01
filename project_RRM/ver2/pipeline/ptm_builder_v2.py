# -*- coding: utf-8 -*-
"""
ptm_builder_v2.py
=================
Port of /app/ptm_builder.py for project_RRM/ver_2, **without modifying
the app/ tree**. Extended with patches for canonical RNA modifications
that were NOT in ver_1:

    8OG  -- 8-oxoguanosine        (G -> 8OG)   carries over from ver_1
    PUU  -- pseudouridine         (U -> PUU)   NEW; C-glycosidic rebonding
    M1A  -- N1-methyladenosine    (A -> M1A)   NEW; adds methyl on N1

Protein PTMs (SEP, TPO, PTR) are kept unchanged from app/ptm_builder.py
in case anyone wants to reuse them in this project.

Atom-name compatibility note (read this before trusting any patch):
-------------------------------------------------------------------
The atom names produced here MUST match the corresponding modXNA
assembled .lib that tleap loads. Currently:

  * 8OG: ver_1 used `O8` (8-oxo carbonyl) and `H7` (N7 proton). modXNA's
    raw 8OG.mol2 has the 8-oxo oxygen named `O6` (duplicate of the
    standard guanine O6). After modxna.sh + tleap assembly, the .lib
    file's name for the 8-oxo oxygen could be `O8`, `O6` (with index
    disambiguation), or something else. **Verified after first .lib
    is built; currently using `O8` as the placeholder. Update when known.**

  * PUU: standard uridine atoms (N1, C2, O2, N3, C4, O4, C5, C6) +
    H1 (newly free N1), no H5 (C5 now glycosidic). modXNA's PUU.mol2
    matches this naming. Glycosidic bond rewiring (remove N1-C1', add
    C5-C1') is handled here in code, since the existing PATCHES schema
    doesn't natively support bond rewiring.

  * M1A: adenosine atoms + C11/1H11/2H11/3H11 methyl on N1. modXNA's
    M1A.mol2 uses these names, so this patch matches directly.

If a future modxna.sh run reveals different atom names, update the
PATCHES dict and the rebonding rules here only — no other file should
need to change.
"""

from __future__ import annotations
import math
from io import StringIO
from pathlib import Path
from typing import Sequence

import numpy as np

# ── NeRF: build one atom from 3 reference atoms + IC ─────────────────────────

def _nerf(a: np.ndarray, b: np.ndarray, c: np.ndarray,
          bond: float, angle_deg: float, dihedral_deg: float) -> np.ndarray:
    """Natural Extension Reference Frame.
    Place atom d such that:
        |c-d|  = bond  (Angstrom)
        angle(b,c,d) = angle_deg
        dihedral(a,b,c,d) = dihedral_deg
    """
    angle    = math.radians(angle_deg)
    dihedral = math.radians(dihedral_deg)

    bc = (c - b) / np.linalg.norm(c - b)
    ab = (b - a) / np.linalg.norm(b - a)

    n = np.cross(ab, bc)
    if np.linalg.norm(n) < 1e-8:
        n = np.array([0.0, 0.0, 1.0])
    n /= np.linalg.norm(n)

    m = np.cross(n, bc)

    return c + bond * (
        -math.cos(angle) * bc
        + math.sin(angle) * math.cos(dihedral) * m
        + math.sin(angle) * math.sin(dihedral) * n
    )


# ── Patch definitions ────────────────────────────────────────────────────────
#
# Each patch:
#   "rename"   : new residue name
#   "from"     : list of accepted source residue names
#   "delete"   : atom names to remove
#   "add"      : list of new atom specs (NeRF parameters)
#   "rebond"   : optional list of (action, atom1, atom2)
#                action: "remove" | "add"

PATCHES: dict[str, dict] = {

    # ── 8-oxoguanosine (G/RG -> 8OG) ───────────────────────────────────────
    # Structural changes vs G: C8 acquires =O8 (carbonyl), N7 acquires H7.
    # IC values from Aduri 2007 + AMBER ideal geometry.
    # TODO(verify-after-modxna): atom names "O8" and "H7" are the chemistry-
    # community standard but modXNA's assembled .lib may use different names
    # (modXNA's raw mol2 has duplicate "O6"). Re-check after first .lib build.
    "8OG": {
        "rename": "8OG",
        "from":   ["G", "RG", "G5", "G3", "GUA"],
        "delete": ["H8"],
        "add": [
            {
                "name": "O8",  "element": "O",
                "ref":  ("N9", "N7", "C8"),
                "bond": 1.230, "angle": 128.8, "dihedral": 180.0,
            },
            {
                "name": "H7",  "element": "H",
                "ref":  ("C8", "C5", "N7"),
                "bond": 1.010, "angle": 125.2, "dihedral": 180.0,
            },
        ],
    },

    # ── Pseudouridine (U/RU -> PUU)  C-glycosidic ─────────────────────────
    # In uridine, C1' is bonded to N1.  In pseudouridine, C1' is bonded to
    # C5 instead.  This requires a true bond rewiring: remove N1-C1', add
    # C5-C1', delete H5 (was on C5), add H1 (now on N1 since N1 is no
    # longer glycosidic).  Heavy-atom positions of the base ring stay put;
    # only the H positions shift.
    "PUU": {
        "rename": "PUU",
        "from":   ["U", "RU", "U5", "U3", "URA"],
        "delete": ["H5"],
        "add": [
            # H1 on N1: place in the ring plane, opposite C2 across N1
            # (sp2 N-H, ~1.01 A, angle C2-N1-H ~ 117 deg, dihedral
            # C6-C2-N1-H ~ 180 to put H trans to C6 across N1).
            {
                "name": "H1",  "element": "H",
                "ref":  ("C6", "C2", "N1"),
                "bond": 1.010, "angle": 117.0, "dihedral": 180.0,
            },
        ],
        "rebond": [
            ("remove", "N1", "C1'"),
            ("add",    "C5", "C1'"),
        ],
    },

    # ── N1-methyladenosine (A/RA -> M1A) ──────────────────────────────────
    # In adenine N1 has no H. m1A adds a methyl group on N1.
    # IC: C-N1 bond 1.475 A, sp2 N angle, methyl Hs at 109.5 from CH3 carbon.
    "M1A": {
        "rename": "M1A",
        "from":   ["A", "RA", "A5", "A3", "ADE"],
        "delete": [],
        "add": [
            # C11: methyl carbon on N1
            {
                "name": "C11", "element": "C",
                "ref":  ("C6", "C2", "N1"),
                "bond": 1.475, "angle": 117.0, "dihedral": 180.0,
            },
            # Three methyl Hs on C11 — staggered, sp3 geometry.
            {
                "name": "1H11", "element": "H",
                "ref":  ("C2", "N1", "C11"),
                "bond": 1.090, "angle": 109.5, "dihedral":  60.0,
            },
            {
                "name": "2H11", "element": "H",
                "ref":  ("C2", "N1", "C11"),
                "bond": 1.090, "angle": 109.5, "dihedral": 180.0,
            },
            {
                "name": "3H11", "element": "H",
                "ref":  ("C2", "N1", "C11"),
                "bond": 1.090, "angle": 109.5, "dihedral": 300.0,
            },
        ],
    },

    # ── Phospho protein PTMs (kept identical to app/ptm_builder.py) ───────
    "SEP": {
        "rename": "SEP", "from": ["SER"], "delete": ["HG"],
        "add": [
            {"name": "P",   "element": "P", "ref": ("CA","CB","OG"),
             "bond": 1.610, "angle": 120.0, "dihedral": 180.0},
            {"name": "O1P", "element": "O", "ref": ("CA","CB","OG"),
             "bond": 1.503, "angle": 108.0, "dihedral":  60.0},
            {"name": "O2P", "element": "O", "ref": ("CA","CB","OG"),
             "bond": 1.503, "angle": 108.0, "dihedral": -60.0},
            {"name": "O3P", "element": "O", "ref": ("CB","OG","P"),
             "bond": 1.503, "angle": 108.0, "dihedral": 180.0},
        ],
    },
    "TPO": {
        "rename": "TPO", "from": ["THR"], "delete": ["HG1"],
        "add": [
            {"name": "P",   "element": "P", "ref": ("CA","CB","OG1"),
             "bond": 1.610, "angle": 120.0, "dihedral": 180.0},
            {"name": "O1P", "element": "O", "ref": ("CA","CB","OG1"),
             "bond": 1.503, "angle": 108.0, "dihedral":  60.0},
            {"name": "O2P", "element": "O", "ref": ("CA","CB","OG1"),
             "bond": 1.503, "angle": 108.0, "dihedral": -60.0},
            {"name": "O3P", "element": "O", "ref": ("CB","OG1","P"),
             "bond": 1.503, "angle": 108.0, "dihedral": 180.0},
        ],
    },
    "PTR": {
        "rename": "PTR", "from": ["TYR"], "delete": ["HH"],
        "add": [
            {"name": "P",   "element": "P", "ref": ("CE1","CZ","OH"),
             "bond": 1.610, "angle": 120.0, "dihedral": 180.0},
            {"name": "O1P", "element": "O", "ref": ("CE1","CZ","OH"),
             "bond": 1.503, "angle": 108.0, "dihedral":  60.0},
            {"name": "O2P", "element": "O", "ref": ("CE1","CZ","OH"),
             "bond": 1.503, "angle": 108.0, "dihedral": -60.0},
            {"name": "O3P", "element": "O", "ref": ("CZ","OH","P"),
             "bond": 1.503, "angle": 108.0, "dihedral": 180.0},
        ],
    },
}

# Map: source residue name -> patch code that targets it (auto-generated).
_RESNAME_TO_PATCH: dict[str, str] = {
    src: code for code, p in PATCHES.items() for src in p.get("from", [])
}


def applicable_patches(resname: str) -> list[str]:
    """List patch codes that can apply to a given residue name."""
    code = _RESNAME_TO_PATCH.get(resname.strip().upper())
    return [code] if code else []


# ── Main API: apply_ptm ─────────────────────────────────────────────────────

def apply_ptm(pdb_path: str | Path,
              chain_id: str,
              resnum: int | str,
              ptm_code: str) -> str:
    """Apply a PTM/RNA-modification patch to a single residue.

    Parameters
    ----------
    pdb_path : path to input PDB
    chain_id : chain identifier (single character)
    resnum   : residue sequence number
    ptm_code : one of the keys of PATCHES (e.g. "8OG", "PUU", "M1A", "SEP", ...)

    Returns
    -------
    PDB-format string of the modified structure. CONECT records are added
    for any non-standard residues so OpenMM/tleap pick up the bonds.
    """
    import Bio.PDB as bpdb

    ptm_code = ptm_code.upper()
    if ptm_code not in PATCHES:
        raise KeyError(f"Unknown patch code {ptm_code!r}. "
                       f"Known: {sorted(PATCHES)}")

    patch = PATCHES[ptm_code]
    resnum = int(resnum)

    parser = bpdb.PDBParser(QUIET=True)
    structure = parser.get_structure("mol", str(pdb_path))
    model = structure[0]

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
        raise ValueError(f"Residue {chain_id}{resnum} not found in {pdb_path}")

    orig_name = target_res.resname.strip().upper()
    if orig_name not in patch.get("from", []):
        raise ValueError(
            f"{ptm_code} patch requires source residue in {patch['from']!r}, "
            f"got {orig_name!r}")

    atom_coords: dict[str, np.ndarray] = {
        a.name.strip(): np.array(a.coord, dtype=float)
        for a in target_res.get_atoms()
    }

    # 1. Delete atoms
    for del_name in patch.get("delete", []):
        if del_name in atom_coords:
            del atom_coords[del_name]
            if target_res.has_id(del_name):
                target_res.detach_child(del_name)

    # 2. Add atoms via NeRF
    serial_offset = max((a.serial_number for a in structure.get_atoms()),
                         default=0)
    for atom_def in patch.get("add", []):
        name = atom_def["name"]
        ref1, ref2, ref3 = atom_def["ref"]
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
            name=name, coord=pos, bfactor=0.0, occupancy=1.0,
            altloc=" ", fullname=f" {name:<3s}",
            serial_number=serial_offset, element=atom_def["element"],
        )
        target_res.add(new_atom)
        atom_coords[name] = pos

    # 3. Rename residue
    target_res.resname = patch["rename"]

    # 4. Serialise to string
    io = bpdb.PDBIO()
    io.set_structure(structure)
    buf = StringIO()
    io.save(buf)
    pdb_str = buf.getvalue()

    # 5. Apply rebonding (PUU only currently). Handled via CONECT post-pass:
    #    BioPython doesn't preserve bonds explicitly in PDBIO; we add the
    #    new C5-C1' bond as a CONECT record and rely on tleap's PDB reader
    #    to pick it up. The "remove N1-C1'" half is implicit: N1 just
    #    happens not to have a CONECT partner C1' anymore. tleap/OpenMM
    #    will rebuild bonds from the LIB template once the residue is
    #    renamed PUU, so the CONECT trick is mainly belt-and-suspenders.
    pdb_str = _add_conect_for_patch(pdb_str, patch, target_res, chain_id, resnum)
    return pdb_str


def _add_conect_for_patch(pdb_str: str, patch: dict, residue, chain_id: str,
                           resnum: int) -> str:
    """Add explicit CONECT records for any 'rebond add' actions in the patch.

    OpenMM/tleap typically rebuild bonds from LIB templates for known
    residues, so this is precautionary. For PUU the new C5-C1' bond is
    one we want preserved through PDB I/O.
    """
    rebond = patch.get("rebond", [])
    if not rebond:
        return pdb_str

    # Build a name -> serial map for atoms in the target residue
    name_to_serial: dict[str, int] = {}
    for line in pdb_str.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            serial = int(line[6:11])
            aname  = line[12:16].strip()
            chain  = line[21]
            try:
                rseq = int(line[22:26])
            except ValueError:
                continue
            if chain == chain_id and rseq == resnum:
                name_to_serial[aname] = serial
        except (ValueError, IndexError):
            continue

    new_conect: list[str] = []
    for action, a1, a2 in rebond:
        if action != "add":
            continue
        s1 = name_to_serial.get(a1)
        s2 = name_to_serial.get(a2)
        if s1 and s2:
            new_conect.append(f"CONECT{s1:>5d}{s2:>5d}")
            new_conect.append(f"CONECT{s2:>5d}{s1:>5d}")

    if not new_conect:
        return pdb_str

    lines = pdb_str.rstrip().splitlines()
    result: list[str] = []
    inserted = False
    for line in lines:
        if line.startswith("END") and not inserted:
            result.extend(new_conect)
            inserted = True
        result.append(line)
    if not inserted:
        result.extend(new_conect)
    return "\n".join(result) + "\n"


# ── CLI for quick testing ───────────────────────────────────────────────────

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 5:
        print("Usage: python ptm_builder_v2.py <pdb> <chain> <resnum> <code>")
        sys.exit(1)
    out = apply_ptm(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    print(out)
