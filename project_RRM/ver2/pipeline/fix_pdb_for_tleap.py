# -*- coding: utf-8 -*-
"""
fix_pdb_for_tleap.py
====================
Post-process a PDBFixer-produced PDB so tleap can load it without the
two most common pitfalls:

  1. **Histidine protonation state**: PDBFixer keeps the residue name
     as plain ``HIS`` while sometimes adding HD1 alone, HE2 alone, or
     both. tleap (ff14SB) requires the residue name to indicate the
     protonation:
        HID  -- delta-protonated  (HD1 present, HE2 absent)
        HIE  -- epsilon-protonated (HE2 present, HD1 absent)  [default at pH 7]
        HIP  -- doubly protonated (both)  [+1 charge]
     We rename HIS -> HID/HIE/HIP based on which H atoms are present.

  2. **N-terminal amine H naming**: PDBFixer writes the backbone amide
     H as just ``H`` plus optional ``H2``/``H3``. tleap's NXXX templates
     (e.g. NGLY) require ``H1``/``H2``/``H3``. Rename ``H`` -> ``H1`` on
     N-terminal residues (residues whose name is ``NXXX`` style or that
     are first in their chain and have H2/H3 already).

This is invoked by prepare_system.py after PDBFixer; can also be used
standalone:

    python fix_pdb_for_tleap.py <input.pdb> <output.pdb>

Idempotent: re-running on a fixed PDB makes no changes.
"""

from __future__ import annotations
import sys
from collections import defaultdict
from pathlib import Path


HIS_NAMES = {"HIS", "HID", "HIE", "HIP"}


def fix_pdb(in_path: Path, out_path: Path) -> dict:
    lines = in_path.read_text(encoding="utf-8").splitlines(keepends=True)

    # Pass 1: collect atom names per residue (chain, resseq) for HIS-class
    his_atoms: dict[tuple, set[str]] = defaultdict(set)
    # Also: collect whether residue has H/H2/H3 patterns (for N-term detection)
    h_atoms_per_res: dict[tuple, set[str]] = defaultdict(set)
    first_res_per_chain: dict[str, tuple] = {}

    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            atom_name = line[12:16].strip()
            resname   = line[17:20].strip()
            chain     = line[21]
            resseq    = int(line[22:26])
        except (ValueError, IndexError):
            continue

        key = (chain, resseq)
        if resname in HIS_NAMES:
            his_atoms[key].add(atom_name)
        h_atoms_per_res[key].add(atom_name)
        if chain not in first_res_per_chain:
            first_res_per_chain[chain] = key

    # Decide HIS rename
    his_rename: dict[tuple, str] = {}
    for key, atoms in his_atoms.items():
        has_hd1 = "HD1" in atoms
        has_he2 = "HE2" in atoms
        if has_hd1 and has_he2:
            his_rename[key] = "HIP"
        elif has_hd1:
            his_rename[key] = "HID"
        else:                # default / has_he2
            his_rename[key] = "HIE"

    # Decide which residues to apply N-terminal H -> H1 rename:
    # - First residue of each chain
    # - AND has both H and (H2 or H3) — indicates PDBFixer-added N-term protons
    nterm_residues: set[tuple] = set()
    for chain, first_key in first_res_per_chain.items():
        atoms = h_atoms_per_res.get(first_key, set())
        if "H" in atoms and ("H2" in atoms or "H3" in atoms):
            nterm_residues.add(first_key)

    # Pass 2: rewrite lines
    n_his_renamed = 0
    n_h_renamed = 0
    new_lines = []
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            new_lines.append(line)
            continue
        try:
            atom_name = line[12:16].strip()
            resname   = line[17:20].strip()
            chain     = line[21]
            resseq    = int(line[22:26])
        except (ValueError, IndexError):
            new_lines.append(line)
            continue
        key = (chain, resseq)

        # HIS rename
        new_resname = resname
        if key in his_rename and resname in HIS_NAMES:
            new_resname = his_rename[key]
            if new_resname != resname:
                n_his_renamed += 1

        # N-term H -> H1
        new_atom_name = atom_name
        if key in nterm_residues and atom_name == "H":
            new_atom_name = "H1"
            n_h_renamed += 1

        if new_resname == resname and new_atom_name == atom_name:
            new_lines.append(line)
            continue

        # Rebuild the line keeping fixed-column layout
        # PDB cols: 1-6 record, 7-11 serial, 13-16 name, 17 altloc,
        # 18-20 resname, 22 chain, 23-26 resseq.
        # Atom-name column is 13-16 (4 chars), right-justified to col 14
        # for 1-2 char names (with leading space) — but let's just use
        # left-justified 4-char field starting at col 13.
        formatted_name = new_atom_name.ljust(4)[:4]
        # PDB convention: if first char of element is a digit, name starts at col 13
        # else name starts at col 14 (right-shifted). For "H1" we want " H1 "
        # not "H1  ". Use a simple rule: 1-3 char names get a leading space.
        if len(new_atom_name) <= 3 and not new_atom_name[0].isdigit():
            formatted_name = (" " + new_atom_name).ljust(4)[:4]
        new_line = (line[:12] + formatted_name + line[16:17]
                    + new_resname.ljust(3) + line[20:])
        new_lines.append(new_line)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("".join(new_lines), encoding="utf-8")

    return {
        "his_renames": dict(his_rename),
        "n_his_renamed": n_his_renamed,
        "nterm_residues": sorted(nterm_residues),
        "n_h_renamed": n_h_renamed,
    }


def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: python fix_pdb_for_tleap.py <input.pdb> <output.pdb>")
    in_path  = Path(sys.argv[1]).resolve()
    out_path = Path(sys.argv[2]).resolve()
    if not in_path.exists():
        sys.exit(f"ERROR: {in_path} not found")

    info = fix_pdb(in_path, out_path)
    print(f"Wrote {out_path}")
    print(f"  HIS -> protonation-aware rename: {info['n_his_renamed']} atoms changed")
    if info['his_renames']:
        for (chain, resseq), tgt in sorted(info['his_renames'].items()):
            print(f"    {chain}{resseq}: -> {tgt}")
    print(f"  N-terminal H -> H1: {info['n_h_renamed']} atoms changed")
    if info['nterm_residues']:
        print(f"    N-term residues: {info['nterm_residues']}")


if __name__ == "__main__":
    main()
