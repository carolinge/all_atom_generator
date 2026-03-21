# -*- coding: utf-8 -*-
"""
build_8oxoG_rna.py
==================
Derives CHARMM36 parameters for RNA 8-oxoguanosine (8OG) from:
  - toppar_all36_na_modifications.str  (DNA 8OG = 8-oxo-dG)
  - top_all36_na.rtf                   (standard RNA GUA residue, for 2'-OH)

Usage:
    conda activate allatom
    python scripts/build_8oxoG_rna.py \
        --modifications toppar/stream/na/toppar_all36_na_modifications.str \
        --na_rtf        toppar/top_all36_na.rtf \
        --output        app/ff_custom/toppar_8OG_RNA.str

Required files (download once from Mackerell lab):
    https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/toppar_c36_jul22.tgz

Scientific basis:
    8-oxoguanosine (RNA) differs from 8-oxo-dG (DNA) only in the sugar:
      - C2': type CN8 (dG, -CH2-) -> CN7 (rG, -CHOH-)
      - Remove H2'' atom (HN8 type)
      - Add O2' atom  (ON5 type, q = -0.61) matching standard GUA
      - Add HO2' atom (HN5 type, q = +0.43) matching standard GUA
      - Adjust C2' charge: -0.18 (dG) -> +0.14 (rG)
      - Adjust H2' charge: +0.09 (dG) -> +0.09 (rG, unchanged)
    Total residue charge remains 0 (nucleoside, before phosphate patch).

    Reference atom types from CHARMM36 NA force field:
      ON5 / HN5 -- 2'-hydroxyl (from top_all36_na.rtf GUA)
      CN7       -- C2' in ribose (CHOH carbon)
      HN7       -- H2' in ribose
"""

import argparse
import re
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Sugar atom changes: DNA 8-oxo-dG  ->  RNA 8-oxoguanosine
# ---------------------------------------------------------------------------
# Atoms to REMOVE from DNA 8OG topology
REMOVE_ATOMS = {"H2''"}           # second H on C2' (deoxyribose only)

# Atoms to ADD (name, type, charge) -- copied from CHARMM36 GUA (RNA)
ADD_ATOMS = [
    ("O2'",  "ON5",  -0.61),
    ("HO2'", "HN5",  +0.43),
]

# Atom-by-atom charge/type overrides for the sugar atoms that change
PATCH_ATOMS = {
    "C2'": ("CN7",  +0.14),   # dG: CN8 -0.18  ->  rG: CN7 +0.14
    # H2' charge stays the same (+0.09), type HN7 stays the same
}

# Bonds to ADD (involving new atoms)
ADD_BONDS = [
    ("C2'", "O2'"),
    ("O2'", "HO2'"),
]

# Angles to ADD (involving new atoms -- CHARMM will look up parameters by atom type)
ADD_ANGLES = [
    ("C1'", "C2'", "O2'"),
    ("C3'", "C2'", "O2'"),
    ("H2'", "C2'", "O2'"),
    ("C2'", "O2'", "HO2'"),
]


# ---------------------------------------------------------------------------
# Parser helpers
# ---------------------------------------------------------------------------

def extract_residue_block(text: str, resname: str) -> str | None:
    """Return the full RESI block for *resname* from a CHARMM RTF/stream text."""
    # Match from 'RESI 8OG' up to next RESI/PRES/END line
    pattern = rf"(?m)^(RESI\s+{re.escape(resname)}\b.*?)(?=^RESI\s|^PRES\s|^END\b)"
    m = re.search(pattern, text, re.DOTALL | re.IGNORECASE)
    return m.group(1).strip() if m else None


def parse_atom_lines(block: str) -> list[dict]:
    """Parse ATOM lines from a RESI block."""
    atoms = []
    for line in block.splitlines():
        m = re.match(r"\s*ATOM\s+(\S+)\s+(\S+)\s+([+-]?\d+\.\d+)", line, re.IGNORECASE)
        if m:
            atoms.append({"name": m.group(1), "type": m.group(2), "charge": float(m.group(3))})
    return atoms


def parse_bond_lines(block: str) -> list[tuple[str, str]]:
    """Parse BOND lines, returning list of (a1, a2) pairs."""
    bonds = []
    for line in block.splitlines():
        if re.match(r"\s*BOND\s", line, re.IGNORECASE):
            parts = line.split()[1:]
            for i in range(0, len(parts) - 1, 2):
                bonds.append((parts[i], parts[i+1]))
    return bonds


# ---------------------------------------------------------------------------
# Main build function
# ---------------------------------------------------------------------------

def build_rna_8og(mod_str_text: str, na_rtf_text: str) -> str:
    """Return the text of a CHARMM stream file defining RNA 8OG."""

    # 1. Extract DNA 8OG block
    dna_block = extract_residue_block(mod_str_text, "8OG")
    if dna_block is None:
        sys.exit("ERROR: Could not find 'RESI 8OG' in modifications stream file.")

    dna_atoms = parse_atom_lines(dna_block)
    dna_bonds = parse_bond_lines(dna_block)

    # 2. Build RNA atom list
    rna_atoms = []
    for atom in dna_atoms:
        name = atom["name"]
        if name in REMOVE_ATOMS:
            continue
        atype, charge = PATCH_ATOMS.get(name, (atom["type"], atom["charge"]))
        rna_atoms.append({"name": name, "type": atype, "charge": charge})

    # Insert O2' and HO2' after C2'
    insert_idx = next((i for i, a in enumerate(rna_atoms) if a["name"] == "C2'"), len(rna_atoms)) + 1
    for name, atype, charge in ADD_ATOMS:
        rna_atoms.insert(insert_idx, {"name": name, "type": atype, "charge": charge})
        insert_idx += 1

    # Verify total charge is 0
    total_q = sum(a["charge"] for a in rna_atoms)
    if abs(total_q) > 0.01:
        print(f"WARNING: RNA 8OG total charge = {total_q:.4f} (expected 0). "
              "Manual charge adjustment may be needed.", file=sys.stderr)

    # 3. Build RNA bond list (remove bonds involving H2'', add new ones)
    rna_bonds = [(a, b) for a, b in dna_bonds if "H2''" not in (a, b)]
    rna_bonds.extend(ADD_BONDS)

    # 4. Format output stream file
    lines = [
        "* CHARMM36 parameters for RNA 8-oxoguanosine (8OG)",
        "* Derived from toppar_all36_na_modifications.str (DNA 8OG)",
        "* by build_8oxoG_rna.py -- 2'-OH added from GUA (top_all36_na.rtf)",
        "* Validate with short vacuum MD before production use.",
        "*",
        "",
        "read rtf card append",
        "* RNA 8-oxoguanosine topology",
        "*",
        "36  1",
        "",
        "RESI 8OG      0.00  ! RNA 8-oxoguanosine; same base as DNA 8OG, ribose sugar",
        "GROUP",
    ]

    for atom in rna_atoms:
        lines.append(f"ATOM {atom['name']:<6s} {atom['type']:<6s} {atom['charge']:>8.4f}")

    lines.append("")
    # Write bonds in pairs per line (CHARMM convention)
    bond_pairs = rna_bonds
    bond_lines_out = []
    for i in range(0, len(bond_pairs), 4):
        chunk = bond_pairs[i:i+4]
        bond_str = "  ".join(f"{a} {b}" for a, b in chunk)
        bond_lines_out.append(f"BOND {bond_str}")
    lines.extend(bond_lines_out)

    lines.append("")
    for angle in ADD_ANGLES:
        lines.append(f"ANGLE {angle[0]} {angle[1]} {angle[2]}")

    lines.extend([
        "",
        "! Dihedral and improper parameters are inherited from DNA 8OG via atom types.",
        "! No new DIHE/IMPR entries needed unless validation reveals discrepancies.",
        "",
        "END",
        "",
        "read para card flex append",
        "* No new parameters needed -- all atom types (ON5, HN5, CN7, NN2U, CN1, ON1)",
        "* are already defined in par_all36_na.prm",
        "*",
        "",
        "END",
        "",
    ])

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Build RNA 8-oxoguanosine CHARMM36 stream file")
    parser.add_argument("--modifications", required=True,
                        help="Path to toppar_all36_na_modifications.str")
    parser.add_argument("--na_rtf",        required=True,
                        help="Path to top_all36_na.rtf (for reference; 2'-OH types are hardcoded)")
    parser.add_argument("--output",        required=True,
                        help="Output .str file path")
    args = parser.parse_args()

    mod_text = Path(args.modifications).read_text(encoding="utf-8", errors="replace")
    rtf_text = Path(args.na_rtf).read_text(encoding="utf-8", errors="replace")

    result = build_rna_8og(mod_text, rtf_text)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(result, encoding="utf-8")
    print(f"Written: {out}")
    print("Next step: load alongside charmm36.xml in server.py FF_RNA entry.")


if __name__ == "__main__":
    main()
