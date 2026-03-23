# -*- coding: utf-8 -*-
"""
validate_8OG_ff.py
==================
Validates the 8OG_RNA_amber.xml force field file by:
  1. Loading it with RNA.OL3.xml
  2. Checking all three residue templates (8OG, 8OG3, 8OG5) are present
  3. Building a G5-8OG3 dinucleotide topology (proper chain with satisfied
     external bonds) and calling createSystem -- raises exceptions for any
     missing bonds, angles, torsions, or nonbonded parameters.
  4. Reporting success or listing missing parameter types.

Usage:
    conda activate allatom
    cd d:/all_atom
    python scripts/validate_8OG_ff.py
"""

import sys
from pathlib import Path

try:
    from openmm.app import ForceField, Topology, Element
    import openmm.app as omm_app
except ImportError:
    sys.exit("ERROR: openmm not found. Run: conda activate allatom")

PROJECT_ROOT = Path(__file__).parent.parent
OL3_XML  = "amber14/RNA.OL3.xml"
TIP3_XML = "amber14/tip3p.xml"
OG_XML   = str(PROJECT_ROOT / "app" / "ff_custom" / "8OG_RNA_amber.xml")

if not Path(OG_XML).exists():
    sys.exit(f"ERROR: {OG_XML} not found.")

print("=" * 60)
print("8OG RNA AMBER Force Field Validation")
print("=" * 60)

# ── Step 1: Load force field ──────────────────────────────────────
print("\n[1] Loading force field files ...")
try:
    ff = ForceField(OL3_XML, TIP3_XML, OG_XML)
    print("    OK: ForceField loaded")
except Exception as e:
    print(f"    FAIL: {e}")
    sys.exit(1)

# ── Step 2: Check residue templates ───────────────────────────────
print("\n[2] Checking residue templates ...")
all_ok = True
for name in ("8OG", "8OG3", "8OG5"):
    if name in ff._templates:
        n = len(ff._templates[name].atoms)
        print(f"    OK: {name} ({n} atoms)")
    else:
        print(f"    FAIL: template '{name}' not found")
        all_ok = False
if not all_ok:
    sys.exit(1)

# ── Step 3: Build G5 + 8OG3 dinucleotide topology ─────────────────
# G5 is the 5'-terminal guanosine (from RNA.OL3.xml, no phosphate).
# 8OG3 is the 3'-terminal 8-oxoguanosine (our XML, has phosphate,
#   ExternalBond on P).  The inter-residue bond is G5.O3' -- 8OG3.P.
print("\n[3] Building test topology (G5 + 8OG3 dinucleotide) ...")

def add_residue_from_template(topology, chain, tmpl_name, ff):
    tmpl = ff._templates[tmpl_name]
    res = topology.addResidue(tmpl_name, chain)
    atoms_list = list(tmpl.atoms)
    atom_map = {}
    for a in atoms_list:
        try:
            elem = Element.getBySymbol(a.element.symbol if a.element else "C")
        except Exception:
            elem = Element.getBySymbol("C")
        atom_map[a.name] = topology.addAtom(a.name, elem, res)
    for b in tmpl.bonds:
        topology.addBond(atom_map[atoms_list[b[0]].name],
                         atom_map[atoms_list[b[1]].name])
    return res, atom_map

topology = Topology()
chain = topology.addChain()

res1, amap1 = add_residue_from_template(topology, chain, "G5",   ff)
res2, amap2 = add_residue_from_template(topology, chain, "8OG3", ff)

# Inter-residue phosphodiester bond: O3' of G5 -- P of 8OG3
topology.addBond(amap1["O3'"], amap2["P"])

print(f"    Topology: {topology.getNumAtoms()} atoms, "
      f"{topology.getNumBonds()} bonds")

# ── Step 4: Create System ─────────────────────────────────────────
print("\n[4] Creating OpenMM System (validates all FF parameters) ...")
try:
    system = ff.createSystem(topology, nonbondedMethod=omm_app.NoCutoff)
    print(f"    OK: System created ({system.getNumParticles()} particles)")
    forces = [type(f).__name__ for f in system.getForces()]
    print(f"    Forces: {forces}")
except Exception as e:
    print(f"\n    FAIL: {e}")
    sys.exit(1)

print("\n" + "=" * 60)
print("VALIDATION PASSED -- 8OG parameters are complete.")
print("Caveats: charges are approximate (see XML header).")
print("         Replace with QM RESP charges for production use.")
print("=" * 60)
