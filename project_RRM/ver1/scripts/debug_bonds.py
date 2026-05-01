# -*- coding: utf-8 -*-
"""Debug: compare WT vs 8OG PDB bond topology for chain B."""
from openmm.app import PDBFile
import numpy as np

def show_chain_b(pdb_path, label):
    pdb = PDBFile(str(pdb_path))
    top = pdb.topology
    from openmm.unit import nanometer
    pos = pdb.positions.value_in_unit(nanometer)

    print(f"\n=== {label}: {pdb_path} ===")
    for c in top.chains():
        if c.index != 1:
            continue
        for r in c.residues():
            cross = []
            for b in top.bonds():
                a1, a2 = b[0], b[1]
                if (a1.residue == r) != (a2.residue == r):
                    if a1.residue == r:
                        cross.append(f"{a1.name}->{a2.residue.name}{a2.residue.id}.{a2.name}")
                    else:
                        cross.append(f"{a1.residue.name}{a1.residue.id}.{a1.name}->{a2.name}")
            print(f"  {r.name:>4s} {r.id:>3s} idx={r.index:3d}  cross_bonds={cross}")

    # Check O3'-P distances between consecutive residues in chain B
    print(f"\n  Backbone O3'->P distances:")
    residues = []
    for c in top.chains():
        if c.index != 1:
            continue
        residues = list(c.residues())

    for i in range(len(residues) - 1):
        r1, r2 = residues[i], residues[i+1]
        o3p = None
        p_next = None
        for a in r1.atoms():
            if a.name == "O3'":
                o3p = a
        for a in r2.atoms():
            if a.name == "P":
                p_next = a
        if o3p and p_next:
            p1 = np.array(pos[o3p.index]) * 10  # nm -> A
            p2 = np.array(pos[p_next.index]) * 10
            dist = np.linalg.norm(p1 - p2)
            print(f"    {r1.name}{r1.id}.O3' -> {r2.name}{r2.id}.P  = {dist:.2f} A")

show_chain_b("D:/all_atom/output/4BS2_WT_MD/4BS2_prepared.pdb", "WT")
show_chain_b("D:/all_atom/output/4BS2_8OG_G3_MD/4BS2_8OG_G3.pdb", "8OG")
