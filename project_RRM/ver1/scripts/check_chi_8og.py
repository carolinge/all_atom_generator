"""
Check χ glycosidic torsion of G3 (chain B) in WT and 8OG starting structures.

χ = O4'-C1'-N9-C4 (for purines)
  anti:  χ ≈ -180° to -90°  (canonical, base over sugar)
  syn:   χ ≈ -90° to +90°   (base flipped over ribose; 8OG prefers this in
                              DNA and free RNA due to O8-H1' clash)

If both WT and 8OG start in anti AND the 8OG MD never reaches syn within
100 ns, the cause is initial-state trapping, not the force field.
"""
import math
import numpy as np
import sys
from pathlib import Path

try:
    import Bio.PDB as bpdb
except ImportError:
    sys.stderr.write("BioPython not available; using minimal parser\n")
    bpdb = None


def parse_pdb_min(path: Path):
    """Minimal PDB parser — returns dict[(chain, resseq)] -> dict[atom_name -> xyz]."""
    res = {}
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        aname  = line[12:16].strip()
        chain  = line[21]
        try:
            resseq = int(line[22:26])
        except ValueError:
            continue
        x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
        key = (chain, resseq)
        res.setdefault(key, {"resname": line[17:20].strip(), "atoms": {}})
        res[key]["atoms"][aname] = np.array([x, y, z])
    return res


def dihedral(p0, p1, p2, p3):
    b0 = p1 - p0
    b1 = p2 - p1
    b2 = p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    x = np.dot(v, w)
    y = np.dot(np.cross(b1n, v), w)
    return math.degrees(math.atan2(y, x))


def classify_chi(chi):
    """Return rotamer label for χ (in -180..180)."""
    c = chi
    if c > 180:  c -= 360
    if c < -180: c += 360
    if -180 <= c <= -90:   return "anti"
    if -90  <  c <  -45:   return "high-anti"
    if -45  <= c <= 90:    return "SYN"
    if 90   <  c <= 180:   return "anti"
    return "?"


def chi_for_residue(atoms):
    needed = ["O4'", "C1'", "N9", "C4"]
    if not all(a in atoms for a in needed):
        return None
    return dihedral(atoms["O4'"], atoms["C1'"], atoms["N9"], atoms["C4"])


def report(label, pdb_path, chain, resseq):
    if not pdb_path.exists():
        print(f"  [{label}] file not found: {pdb_path}")
        return
    res = parse_pdb_min(pdb_path)
    key = (chain, resseq)
    if key not in res:
        print(f"  [{label}] residue {chain}{resseq} not found")
        # Show what's in chain B
        avail = sorted(k for k in res if k[0] == chain)
        print(f"           chain {chain} has residues: {avail[:20]}")
        return
    rname = res[key]["resname"]
    chi = chi_for_residue(res[key]["atoms"])
    if chi is None:
        print(f"  [{label}] {rname} {chain}{resseq}: missing atoms for χ "
              f"(have: {sorted(res[key]['atoms'])})")
        return
    rot = classify_chi(chi)
    flag = "  <-- SYN!" if rot == "SYN" else ""
    print(f"  [{label}] {rname:>4s} {chain}{resseq}: χ = {chi:+7.1f}°  →  {rot}{flag}")
    # also show O8-H1' distance for 8OG (the steric driver of syn preference)
    atoms = res[key]["atoms"]
    if "O8" in atoms and "H1'" in atoms:
        d = np.linalg.norm(atoms["O8"] - atoms["H1'"])
        warn = " (CLASH!)" if d < 2.5 else ""
        print(f"          O8 - H1' distance = {d:.2f} A{warn}")
    if "O8" in atoms and "C1'" in atoms:
        d = np.linalg.norm(atoms["O8"] - atoms["C1'"])
        print(f"          O8 - C1' distance = {d:.2f} A (anti=ok if >3.0)")


def main():
    base = Path("d:/all_atom")
    print("=" * 70)
    print("χ glycosidic torsion check for 4BS2 G3/8OG (chain B residue 3)")
    print("=" * 70)
    print("\nχ rotamer convention:")
    print("  anti:  -180° ≤ χ ≤ -90°  (or +90° to +180°)  — canonical bound state")
    print("  syn:    -45° ≤ χ ≤  +90°                     — 8OG's preferred state")
    print()

    # WT prepared
    report("WT prepared",
           base / "output/4BS2_WT_MD/4BS2_prepared.pdb",
           "B", 3)

    # 8OG starting (NeRF-built)
    report("8OG starting",
           base / "output/4BS2_8OG_G3_MD/4BS2_8OG_G3.pdb",
           "B", 3)

    # 8OG in each replica
    for r in (1, 2, 3):
        report(f"8OG_r{r} input",
               base / f"output/4BS2_8OG_r{r}/4BS2_8OG_G3.pdb",
               "B", 3)

    print()
    print("Diagnostic:")
    print("  If all reports show anti AND your MD trajectory also stays anti")
    print("  → it's an INITIAL STATE TRAPPING problem, not the force field.")
    print("  → 8OG can't flip in 100 ns when locked in the protein binding pocket.")
    print("  Fix: build a syn starting structure, or run BOTH anti and syn replicas")
    print("       and compare populations / binding stability.")


if __name__ == "__main__":
    main()
