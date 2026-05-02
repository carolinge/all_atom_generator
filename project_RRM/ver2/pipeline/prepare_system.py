# -*- coding: utf-8 -*-
"""
prepare_system.py
=================
End-to-end structure preparation for ver_2:

    CIF/PDB  -->  PDBFixer (clean + add H)  -->  apply 1+ patches
              -->  PDB ready for tleap (build_amber_system.sh)

Designed to live next to ptm_builder_v2.py. Output goes to
project_RRM/ver2/structures/modified/<NAME>.pdb.

Usage:
    python prepare_system.py <input.cif|input.pdb> <output_name> \\
        --patch <chain>:<resnum>:<code> [--patch ...]

Examples:
    # WT 4BS2 (no patches), just clean
    python prepare_system.py structures/original/4BS2.cif 4BS2_WT

    # 4BS2 with G3 -> 8OG on chain B
    python prepare_system.py structures/original/4BS2.cif 4BS2_8OG_G3 \\
        --patch B:3:8OG

    # 4BS2 with two PUU substitutions
    python prepare_system.py structures/original/4BS2.cif 4BS2_PUU_U2_U4 \\
        --patch B:2:PUU --patch B:4:PUU

    # adenine-context tetraloop with M1A in middle (free-RNA, no protein)
    python prepare_system.py structures/original/GAAA.pdb GAAA_M1A \\
        --patch A:2:M1A

The output PDB is what build_amber_system.sh consumes on the cluster.
"""

from __future__ import annotations
import argparse
import sys
import tempfile
from pathlib import Path

from ptm_builder_v2 import apply_ptm, PATCHES   # local import
from fix_pdb_for_tleap import fix_pdb            # local import


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT_DIR = REPO_ROOT / "project_RRM" / "ver2" / "structures" / "modified"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("input", help="input structure file (.cif or .pdb)")
    p.add_argument("name", help="output system name (no extension)")
    p.add_argument("--patch", action="append", default=[],
                   help="<chain>:<resnum>:<code> e.g. B:3:8OG. Repeat per patch.")
    p.add_argument("--ph", type=float, default=7.0,
                   help="pH for hydrogen placement (default 7.0)")
    p.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR,
                   help=f"output directory (default {DEFAULT_OUT_DIR})")
    p.add_argument("--keep-water", action="store_true",
                   help="keep crystallographic waters (default: remove)")
    return p.parse_args()


def parse_patch_spec(spec: str) -> tuple[str, int, str]:
    parts = spec.split(":")
    if len(parts) != 3:
        raise ValueError(f"Bad --patch spec {spec!r}; expected chain:resnum:code")
    chain, resnum, code = parts
    if code.upper() not in PATCHES:
        raise ValueError(f"Unknown patch code {code!r}. "
                         f"Known: {sorted(PATCHES)}")
    return chain, int(resnum), code.upper()


def run_pdbfixer(input_path: Path, ph: float, keep_water: bool,
                 out_path: Path) -> Path:
    """Clean structure with PDBFixer; write to out_path."""
    from pdbfixer import PDBFixer
    import openmm.app as mm

    print(f"[pdbfixer] loading {input_path}")
    fixer = PDBFixer(filename=str(input_path))

    print("[pdbfixer] findMissingResidues")
    fixer.findMissingResidues()

    print("[pdbfixer] removeHeterogens "
          f"(keep_water={keep_water})")
    fixer.removeHeterogens(keepWater=keep_water)

    print("[pdbfixer] findMissingAtoms / addMissingAtoms")
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    print(f"[pdbfixer] addMissingHydrogens(pH={ph})")
    fixer.addMissingHydrogens(pH=ph)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        mm.PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
    print(f"[pdbfixer] -> {out_path}  ({fixer.topology.getNumAtoms()} atoms)")
    return out_path


def apply_patches(start_pdb: Path, patches: list[tuple[str, int, str]],
                  out_path: Path) -> Path:
    """Apply patches sequentially (each one re-reads the previous output)."""
    if not patches:
        # No patches — just copy.
        out_path.write_text(start_pdb.read_text(encoding="utf-8"),
                            encoding="utf-8")
        print(f"[patches] no patches; copied -> {out_path}")
        return out_path

    current = start_pdb
    tmp_files: list[Path] = []
    for i, (chain, resnum, code) in enumerate(patches, start=1):
        print(f"[patches] {i}/{len(patches)}: {code} on {chain}{resnum}")
        new_str = apply_ptm(str(current), chain, resnum, code)
        if i == len(patches):
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_text(new_str, encoding="utf-8")
            current = out_path
        else:
            tmp = Path(tempfile.NamedTemporaryFile(
                mode="w", suffix=".pdb", delete=False).name)
            tmp.write_text(new_str, encoding="utf-8")
            tmp_files.append(tmp)
            current = tmp

    for t in tmp_files:
        try: t.unlink()
        except OSError: pass

    print(f"[patches] all applied -> {out_path}")
    return out_path


def main():
    args = parse_args()

    input_path = Path(args.input).resolve()
    if not input_path.exists():
        # Try resolving relative to repo root
        rel = REPO_ROOT / args.input
        if rel.exists():
            input_path = rel.resolve()
        else:
            sys.exit(f"ERROR: input not found: {args.input}")

    patches = [parse_patch_spec(s) for s in args.patch]

    out_dir = Path(args.out_dir).resolve()
    final_path = out_dir / f"{args.name}.pdb"
    cleaned_path = out_dir / f"{args.name}_cleaned.pdb"

    # Stage 1: PDBFixer
    run_pdbfixer(input_path, args.ph, args.keep_water, cleaned_path)

    # Stage 1b: tleap-friendly fixes (HIS protonation, N-term H1 naming)
    fixed_path = out_dir / f"{args.name}_fixed.pdb"
    fix_info = fix_pdb(cleaned_path, fixed_path)
    print(f"[fix_pdb] HIS renames: {fix_info['n_his_renamed']}; "
          f"N-term H->H1: {fix_info['n_h_renamed']}")
    cleaned_path.unlink(missing_ok=True)

    # Stage 2: patches
    apply_patches(fixed_path, patches, final_path)

    # Cleanup intermediate
    if fixed_path != final_path:
        fixed_path.unlink(missing_ok=True)

    print("")
    print(f"=== prepare_system.py DONE ===")
    print(f"Output: {final_path}")
    if patches:
        codes = "+".join(f"{c}{r}->{p}" for c, r, p in patches)
        print(f"Patches applied: {codes}")
    print("")
    print("Next on the cluster:")
    print(f"  bash project_RRM/ver2/pipeline/build_amber_system.sh \\")
    print(f"      {args.name} \\")
    print(f"      project_RRM/ver2/structures/modified/{args.name}.pdb")


if __name__ == "__main__":
    main()
