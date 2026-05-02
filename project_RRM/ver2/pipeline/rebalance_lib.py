# -*- coding: utf-8 -*-
"""
rebalance_lib.py
================
Post-process AMBER .lib files produced by our patched modxna.sh to:

  1. Adjust the residue net charge to the expected integer/cap target
     by distributing the residual to a designated **anchor atom**. The
     residual exists because we removed the malformed
     `strip MASK charge VAL` modxna.sh syntax to get past the cpptraj
     bug; this restores what those corrections were meant to do.

  2. Fix the duplicate atom-name issue in 8OG: modXNA's 8OG.mol2 has
     two atoms named "O6" (the 8-oxo carbonyl O and the standard
     guanine 6-O). Rename the FIRST one (bonded to C8) to "O8".

Inputs:  $PROJECT/force_fields/lib_amber/<CODE>{I,3,5}.lib
Outputs: same files, in-place modified. Backups saved as *.lib.preBalance.

Targets (per modxna's CHARGE constants):
    internal: net charge = -1.000000
    3cap   : net charge = -0.679652   (CHARGE_3CAP)
    5cap   : net charge = -0.320348   (CHARGE_5CAP)

Anchor atom (where to deposit the correction):
    8OG / M1A / M6A : N9   (purines: glycosidic atom)
    PUU / PSU       : C5   (canonical psi: C-glycosidic at C5)

Usage on the cluster:
    bash project_RRM/ver2/pipeline/server_setup.sh   # ensures allatom_v2 active
    python project_RRM/ver2/pipeline/rebalance_lib.py
"""

from __future__ import annotations
import os
import re
import shutil
import sys
from pathlib import Path

# ── Configuration ────────────────────────────────────────────────────────────
LIB_DIR = Path(os.environ.get(
    "MODXNA_LIB_DIR",
    "/data/biophys/carolinge/clawork/37_OXR/force_fields/lib_amber"))

CHARGE_INTERNAL = -1.000000
CHARGE_3CAP     = -0.679652
CHARGE_5CAP     = -0.320348

# code (3-letter base) -> anchor atom name
ANCHOR = {
    "8OG": "N9",
    "M1A": "N9",
    "M6A": "N9",
    "PUU": "C5",
    "PSU": "C5",
}

# Variant suffix -> target charge
TARGET_BY_VARIANT = {
    "I": CHARGE_INTERNAL,   # internal (e.g. 8OGI.lib)
    "3": CHARGE_3CAP,       # 3'-cap
    "5": CHARGE_5CAP,       # 5'-cap
}

# 8OG-specific atom rename: which "O6" needs to become "O8"?
# In modXNA's 8OG.mol2 ATOM section, atom 3 is the 8-oxo O (bonded to C8 / type CK).
# In the assembled lib, this atom appears as the FIRST "O6" entry (it's listed
# right after the C8 it bonds to, before the standard guanine O6 atom).
#
# We use a heuristic: the first occurrence of name="O6" with type=O appearing
# *before* atom name "C6" in the entry is the 8-oxo. Rename it to "O8".


# ── lib parsing ─────────────────────────────────────────────────────────────

ATOM_LINE_RE = re.compile(
    r'^( +)"([^"]+)"\s+"([^"]+)"\s+(\S+\s+\S+\s+\S+\s+\S+\s+\S+)\s+(-?\d+\.\d+)\s*$'
)


def find_residue_blocks(lines: list[str]) -> dict:
    """Return {residue_name: (atoms_start, atoms_end_exclusive)} index ranges
    where each .lib entry's atoms section is."""
    blocks: dict = {}
    cur_res = None
    cur_start = None
    for i, line in enumerate(lines):
        m = re.match(r'^!entry\.([A-Za-z0-9]+)\.unit\.atoms ', line)
        if m:
            cur_res = m.group(1)
            cur_start = i + 1
            continue
        if cur_res and line.startswith("!entry.") and cur_start:
            # Hit the next section (atomspertinfo etc.)
            blocks.setdefault(cur_res, (cur_start, i))
            cur_res = None
            cur_start = None
    return blocks


def parse_atom_line(line: str) -> dict | None:
    m = ATOM_LINE_RE.match(line)
    if not m:
        return None
    return {
        "indent":   m.group(1),
        "name":     m.group(2),
        "type":     m.group(3),
        "middle":   m.group(4),
        "charge":   float(m.group(5)),
    }


def reformat_atom_line(parts: dict) -> str:
    return (f'{parts["indent"]}"{parts["name"]}" "{parts["type"]}" '
            f'{parts["middle"]} {parts["charge"]:.6f}\n')


# ── Per-lib processing ─────────────────────────────────────────────────────

def rebalance_one_lib(lib_path: Path) -> dict:
    """Read, rebalance, and write back a single .lib file.

    Returns a dict with diagnostics.
    """
    name = lib_path.stem            # e.g. "8OGI"
    if len(name) != 4:
        print(f"  [skip] {lib_path.name}: residue code not 4 chars")
        return {}
    base, variant = name[:3], name[3]
    if variant not in TARGET_BY_VARIANT:
        print(f"  [skip] {lib_path.name}: variant '{variant}' not in (I,3,5)")
        return {}
    if base not in ANCHOR:
        print(f"  [skip] {lib_path.name}: no anchor for base {base!r}")
        return {}
    target = TARGET_BY_VARIANT[variant]
    anchor = ANCHOR[base]

    text = lib_path.read_text(encoding="utf-8").splitlines(keepends=True)
    blocks = find_residue_blocks(text)
    if name not in blocks:
        print(f"  [skip] {lib_path.name}: no entry for residue {name}")
        return {}

    a_start, a_end = blocks[name]

    # Parse atoms in the block
    parsed = []
    seen_c6 = False
    o8_renamed = False
    for i in range(a_start, a_end):
        p = parse_atom_line(text[i])
        if p is None:
            continue
        # 8OG: rename the FIRST O6 (which appears before C6) to O8
        if base == "8OG" and not o8_renamed and p["name"] == "O6" and not seen_c6:
            p["name"] = "O8"
            o8_renamed = True
            text[i] = reformat_atom_line(p)
        if p["name"] == "C6":
            seen_c6 = True
        parsed.append((i, p))

    sum_q = sum(p["charge"] for _, p in parsed)
    residual = sum_q - target

    # Find the anchor atom in this block
    anchor_idx = None
    for line_i, p in parsed:
        if p["name"] == anchor:
            anchor_idx = line_i
            break
    if anchor_idx is None:
        print(f"  [warn] {lib_path.name}: anchor atom {anchor!r} not found "
              f"in residue {name}; skipping rebalance")
        return {"residue": name, "atoms": len(parsed), "sum": sum_q,
                "target": target, "residual": residual, "renamed_O8": o8_renamed,
                "rebalanced": False}

    # Adjust the anchor atom's charge by -residual
    p_anchor = parse_atom_line(text[anchor_idx])
    old_anchor_q = p_anchor["charge"]
    p_anchor["charge"] = old_anchor_q - residual
    text[anchor_idx] = reformat_atom_line(p_anchor)

    # Verify
    new_sum = sum(parse_atom_line(text[line_i])["charge"]
                  for line_i, _ in parsed)

    # Backup + write
    backup = lib_path.with_suffix(lib_path.suffix + ".preBalance")
    if not backup.exists():
        shutil.copy2(lib_path, backup)
    lib_path.write_text("".join(text), encoding="utf-8")

    return {
        "residue": name,
        "atoms": len(parsed),
        "anchor": anchor,
        "anchor_old_q": old_anchor_q,
        "anchor_new_q": p_anchor["charge"],
        "sum_old": sum_q,
        "sum_new": new_sum,
        "target": target,
        "residual_subtracted": residual,
        "renamed_O8": o8_renamed,
        "rebalanced": True,
    }


def main():
    if not LIB_DIR.is_dir():
        sys.exit(f"ERROR: lib dir not found: {LIB_DIR}")

    libs = sorted(LIB_DIR.glob("*.lib"))
    if not libs:
        sys.exit(f"ERROR: no .lib in {LIB_DIR}")

    print(f"Processing {len(libs)} libs in {LIB_DIR}\n")
    print(f"{'lib':<10} {'atoms':>5} {'anchor':>6} "
          f"{'old anchor q':>12} {'new anchor q':>12} "
          f"{'sum_old':>10} {'target':>10} {'sum_new':>10} {'O8 fix':>7}")
    print("-" * 100)

    all_results = []
    for lib in libs:
        r = rebalance_one_lib(lib)
        if not r or not r.get("rebalanced"):
            continue
        all_results.append(r)
        print(f"{r['residue']:<10} {r['atoms']:>5} {r['anchor']:>6} "
              f"{r['anchor_old_q']:>+12.6f} {r['anchor_new_q']:>+12.6f} "
              f"{r['sum_old']:>+10.6f} {r['target']:>+10.6f} {r['sum_new']:>+10.6f} "
              f"{'YES' if r['renamed_O8'] else '-':>7}")

    # Sanity flags
    print("")
    print("Sanity checks:")
    for r in all_results:
        if abs(r["sum_new"] - r["target"]) > 1e-4:
            print(f"  ⚠️  {r['residue']}: sum {r['sum_new']:+.6f} "
                  f"!= target {r['target']:+.6f}")
        if abs(r["residual_subtracted"]) > 0.5:
            print(f"  ⚠️  {r['residue']}: large adjustment "
                  f"{r['residual_subtracted']:+.4f} e on anchor {r['anchor']} "
                  f"({r['anchor_old_q']:+.4f} -> {r['anchor_new_q']:+.4f})")
    print("Done.")


if __name__ == "__main__":
    main()
