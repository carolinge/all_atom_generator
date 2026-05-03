# -*- coding: utf-8 -*-
"""
audit_charges.py
================
Diagnose the +0.32 e charge residual that rebalance_lib.py is dumping
on N9 of internal 8OG. The expert reviewer flagged this as
"the worst possible atom for chi torsion fidelity" and noted that
typical RESP residuals are <0.05 e, so 0.32 e is ~6x too large.

This script audits the modxna fragment-stitching arithmetic to find
where the residual comes from.

For each modification (8OG, PUU, M1A, also reference standard G/U/A
from modrna08), we compute:

  1. Raw base mol2 total charge (24 atoms incl. HEAD01STRIP context)
  2. After-strip charge (the atoms that survive into the assembled
     residue, per HEAD01STRIP annotation)
  3. Add backbone (RPO) and sugar (RC3) base charge contributions
     for an INTERNAL residue (target -1.0 e total)
  4. Show the residual + per-atom breakdown

Goal: identify whether 0.32 e is
  (a) a bug in our patched modxna.sh (wrong atoms stripped)
  (b) a feature of modxna's modular charge protocol that requires
      multi-atom redistribution (not just N9)
  (c) a known property documented in Bergonzo 2024 that we missed

Outputs: results/charge_audit.md with per-modification tables.
"""

from __future__ import annotations
import re
import sys
from pathlib import Path
from typing import Dict, List

REPO_ROOT = Path(__file__).resolve().parents[3]
MODXNA    = REPO_ROOT / "project_RRM" / "ver2" / "force_fields" / "modxna"
LIB_BASE  = MODXNA / "dat" / "lib_base"
LIB_BACKBONE = MODXNA / "dat" / "lib_backbone"
LIB_SUGAR    = MODXNA / "dat" / "lib_sugar"

# Targets per modxna.sh:
TARGET_INTERNAL = -1.000000
TARGET_3CAP     = -0.679652
TARGET_5CAP     = -0.320348


def parse_mol2(path: Path) -> dict:
    text = path.read_text(encoding="utf-8").splitlines()
    section = None
    atoms = []
    header = ""
    for raw in text:
        line = raw.rstrip()
        if line.startswith("@<TRIPOS>"):
            section = line.split(">", 1)[1].strip()
            continue
        if section == "MOLECULE":
            if line.strip() and not header:
                header = line.strip()
        elif section == "ATOM":
            parts = line.split()
            if len(parts) >= 9:
                atoms.append({
                    "idx":  int(parts[0]),
                    "name": parts[1],
                    "type": parts[5],
                    "charge": float(parts[8]),
                })
    return {"path": path, "header": header, "atoms": atoms}


def parse_strip_list(header: str) -> dict:
    """Extract HEAD01STRIP, HEAD03STRIP atom name lists from modxna header."""
    # modXNA-fragment:8OG:HEAD01:@N9:HEAD01STRIP:@C1',H1',O4',HHO4,C2',H2'1,H2'3,H2'2
    out = {"head01": None, "strip01": [], "strip03": [],
           "tail01": None, "strip_tail01": []}
    parts = header.split(":")
    i = 2
    while i + 1 < len(parts):
        key = parts[i]; val = parts[i + 1].lstrip("@")
        if key == "HEAD01":          out["head01"] = val
        elif key == "HEAD01STRIP":
            out["strip01"] = [a.strip() for a in val.split(",") if a.strip()]
        elif key == "HEAD03STRIP":
            out["strip03"] = [a.strip() for a in val.split(",") if a.strip()]
        elif key == "TAIL01":        out["tail01"] = val
        elif key == "TAIL01STRIP":
            out["strip_tail01"] = [a.strip() for a in val.split(",") if a.strip()]
        i += 2
    return out


def audit_base(code: str) -> dict:
    p = LIB_BASE / f"{code}.mol2"
    if not p.exists():
        return None
    mol = parse_mol2(p)
    info = parse_strip_list(mol["header"])
    strip_set = set(info["strip01"]) | set(info["strip03"])
    raw_q  = sum(a["charge"] for a in mol["atoms"])
    kept   = [a for a in mol["atoms"] if a["name"] not in strip_set]
    stripped = [a for a in mol["atoms"] if a["name"] in strip_set]
    kept_q = sum(a["charge"] for a in kept)
    stripped_q = sum(a["charge"] for a in stripped)
    return {
        "code": code, "head01": info["head01"],
        "strip_atoms": info["strip01"], "n_atoms": len(mol["atoms"]),
        "n_kept": len(kept), "n_stripped": len(stripped),
        "raw_q": raw_q, "kept_q": kept_q, "stripped_q": stripped_q,
        "kept_atoms": kept,
    }


def audit_backbone(code: str = "RPO") -> dict:
    p = LIB_BACKBONE / f"{code}.mol2"
    mol = parse_mol2(p)
    info = parse_strip_list(mol["header"])
    strip_set = (set(info["strip01"]) | set(info["strip_tail01"]))
    raw_q  = sum(a["charge"] for a in mol["atoms"])
    kept   = [a for a in mol["atoms"] if a["name"] not in strip_set]
    kept_q = sum(a["charge"] for a in kept)
    return {"code": code, "raw_q": raw_q, "kept_q": kept_q,
            "n_kept": len(kept), "n_stripped": len(mol["atoms"]) - len(kept),
            "head01": info["head01"], "tail01": info["tail01"],
            "strip01": info["strip01"], "strip_tail01": info["strip_tail01"]}


def audit_sugar(code: str = "RC3") -> dict:
    p = LIB_SUGAR / f"{code}.mol2"
    mol = parse_mol2(p)
    info = parse_strip_list(mol["header"])
    strip_set = (set(info["strip01"]) | set(info["strip03"])
                 | set(info["strip_tail01"]))
    raw_q  = sum(a["charge"] for a in mol["atoms"])
    kept   = [a for a in mol["atoms"] if a["name"] not in strip_set]
    kept_q = sum(a["charge"] for a in kept)
    return {"code": code, "raw_q": raw_q, "kept_q": kept_q,
            "n_kept": len(kept), "n_stripped": len(mol["atoms"]) - len(kept),
            "head01": info["head01"], "tail01": info["tail01"],
            "strip01": info["strip01"], "strip03": info["strip03"],
            "strip_tail01": info["strip_tail01"]}


def main():
    out_lines = []
    out_lines.append("# Charge audit — modXNA fragment stitching residuals")
    out_lines.append("")
    out_lines.append("Goal: explain the +0.32 e residual on internal 8OG that "
                     "rebalance_lib.py was dumping onto N9. Per expert review:")
    out_lines.append("")
    out_lines.append("> Typical RESP integer-rounding residuals are <0.05 e; "
                     "0.32 e suggests something was off in the modXNA "
                     "fragment-stitching or capping scheme")
    out_lines.append("")

    # ── Backbone & Sugar baseline ─────────────────────────────────────────
    bb = audit_backbone("RPO")
    su = audit_sugar("RC3")
    out_lines.append("## Backbone (RPO) and sugar (RC3) baseline")
    out_lines.append("")
    out_lines.append(f"### RPO  (`{bb['head01']}` head, `{bb['tail01']}` tail)")
    out_lines.append(f"  raw_q = {bb['raw_q']:+.6f} e   "
                     f"kept_q = {bb['kept_q']:+.6f} e   "
                     f"({bb['n_kept']} kept / {bb['n_stripped']} stripped)")
    out_lines.append(f"  HEAD01STRIP atoms: {bb['strip01']}")
    out_lines.append(f"  TAIL01STRIP atoms: {bb['strip_tail01']}")
    out_lines.append("")
    out_lines.append(f"### RC3  (`{su['head01']}` head, `{su['tail01']}` tail)")
    out_lines.append(f"  raw_q = {su['raw_q']:+.6f} e   "
                     f"kept_q = {su['kept_q']:+.6f} e   "
                     f"({su['n_kept']} kept / {su['n_stripped']} stripped)")
    out_lines.append(f"  HEAD01STRIP: {su['strip01']}")
    out_lines.append(f"  HEAD03STRIP: {su['strip03']}")
    out_lines.append(f"  TAIL01STRIP: {su['strip_tail01']}")
    out_lines.append("")

    # ── Each base + arithmetic ────────────────────────────────────────────
    out_lines.append("## Base fragments + assembly arithmetic")
    out_lines.append("")
    out_lines.append("If the modXNA model is `internal residue charge = "
                     "kept_backbone + kept_sugar + kept_base`, then the "
                     "expected internal residue charge depends only on these "
                     "three numbers. Target for ANY internal nucleotide = "
                     f"**{TARGET_INTERNAL:+.6f} e**.")
    out_lines.append("")
    out_lines.append("| code | base raw_q | base kept_q | "
                     "BB.kept + SU.kept + base.kept | residual vs −1.0 |")
    out_lines.append("|------|------:|------:|------:|------:|")
    bb_su = bb['kept_q'] + su['kept_q']
    for code in ["8OG", "PUU", "M1A"]:
        b = audit_base(code)
        if not b: continue
        total_internal = bb_su + b["kept_q"]
        residual = total_internal - TARGET_INTERNAL
        out_lines.append(
            f"| `{code}` | {b['raw_q']:+.4f} | {b['kept_q']:+.4f} | "
            f"{total_internal:+.4f} | {residual:+.4f} |"
        )
    out_lines.append("")

    # ── 8OG detail ────────────────────────────────────────────────────────
    out_lines.append("## 8OG detail")
    out_lines.append("")
    b = audit_base("8OG")
    out_lines.append(f"- Raw 8OG.mol2: {b['n_atoms']} atoms, "
                     f"sum charge = {b['raw_q']:+.6f} e")
    out_lines.append(f"- HEAD01STRIP removes {b['n_stripped']} atoms, "
                     f"removed charge {b['stripped_q']:+.6f} e")
    out_lines.append(f"- Kept (base of 16 atoms): "
                     f"{b['kept_q']:+.6f} e")
    out_lines.append("")
    out_lines.append("### Per-atom (kept) charges in 8OG base fragment:")
    out_lines.append("")
    out_lines.append("| name | type | charge |")
    out_lines.append("|------|------|------:|")
    for a in b["kept_atoms"]:
        out_lines.append(f"| `{a['name']}` | `{a['type']}` | {a['charge']:+.6f} |")
    out_lines.append("")

    out_lines.append("## Interpretation")
    out_lines.append("")
    out_lines.append("If `BB.kept + SU.kept + base.kept` ≈ −1.0 e for 8OG, "
                     "then modxna.sh's stripped charge corrections "
                     "(`charge -0.8832`, `charge -0.01191`, `charge -0.10489`) "
                     "would NOT be needed — the kept-charge sum is already "
                     "close to integer. Our patch removed those corrections, "
                     "which is correct IF the corrections were stale code.")
    out_lines.append("")
    out_lines.append("If the residual is ~0.3 e for 8OG but ~0 for standard G "
                     "(modrna08), then the 8OG base mol2's raw charges have "
                     "been adjusted by modxna's authors to compensate for the "
                     "removed `charge VAL` lines — but the compensation is on "
                     "the WRONG atom, OR the mol2 has not been updated to "
                     "match the broken modxna.sh.")
    out_lines.append("")
    out_lines.append("Check this hypothesis: compare 8OG residual vs PUU and "
                     "M1A residuals above. If all 3 are ~0.3 e, it's a "
                     "consistent modxna issue. If only 8OG is off, it's "
                     "specific to that fragment.")

    out_path = REPO_ROOT / "project_RRM" / "ver2" / "results" / "charge_audit.md"
    out_path.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    print(f"Wrote {out_path}")
    # Also print the table to stdout for quick inspection
    for ln in out_lines:
        print(ln)


if __name__ == "__main__":
    main()
