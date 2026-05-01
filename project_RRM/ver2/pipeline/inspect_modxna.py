# -*- coding: utf-8 -*-
"""
inspect_modxna.py
=================
Read the four base-fragment .mol2 files we need (8OG, PSU, M1A, M6A) and
report:

  - Total atom count, total charge
  - Annotated HEAD atom (connection to sugar) and STRIP atoms
  - Charge sum after stripping (= what the base contributes when assembled
    onto the standard sugar + phosphate backbone)
  - Atom type list, flagging any types not in standard amber14 OL3
  - Atom name duplicates (e.g. the suspected `O6` × 2 in 8OG)

This is purely a verification step. It does not run AmberTools / cpptraj /
tleap. It writes a markdown report to results/inspect_modxna.md.

Usage (from repo root):
    python project_RRM/ver2/pipeline/inspect_modxna.py
"""

from __future__ import annotations
import re
import sys
from collections import Counter
from pathlib import Path

REPO     = Path(__file__).resolve().parents[3]
MODXNA   = REPO / "project_RRM" / "ver2" / "force_fields" / "modxna"
LIB_BASE = MODXNA / "dat" / "lib_base"

# Canonical natural-RNA-modification codes inside modXNA:
#   8OG = 8-oxoguanosine                                              ✅ canonical
#   PUU = pseudouridine (HEAD01 = C5, C-glycosidic — defining ψ trait) ✅ canonical
#         NOTE: modXNA's PSU code is a 2-thio-5-isobutyl variant, NOT
#         canonical pseudouridine — designer chemistry for ASOs.
#   M1A = N1-methyladenosine                                           ✅ canonical
#   m6A = N6-methyladenosine                                           ❌ NOT in modXNA.
#         modXNA's M6A is N6,N6-dimethyladenosine (two methyls on N6).
#         modXNA's DMA is 2,8-dimethyladenosine (C2-Me + C8-Me).
#         Use Bussi 2022 (github.com/bussilab/m6a-charge-fitting) OR
#         derive from Aduri 2007 SI. Tracked separately in force_field_sources.md.
TARGETS  = ["8OG", "PUU", "M1A"]
ALSO_INSPECT_FOR_REFERENCE = ["PSU", "M6A", "DMA"]   # diagnose what these actually are

# Standard atom types seen in amber14/RNA.OL3.xml. Anything outside this set
# is "new" and requires explicit parameters in frcmod.modxna.
OL3_TYPES = {
    "N*", "NA", "NB", "NC", "N2",
    "C", "CA", "CB", "CK", "CM", "CQ", "CR", "CT", "CJ", "CI", "C7", "C2",
    "O", "OS", "OH", "O2",
    "P",
    "H", "HC", "HA", "H1", "H2", "H4", "H5", "HO", "HS",
    "S",
    "F", "Cl", "Br", "I",
}


def parse_mol2_header(line: str) -> dict:
    """Parse a modXNA fragment header.

    Header format example:
      modXNA-fragment:8OG:HEAD01:@N9:HEAD01STRIP:@C1',H1',O4',...

    Returns dict with keys: code, head01, strip01, head03 (optional),
    strip03 (optional), tail01 (optional), strip_tail01 (optional).
    """
    out = {"raw": line.strip()}
    parts = line.strip().split(":")
    # parts[0] = "modXNA-fragment"  or  "modXNA-5'-fragment"
    if len(parts) < 2:
        return out
    out["code"] = parts[1]
    i = 2
    while i + 1 < len(parts):
        key = parts[i]
        val = parts[i + 1].lstrip("@")
        if key == "HEAD01":          out["head01"] = val
        elif key == "HEAD01STRIP":   out["strip01"] = val
        elif key == "HEAD03":        out["head03"] = val
        elif key == "HEAD03STRIP":   out["strip03"] = val
        elif key == "TAIL01":        out["tail01"] = val
        elif key == "TAIL01STRIP":   out["strip_tail01"] = val
        i += 2
    return out


def parse_mol2(path: Path) -> dict:
    text = path.read_text(encoding="utf-8").splitlines()
    section = None
    atoms: list[dict] = []
    bonds: list[tuple[int, int]] = []
    header: dict = {}

    for raw in text:
        line = raw.rstrip()
        if line.startswith("@<TRIPOS>"):
            section = line.split(">", 1)[1].strip()
            continue
        if section == "MOLECULE":
            # First non-empty line under MOLECULE is the name; modXNA
            # encodes fragment metadata in this name string.
            if line.strip() and not header:
                header = parse_mol2_header(line)
        elif section == "ATOM":
            # AtomNo  Name  X Y Z  Type  ResNo  ResName  Charge
            parts = line.split()
            if len(parts) < 9:
                continue
            atoms.append({
                "idx":  int(parts[0]),
                "name": parts[1],
                "type": parts[5],
                "charge": float(parts[8]),
            })
        elif section == "BOND":
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                bonds.append((int(parts[1]), int(parts[2])))
            except ValueError:
                continue

    return {"path": path, "header": header, "atoms": atoms, "bonds": bonds}


def comma_split(s: str | None) -> list[str]:
    if not s:
        return []
    return [x.strip() for x in s.split(",") if x.strip()]


def report_one(mol: dict) -> list[str]:
    out: list[str] = []
    code = mol["header"].get("code", "?")
    atoms = mol["atoms"]
    out.append(f"## {code}  —  `{mol['path'].name}`")
    out.append("")
    out.append(f"- Header: `{mol['header']['raw']}`")
    out.append(f"- Total atoms: **{len(atoms)}**")
    out.append(f"- Total charge (all atoms): **{sum(a['charge'] for a in atoms):+.6f} e**")

    # Strip lists
    strip01 = comma_split(mol["header"].get("strip01"))
    strip03 = comma_split(mol["header"].get("strip03"))
    strip_all = set(strip01) | set(strip03)
    out.append(f"- HEAD01 atom (connects to sugar): `{mol['header'].get('head01', '?')}`")
    out.append(f"- HEAD01STRIP atoms (dropped on assembly): `{', '.join(strip01) or 'none'}`")
    if strip03:
        out.append(f"- HEAD03STRIP atoms: `{', '.join(strip03)}`")

    # Net base contribution = total - stripped
    kept = [a for a in atoms if a["name"] not in strip_all]
    stripped = [a for a in atoms if a["name"] in strip_all]
    out.append(
        f"- After stripping ({len(stripped)} atoms removed): "
        f"**{len(kept)} atoms** remain, "
        f"sum charge = **{sum(a['charge'] for a in kept):+.6f} e**"
    )

    # Atom name duplicates (illegal in PDB; need rename before tleap/OpenMM)
    name_counts = Counter(a["name"] for a in atoms)
    dup = {n: c for n, c in name_counts.items() if c > 1}
    if dup:
        out.append(f"- ⚠️ **Duplicate atom names**: {dup}")
        for n in dup:
            occurrences = [(a["idx"], a["type"], a["charge"]) for a in atoms if a["name"] == n]
            out.append(f"    - `{n}`: " + "; ".join(
                f"idx={i} type={t} q={q:+.4f}" for i, t, q in occurrences))
    else:
        out.append(f"- Atom name uniqueness: ✅ all {len(atoms)} unique")

    # Atom types — flag any not in standard OL3 set
    types_used = sorted({a["type"] for a in atoms})
    new_types = [t for t in types_used if t not in OL3_TYPES]
    out.append(f"- Atom types used: `{', '.join(types_used)}`")
    if new_types:
        out.append(f"- ⚠️ Types **outside standard OL3**: `{', '.join(new_types)}` "
                   f"(must be defined in frcmod.modxna)")
    else:
        out.append(f"- All types are in standard OL3 set ✅")

    # Detailed atom table
    out.append("")
    out.append("| idx | name | type | charge | strip? |")
    out.append("|----:|------|------|-------:|--------|")
    for a in atoms:
        flag = "STRIP" if a["name"] in strip_all else ""
        out.append(f"| {a['idx']:>3d} | `{a['name']}` | `{a['type']}` | "
                   f"{a['charge']:+.6f} | {flag} |")
    out.append("")
    return out


def main():
    report: list[str] = []
    report.append("# modXNA Base-Fragment Inspection Report")
    report.append("")
    report.append(f"Source: `{MODXNA.relative_to(REPO)}` "
                  f"(commit recorded in [`VERSIONS.md`]"
                  f"(../force_fields/VERSIONS.md))")
    report.append("")
    report.append(f"Targets: {', '.join(TARGETS)}")
    report.append("")
    report.append("## Verification questions")
    report.append("")
    report.append("1. Do the 4 base mol2 files exist and parse cleanly?")
    report.append("2. What's the **net charge contribution** of the base "
                  "(after STRIP atoms removed)? In a properly assembled "
                  "OL3-RNA nucleotide the base should contribute "
                  "**~+0.00 to +0.40 e** (the rest of the −1.0 sits on "
                  "the phosphate + sugar).")
    report.append("3. Are there any **duplicate atom names** that would "
                  "break tleap or OpenMM template matching?")
    report.append("4. Are there any **atom types outside standard OL3** "
                  "that would need new definitions in the OpenMM XML?")
    report.append("")

    found = []
    for code in TARGETS:
        path = LIB_BASE / f"{code}.mol2"
        if not path.exists():
            report.append(f"## {code}  —  **MISSING**: {path}")
            report.append("")
            continue
        mol = parse_mol2(path)
        report.extend(report_one(mol))
        found.append(code)

    report.append("---")
    report.append("")
    report.append("# Reference: misleadingly-named modXNA codes")
    report.append("")
    report.append("These are **NOT** what their name suggests in modXNA's library. "
                  "Inspecting them here so the gotchas are documented in one place.")
    report.append("")
    for code in ALSO_INSPECT_FOR_REFERENCE:
        path = LIB_BASE / f"{code}.mol2"
        if not path.exists():
            continue
        mol = parse_mol2(path)
        report.extend(report_one(mol))

    out_dir = REPO / "project_RRM" / "ver2" / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    report_path = out_dir / "inspect_modxna.md"
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"Found: {found}")
    print(f"Report written: {report_path.relative_to(REPO)}")


if __name__ == "__main__":
    main()
