# -*- coding: utf-8 -*-
"""Plot RMSF comparison WT vs 8OG, highlight contact-region residues."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
import json

ANAL = Path("d:/all_atom/output/analysis")

wt_p = np.loadtxt(ANAL / "WT_rmsf_protein.csv",  delimiter=",", skiprows=1)
og_p = np.loadtxt(ANAL / "8OG_rmsf_protein.csv", delimiter=",", skiprows=1)
wt_r = np.loadtxt(ANAL / "WT_rmsf_rna.csv",      delimiter=",", skiprows=1)
og_r = np.loadtxt(ANAL / "8OG_rmsf_rna.csv",     delimiter=",", skiprows=1)

# Load contact map metadata to identify contact-region residues
import json
meta = json.load(open(ANAL / "WT_meta.json"))
cm   = np.loadtxt(ANAL / "WT_contact_map.csv", delimiter=",")
prot_residues = meta["prot_residues"]
# A residue is "contact-region" if it has > 1 contact/frame on average for any
# RNA residue (i.e., its column max > 1)
contact_strength = cm.max(axis=0)
contact_protein  = [prot_residues[i] for i, c in enumerate(contact_strength) if c > 1.0]
print(f"Contact protein residues (>1 contact/frame to any RNA): {contact_protein}")

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 1, figsize=(12, 7))

# Protein
ax = axes[0]
ax.plot(wt_p[:, 0], wt_p[:, 1], color="#2563eb", lw=1.2, label="WT")
ax.plot(og_p[:, 0], og_p[:, 1], color="#dc2626", lw=1.2, label="8OG (G3)")
# Highlight contact residues with vertical bands
for r in contact_protein:
    ax.axvspan(r - 0.5, r + 0.5, color="gold", alpha=0.18, zorder=0)
ax.set_xlabel("Protein residue (chain A)")
ax.set_ylabel("RMSF (Å)")
ax.set_title("Per-residue RMSF — protein (chain A).  "
             "Gold bands = RNA-contact residues (>1 contact/frame)")
ax.legend(loc="best", frameon=False)
ax.grid(alpha=0.3)
ax.set_xlim(prot_residues[0], prot_residues[-1])

# RNA
ax = axes[1]
ax.plot(wt_r[:, 0], wt_r[:, 1], "o-", color="#2563eb", lw=2, ms=7, label="WT")
ax.plot(og_r[:, 0], og_r[:, 1], "s-", color="#dc2626", lw=2, ms=7, label="8OG (G3)")
# Highlight modified residue
ax.axvspan(2.5, 3.5, color="gold", alpha=0.25, zorder=0,
           label="Modified residue (G3)")
ax.set_xlabel("RNA residue (chain B)")
ax.set_ylabel("RMSF (Å)")
ax.set_title("Per-residue RMSF — RNA (chain B)")
ax.legend(loc="best", frameon=False)
ax.grid(alpha=0.3)
ax.set_xticks(range(1, 13))

plt.tight_layout()
plt.savefig(ANAL / "compare_rmsf.png", dpi=150, bbox_inches="tight")
plt.close()
print("Wrote compare_rmsf.png")

# ── Quantitative summary: contact-region vs non-contact ───────────────────────
print("\n=== Protein RMSF averages ===")
contact_set = set(contact_protein)
wt_p_contact = [r for rid, r in wt_p if int(rid) in contact_set]
og_p_contact = [r for rid, r in og_p if int(rid) in contact_set]
wt_p_other   = [r for rid, r in wt_p if int(rid) not in contact_set]
og_p_other   = [r for rid, r in og_p if int(rid) not in contact_set]
print(f"  Contact residues  WT: {np.mean(wt_p_contact):.2f} ± {np.std(wt_p_contact):.2f}  "
      f"8OG: {np.mean(og_p_contact):.2f} ± {np.std(og_p_contact):.2f}")
print(f"  Non-contact      WT: {np.mean(wt_p_other):.2f} ± {np.std(wt_p_other):.2f}  "
      f"8OG: {np.mean(og_p_other):.2f} ± {np.std(og_p_other):.2f}")

lines = []
lines.append("\n=== RNA RMSF (residue-level) ===")
lines.append(f"  resid   WT (A)   8OG (A)   delta (8OG-WT)")
for i in range(len(wt_r)):
    rid = int(wt_r[i, 0])
    delta = og_r[i, 1] - wt_r[i, 1]
    marker = " <-- G3 modified" if rid == 3 else ""
    lines.append(f"   {rid:3d}    {wt_r[i,1]:5.2f}    {og_r[i,1]:5.2f}    {delta:+5.2f}{marker}")

lines.append("\n=== Core binding region (residues 3-9) ===")
core_mask_wt = (wt_r[:, 0] >= 3) & (wt_r[:, 0] <= 9)
core_mask_og = (og_r[:, 0] >= 3) & (og_r[:, 0] <= 9)
lines.append(f"  WT core RMSF mean : {wt_r[core_mask_wt, 1].mean():.2f} A")
lines.append(f"  8OG core RMSF mean: {og_r[core_mask_og, 1].mean():.2f} A")

txt = "\n".join(lines)
print(txt)
(ANAL / "rmsf_comparison.txt").write_text(txt, encoding="utf-8")
