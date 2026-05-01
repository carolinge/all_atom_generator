# -*- coding: utf-8 -*-
"""Generate side-by-side comparison plots for WT vs 8OG binding metrics."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

ANAL = Path("d:/all_atom/output/analysis")

wt  = np.loadtxt(ANAL / "WT_cv.csv",  delimiter=",", skiprows=1)
og  = np.loadtxt(ANAL / "8OG_cv.csv", delimiter=",", skiprows=1)
# columns: time_ps, rna_rmsd_A, n_contacts, com_dist_A, g3_contacts, g3_hbond_polar

t_wt = wt[:, 0] / 1000.0  # ns
t_og = og[:, 0] / 1000.0

fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
labels = ["WT", "8OG (G3)"]
colors = ["#2563eb", "#dc2626"]

panels = [
    (1, "RNA RMSD", "RMSD (Å)"),
    (2, "Total protein-RNA contacts (heavy atoms <4 Å)", "contacts"),
    (4, "G3 contacts with protein", "contacts"),
    (5, "G3 base polar contacts (<3.5 Å)", "polar contacts"),
]

for ax, (col, title, ylab) in zip(axes.flat, panels):
    ax.plot(t_wt, wt[:, col], color=colors[0], alpha=0.4, lw=0.6)
    ax.plot(t_og, og[:, col], color=colors[1], alpha=0.4, lw=0.6)
    # rolling mean (window of 25 frames = 0.5 ns)
    w = 25
    if len(wt) > w:
        rm_wt = np.convolve(wt[:, col], np.ones(w)/w, mode="valid")
        rm_og = np.convolve(og[:, col], np.ones(w)/w, mode="valid")
        ax.plot(t_wt[w-1:], rm_wt, color=colors[0], lw=2.0, label=labels[0])
        ax.plot(t_og[w-1:], rm_og, color=colors[1], lw=2.0, label=labels[1])
    ax.set_title(title, fontsize=11)
    ax.set_ylabel(ylab)
    ax.grid(alpha=0.3)
    ax.legend(loc="best", frameon=False, fontsize=9)

axes[1, 0].set_xlabel("Time (ns)")
axes[1, 1].set_xlabel("Time (ns)")
plt.suptitle("4BS2: Wild-type vs 8-oxoG (G3) — protein-RNA binding (10 ns)",
             fontsize=13, y=0.995)
plt.tight_layout()
plt.savefig(ANAL / "compare_timeseries.png", dpi=150, bbox_inches="tight")
plt.close()
print("Wrote compare_timeseries.png")

# ── Difference contact map ────────────────────────────────────────────────────
cm_wt = np.loadtxt(ANAL / "WT_contact_map.csv",  delimiter=",")
cm_og = np.loadtxt(ANAL / "8OG_contact_map.csv", delimiter=",")
import json
meta_wt = json.load(open(ANAL / "WT_meta.json"))
meta_og = json.load(open(ANAL / "8OG_meta.json"))
prot_wt = meta_wt["prot_residues"]
rna_wt  = meta_wt["rna_residues"]

# Difference: 8OG minus WT
diff = cm_og - cm_wt

fig, axes = plt.subplots(1, 3, figsize=(20, 4), gridspec_kw={"width_ratios": [1, 1, 1.05]})
vmax = max(cm_wt.max(), cm_og.max())
for ax, cm, ttl in zip(axes[:2], [cm_wt, cm_og], ["WT", "8OG (G3)"]):
    im = ax.imshow(cm, aspect="auto", cmap="hot", vmin=0, vmax=vmax,
                   extent=[prot_wt[0], prot_wt[-1], rna_wt[-1], rna_wt[0]])
    ax.set_title(f"{ttl} contact map")
    ax.set_xlabel("Protein residue (chain A)")
    ax.set_ylabel("RNA residue (chain B)")
    plt.colorbar(im, ax=ax, label="contacts/frame")

# Diff with diverging colormap
vmax_d = max(abs(diff.min()), abs(diff.max()))
im = axes[2].imshow(diff, aspect="auto", cmap="RdBu_r", vmin=-vmax_d, vmax=vmax_d,
                    extent=[prot_wt[0], prot_wt[-1], rna_wt[-1], rna_wt[0]])
axes[2].set_title("Δ contacts (8OG − WT)\nblue = lost, red = gained")
axes[2].set_xlabel("Protein residue (chain A)")
axes[2].set_ylabel("RNA residue (chain B)")
plt.colorbar(im, ax=axes[2], label="Δ contacts/frame")

plt.tight_layout()
plt.savefig(ANAL / "compare_contactmap.png", dpi=150, bbox_inches="tight")
plt.close()
print("Wrote compare_contactmap.png")
