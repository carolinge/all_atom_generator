# -*- coding: utf-8 -*-
"""
plot_per_replica.py
===================
Generate per-replica time-series figures (one row per replica) so we can
see whether any specific replica behaves anomalously (e.g. partial
unbinding) — information that gets averaged out in the cross-replica
mean ± SEM aggregate plot.

Usage:
    python plot_per_replica.py \\
        --in <results_dir> \\
        --out <results_dir>/_aggregate/figures
"""

from __future__ import annotations
import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


COLORS = {"WT": "#2563eb", "8OG": "#dc2626"}


def load_cv(rep_dir: Path) -> dict | None:
    cv_p = rep_dir / "cv.csv"
    if not cv_p.exists():
        return None
    arr = np.loadtxt(cv_p, delimiter=",", skiprows=1)
    return {
        "t_ns":         arr[:, 0] / 1000.0,
        "prot_rmsd":    arr[:, 1],
        "rna_rmsd":     arr[:, 2],
        "n_contacts":   arr[:, 3],
        "com_dist":     arr[:, 4],
        "g3_contacts":  arr[:, 5],
        "g3_hbonds":    arr[:, 6],
    }


def pbc_mask(d):
    return ((d["rna_rmsd"] < 25.0)
            & (d["n_contacts"] < 1500)
            & (d["com_dist"] < 30.0))


def plot_systems(in_root: Path, out_dir: Path,
                 systems: dict[str, list[Path]]):
    out_dir.mkdir(parents=True, exist_ok=True)

    METRICS = [
        ("rna_rmsd",     "RNA heavy-atom RMSD (Å)",     5,  None),
        ("n_contacts",   "Total prot-RNA contacts",     0,  None),
        ("g3_contacts",  "G3 contacts with protein",    0,  None),
        ("g3_hbonds",    "G3 polar contacts with protein", 0, None),
    ]

    n_metrics = len(METRICS)
    n_replicas = max(len(v) for v in systems.values())
    fig, axes = plt.subplots(n_replicas, n_metrics,
                              figsize=(4 * n_metrics, 3 * n_replicas),
                              sharex=True, squeeze=False)

    for col, (key, ylab, ymin, ymax) in enumerate(METRICS):
        for row in range(n_replicas):
            ax = axes[row, col]
            for sys_label, replicas in systems.items():
                if row >= len(replicas):
                    continue
                rep_dir = replicas[row]
                d = load_cv(rep_dir)
                if d is None:
                    continue
                mask = pbc_mask(d)
                col_color = COLORS.get(sys_label, "k")
                # Plot all data (faded) and PBC-good data (bold)
                ax.plot(d["t_ns"], d[key], color=col_color, lw=0.5, alpha=0.25)
                # Mask out outliers as gaps
                masked = d[key].copy()
                masked[~mask] = np.nan
                ax.plot(d["t_ns"], masked, color=col_color, lw=1.4,
                        label=f"{sys_label} ({rep_dir.name})")
            if row == 0:
                ax.set_title(ylab.split(" (")[0])
            if col == 0:
                ax.set_ylabel(f"r{row+1}\n{ylab}")
            if row == n_replicas - 1:
                ax.set_xlabel("Time (ns)")
            ax.legend(loc="best", fontsize=8, frameon=False)
            ax.grid(alpha=0.3)
            if ymin is not None:
                ax.set_ylim(bottom=ymin)
    plt.suptitle("Per-replica time series — solid = PBC-good frames, "
                 "faded thin = raw (with PBC artefacts)",
                 fontsize=11, y=0.995)
    plt.tight_layout()
    plt.savefig(out_dir / "per_replica_timeseries.png",
                dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Wrote {out_dir / 'per_replica_timeseries.png'}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in",  dest="in_root",  type=Path, required=True)
    ap.add_argument("--out", dest="out_dir",  type=Path, required=True)
    args = ap.parse_args()

    systems = {
        "WT":  sorted((args.in_root / "4BS2_WT").glob("r*")),
        "8OG": sorted((args.in_root / "4BS2_8OG_G3").glob("r*")),
    }
    plot_systems(args.in_root, args.out_dir, systems)


if __name__ == "__main__":
    main()
