# -*- coding: utf-8 -*-
"""
Multi-replica analysis with smart handling of partial trajectories.

For each replica:
  - load solvated.pdb + traj.dcd (whatever frames exist)
  - compute per-frame CVs and per-residue RMSF

Across replicas:
  - at each time point t, average only those replicas with len(traj) >= t
  - report mean ± SEM(n_avail) at each t
  - report n_replicas covering each time point

Usage:
  python analyze_replicas.py \\
      --label WT  --replicas <dir1> <dir2> <dir3> \\
      --label 8OG --replicas <dir4> <dir5> <dir6> \\
      --out anal_gyrate/multi
"""

import sys, json, argparse, os
from pathlib import Path
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
from MDAnalysis.analysis.distances import distance_array


def parse_args():
    ap = argparse.ArgumentParser()
    # Repeated --label / --replicas pairs
    ap.add_argument("--label", action="append", required=True,
                    help="system label (e.g., WT, 8OG); repeat per system")
    ap.add_argument("--replicas", action="append", nargs="+", required=True,
                    help="replica directories for the preceding --label")
    ap.add_argument("--out", required=True)
    ap.add_argument("--mod-resid", type=int, default=3)
    ap.add_argument("--contact-cutoff", type=float, default=4.0)
    ap.add_argument("--stride", type=int, default=1)
    return ap.parse_args()


def analyze_one_replica(rep_dir: Path, mod_resid: int, contact_cutoff: float,
                         stride: int):
    """Return dict with per-frame CVs and per-residue RMSF for one replica."""
    pdb_path  = rep_dir / "solvated.pdb"
    traj_path = rep_dir / "traj.dcd"
    if not pdb_path.exists() or not traj_path.exists():
        print(f"  [{rep_dir.name}] missing solvated.pdb or traj.dcd, skip")
        return None

    print(f"  [{rep_dir.name}] loading ...")
    u = mda.Universe(str(pdb_path), str(traj_path))
    n_frames = len(u.trajectory) // stride
    if n_frames == 0:
        print(f"  [{rep_dir.name}] empty trajectory, skip")
        return None

    # Selections (try segid then chainID)
    prot = u.select_atoms("segid A and not name H*")
    rna  = u.select_atoms("segid B and not name H*")
    if len(prot) == 0:
        prot = u.select_atoms("chainID A and not name H*")
        rna  = u.select_atoms("chainID B and not name H*")
    g3       = u.select_atoms(f"resid {mod_resid} and not name H*") & rna
    g3_base  = u.select_atoms(
        f"resid {mod_resid} and "
        "(name N1 N2 N3 N7 N9 C2 C4 C5 C6 C8 O6 O8)") & rna

    align_sel = "segid A and name CA" if len(u.select_atoms("segid A and name CA")) > 0 \
                else "chainID A and name CA"
    ref = mda.Universe(str(pdb_path))
    align.AlignTraj(u, ref, select=align_sel, in_memory=True).run()

    rna_residues  = sorted({a.resid for a in rna})
    prot_residues = sorted({a.resid for a in prot})
    rna_groups  = [rna.select_atoms(f"resid {r}")  for r in rna_residues]
    prot_groups = [prot.select_atoms(f"resid {r}") for r in prot_residues]
    for ag in rna_groups + prot_groups:
        ag._buf = np.zeros((n_frames, len(ag), 3), dtype=np.float32)

    # Reference RNA positions for RMSD (PBC-corrected at frame 0)
    u.trajectory[0]
    box0 = u.trajectory.ts.dimensions[:3].astype(float)
    prot_com0 = prot.center_of_mass()
    rna_ref = rna.positions.copy()
    rna_ref -= np.round((rna_ref - prot_com0) / box0) * box0

    # Polar protein atoms for H-bond proxy
    prot_polar = u.select_atoms(
        "(segid A or chainID A) and (name N* or name O* or name S*) and not name H*")

    times = np.zeros(n_frames)
    rna_rmsd = np.zeros(n_frames)
    n_contacts = np.zeros(n_frames, dtype=np.int32)
    com_dist = np.zeros(n_frames)
    g3_contacts = np.zeros(n_frames, dtype=np.int32)
    g3_polar    = np.zeros(n_frames, dtype=np.int32)

    f = 0
    for ts in u.trajectory[::stride]:
        d = distance_array(prot.positions, rna.positions, box=ts.dimensions)
        d_g3 = distance_array(g3.positions, prot.positions, box=ts.dimensions)
        d_hb = (distance_array(g3_base.positions, prot_polar.positions,
                               box=ts.dimensions)
                if len(g3_base) > 0 else np.empty(0))

        # PBC-aware COM distance: minimum-image
        box3 = ts.dimensions[:3].astype(float)
        prot_com = prot.center_of_mass()
        rna_com  = rna.center_of_mass()
        com_vec  = rna_com - prot_com
        com_vec -= np.round(com_vec / box3) * box3   # min-image
        com_dist[f]    = float(np.linalg.norm(com_vec))

        # PBC-aware RNA RMSD: shift RNA positions to nearest image of protein COM
        rna_pos = rna.positions.copy()
        shift   = np.round((rna_pos - prot_com) / box3) * box3
        rna_pos -= shift
        # rna_ref was captured at frame 0 (already correct for that frame)
        rna_rmsd[f]    = rms.rmsd(rna_pos, rna_ref, superposition=False)

        times[f]       = ts.time
        n_contacts[f]  = int(np.sum(d < contact_cutoff))
        g3_contacts[f] = int(np.sum(d_g3 < contact_cutoff))
        g3_polar[f]    = int(np.sum(d_hb < 3.5)) if d_hb.size else 0

        # For RMSF: store PBC-corrected positions for RNA; protein is fine after
        # CA alignment, but apply same correction defensively.
        for ag in rna_groups:
            pos = ag.positions.copy()
            pos -= np.round((pos - prot_com) / box3) * box3
            ag._buf[f] = pos
        for ag in prot_groups:
            ag._buf[f] = ag.positions
        f += 1
        if ts.frame % 200 == 0:
            print(f"    frame {ts.frame}  t={ts.time/1000:.1f} ns")

    # RMSF per residue
    rna_rmsf = []
    for ag in rna_groups:
        coords = ag._buf
        mean = coords.mean(axis=0)
        sq   = ((coords - mean) ** 2).sum(axis=2)
        rna_rmsf.append(np.sqrt(sq.mean(axis=1).mean()))
    prot_rmsf = []
    for ag in prot_groups:
        coords = ag._buf
        mean = coords.mean(axis=0)
        sq   = ((coords - mean) ** 2).sum(axis=2)
        prot_rmsf.append(np.sqrt(sq.mean(axis=1).mean()))

    # Free buffers
    for ag in rna_groups + prot_groups:
        del ag._buf

    return {
        "n_frames": n_frames,
        "times":      times,
        "rna_rmsd":   rna_rmsd,
        "n_contacts": n_contacts,
        "com_dist":   com_dist,
        "g3_contacts": g3_contacts,
        "g3_polar":    g3_polar,
        "rna_residues":  rna_residues,
        "prot_residues": prot_residues,
        "rna_rmsf":      np.array(rna_rmsf),
        "prot_rmsf":     np.array(prot_rmsf),
    }


def stack_with_padding(arrs):
    """Stack 1-D arrays of different lengths into a 2-D array padded with NaN.
    Shape: (n_arr, max_len). NaN means 'this replica did not reach this point'."""
    if not arrs:
        return np.empty((0, 0))
    max_len = max(len(a) for a in arrs)
    out = np.full((len(arrs), max_len), np.nan, dtype=float)
    for i, a in enumerate(arrs):
        out[i, :len(a)] = a
    return out


def mean_sem_with_n(stacked):
    """Return mean, SEM and n_avail at each column ignoring NaN."""
    n_avail = np.sum(~np.isnan(stacked), axis=0)
    with np.errstate(invalid="ignore"):
        mean = np.nanmean(stacked, axis=0)
        sd   = np.nanstd(stacked, axis=0, ddof=1)
    sem = sd / np.sqrt(np.maximum(n_avail, 1))
    return mean, sem, n_avail


def main():
    args = parse_args()
    out_prefix = Path(args.out)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    if len(args.label) != len(args.replicas):
        sys.exit("--label and --replicas must alternate (one --label per group)")

    systems = {}   # label -> list of replica analyses
    for lab, reps in zip(args.label, args.replicas):
        systems[lab] = []
        print(f"\n=== System: {lab} ===")
        for rd in reps:
            res = analyze_one_replica(Path(rd), args.mod_resid,
                                      args.contact_cutoff, args.stride)
            if res is not None:
                systems[lab].append((Path(rd).name, res))

    # Save raw per-replica results as a JSON-friendly summary
    summary = {}
    for lab, reps in systems.items():
        summary[lab] = []
        for name, res in reps:
            summary[lab].append({
                "replica": name,
                "n_frames": int(res["n_frames"]),
                "t_max_ns": float(res["times"][-1] / 1000.0),
                "mean_rna_rmsd_A": float(np.mean(res["rna_rmsd"])),
                "mean_n_contacts": float(np.mean(res["n_contacts"])),
                "mean_g3_contacts": float(np.mean(res["g3_contacts"])),
                "mean_g3_polar":    float(np.mean(res["g3_polar"])),
            })
    json.dump(summary, open(f"{out_prefix}_per_replica.json", "w"), indent=2)

    # ── Time-series with cross-replica averaging ──────────────────────────────
    # We assume all replicas use the same dt (they do because run_openmm.py
    # writes every 5000 steps = 20 ps).  Use the time axis of the longest
    # replica from each system.
    ts_data = {}   # label -> dict of stacked arrays
    for lab, reps in systems.items():
        if not reps:
            continue
        # Stack each metric
        stacks = {k: stack_with_padding([r[k] for _, r in reps])
                  for k in ("times", "rna_rmsd", "n_contacts", "com_dist",
                            "g3_contacts", "g3_polar")}
        # Collapse times: take max for each column
        times = np.nanmax(stacks["times"], axis=0) / 1000.0   # ns
        ts_data[lab] = {"times_ns": times}
        for k in ("rna_rmsd", "n_contacts", "com_dist", "g3_contacts", "g3_polar"):
            mean, sem, n_avail = mean_sem_with_n(stacks[k])
            ts_data[lab][k] = {"mean": mean, "sem": sem, "n": n_avail}

    # Save numpy archive
    np.savez(f"{out_prefix}_timeseries.npz", **{
        f"{lab}_{k}_{stat}": arr
        for lab, dat in ts_data.items()
        for k, v in dat.items()
        if k != "times_ns"
        for stat, arr in v.items()
    }, **{f"{lab}_times_ns": dat["times_ns"] for lab, dat in ts_data.items()})

    # ── Plot ──────────────────────────────────────────────────────────────────
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    metrics = [
        ("rna_rmsd",    "RNA RMSD (Å)",                     1.0),
        ("n_contacts",  "Total contacts (<4 Å)",            1.0),
        ("g3_contacts", "G3 contacts with protein",         1.0),
        ("g3_polar",    "G3 base polar contacts (<3.5 Å)",  1.0),
    ]
    colors = {"WT": "#2563eb", "8OG": "#dc2626"}
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    for ax, (key, label, _) in zip(axes.flat, metrics):
        for lab, dat in ts_data.items():
            t = dat["times_ns"]
            m = dat[key]["mean"]
            s = dat[key]["sem"]
            n = dat[key]["n"]
            col = colors.get(lab, "gray")
            ax.plot(t, m, color=col, lw=1.8, label=f"{lab} (n={int(n.max())})")
            ax.fill_between(t, m - s, m + s, color=col, alpha=0.18)
            # Mark where n drops below max (i.e., partial coverage)
            n_max = n.max()
            if any(n < n_max):
                # Vertical line where coverage drops
                drop_idx = np.where(n < n_max)[0]
                if drop_idx.size:
                    ax.axvline(t[drop_idx[0]], color=col, ls=":", alpha=0.6)
        ax.set_title(label)
        ax.set_ylabel(label.split("(")[0].strip())
        ax.grid(alpha=0.3)
        ax.legend(loc="best", fontsize=9, frameon=False)
    axes[1, 0].set_xlabel("Time (ns)")
    axes[1, 1].set_xlabel("Time (ns)")
    plt.suptitle("4BS2 multi-replica binding metrics — mean ± SEM across replicas\n"
                 "(dotted vertical = coverage drop; band = SEM over n available replicas)",
                 fontsize=11, y=0.995)
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_timeseries.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\nWrote {out_prefix}_timeseries.png")

    # ── RMSF mean/SEM across replicas ─────────────────────────────────────────
    rmsf_summary = {}
    for lab, reps in systems.items():
        if not reps:
            continue
        prot_residues = reps[0][1]["prot_residues"]
        rna_residues  = reps[0][1]["rna_residues"]
        prot_stack = np.stack([r["prot_rmsf"] for _, r in reps])  # (n_rep, n_res)
        rna_stack  = np.stack([r["rna_rmsf"]  for _, r in reps])
        rmsf_summary[lab] = {
            "prot_residues": prot_residues,
            "rna_residues":  rna_residues,
            "prot_mean": prot_stack.mean(axis=0),
            "prot_sem":  prot_stack.std(axis=0, ddof=1) / np.sqrt(len(reps)),
            "rna_mean":  rna_stack.mean(axis=0),
            "rna_sem":   rna_stack.std(axis=0, ddof=1) / np.sqrt(len(reps)),
            "n_rep":     len(reps),
        }

    np.savez(f"{out_prefix}_rmsf.npz", **{
        f"{lab}_{k}": v for lab, d in rmsf_summary.items() for k, v in d.items()
        if k != "n_rep" and not isinstance(v, list)
    })

    fig, axes = plt.subplots(2, 1, figsize=(12, 7))
    ax = axes[0]
    for lab, d in rmsf_summary.items():
        col = colors.get(lab, "gray")
        ax.plot(d["prot_residues"], d["prot_mean"], color=col, lw=1.2,
                label=f"{lab} (n={d['n_rep']})")
        ax.fill_between(d["prot_residues"],
                        d["prot_mean"] - d["prot_sem"],
                        d["prot_mean"] + d["prot_sem"],
                        color=col, alpha=0.18)
    ax.set_xlabel("Protein residue (chain A)")
    ax.set_ylabel("RMSF (Å)")
    ax.set_title("Per-residue RMSF — protein chain A (mean ± SEM across replicas)")
    ax.legend(loc="best", frameon=False)
    ax.grid(alpha=0.3)

    ax = axes[1]
    for lab, d in rmsf_summary.items():
        col = colors.get(lab, "gray")
        ax.errorbar(d["rna_residues"], d["rna_mean"], yerr=d["rna_sem"],
                    fmt="o-", color=col, lw=2, ms=7, capsize=3,
                    label=f"{lab} (n={d['n_rep']})")
    ax.axvspan(2.5, 3.5, color="gold", alpha=0.25, label="Modified G3")
    ax.set_xlabel("RNA residue (chain B)")
    ax.set_ylabel("RMSF (Å)")
    ax.set_title("Per-residue RMSF — RNA chain B")
    ax.legend(loc="best", frameon=False)
    ax.grid(alpha=0.3)
    ax.set_xticks(range(1, 13))
    plt.tight_layout()
    plt.savefig(f"{out_prefix}_rmsf.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Wrote {out_prefix}_rmsf.png")

    # ── Print summary ─────────────────────────────────────────────────────────
    lines = ["\n=== Per-replica summary ==="]
    for lab, items in summary.items():
        lines.append(f"\n[{lab}]")
        for it in items:
            lines.append(f"  {it['replica']}  t_max={it['t_max_ns']:5.1f}ns "
                         f"contacts={it['mean_n_contacts']:.0f} "
                         f"g3_polar={it['mean_g3_polar']:.1f}")
    txt = "\n".join(lines)
    print(txt)
    Path(f"{out_prefix}_summary.txt").write_text(txt, encoding="utf-8")


if __name__ == "__main__":
    main()
