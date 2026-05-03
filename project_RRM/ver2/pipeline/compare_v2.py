# -*- coding: utf-8 -*-
"""
compare_v2.py
=============
Cross-replica aggregation + WT vs 8OG comparison + statistical tests +
publication-quality plots.

Reads per-replica CSVs/npz produced by analyze_v2.py and writes a single
results bundle in <out_root>/_aggregate/:

    timeseries.npz     mean +/- SEM at each time point per system, per metric
    rmsf.npz           per-residue RMSF mean +/- SEM per system
    contact_diff.npz   8OG mean - WT mean per (RNA-residue, prot-residue)
    chi_distributions.npz histograms of chi for the modified base
    pucker_distributions.npz histograms of P (pseudorotation) for modified base
    stats.csv          Welch's t / KS / Mann-Whitney p-values for key metrics
    figures/*.png      publication-quality figures
    report.md          markdown report tying it all together

Usage:
    python compare_v2.py \\
        --label-replicas WT  4BS2_WT/r1 4BS2_WT/r2 4BS2_WT/r3 \\
        --label-replicas 8OG 4BS2_8OG_G3/r1 4BS2_8OG_G3/r2 4BS2_8OG_G3/r3 \\
        --mod-resid 272 \\
        --in  /data/.../results \\
        --out /data/.../results/_aggregate
"""

from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
from typing import Sequence

import numpy as np

# Stats — use scipy if available, else a tiny fallback
try:
    from scipy import stats as scstats
    SCIPY_OK = True
except ImportError:
    SCIPY_OK = False


# ── Loaders ─────────────────────────────────────────────────────────────────

def load_replica(in_root: Path, rel_path: str) -> dict:
    """Load all artifacts for one replica."""
    rep_dir = in_root / rel_path
    if not rep_dir.is_dir():
        sys.exit(f"ERROR: {rep_dir} not found (expected analyze_v2 output)")
    cv_path = rep_dir / "cv.csv"
    if not cv_path.exists():
        sys.exit(f"ERROR: {cv_path} missing — analyze_v2 not run for this replica")
    cv = np.loadtxt(cv_path, delimiter=",", skiprows=1)
    rmsf_p = np.loadtxt(rep_dir / "rmsf_protein.csv", delimiter=",", skiprows=1)
    rmsf_r = np.loadtxt(rep_dir / "rmsf_rna.csv",     delimiter=",", skiprows=1)
    cmap   = np.load(rep_dir / "contact_map.npy")
    chi    = np.loadtxt(rep_dir / "chi_dihedral.csv", delimiter=",", skiprows=1)
    pucker = np.loadtxt(rep_dir / "pucker.csv",       delimiter=",", skiprows=1)
    summary = json.loads((rep_dir / "summary.json").read_text())

    # chi/pucker headers (resids)
    with (rep_dir / "chi_dihedral.csv").open() as f:
        head = f.readline().strip().lstrip("#").strip().split(",")
        chi_resids = [int(x) for x in head]

    return {
        "rep_dir": rep_dir,
        "times_ps":     cv[:, 0],
        "prot_rmsd":    cv[:, 1],
        "rna_rmsd":     cv[:, 2],
        "n_contacts":   cv[:, 3],
        "com_dist":     cv[:, 4],
        "g3_contacts":  cv[:, 5],
        "g3_hbonds":    cv[:, 6],
        "rmsf_protein": rmsf_p,        # cols: resid, rmsf
        "rmsf_rna":     rmsf_r,
        "contact_map":  cmap,
        "chi_arr":      chi,            # shape (n_frames, n_RNA_res)
        "pucker_arr":   pucker,
        "chi_resids":   chi_resids,
        "summary":      summary,
    }


def stack_with_padding(arrs: Sequence[np.ndarray]) -> np.ndarray:
    """Stack 1D arrays of possibly different lengths into 2D with NaN."""
    max_len = max(len(a) for a in arrs)
    out = np.full((len(arrs), max_len), np.nan)
    for i, a in enumerate(arrs):
        out[i, :len(a)] = a
    return out


def mean_sem_n(stack: np.ndarray) -> tuple:
    """Returns (mean, sem, n_avail) ignoring NaN per column."""
    n_avail = np.sum(~np.isnan(stack), axis=0)
    with np.errstate(invalid="ignore"):
        m = np.nanmean(stack, axis=0)
        s = np.nanstd(stack, axis=0, ddof=1)
    sem = s / np.sqrt(np.maximum(n_avail, 1))
    return m, sem, n_avail


# ── Statistical tests ──────────────────────────────────────────────────────

def stats_for_pair(values_a: np.ndarray, values_b: np.ndarray) -> dict:
    """Return Welch t-test on means + KS + Mann-Whitney on distributions."""
    out = {"n_a": int(len(values_a)), "n_b": int(len(values_b)),
           "mean_a": float(np.nanmean(values_a)),
           "mean_b": float(np.nanmean(values_b)),
           "mean_diff": float(np.nanmean(values_b) - np.nanmean(values_a))}
    if SCIPY_OK and len(values_a) > 2 and len(values_b) > 2:
        a = values_a[~np.isnan(values_a)]
        b = values_b[~np.isnan(values_b)]
        if len(a) > 2 and len(b) > 2:
            t, p_t = scstats.ttest_ind(a, b, equal_var=False)
            out["welch_t"] = float(t); out["welch_p"] = float(p_t)
            ks, p_ks = scstats.ks_2samp(a, b)
            out["ks"] = float(ks); out["ks_p"] = float(p_ks)
            U, p_mw = scstats.mannwhitneyu(a, b, alternative="two-sided")
            out["mw_U"] = float(U); out["mw_p"] = float(p_mw)
    return out


# ── Plot helpers ────────────────────────────────────────────────────────────

def _setup_mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.size": 10, "axes.labelsize": 11, "axes.titlesize": 12,
        "figure.dpi": 150, "savefig.dpi": 200,
        "axes.spines.top": False, "axes.spines.right": False,
    })
    return plt


COLORS = {"WT": "#2563eb", "8OG": "#dc2626"}


def plot_timeseries(ts_data: dict, fig_dir: Path):
    plt = _setup_mpl()
    metrics = [
        ("prot_rmsd",   "Protein backbone RMSD", "RMSD (Å)"),
        ("rna_rmsd",    "RNA heavy-atom RMSD",   "RMSD (Å)"),
        ("n_contacts",  "Total protein-RNA contacts (<4 Å)",   "contacts"),
        ("g3_contacts", "Modified residue contacts with protein", "contacts"),
        ("com_dist",    "Protein-RNA COM distance",  "distance (Å)"),
        ("g3_hbonds",   "Modified residue polar contacts with protein",
                        "polar contacts"),
    ]
    fig, axes = plt.subplots(3, 2, figsize=(11, 9), sharex=True)
    for ax, (k, title, ylab) in zip(axes.flat, metrics):
        for lab, col in COLORS.items():
            if lab not in ts_data: continue
            d = ts_data[lab]
            t = d["times_ps"] / 1000.0
            ax.plot(t, d[k]["mean"], color=col, lw=1.6,
                    label=f"{lab} (n={int(d[k]['n'].max())})")
            ax.fill_between(t, d[k]["mean"] - d[k]["sem"],
                            d[k]["mean"] + d[k]["sem"], color=col, alpha=0.18)
        ax.set_title(title); ax.set_ylabel(ylab); ax.grid(alpha=0.3)
        ax.legend(loc="best", fontsize=8, frameon=False)
    for ax in axes[-1]:
        ax.set_xlabel("Time (ns)")
    plt.suptitle("4BS2 ver_2: WT vs 8oxo-G3 — time series across replicas\n"
                 "(line = mean, band = SEM)",
                 y=0.995, fontsize=11)
    plt.tight_layout()
    plt.savefig(fig_dir / "timeseries.png", bbox_inches="tight")
    plt.close()


def plot_rmsf(rmsf_data: dict, fig_dir: Path, mod_resid: int):
    plt = _setup_mpl()
    fig, axes = plt.subplots(2, 1, figsize=(12, 7))

    # Protein
    ax = axes[0]
    for lab, col in COLORS.items():
        if lab not in rmsf_data: continue
        d = rmsf_data[lab]
        ax.plot(d["prot_residues"], d["prot_mean"], color=col, lw=1.0,
                label=f"{lab} (n={d['n_rep']})")
        ax.fill_between(d["prot_residues"],
                        d["prot_mean"] - d["prot_sem"],
                        d["prot_mean"] + d["prot_sem"],
                        color=col, alpha=0.18)
    ax.set_title("Per-residue RMSF — protein chain A (mean ± SEM across replicas)")
    ax.set_xlabel("Protein residue"); ax.set_ylabel("RMSF (Å)")
    ax.legend(loc="best", frameon=False); ax.grid(alpha=0.3)

    # RNA
    ax = axes[1]
    for lab, col in COLORS.items():
        if lab not in rmsf_data: continue
        d = rmsf_data[lab]
        ax.errorbar(d["rna_residues"], d["rna_mean"], yerr=d["rna_sem"],
                    fmt="o-", color=col, lw=2, ms=7, capsize=3,
                    label=f"{lab} (n={d['n_rep']})")
    ax.axvspan(mod_resid - 0.5, mod_resid + 0.5, color="gold", alpha=0.25,
               label=f"Modified residue (resid {mod_resid})")
    ax.set_title("Per-residue RMSF — RNA chain B")
    ax.set_xlabel("RNA residue"); ax.set_ylabel("RMSF (Å)")
    ax.legend(loc="best", frameon=False); ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(fig_dir / "rmsf.png", bbox_inches="tight")
    plt.close()


def plot_contact_diff(cmap_data: dict, fig_dir: Path, prot_residues, rna_residues):
    plt = _setup_mpl()
    if not ("WT" in cmap_data and "8OG" in cmap_data):
        return
    cm_wt = cmap_data["WT"]
    cm_og = cmap_data["8OG"]
    diff  = cm_og - cm_wt
    vmax = max(cm_wt.max(), cm_og.max())
    vmax_d = max(abs(diff.min()), abs(diff.max())) or 1e-6

    fig, axes = plt.subplots(1, 3, figsize=(20, 4),
                             gridspec_kw={"width_ratios": [1, 1, 1.05]})
    ext = [prot_residues[0], prot_residues[-1], rna_residues[-1], rna_residues[0]]
    for ax, cm, ttl in zip(axes[:2], [cm_wt, cm_og], ["WT mean", "8OG mean"]):
        im = ax.imshow(cm, aspect="auto", cmap="hot", vmin=0, vmax=vmax,
                       extent=ext)
        ax.set_title(f"{ttl} contact map")
        ax.set_xlabel("Protein residue"); ax.set_ylabel("RNA residue")
        plt.colorbar(im, ax=ax, label="contacts/frame")
    im = axes[2].imshow(diff, aspect="auto", cmap="RdBu_r",
                       vmin=-vmax_d, vmax=vmax_d, extent=ext)
    axes[2].set_title("Δ (8OG − WT) — blue=lost, red=gained")
    axes[2].set_xlabel("Protein residue"); axes[2].set_ylabel("RNA residue")
    plt.colorbar(im, ax=axes[2], label="Δ contacts/frame")
    plt.tight_layout()
    plt.savefig(fig_dir / "contact_diff.png", bbox_inches="tight")
    plt.close()


def plot_chi_pucker(chi_data: dict, pucker_data: dict, fig_dir: Path,
                    mod_resid: int):
    plt = _setup_mpl()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    # χ histogram
    ax = axes[0]
    bins = np.arange(-180, 181, 10)
    for lab, col in COLORS.items():
        if lab not in chi_data: continue
        vals = chi_data[lab]
        vals = vals[~np.isnan(vals)]
        ax.hist(vals, bins=bins, density=True, alpha=0.45, color=col,
                label=f"{lab} (n={len(vals):,} frames-x-replicas)")
    ax.axvspan(-130, -90, color="lightgray", alpha=0.5,
               label="anti band (canonical)")
    ax.axvspan( 30,   90, color="lightgray", alpha=0.3,
               label="syn band")
    ax.set_xlabel("Glycosidic χ (°)"); ax.set_ylabel("density")
    ax.set_title(f"χ distribution at modified residue {mod_resid}")
    ax.legend(loc="best", fontsize=8, frameon=False)

    # Pucker histogram
    ax = axes[1]
    bins = np.arange(0, 361, 10)
    for lab, col in COLORS.items():
        if lab not in pucker_data: continue
        vals = pucker_data[lab]
        vals = vals[~np.isnan(vals)]
        ax.hist(vals, bins=bins, density=True, alpha=0.45, color=col,
                label=lab)
    ax.axvspan(0, 36, color="lightgreen", alpha=0.3,
               label="C3'-endo (north)")
    ax.axvspan(144, 180, color="lightyellow", alpha=0.3,
               label="C2'-endo (south)")
    ax.set_xlabel("Pseudorotation P (°)"); ax.set_ylabel("density")
    ax.set_title(f"Sugar pucker at modified residue {mod_resid}")
    ax.legend(loc="best", fontsize=8, frameon=False)
    plt.tight_layout()
    plt.savefig(fig_dir / "chi_pucker.png", bbox_inches="tight")
    plt.close()


# ── Main ───────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label-replicas", nargs="+", action="append", required=True,
                    help='LABEL replica/path... — repeat per system. '
                         'Example: --label-replicas WT 4BS2_WT/r1 ...')
    ap.add_argument("--mod-resid", type=int, required=True,
                    help="resid of the modified RNA residue (in solvated prmtop)")
    ap.add_argument("--in",  dest="in_root",  type=Path, required=True)
    ap.add_argument("--out", dest="out_root", type=Path, required=True)
    args = ap.parse_args()

    # Parse label-replicas args
    systems: dict = {}
    for chunk in args.label_replicas:
        if len(chunk) < 2:
            sys.exit("each --label-replicas needs LABEL plus >=1 path")
        lab = chunk[0]
        systems[lab] = chunk[1:]

    out_root = args.out_root
    out_root.mkdir(parents=True, exist_ok=True)
    fig_dir = out_root / "figures"; fig_dir.mkdir(exist_ok=True)

    # Load all replicas
    all_data: dict = {}
    for lab, paths in systems.items():
        all_data[lab] = []
        for p in paths:
            r = load_replica(args.in_root, p)
            all_data[lab].append(r)
            print(f"loaded {lab} / {p}")

    # ── Aggregate time series ──────────────────────────────────────────────
    ts_keys = ["prot_rmsd", "rna_rmsd", "n_contacts", "com_dist",
               "g3_contacts", "g3_hbonds"]
    ts_data: dict = {}
    for lab, reps in all_data.items():
        # Time axis: take longest available
        times = max((r["times_ps"] for r in reps), key=len)
        ts_data[lab] = {"times_ps": times}
        for k in ts_keys:
            stack = stack_with_padding([r[k] for r in reps])
            m, sem, n = mean_sem_n(stack)
            ts_data[lab][k] = {"mean": m, "sem": sem, "n": n}

    # ── Aggregate RMSF ─────────────────────────────────────────────────────
    rmsf_data: dict = {}
    for lab, reps in all_data.items():
        prot_residues = reps[0]["rmsf_protein"][:, 0].astype(int)
        rna_residues  = reps[0]["rmsf_rna"][:, 0].astype(int)
        prot_stack = np.stack([r["rmsf_protein"][:, 1] for r in reps])
        rna_stack  = np.stack([r["rmsf_rna"][:, 1]     for r in reps])
        rmsf_data[lab] = {
            "prot_residues": prot_residues, "rna_residues":  rna_residues,
            "prot_mean": prot_stack.mean(axis=0),
            "prot_sem":  prot_stack.std(axis=0, ddof=1) / np.sqrt(len(reps)),
            "rna_mean":  rna_stack.mean(axis=0),
            "rna_sem":   rna_stack.std(axis=0, ddof=1) / np.sqrt(len(reps)),
            "n_rep": len(reps),
        }

    # ── Aggregate contact map ──────────────────────────────────────────────
    cmap_data: dict = {}
    for lab, reps in all_data.items():
        cmap_stack = np.stack([r["contact_map"] for r in reps])
        cmap_data[lab] = cmap_stack.mean(axis=0)

    # ── Modified-residue χ + pucker concatenated across replicas ──────────
    mod_resid = args.mod_resid
    chi_concat: dict = {}; pucker_concat: dict = {}
    for lab, reps in all_data.items():
        chi_vals = []; pkr_vals = []
        for r in reps:
            if mod_resid in r["chi_resids"]:
                col = r["chi_resids"].index(mod_resid)
                chi_vals.append(r["chi_arr"][:, col])
                pkr_vals.append(r["pucker_arr"][:, col])
        if chi_vals:
            chi_concat[lab]    = np.concatenate(chi_vals)
            pucker_concat[lab] = np.concatenate(pkr_vals)

    # ── Statistical tests on key scalar metrics ────────────────────────────
    stats_rows = []
    for k in ts_keys:
        if "WT" in ts_data and "8OG" in ts_data:
            # Aggregate across all (frame, replica) → distribution
            wt_arr = np.concatenate([r[k] for r in all_data["WT"]])
            og_arr = np.concatenate([r[k] for r in all_data["8OG"]])
            row = {"metric": k}
            row.update(stats_for_pair(wt_arr, og_arr))
            stats_rows.append(row)
    if "WT" in chi_concat and "8OG" in chi_concat:
        row = {"metric": f"chi_at_resid_{mod_resid}"}
        row.update(stats_for_pair(chi_concat["WT"], chi_concat["8OG"]))
        stats_rows.append(row)
    if "WT" in pucker_concat and "8OG" in pucker_concat:
        row = {"metric": f"pucker_at_resid_{mod_resid}"}
        row.update(stats_for_pair(pucker_concat["WT"], pucker_concat["8OG"]))
        stats_rows.append(row)

    # Write stats.csv
    if stats_rows:
        keys = sorted({k for r in stats_rows for k in r.keys()})
        with (out_root / "stats.csv").open("w") as f:
            f.write(",".join(keys) + "\n")
            for r in stats_rows:
                f.write(",".join(str(r.get(k, "")) for k in keys) + "\n")
        print(f"Wrote {out_root/'stats.csv'}")

    # ── Plots ──────────────────────────────────────────────────────────────
    plot_timeseries(ts_data, fig_dir)
    plot_rmsf(rmsf_data, fig_dir, mod_resid)
    plot_contact_diff(cmap_data, fig_dir,
                      prot_residues=rmsf_data[next(iter(rmsf_data))]["prot_residues"],
                      rna_residues =rmsf_data[next(iter(rmsf_data))]["rna_residues"])
    if chi_concat or pucker_concat:
        plot_chi_pucker(chi_concat, pucker_concat, fig_dir, mod_resid)
    print(f"Figures: {fig_dir}")

    # ── Save aggregated arrays ─────────────────────────────────────────────
    save_dict = {}
    for lab, dat in ts_data.items():
        save_dict[f"{lab}_times_ps"] = dat["times_ps"]
        for k in ts_keys:
            for stat, arr in dat[k].items():
                save_dict[f"{lab}_{k}_{stat}"] = arr
    np.savez(out_root / "timeseries.npz", **save_dict)
    np.savez(out_root / "rmsf.npz",
             **{f"{lab}_{k}": (np.asarray(v) if not isinstance(v, int) else v)
                 for lab, d in rmsf_data.items() for k, v in d.items()
                 if k != "n_rep"})
    if "WT" in cmap_data and "8OG" in cmap_data:
        np.savez(out_root / "contact_diff.npz",
                 wt_mean=cmap_data["WT"], og_mean=cmap_data["8OG"],
                 diff=cmap_data["8OG"] - cmap_data["WT"])
    np.savez(out_root / "chi_distributions.npz",
             **{lab: arr for lab, arr in chi_concat.items()})
    np.savez(out_root / "pucker_distributions.npz",
             **{lab: arr for lab, arr in pucker_concat.items()})

    # ── Markdown report ────────────────────────────────────────────────────
    write_report(out_root, ts_data, rmsf_data, stats_rows,
                 chi_concat, pucker_concat, mod_resid, all_data)
    print(f"Done. Output: {out_root}")


def write_report(out, ts_data, rmsf_data, stats_rows,
                 chi_concat, pucker_concat, mod_resid, all_data):
    L = []
    L.append("# ver_2 4BS2 / WT vs 8oxo-G3 — comparative MD analysis")
    L.append("")
    L.append("## Replicas analysed")
    for lab, reps in all_data.items():
        L.append(f"### {lab}")
        for r in reps:
            s = r["summary"]
            L.append(f"- `{r['rep_dir'].name}`  "
                     f"frames {s['n_frames_analysed']} / t_max {s['t_max_ns']:.0f} ns")
        L.append("")

    # Time-series scalar means
    L.append("## Time-series scalar metrics (means across all frames+replicas)")
    L.append("")
    L.append("| Metric | WT mean | 8OG mean | Δ |")
    L.append("|------|------:|------:|------:|")
    for r in stats_rows:
        L.append(f"| `{r['metric']}` | {r['mean_a']:.3f} | {r['mean_b']:.3f} | "
                 f"{r['mean_diff']:+.3f} |")
    L.append("")

    # Stats table
    if SCIPY_OK and stats_rows and any("welch_p" in r for r in stats_rows):
        L.append("## Statistical tests (WT vs 8OG)")
        L.append("")
        L.append("| Metric | n_WT | n_8OG | Welch t (p) | KS (p) | MW U (p) |")
        L.append("|------|------:|------:|------:|------:|------:|")
        for r in stats_rows:
            wp = r.get("welch_p", float("nan"))
            kp = r.get("ks_p",   float("nan"))
            mp = r.get("mw_p",   float("nan"))
            L.append(f"| `{r['metric']}` | {r['n_a']} | {r['n_b']} | "
                     f"{r.get('welch_t', 0):.2f} ({wp:.2e}) | "
                     f"{r.get('ks', 0):.3f} ({kp:.2e}) | "
                     f"{r.get('mw_U', 0):.0f} ({mp:.2e}) |")
        L.append("")

    L.append("## Figures")
    L.append("")
    L.append("![Time series](figures/timeseries.png)")
    L.append("")
    L.append("![Per-residue RMSF](figures/rmsf.png)")
    L.append("")
    L.append("![Contact map difference](figures/contact_diff.png)")
    L.append("")
    L.append(f"![χ + sugar pucker at resid {mod_resid}](figures/chi_pucker.png)")
    L.append("")

    L.append("## Methodology notes")
    L.append("")
    L.append("- **Force field**: Amber14 ff14SB (protein) + RNA.OL3 + TIP3P + "
             "modxna (Bergonzo 2024) for 8OG residue. modxna's strip+charge "
             "bug patched + post-build charge rebalancing applied "
             "(see `pipeline/rebalance_lib.py`). Net residue charge of "
             "internal 8OG = -1.000000 e exactly.")
    L.append("- **MD protocol**: minimize 1000 → NVT 100 ps → NPT 100 ps → "
             "Production 100 ns @ 4 fs HMR Langevin (300 K, 1/ps). PME "
             "1.2 nm cutoff, HBonds constraints. n=3 replicas per system, "
             "independent seeds.")
    L.append("- **Alignment**: protein C-alpha to frame 0 in memory before "
             "RNA RMSD (so RMSD captures motion *relative to protein*).")
    L.append("- **PBC**: minimum-image for RNA RMSD, COM distance, and per-"
             "residue position storage.")
    L.append("- **Statistical tests**: Welch's t-test (mean), Kolmogorov-"
             "Smirnov (distribution), Mann-Whitney U (rank). With n=3 "
             "replicas pooled across frames, p-values are over the empirical "
             "distribution of frames; treat as descriptive, not formal "
             "hypothesis tests, given correlated samples.")
    L.append("")

    (out / "report.md").write_text("\n".join(L), encoding="utf-8")
    print(f"Report: {out/'report.md'}")


if __name__ == "__main__":
    main()
