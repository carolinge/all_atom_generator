# -*- coding: utf-8 -*-
"""
Analyze protein-RNA binding stability for 4BS2 trajectories.

CVs computed:
  1. RNA RMSD (after protein alignment)
  2. Protein-RNA heavy-atom contact count (< 4 A)
  3. Protein-RNA COM-COM distance
  4. G3-specific contacts with protein
  5. G3 base hydrogen bonds (key: N7 acceptor in WT, donor in 8OG)
  6. Per-residue contact map (RNA nt x protein residue)

Inputs : pdb, traj.dcd, optional residue id of modified base (default 3)
Outputs: <prefix>_cv.csv, <prefix>_per_res.csv, <prefix>_summary.txt
         <prefix>_*.png plots
"""

import sys, os, argparse, json
from pathlib import Path
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align, contacts
from MDAnalysis.analysis.distances import distance_array


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb",   required=True)
    ap.add_argument("--traj",  required=True)
    ap.add_argument("--out",   required=True, help="output prefix (no ext)")
    ap.add_argument("--mod-resid", type=int, default=3,
                    help="residue id of modified base on chain B (default 3)")
    ap.add_argument("--mod-resname", default="G",
                    help="residue name of modified base (G or 8OG)")
    ap.add_argument("--contact-cutoff", type=float, default=4.0,
                    help="heavy-atom contact cutoff in A")
    ap.add_argument("--stride", type=int, default=1)
    return ap.parse_args()


def main():
    args = parse_args()
    out_prefix = Path(args.out)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    print(f"Loading {args.pdb} + {args.traj} ...")
    u = mda.Universe(args.pdb, args.traj)
    n_frames = len(u.trajectory) // args.stride
    print(f"  Atoms: {len(u.atoms)}  Frames: {len(u.trajectory)} (stride={args.stride})")

    # ── selections ────────────────────────────────────────────────────────────
    # Protein = chain A (heavy atoms only for contacts/RMSD)
    # RNA     = chain B
    prot      = u.select_atoms("segid A and not name H*")
    rna       = u.select_atoms("segid B and not name H*")
    prot_bb   = u.select_atoms("segid A and name CA")
    g3        = u.select_atoms(f"segid B and resid {args.mod_resid} and not name H*")
    g3_base   = u.select_atoms(
        f"segid B and resid {args.mod_resid} and "
        "(name N1 N2 N3 N7 N9 C2 C4 C5 C6 C8 O6 O8)")

    print(f"  Protein heavy: {len(prot)}  RNA heavy: {len(rna)}  "
          f"G{args.mod_resid} heavy: {len(g3)}")

    # Try alternative chain selectors if 'segid' empty
    if len(prot) == 0:
        prot      = u.select_atoms("chainID A and not name H*")
        rna       = u.select_atoms("chainID B and not name H*")
        prot_bb   = u.select_atoms("chainID A and name CA")
        g3        = u.select_atoms(
            f"chainID B and resid {args.mod_resid} and not name H*")
        g3_base   = u.select_atoms(
            f"chainID B and resid {args.mod_resid} and "
            "(name N1 N2 N3 N7 N9 C2 C4 C5 C6 C8 O6 O8)")
        print(f"  (using chainID) Protein heavy: {len(prot)}  RNA heavy: {len(rna)}")

    if len(prot) == 0 or len(rna) == 0:
        print("ERROR: empty selections")
        sys.exit(1)

    # ── alignment to first frame on protein backbone ──────────────────────────
    print("Aligning trajectory on protein CA (frame 0 reference) ...")
    ref = mda.Universe(args.pdb)
    align.AlignTraj(u, ref, select="segid A and name CA" if len(prot_bb) > 0
                    else "chainID A and name CA",
                    in_memory=True).run()

    # ── per-frame computations ────────────────────────────────────────────────
    times, rna_rmsd, n_contacts, com_dist, g3_contacts = [], [], [], [], []
    g3_hbond_counts = []  # H-bonds involving G3 base
    rna_residues = sorted({a.resid for a in rna})
    prot_residues = sorted({a.resid for a in prot})
    contact_matrix = np.zeros((len(rna_residues), len(prot_residues)),
                              dtype=np.float32)

    rid_to_i = {r: i for i, r in enumerate(rna_residues)}
    pid_to_j = {r: j for j, r in enumerate(prot_residues)}

    # RNA reference positions (after alignment)
    rna_ref = rna.positions.copy() if hasattr(rna, "positions") else None
    u.trajectory[0]
    rna_ref = rna.positions.copy()

    print("Iterating frames ...")
    for ts in u.trajectory[::args.stride]:
        # 1. RNA RMSD vs frame 0
        rmsd = float(rms.rmsd(rna.positions, rna_ref, superposition=False))

        # 2. Heavy-atom contacts < cutoff (counted)
        d = distance_array(prot.positions, rna.positions, box=ts.dimensions)
        ncon = int(np.sum(d < args.contact_cutoff))

        # 3. COM-COM distance
        com = float(np.linalg.norm(prot.center_of_mass() - rna.center_of_mass()))

        # 4. G3-specific contacts (G3 heavy x protein heavy)
        d_g3 = distance_array(g3.positions, prot.positions, box=ts.dimensions)
        ng3 = int(np.sum(d_g3 < args.contact_cutoff))

        # 5. G3 base H-bond candidates: heavy donor-acceptor pairs <3.5 A
        #    (geometric check only; counts polar contacts as a proxy for H-bonds)
        if len(g3_base) > 0:
            prot_polar = u.select_atoms(
                "(segid A or chainID A) and (name N* or name O* or name S*) "
                "and not name H*")
            d_hb = distance_array(g3_base.positions, prot_polar.positions,
                                  box=ts.dimensions)
            nhb = int(np.sum(d_hb < 3.5))
        else:
            nhb = 0

        # 6. Per-residue contact map
        # For each RNA residue, count contacts with each protein residue
        # (use already-computed d for protein x RNA matrix)
        for i_rna_atom, atom in enumerate(rna):
            close_prot_idx = np.where(d[:, i_rna_atom] < args.contact_cutoff)[0]
            if len(close_prot_idx) == 0:
                continue
            ri = rid_to_i[atom.resid]
            for k in close_prot_idx:
                pj = pid_to_j[prot[k].resid]
                contact_matrix[ri, pj] += 1.0

        times.append(float(ts.time))
        rna_rmsd.append(rmsd)
        n_contacts.append(ncon)
        com_dist.append(com)
        g3_contacts.append(ng3)
        g3_hbond_counts.append(nhb)

        if ts.frame % 50 == 0:
            print(f"  frame {ts.frame:5d}  t={ts.time/1000:.1f} ns  "
                  f"RMSD={rmsd:.2f} contacts={ncon} G3con={ng3} G3hb={nhb}")

    # Normalize contact matrix to mean per frame
    contact_matrix /= len(times)

    # ── save CSVs ─────────────────────────────────────────────────────────────
    np.savetxt(f"{out_prefix}_cv.csv",
               np.column_stack([times, rna_rmsd, n_contacts,
                                com_dist, g3_contacts, g3_hbond_counts]),
               delimiter=",",
               header="time_ps,rna_rmsd_A,n_contacts,com_dist_A,g3_contacts,g3_hbond_polar",
               comments="")
    np.savetxt(f"{out_prefix}_contact_map.csv", contact_matrix, delimiter=",")
    json.dump({
        "rna_residues":  list(map(int, rna_residues)),
        "prot_residues": list(map(int, prot_residues)),
        "n_frames": len(times),
        "stride": args.stride,
        "mod_resid": args.mod_resid,
        "mod_resname": args.mod_resname,
    }, open(f"{out_prefix}_meta.json", "w"), indent=2)

    # ── summary stats ─────────────────────────────────────────────────────────
    rmsd_arr = np.array(rna_rmsd)
    cont_arr = np.array(n_contacts)
    g3_arr   = np.array(g3_contacts)
    hb_arr   = np.array(g3_hbond_counts)
    com_arr  = np.array(com_dist)

    summary = (
        f"=== {out_prefix.name} ===\n"
        f"Frames analyzed     : {len(times)}\n"
        f"Time range          : {times[0]/1000:.1f} - {times[-1]/1000:.1f} ns\n"
        f"\n"
        f"RNA RMSD (A)        : mean={rmsd_arr.mean():.2f}  "
        f"std={rmsd_arr.std():.2f}  max={rmsd_arr.max():.2f}\n"
        f"Total contacts (<{args.contact_cutoff}A): "
        f"mean={cont_arr.mean():.0f}  std={cont_arr.std():.0f}\n"
        f"COM-COM dist (A)    : mean={com_arr.mean():.2f}  "
        f"std={com_arr.std():.2f}\n"
        f"G{args.mod_resid} contacts        : mean={g3_arr.mean():.1f}  "
        f"std={g3_arr.std():.1f}\n"
        f"G{args.mod_resid} base polar (<3.5A): mean={hb_arr.mean():.1f}  "
        f"std={hb_arr.std():.1f}\n"
        f"\n"
        f"Top contacted protein residues for RNA (mean contacts/frame >= 1):\n"
    )
    # Top contact partners for the modified residue
    if args.mod_resid in rid_to_i:
        mod_row = contact_matrix[rid_to_i[args.mod_resid]]
        order = np.argsort(mod_row)[::-1]
        summary += f"\nTop protein partners of G{args.mod_resid}:\n"
        for k in order[:15]:
            if mod_row[k] < 0.5:
                break
            summary += f"  prot resid {prot_residues[k]:>4d}: {mod_row[k]:.1f} contacts/frame\n"

    Path(f"{out_prefix}_summary.txt").write_text(summary, encoding="utf-8")
    print("\n" + summary)
    print(f"\nWrote: {out_prefix}_cv.csv, _contact_map.csv, _meta.json, _summary.txt")

    # ── plots (matplotlib if available) ───────────────────────────────────────
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        t_ns = np.array(times) / 1000.0
        axes[0, 0].plot(t_ns, rmsd_arr); axes[0, 0].set_title("RNA RMSD")
        axes[0, 0].set_xlabel("Time (ns)"); axes[0, 0].set_ylabel("RMSD (A)")
        axes[0, 1].plot(t_ns, cont_arr); axes[0, 1].set_title("Total contacts")
        axes[0, 1].set_xlabel("Time (ns)"); axes[0, 1].set_ylabel(f"Contacts (<{args.contact_cutoff}A)")
        axes[1, 0].plot(t_ns, g3_arr); axes[1, 0].set_title(f"G{args.mod_resid} contacts")
        axes[1, 0].set_xlabel("Time (ns)"); axes[1, 0].set_ylabel("Contacts")
        axes[1, 1].plot(t_ns, com_arr); axes[1, 1].set_title("COM-COM distance")
        axes[1, 1].set_xlabel("Time (ns)"); axes[1, 1].set_ylabel("Distance (A)")
        plt.tight_layout()
        plt.savefig(f"{out_prefix}_cv.png", dpi=120)
        plt.close()

        fig, ax = plt.subplots(figsize=(14, 4))
        im = ax.imshow(contact_matrix, aspect="auto", cmap="hot",
                       extent=[prot_residues[0], prot_residues[-1],
                               rna_residues[-1], rna_residues[0]])
        ax.set_xlabel("Protein residue (chain A)")
        ax.set_ylabel("RNA residue (chain B)")
        ax.set_title(f"Per-residue contact map (mean contacts/frame)")
        plt.colorbar(im, ax=ax)
        plt.tight_layout()
        plt.savefig(f"{out_prefix}_contact_map.png", dpi=120)
        plt.close()

        print(f"Wrote: {out_prefix}_cv.png, _contact_map.png")
    except ImportError:
        print("matplotlib not available, skipping plots")


if __name__ == "__main__":
    main()
