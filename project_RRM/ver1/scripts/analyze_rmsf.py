# -*- coding: utf-8 -*-
"""
Per-residue RMSF analysis for protein-RNA binding stability.

Computes RMSF (root-mean-square fluctuation) of each residue in:
  - chain B (RNA): all 12 nucleotides
  - chain A (protein): all residues, but highlight contact-region residues

Outputs:
  - <out>_rmsf.csv : per-residue RMSF for both chains
  - <out>_rmsf.png : RMSF plot, contact residues marked

Trajectory is aligned on protein C-alpha first, then RMSF is computed
on heavy atoms per residue.
"""

import sys, json, argparse
from pathlib import Path
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb",  required=True)
    ap.add_argument("--traj", required=True)
    ap.add_argument("--out",  required=True, help="output prefix")
    ap.add_argument("--mod-resid", type=int, default=3,
                    help="modified RNA residue id (default 3)")
    ap.add_argument("--stride", type=int, default=1)
    return ap.parse_args()


def per_residue_rmsf(atomgroup_residues):
    """Given a list of AtomGroups (one per residue), return mean-heavy-atom
    RMSF for each residue (Å). The Universe must already be aligned, with
    positions accumulated externally."""
    rmsf_list = []
    for ag in atomgroup_residues:
        if len(ag) == 0:
            rmsf_list.append(np.nan)
            continue
        # ag.positions has shape (n_atoms, 3); we accumulated coords externally
        coords = ag._coord_buffer  # set below
        mean_pos = coords.mean(axis=0)              # (n_atoms, 3)
        sq_disp  = ((coords - mean_pos)**2).sum(axis=2)  # (n_frames, n_atoms)
        # average over atoms within residue, then over frames
        residue_msf = sq_disp.mean(axis=1)          # (n_frames,)
        rmsf = np.sqrt(residue_msf.mean())          # scalar
        rmsf_list.append(rmsf)
    return rmsf_list


def main():
    args = parse_args()
    out_prefix = Path(args.out)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    print(f"Loading {args.pdb} + {args.traj} ...")
    u = mda.Universe(args.pdb, args.traj)
    print(f"  Frames: {len(u.trajectory)}")

    # Try segid then chainID
    prot = u.select_atoms("segid A and not name H*")
    rna  = u.select_atoms("segid B and not name H*")
    if len(prot) == 0:
        prot = u.select_atoms("chainID A and not name H*")
        rna  = u.select_atoms("chainID B and not name H*")
    print(f"  Protein heavy: {len(prot)}  RNA heavy: {len(rna)}")

    # Align on protein CA to first frame (in memory)
    print("Aligning on protein CA ...")
    ref = mda.Universe(args.pdb)
    align_sel = "segid A and name CA" if len(u.select_atoms("segid A and name CA")) > 0 \
                else "chainID A and name CA"
    align.AlignTraj(u, ref, select=align_sel, in_memory=True).run()

    # Iterate frames, accumulate per-residue heavy-atom positions
    prot_residues = sorted({a.resid for a in prot})
    rna_residues  = sorted({a.resid for a in rna})

    # Build residue-grouped AtomGroups
    prot_groups = [prot.select_atoms(f"resid {rid}") for rid in prot_residues]
    rna_groups  = [rna.select_atoms(f"resid {rid}")  for rid in rna_residues]

    # Pre-allocate coordinate buffers
    n_frames = len(u.trajectory) // args.stride
    for ag in prot_groups + rna_groups:
        ag._coord_buffer = np.zeros((n_frames, len(ag), 3), dtype=np.float32)

    print(f"Iterating {n_frames} frames ...")
    f = 0
    for ts in u.trajectory[::args.stride]:
        for ag in prot_groups:
            ag._coord_buffer[f] = ag.positions
        for ag in rna_groups:
            ag._coord_buffer[f] = ag.positions
        f += 1
        if ts.frame % 100 == 0:
            print(f"  frame {ts.frame}")

    # Compute RMSF per residue
    print("Computing per-residue RMSF ...")
    prot_rmsf = per_residue_rmsf(prot_groups)
    rna_rmsf  = per_residue_rmsf(rna_groups)

    # Free memory
    for ag in prot_groups + rna_groups:
        del ag._coord_buffer

    # Save CSV
    np.savetxt(f"{out_prefix}_rmsf_protein.csv",
               np.column_stack([prot_residues, prot_rmsf]),
               delimiter=",", header="resid,rmsf_A", comments="")
    np.savetxt(f"{out_prefix}_rmsf_rna.csv",
               np.column_stack([rna_residues, rna_rmsf]),
               delimiter=",", header="resid,rmsf_A", comments="")

    print(f"\nProtein RMSF:")
    for rid, r in zip(prot_residues, prot_rmsf):
        if r > 1.5:
            print(f"  resid {rid:>4d}: {r:.2f} Å")
    print(f"\nRNA RMSF:")
    for rid, r in zip(rna_residues, rna_rmsf):
        print(f"  resid {rid:>4d}: {r:.2f} Å")

    # Save metadata
    json.dump({
        "prot_residues": list(map(int, prot_residues)),
        "prot_rmsf":     list(map(float, prot_rmsf)),
        "rna_residues":  list(map(int, rna_residues)),
        "rna_rmsf":      list(map(float, rna_rmsf)),
        "n_frames":      n_frames,
        "mod_resid":     args.mod_resid,
    }, open(f"{out_prefix}_rmsf.json", "w"), indent=2)

    print(f"\nWrote: {out_prefix}_rmsf_protein.csv, _rmsf_rna.csv, _rmsf.json")


if __name__ == "__main__":
    main()
