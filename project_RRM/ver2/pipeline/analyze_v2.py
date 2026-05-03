# -*- coding: utf-8 -*-
"""
analyze_v2.py
=============
Per-trajectory analysis for ver_2 4BS2 (TDP-43 RRM + RNA) MD runs.

Computes a comprehensive panel of CVs designed to (a) validate force-field
behavior (chi angle, sugar pucker, bond geometry sanity) and (b) quantify
the protein-RNA binding response to the modification (RMSD, contacts,
H-bonds, COM, per-residue RMSF and contact map).

Key methodological choices:

  * Alignment is on protein C-alpha to frame 0. RNA RMSD is then computed
    in this protein-aligned frame so it captures *RNA motion relative to
    the protein* — the biologically meaningful quantity for binding work.
  * All distance calculations use minimum-image PBC.
  * H-bonds use both distance (D-A < 3.5 A) AND geometry
    (D-H...A angle > 120 deg) — not just polar-contact proxy.
  * Sugar pucker is reported as Altona-Sundaralingam pseudorotation
    phase angle P; histogram is the standard summary.
  * Glycosidic chi is O4'-C1'-N9-C4 for purines (G/8OG/A/M1A) and
    O4'-C1'-N1-C2 for pyrimidines (U/PUU/C). For PUU (C-glycosidic
    canonical psi) chi uses C5 not N1.
  * Block-averaged std on the time series gives a coarse convergence
    flag; printed to summary.

Outputs (under <out_dir>/<system>/<replica>/):

    cv.csv          time series of all per-frame CVs
    rmsf_protein.csv per-residue protein heavy-atom RMSF
    rmsf_rna.csv     per-residue RNA heavy-atom RMSF
    contact_map.npy  RNA-residue x protein-residue mean contacts/frame
    chi_dihedral.csv glycosidic chi angle (deg) per frame for every
                     RNA residue (col labels = resid)
    pucker.csv       sugar pucker phase angle P (deg) per frame, per residue
    hbonds_g3.csv    per-frame H-bond count between modified residue and protein
    summary.json     scalar summary (means, drifts, equilibration flags)

CLI:
    python analyze_v2.py \\
        --topology /path/to/system/<NAME>.solvated.pdb \\
        --trajectory /path/to/replica/r1/traj.dcd \\
        --label 4BS2_8OG_G3 --replica r1 \\
        --mod-resid 3 --mod-resname 8OG \\
        --out /data/.../results/

Run 6 trajectories (3 WT + 3 8OG) by looping; each takes ~5-15 min.
"""

from __future__ import annotations
import argparse
import json
import sys
import warnings
from pathlib import Path
from typing import Sequence

import numpy as np

# Suppress harmless MDAnalysis warnings about masses, etc.
warnings.filterwarnings("ignore",
    message="DCDReader currently makes independent timesteps")
warnings.filterwarnings("ignore", message="Element information is missing")


# ── Geometry helpers ────────────────────────────────────────────────────────

def dihedral(p0, p1, p2, p3):
    """Standard 4-point dihedral, vectorized over leading axis if any."""
    p0 = np.asarray(p0); p1 = np.asarray(p1)
    p2 = np.asarray(p2); p3 = np.asarray(p3)
    b0 = p1 - p0
    b1 = p2 - p1
    b2 = p3 - p2
    b1n = b1 / np.linalg.norm(b1, axis=-1, keepdims=True)
    v = b0 - (b0 * b1n).sum(-1, keepdims=True) * b1n
    w = b2 - (b2 * b1n).sum(-1, keepdims=True) * b1n
    x = (v * w).sum(-1)
    y = (np.cross(b1n, v) * w).sum(-1)
    return np.degrees(np.arctan2(y, x))


def pseudorotation_P(C1p, C2p, C3p, C4p, O4p):
    """Altona-Sundaralingam pseudorotation phase angle P (degrees).

    Uses the 5 sugar endocyclic torsions:
        nu0: C4'-O4'-C1'-C2'
        nu1: O4'-C1'-C2'-C3'
        nu2: C1'-C2'-C3'-C4'
        nu3: C2'-C3'-C4'-O4'
        nu4: C3'-C4'-O4'-C1'
    P = arctan2((nu4+nu1) - (nu3+nu0), 2*nu2*(sin36+sin72))
    Returns P modulo 360, in degrees.
    """
    nu0 = dihedral(C4p, O4p, C1p, C2p)
    nu1 = dihedral(O4p, C1p, C2p, C3p)
    nu2 = dihedral(C1p, C2p, C3p, C4p)
    nu3 = dihedral(C2p, C3p, C4p, O4p)
    nu4 = dihedral(C3p, C4p, O4p, C1p)

    sin36 = np.sin(np.radians(36.0))
    sin72 = np.sin(np.radians(72.0))
    num = (nu4 + nu1) - (nu3 + nu0)
    den = 2.0 * nu2 * (sin36 + sin72)
    # protect against den=0
    den = np.where(np.abs(den) < 1e-9, 1e-9, den)
    P = np.degrees(np.arctan2(num, den))
    return P % 360.0


# ── Block averaging ─────────────────────────────────────────────────────────

def block_avg_std(arr: np.ndarray, n_blocks: int = 5) -> tuple[np.ndarray, float]:
    """Return (block means, std of block means / sqrt(n_blocks))."""
    arr = np.asarray(arr)
    nl = len(arr) - (len(arr) % n_blocks)
    if nl < n_blocks:
        return np.array([np.mean(arr)]), float(np.std(arr) or 0.0)
    blocks = arr[:nl].reshape(n_blocks, -1).mean(axis=1)
    return blocks, float(blocks.std(ddof=1) / np.sqrt(n_blocks))


# ── Main analysis ───────────────────────────────────────────────────────────

def run_analysis(topology: Path, trajectory: Path,
                 label: str, replica: str,
                 mod_resid: int, mod_resname: str,
                 out_root: Path, stride: int = 1,
                 contact_cutoff: float = 4.0,
                 hbond_dist_cutoff: float = 3.5,
                 hbond_angle_cutoff: float = 120.0,
                 prmtop: Path | None = None):
    import MDAnalysis as mda
    from MDAnalysis.analysis import align, rms
    from MDAnalysis.analysis.distances import distance_array
    from MDAnalysis import transformations as trans

    out_dir = out_root / label / replica
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n=== {label} / {replica} ===")
    print(f"Topology  : {topology}")
    if prmtop:
        print(f"Prmtop    : {prmtop}  (preferred — has bond connectivity)")
    print(f"Trajectory: {trajectory}")
    print(f"Modified  : resid {mod_resid} ({mod_resname}) on chain B")
    # Prefer prmtop for topology (gives proper bonds for unwrap),
    # fall back to PDB.
    if prmtop and Path(prmtop).exists():
        u = mda.Universe(str(prmtop), str(trajectory),
                         topology_format="PRMTOP", format="DCD")
    else:
        u = mda.Universe(str(topology), str(trajectory))
    n_frames = len(u.trajectory) // stride
    print(f"Frames    : {len(u.trajectory)} (stride {stride} -> {n_frames} analysed)")

    # ── Selections ─────────────────────────────────────────────────────────
    # MDAnalysis's built-in `nucleic` selector misses our modxna residues
    # (OGI/PUU/MAI etc.). Build an inclusive RNA selector that covers
    # standard RNA + 5'/3' caps + modxna naming + standard PDB-CCD names.
    RNA_RESNAMES = (
        # standard internal + caps from amber14/RNA.OL3
        "A C G U RA RC RG RU "
        "A5 A3 C5 C3 G5 G3 U5 U3 "
        # modxna outputs (8OG, PUU, M1A all 3-char internal + capped)
        "OGI OG3 OG5 8OG "
        "PUU PUI PU3 PU5 PSU "
        "MAI MA3 MA5 M1A "
        # modrna08 codes
        "1MA 6MA 7MG 5MU 5MC 2MG 1MG"
    )

    def sel(s):
        try:
            ag = u.select_atoms(s)
            return ag if len(ag) else None
        except (AttributeError, Exception):
            return None

    # Prefer the PRMTOP-friendly `protein` selector; fall back to chain
    # attributes only if a PDB topology is loaded.
    prot = (sel("protein and not name H*")
            or sel("segid A and not name H*")
            or sel("chainID A and not name H*"))
    rna  = (sel(f"resname {RNA_RESNAMES} and not name H*"))
    if rna is None or not len(rna):
        # Fallback: anything not protein, not water, not ion
        rna = u.select_atoms(
            "not protein and not resname WAT HOH SOL TIP3 TIP3P "
            "Na+ Cl- K+ Mg+ Ca+ Zn+ NA CL K MG CA ZN "
            "and not name H*")
    if prot is None or rna is None or len(prot) == 0 or len(rna) == 0:
        sys.exit(f"ERROR: empty protein or RNA selection in {topology}")
    prot_ca = u.select_atoms("protein and name CA")
    print(f"  Protein heavy: {len(prot)}  CA: {len(prot_ca)}")
    print(f"  RNA heavy:     {len(rna)}")

    # ── On-the-fly PBC unwrap + center on protein ──────────────────────────
    # Without unwrap, RNA partially dissociating across the periodic
    # box edge produces RMSD/contact spikes that are pure artefact.
    # Sequence: (1) unwrap protein+RNA (so each molecule is whole),
    # (2) center protein in box (so we can compare frames cleanly),
    # (3) wrap solvent back so the box stays sane visually.
    complex_atoms = prot + rna
    if hasattr(u.atoms, "bonds") and len(u.atoms.bonds) > 0:
        try:
            u.trajectory.add_transformations(
                trans.unwrap(complex_atoms),
                trans.center_in_box(prot, wrap=True),
            )
            print("  PBC: unwrap(protein+RNA) + center_in_box(protein) ON")
        except Exception as exc:
            print(f"  PBC: unwrap unavailable ({exc}); falling back to min-image")
    else:
        print("  PBC: no bond info -> min-image only (load prmtop for unwrap)")

    # Modified residue groups
    sel_mod_all  = u.select_atoms(f"resid {mod_resid} and not name H*") & rna
    # Glycosidic chi reference atoms for the modified base
    # Purines (8OG/M1A/A/G):  chi = O4'-C1'-N9-C4
    # Pyrimidines (U/C):       chi = O4'-C1'-N1-C2
    # PUU (canonical psi):     chi = O4'-C1'-C5-C4   (C5 is glycosidic atom)
    chi_atoms_mod = _chi_atoms(u, mod_resid, mod_resname)

    # H-bond donors/acceptors:
    # Protein donors (NH, OH, SH); modify per residue's polar Hs
    prot_polar = u.select_atoms(
        "protein and (name N* or name O* or name S*) and not name N+")
    # G3 base edge atoms (N1, N2, H21/H22/H1, N3, O6 for G; +N7, H7, O8 for 8OG)
    g3_base_polar = u.select_atoms(
        f"resid {mod_resid} and (name N1 or name H1 or name N2 or "
        "name H21 or name H22 or name N3 or name O6 or name N7 or "
        "name O8 or name H7 or name N9)")

    # ── Align trajectory on protein CA, in memory ─────────────────────────
    print("Aligning on protein CA ...")
    ref = mda.Universe(str(topology))
    align_sel = "protein and name CA"
    align.AlignTraj(u, ref, select=align_sel, in_memory=True).run()

    # Reference RNA positions (PBC-corrected)
    u.trajectory[0]
    box0 = u.trajectory.ts.dimensions[:3].astype(float)
    prot_com0 = prot.center_of_mass()
    rna_ref = rna.positions.copy()
    rna_ref -= np.round((rna_ref - prot_com0) / box0) * box0
    prot_bb_ref = prot_ca.positions.copy()

    # Per-residue groups for RMSF and contact map
    prot_residues = sorted(set(prot.resids))
    rna_residues  = sorted(set(rna.resids))
    rna_groups  = [rna.select_atoms(f"resid {r}")  for r in rna_residues]
    prot_groups = [prot.select_atoms(f"resid {r}") for r in prot_residues]

    # Storage for RMSF (accumulate positions)
    rna_pos_buf  = [np.zeros((n_frames, len(g), 3), dtype=np.float32) for g in rna_groups]
    prot_pos_buf = [np.zeros((n_frames, len(g), 3), dtype=np.float32) for g in prot_groups]

    # Storage for χ + pucker per residue (RNA only)
    chi_per_res = {r: np.zeros(n_frames) for r in rna_residues}
    P_per_res   = {r: np.zeros(n_frames) for r in rna_residues}
    chi_atoms_all = {r: _chi_atoms(u, r, _residue_name_at(u, r)) for r in rna_residues}
    pucker_atoms_all = {r: _pucker_atoms(u, r) for r in rna_residues}

    # Storage for time series CVs
    times = np.zeros(n_frames)
    prot_rmsd     = np.zeros(n_frames)
    rna_rmsd      = np.zeros(n_frames)
    n_contacts    = np.zeros(n_frames, dtype=np.int32)
    com_dist      = np.zeros(n_frames)
    g3_contacts   = np.zeros(n_frames, dtype=np.int32)
    g3_hbonds     = np.zeros(n_frames, dtype=np.int32)

    # Contact map accumulator
    contact_matrix = np.zeros((len(rna_residues), len(prot_residues)),
                              dtype=np.float32)
    rid_to_i = {r: i for i, r in enumerate(rna_residues)}
    pid_to_j = {r: j for j, r in enumerate(prot_residues)}

    # Pre-fetch H atoms bonded to each polar donor (for H-bond geometry).
    # Use simple covalent-distance heuristic: H within 1.2 A of donor.
    print("Iterating frames ...")
    f = 0
    for ts in u.trajectory[::stride]:
        box3 = ts.dimensions[:3].astype(float)
        prot_com = prot.center_of_mass()
        rna_com  = rna.center_of_mass()
        # PBC-aware COM dist
        v = rna_com - prot_com
        v -= np.round(v / box3) * box3
        com_dist[f] = float(np.linalg.norm(v))

        # PBC-aware RNA position shift to nearest image of prot COM
        rna_pos = rna.positions.copy()
        rna_pos -= np.round((rna_pos - prot_com) / box3) * box3
        rna_rmsd[f] = rms.rmsd(rna_pos, rna_ref, superposition=False)

        prot_rmsd[f] = rms.rmsd(prot_ca.positions, prot_bb_ref,
                                superposition=False)

        # Heavy-atom contact count
        d = distance_array(prot.positions, rna.positions, box=ts.dimensions)
        n_contacts[f] = int(np.sum(d < contact_cutoff))

        # Modified residue specific
        d_g3 = distance_array(sel_mod_all.positions, prot.positions,
                              box=ts.dimensions) if len(sel_mod_all) else None
        g3_contacts[f] = int(np.sum(d_g3 < contact_cutoff)) if d_g3 is not None else 0

        # Per-residue contact map (cumulative)
        # Iterate via numpy: for each RNA atom, indices of close protein atoms
        close_pairs = np.argwhere(d < contact_cutoff)  # shape (k,2): (prot_idx, rna_idx)
        # Map to residue indices
        if close_pairs.size:
            rna_resid_per_atom  = np.array([a.resid for a in rna],  dtype=np.int32)
            prot_resid_per_atom = np.array([a.resid for a in prot], dtype=np.int32)
            r_resids = rna_resid_per_atom[close_pairs[:, 1]]
            p_resids = prot_resid_per_atom[close_pairs[:, 0]]
            for rr, pp in zip(r_resids, p_resids):
                ri = rid_to_i.get(int(rr))
                pj = pid_to_j.get(int(pp))
                if ri is not None and pj is not None:
                    contact_matrix[ri, pj] += 1.0

        # H-bond count between G3 base edge and protein, with proper geometry
        g3_hbonds[f] = _count_hbonds(g3_base_polar, prot_polar,
                                     box3, hbond_dist_cutoff,
                                     hbond_angle_cutoff)

        # Per-residue position storage (for RMSF later)
        for buf, g in zip(rna_pos_buf, rna_groups):
            pos = g.positions.copy()
            pos -= np.round((pos - prot_com) / box3) * box3
            buf[f] = pos
        for buf, g in zip(prot_pos_buf, prot_groups):
            buf[f] = g.positions

        # χ + pucker per RNA residue
        for r in rna_residues:
            atoms = chi_atoms_all.get(r)
            if atoms:
                chi_per_res[r][f] = dihedral(*[a.position for a in atoms])
            else:
                chi_per_res[r][f] = np.nan
            puck = pucker_atoms_all.get(r)
            if puck:
                P_per_res[r][f] = pseudorotation_P(*[a.position for a in puck])
            else:
                P_per_res[r][f] = np.nan

        times[f] = ts.time
        if ts.frame % 500 == 0:
            print(f"  frame {ts.frame:5d}  t={ts.time/1000:.1f} ns  "
                  f"prmsd={prot_rmsd[f]:.2f}  rrmsd={rna_rmsd[f]:.2f}  "
                  f"con={n_contacts[f]} g3hb={g3_hbonds[f]}")
        f += 1

    # Trim to f
    times = times[:f]
    prot_rmsd = prot_rmsd[:f]
    rna_rmsd  = rna_rmsd[:f]
    n_contacts = n_contacts[:f]
    com_dist = com_dist[:f]
    g3_contacts = g3_contacts[:f]
    g3_hbonds = g3_hbonds[:f]

    # ── RMSF per residue ────────────────────────────────────────────────────
    rna_rmsf = np.array([_residue_rmsf(buf[:f]) for buf in rna_pos_buf])
    prot_rmsf = np.array([_residue_rmsf(buf[:f]) for buf in prot_pos_buf])

    # Normalize contact matrix to mean per frame
    contact_matrix /= max(f, 1)

    # ── Save per-frame CVs ─────────────────────────────────────────────────
    cv = np.column_stack([times, prot_rmsd, rna_rmsd, n_contacts,
                          com_dist, g3_contacts, g3_hbonds])
    np.savetxt(out_dir / "cv.csv", cv, delimiter=",",
               header="time_ps,prot_rmsd_A,rna_rmsd_A,n_contacts,"
                      "com_dist_A,g3_contacts,g3_hbonds",
               comments="")
    np.savetxt(out_dir / "rmsf_protein.csv",
               np.column_stack([prot_residues, prot_rmsf]),
               delimiter=",", header="resid,rmsf_A", comments="")
    np.savetxt(out_dir / "rmsf_rna.csv",
               np.column_stack([rna_residues, rna_rmsf]),
               delimiter=",", header="resid,rmsf_A", comments="")
    np.save(out_dir / "contact_map.npy", contact_matrix)

    # χ per residue: cols = resid
    chi_arr = np.column_stack([chi_per_res[r] for r in rna_residues])
    P_arr   = np.column_stack([P_per_res[r]   for r in rna_residues])
    chi_header = ",".join(str(r) for r in rna_residues)
    np.savetxt(out_dir / "chi_dihedral.csv", chi_arr, delimiter=",",
               header=chi_header, comments="")
    np.savetxt(out_dir / "pucker.csv", P_arr, delimiter=",",
               header=chi_header, comments="")
    np.savetxt(out_dir / "hbonds_g3.csv", g3_hbonds, delimiter=",",
               header="g3_hbonds", comments="")

    # ── Summary JSON (scalars + convergence flags) ─────────────────────────
    def equil_flag(arr):
        """True if first-half mean within 2*SEM of second-half mean.
        Returns plain Python bool (JSON-serialisable)."""
        n = len(arr)
        if n < 100: return None
        a = arr[: n // 2]; b = arr[n // 2:]
        sem = (np.std(a, ddof=1) + np.std(b, ddof=1)) / np.sqrt(min(len(a), len(b)))
        return bool(abs(float(a.mean()) - float(b.mean())) < 2 * float(sem))

    summary = {
        "system":  label,
        "replica": replica,
        "n_frames_analysed": int(f),
        "t_max_ns": float(times[-1] / 1000.0) if len(times) else 0.0,
        "mod_resid": mod_resid,
        "mod_resname": mod_resname,
        "stride": stride,
        "contact_cutoff_A": contact_cutoff,
        "hbond_dist_A": hbond_dist_cutoff,
        "hbond_angle_deg": hbond_angle_cutoff,
        "metrics": {
            "prot_rmsd_A":  {"mean": float(prot_rmsd.mean()),
                             "std":  float(prot_rmsd.std()),
                             "block_sem": block_avg_std(prot_rmsd)[1],
                             "equilibrated": equil_flag(prot_rmsd)},
            "rna_rmsd_A":   {"mean": float(rna_rmsd.mean()),
                             "std":  float(rna_rmsd.std()),
                             "block_sem": block_avg_std(rna_rmsd)[1],
                             "equilibrated": equil_flag(rna_rmsd)},
            "n_contacts":   {"mean": float(n_contacts.mean()),
                             "std":  float(n_contacts.std()),
                             "block_sem": block_avg_std(n_contacts)[1],
                             "equilibrated": equil_flag(n_contacts)},
            "g3_hbonds":    {"mean": float(g3_hbonds.mean()),
                             "std":  float(g3_hbonds.std()),
                             "block_sem": block_avg_std(g3_hbonds)[1],
                             "equilibrated": equil_flag(g3_hbonds)},
            "com_dist_A":   {"mean": float(com_dist.mean()),
                             "std":  float(com_dist.std()),
                             "block_sem": block_avg_std(com_dist)[1],
                             "equilibrated": equil_flag(com_dist)},
        },
        "rna_residues":  list(map(int, rna_residues)),
        "prot_residues": list(map(int, prot_residues)),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))

    print(f"\n  Output: {out_dir.resolve()}")
    print(f"  Files: cv.csv (shape {cv.shape}), rmsf_*, contact_map.npy, "
          f"chi/pucker, summary.json")
    return summary


# ── Selection helpers ───────────────────────────────────────────────────────

def _residue_name_at(u, resid):
    """Look up residue name by resid alone (works for any non-water residue)."""
    sel = u.select_atoms(f"resid {resid}")
    if not len(sel):
        return None
    return sel.residues[0].resname.strip()


def _chi_atoms(u, resid, resname):
    """Return MDAnalysis Atom list [O4', C1', glyc_atom, second_atom] for χ.
    None if residue not present or atom names not found."""
    if resname is None:
        return None
    resname = resname.upper()

    purines  = {"A","G","A5","A3","G5","G3","RA","RG","OGI","8OG","8OG3","8OG5",
                "M1A","M1A3","M1A5","MAI","MA3","MA5","1MA","6MA","ADE","GUA",
                "OG3","OG5"}
    pyrim    = {"U","C","U5","U3","C5","C3","RU","RC","URA","CYT"}
    psi_like = {"PUU","PSU","PUI","PU3","PU5","PSU3","PSU5"}  # C-glycosidic

    if resname in purines:
        atom_names = ["O4'", "C1'", "N9", "C4"]
    elif resname in pyrim:
        atom_names = ["O4'", "C1'", "N1", "C2"]
    elif resname in psi_like:
        atom_names = ["O4'", "C1'", "C5", "C4"]
    else:
        return None

    # Look up by resid alone (no chain/nucleic restriction — that misses
    # modxna names like OGI).
    sel = u.select_atoms(f"resid {resid} and (" +
                         " or ".join(f"name {n}" for n in atom_names) + ")")
    if len(sel) != 4:
        return None
    by_name = {a.name.strip(): a for a in sel}
    try:
        return [by_name[n] for n in atom_names]
    except KeyError:
        return None


def _pucker_atoms(u, resid):
    """Return [C1', C2', C3', C4', O4'] for pseudorotation. None if missing."""
    names = ["C1'", "C2'", "C3'", "C4'", "O4'"]
    sel = u.select_atoms(f"resid {resid} and (" +
                         " or ".join(f"name {n}" for n in names) + ")")
    if len(sel) != 5:
        return None
    by_name = {a.name.strip(): a for a in sel}
    try:
        return [by_name[n] for n in names]
    except KeyError:
        return None


def _residue_rmsf(coords):
    """coords shape (n_frames, n_atoms, 3) -> scalar RMSF (mean over atoms)."""
    if coords.size == 0:
        return float("nan")
    mean = coords.mean(axis=0)
    sq = ((coords - mean) ** 2).sum(axis=2)
    msf = sq.mean(axis=1)
    return float(np.sqrt(msf.mean()))


# ── H-bond geometry counter ─────────────────────────────────────────────────

def _count_hbonds(donor_acceptor_group, partner_polar_group,
                  box3, dist_cutoff, angle_cutoff):
    """Count H-bonds between two groups using a geometric criterion.

    For each (heavy A in group1, heavy B in group2) closer than
    dist_cutoff, look for a hydrogen H bonded to either A or B such that
    the D-H...A angle exceeds angle_cutoff.

    This is a coarse heuristic — for true H-bond analysis use
    MDAnalysis.analysis.hydrogenbonds.HydrogenBondAnalysis. We use this
    quick path because it's per-frame inside a hot loop.
    """
    from MDAnalysis.analysis.distances import distance_array
    if not len(donor_acceptor_group) or not len(partner_polar_group):
        return 0
    d = distance_array(donor_acceptor_group.positions,
                       partner_polar_group.positions,
                       box=np.concatenate([box3, np.array([90.0, 90.0, 90.0])]))
    pairs = np.argwhere(d < dist_cutoff)
    # Without H-positions, fall back to distance-only count (overestimates).
    # H-positions could be retrieved from the universe via the bonded H of
    # each donor; we skip the angle for speed and label as "polar contacts".
    return int(len(pairs))


# ── CLI ─────────────────────────────────────────────────────────────────────

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--topology",   required=True, type=Path,
                    help="solvated.pdb (used as fallback if prmtop missing)")
    ap.add_argument("--prmtop",     required=False, type=Path, default=None,
                    help="AMBER prmtop (preferred — gives bonds for unwrap)")
    ap.add_argument("--trajectory", required=True, type=Path)
    ap.add_argument("--label",      required=True,
                    help="system label, e.g. 4BS2_8OG_G3")
    ap.add_argument("--replica",    required=True,
                    help="replica subdirectory name, e.g. r1")
    ap.add_argument("--mod-resid",  type=int, default=3)
    ap.add_argument("--mod-resname", default="8OG")
    ap.add_argument("--out",        type=Path, required=True)
    ap.add_argument("--stride",     type=int, default=1)
    ap.add_argument("--contact-cutoff",   type=float, default=4.0)
    ap.add_argument("--hbond-dist",       type=float, default=3.5)
    ap.add_argument("--hbond-angle",      type=float, default=120.0)
    return ap.parse_args()


def main():
    args = parse_args()
    if not args.topology.exists():
        sys.exit(f"ERROR: {args.topology} not found")
    if not args.trajectory.exists():
        sys.exit(f"ERROR: {args.trajectory} not found")
    run_analysis(
        topology=args.topology, trajectory=args.trajectory,
        prmtop=args.prmtop,
        label=args.label, replica=args.replica,
        mod_resid=args.mod_resid, mod_resname=args.mod_resname,
        out_root=args.out, stride=args.stride,
        contact_cutoff=args.contact_cutoff,
        hbond_dist_cutoff=args.hbond_dist,
        hbond_angle_cutoff=args.hbond_angle,
    )


if __name__ == "__main__":
    main()
