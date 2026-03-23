"""
all_atom MD Setup — 4-step wizard backend
Step 1: Upload / select PDB
Step 2: Structure preparation (PDBFixer + mutations)
Step 3: Force field + solvation options
Step 4: Generate & download (OpenMM script + GROMACS files)
"""

import json
import queue
import shutil
import threading
import traceback
import zipfile
from datetime import datetime
from pathlib import Path

from flask import Flask, Response, jsonify, render_template, request, send_file

PROJECT_ROOT = Path(__file__).parent.parent
UPLOAD_DIR   = PROJECT_ROOT / "uploaded"
PREPARED_DIR = PROJECT_ROOT / "prepared"
OUTPUT_DIR   = PROJECT_ROOT / "output"

for d in (UPLOAD_DIR, PREPARED_DIR, OUTPUT_DIR):
    d.mkdir(exist_ok=True)

app = Flask(__name__, template_folder="templates", static_folder="static")
app.config["MAX_CONTENT_LENGTH"] = 64 * 1024 * 1024  # 64 MB

# ── Session (single-user local tool) ──────────────────────────────────────
SESSION: dict = {}

# ── Log streaming ─────────────────────────────────────────────────────────
log_queue:  queue.Queue = queue.Queue()
job_running: bool = False

def emit(msg: str, level: str = "info"):
    log_queue.put(json.dumps({"msg": msg, "level": level}))

def clear_queue():
    while not log_queue.empty():
        try: log_queue.get_nowait()
        except queue.Empty: break

# ── Residue classification tables ──────────────────────────────────────────
STANDARD_AA = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY",
    "HIS","HSD","HSE","HSP","HIE","HID","HIP",
    "ILE","LEU","LYS","MET","PHE","PRO",
    "SER","THR","TRP","TYR","VAL","ACE","NME",
}
DNA_RES  = {"DA","DC","DG","DT","DA5","DA3","DC5","DC3","DG5","DG3","DT5","DT3"}
RNA_RES  = {"A","C","G","U","RA","RC","RG","RU","A5","A3","C5","C3","G5","G3","U5","U3"}
WATER    = {"HOH","WAT","TIP","SOL","H2O"}
IONS     = {"NA","CL","K","MG","CA","ZN","FE","MN","CU","LI","RB","CS","BA","SR","F","BR","I","CD"}
KNOWN_PTM = {
    "SEP":"Phosphoserine","TPO":"Phosphothreonine","PTR":"Phosphotyrosine",
    "HYP":"Hydroxyproline","MLY":"N6-methyl-Lys","MSE":"Selenomethionine",
    "CSO":"S-hydroxy-Cys","KCX":"Carbamylated-Lys","OCS":"Cys-sulfinic acid",
}
# Known modified nucleosides (DNA/RNA); NOT counted in protein chain
KNOWN_NA_MOD = {
    "8OG": "8-oxoguanosine (oxidized G)",
    "PSU": "Pseudouridine",
    "M2G": "N2-methylguanosine",
    "M7G": "7-methylguanosine",
    "OMG": "O2'-methylguanosine",
    "5MU": "5-methyluridine (ribothymidine)",
    "OMA": "O2'-methyladenosine",
    "OMC": "O2'-methylcytidine",
    "OMU": "O2'-methyluridine",
    "1MA": "1-methyladenosine",
    "5MC": "5-methylcytidine",
}

# ── Force field tables ─────────────────────────────────────────────────────
FF_PROTEIN = {
    "charmm36m": {"label": "CHARMM36m (recommended)", "eco": "charmm",
                  "files": ["charmm36.xml"]},
    "ff14sb":    {"label": "AMBER ff14SB",             "eco": "amber",
                  "files": ["amber14/protein.ff14SB.xml"]},
}
FF_DNA = {
    "charmm36_na": {"label": "CHARMM36 NA", "eco": "charmm", "files": ["charmm36/na.xml"]},
    "amber_ol15":  {"label": "AMBER OL15",  "eco": "amber",  "files": ["amber14/DNA.OL15.xml"]},
}
_8OG_RNA_XML = PROJECT_ROOT / "app" / "ff_custom" / "8OG_RNA_amber.xml"
FF_RNA = {
    "charmm36_na": {"label": "CHARMM36 NA", "eco": "charmm", "files": ["charmm36/na.xml"]},
    "amber_ol3":   {"label": "AMBER OL3",   "eco": "amber",  "files": ["amber14/RNA.OL3.xml"]},
    **({"amber_ol3_8og": {
        "label": "AMBER OL3 + 8-oxoguanosine",
        "eco":   "amber",
        "files": ["amber14/RNA.OL3.xml"],
        "extra_xml": [str(_8OG_RNA_XML)],
    }} if _8OG_RNA_XML.exists() else {}),
}
FF_WATER = {
    "charmm": {"label": "TIP3P (CHARMM)", "files": ["charmm36/water.xml"]},
    "amber":  {"label": "TIP3P (AMBER)",  "files": ["amber14/tip3p.xml"]},
}

# ── Helpers ────────────────────────────────────────────────────────────────

def find_pdbs() -> list[str]:
    """Only show original input PDBs, not processed outputs."""
    paths = []
    for p in sorted(PROJECT_ROOT.glob("**/*.pdb")):
        rel = str(p.relative_to(PROJECT_ROOT))
        if not any(x in rel for x in ("output", "prepared", "uploaded", "__")):
            paths.append(rel)
    return paths


def _res_type(name: str) -> str:
    if name in STANDARD_AA:   return "protein"
    if name in KNOWN_NA_MOD:  return "ptm"      # modified nucleoside → orange marker
    if name in DNA_RES:       return "dna"
    if name in RNA_RES:       return "rna"
    if name in KNOWN_PTM:     return "ptm"
    if name in WATER:         return "water"
    if name in IONS:          return "ion"
    return "unknown"


def analyze_pdb(pdb_path: Path) -> dict:
    from pdbfixer import PDBFixer
    fixer = PDBFixer(filename=str(pdb_path))
    chains = []
    for chain in fixer.topology.chains():
        residues = list(chain.residues())
        res_list = [
            {"name": r.name, "seqnum": r.id, "type": _res_type(r.name)}
            for r in residues
        ]
        chains.append({
            "id": chain.id,
            "n_residues": len(residues),
            "first_res": residues[0].name if residues else "?",
            "last_res":  residues[-1].name if residues else "?",
            "residues":  res_list,
        })
    fixer.findMissingResidues()
    n_missing = sum(len(v) for v in fixer.missingResidues.values())
    fixer.findMissingAtoms()
    n_missing_atoms = sum(len(v) for v in fixer.missingAtoms.values())
    return {
        "chains": chains,
        "n_missing_residues": n_missing,
        "n_missing_atoms": n_missing_atoms,
        "total_atoms": fixer.topology.getNumAtoms(),
    }


def classify_components(pdb_path: Path) -> dict:
    """Classify every residue in the structure into component types."""
    import openmm.app as mm
    with open(pdb_path) as f:
        pdb = mm.PDBFile(f)

    comp = {
        "protein": {"chains": [], "n_residues": 0},
        "dna":     {"chains": [], "n_residues": 0},
        "rna":     {"chains": [], "n_residues": 0},
        "ptm":     {"items": [], "n_residues": 0},
        "na_mod":  {"items": [], "n_residues": 0},
        "water":   {"n_molecules": 0},
        "ions":    {"counts": {}, "n_total": 0},
        "unknown": {"items": [], "n_residues": 0},
    }

    for chain in pdb.topology.chains():
        residues = list(chain.residues())
        if not residues:
            continue
        chain_protein = chain_dna = chain_rna = 0
        for r in residues:
            n = r.name.strip()
            if n in STANDARD_AA:
                chain_protein += 1
            elif n in KNOWN_NA_MOD:
                # Modified nucleoside: record separately, don't force into protein chain;
                # chain type will be inferred from the surrounding regular residues.
                comp["na_mod"]["items"].append(
                    {"chain": chain.id, "resnum": r.id, "name": n, "desc": KNOWN_NA_MOD[n]}
                )
                comp["na_mod"]["n_residues"] += 1
            elif n in DNA_RES:
                chain_dna += 1
            elif n in RNA_RES:
                chain_rna += 1
            elif n in WATER:
                comp["water"]["n_molecules"] += 1
            elif n in IONS:
                comp["ions"]["counts"][n] = comp["ions"]["counts"].get(n, 0) + 1
                comp["ions"]["n_total"] += 1
            elif n in KNOWN_PTM:
                comp["ptm"]["items"].append(
                    {"chain": chain.id, "resnum": r.id, "name": n, "desc": KNOWN_PTM[n]}
                )
                comp["ptm"]["n_residues"] += 1
                chain_protein += 1  # protein PTMs live in protein chains
            else:
                comp["unknown"]["items"].append(
                    {"chain": chain.id, "resnum": r.id, "name": n}
                )
                comp["unknown"]["n_residues"] += 1

        if chain_protein > 0 and chain.id not in comp["protein"]["chains"]:
            comp["protein"]["chains"].append(chain.id)
            comp["protein"]["n_residues"] += chain_protein
        if chain_dna > 0 and chain.id not in comp["dna"]["chains"]:
            comp["dna"]["chains"].append(chain.id)
            comp["dna"]["n_residues"] += chain_dna
        if chain_rna > 0 and chain.id not in comp["rna"]["chains"]:
            comp["rna"]["chains"].append(chain.id)
            comp["rna"]["n_residues"] += chain_rna

    return comp


def check_ff_compatibility(ff_choices: dict) -> dict:
    errors, warnings = [], []
    ecos = set()
    if ff_choices.get("protein") in FF_PROTEIN:
        ecos.add(FF_PROTEIN[ff_choices["protein"]]["eco"])
    if ff_choices.get("dna") in FF_DNA:
        ecos.add(FF_DNA[ff_choices["dna"]]["eco"])
    if ff_choices.get("rna") in FF_RNA:
        ecos.add(FF_RNA[ff_choices["rna"]]["eco"])
    if len(ecos) > 1:
        errors.append("Mixed ecosystems: CHARMM and AMBER force fields cannot be combined.")
    return {"ok": len(errors) == 0, "errors": errors, "warnings": warnings}


# ── MDP templates ──────────────────────────────────────────────────────────

MDP_EM = """\
; Energy minimization
integrator      = steep
emtol           = 1000.0
emstep          = 0.01
nsteps          = 50000
nstlist         = 1
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.2
rvdw            = 1.2
pbc             = xyz
"""

MDP_NVT = """\
; NVT equilibration  100 ps
define              = -DPOSRES
integrator          = md
dt                  = 0.002
nsteps              = 50000
nstlog              = 500
nstenergy           = 500
nstxout-compressed  = 500
cutoff-scheme       = Verlet
nstlist             = 10
coulombtype         = PME
rcoulomb            = 1.2
rvdw                = 1.2
pbc                 = xyz
tcoupl              = V-rescale
tc-grps             = Protein Non-Protein
tau_t               = 0.1    0.1
ref_t               = 300    300
pcoupl              = no
constraints         = h-bonds
constraint_algorithm= LINCS
gen_vel             = yes
gen_temp            = 300
"""

MDP_NPT = """\
; NPT equilibration  100 ps
define              = -DPOSRES
integrator          = md
dt                  = 0.002
nsteps              = 50000
nstlog              = 500
nstenergy           = 500
nstxout-compressed  = 500
cutoff-scheme       = Verlet
nstlist             = 10
coulombtype         = PME
rcoulomb            = 1.2
rvdw                = 1.2
pbc                 = xyz
tcoupl              = V-rescale
tc-grps             = Protein Non-Protein
tau_t               = 0.1    0.1
ref_t               = 300    300
pcoupl              = C-rescale
pcoupltype          = isotropic
tau_p               = 2.0
ref_p               = 1.0
compressibility     = 4.5e-5
constraints         = h-bonds
constraint_algorithm= LINCS
gen_vel             = no
"""

MDP_MD = """\
; Production MD  10 ns
integrator          = md
dt                  = 0.002
nsteps              = 5000000
nstlog              = 5000
nstenergy           = 5000
nstxout-compressed  = 5000
cutoff-scheme       = Verlet
nstlist             = 10
coulombtype         = PME
rcoulomb            = 1.2
rvdw                = 1.2
pbc                 = xyz
tcoupl              = V-rescale
tc-grps             = Protein Non-Protein
tau_t               = 0.1    0.1
ref_t               = 300    300
pcoupl              = C-rescale
pcoupltype          = isotropic
tau_p               = 2.0
ref_p               = 1.0
compressibility     = 4.5e-5
constraints         = h-bonds
constraint_algorithm= LINCS
"""


def make_openmm_script(pdb_filename: str, ff_files: list[str],
                        water_files: list[str], box: dict,
                        extra_xml: list[str] | None = None) -> str:
    extra_xml = extra_xml or []
    # Extra XML files are copied next to the script; reference by filename only.
    extra_xml_names = [Path(p).name for p in extra_xml]
    ff_list = ", ".join(f"'{f}'" for f in ff_files + water_files + extra_xml_names)
    padding = box.get("padding", 1.0)
    ionic   = box.get("ionic_strength", 0.15)
    pos_ion = box.get("pos_ion", "Na+")
    neg_ion = box.get("neg_ion", "Cl-")
    ts      = datetime.now().strftime("%Y-%m-%d %H:%M")
    all_ff_labels = ff_files + water_files + extra_xml_names
    # NOTE: all comments must stay ASCII-only to avoid encoding issues
    return f"""\
# -*- coding: utf-8 -*-
# Generated by all_atom MD Setup  {ts}
# Force fields : {', '.join(all_ff_labels)}
# IMPORTANT    : ForceField covers standard protein/DNA/RNA + water only.
#                Small molecule ligands require a separate parameterized XML
#                (GAFF2 or CGenFF) passed as an extra argument to ForceField().
from pathlib import Path
from openmm import *
from openmm.app import *
from openmm.unit import *

HERE = Path(__file__).parent

# --- Load prepared structure ---
pdb = PDBFile(str(HERE / '{pdb_filename}'))
ff  = ForceField({ff_list})

# --- Solvation ---
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addSolvent(ff,
    padding={padding}*nanometer,
    ionicStrength={ionic}*molar,
    positiveIon='{pos_ion}', negativeIon='{neg_ion}')
print(f"System: {{modeller.topology.getNumAtoms()}} atoms")

# --- Create system (HMR enables 4 fs timestep) ---
system = ff.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.2*nanometer,
    constraints=HBonds,
    hydrogenMass=1.5*amu,
)

# --- Integrator: Langevin, 300 K, 4 fs ---
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 4*femtoseconds)

# --- Platform: CUDA > OpenCL > CPU ---
for pname in ('CUDA', 'OpenCL', 'CPU'):
    try:
        platform = Platform.getPlatformByName(pname)
        print(f"Using platform: {{pname}}")
        break
    except Exception:
        continue

props = {{'CudaPrecision': 'mixed'}} if pname == 'CUDA' else {{}}
sim = Simulation(modeller.topology, system, integrator, platform, props)
sim.context.setPositions(modeller.positions)

# --- Energy minimization ---
print("Minimizing...")
sim.minimizeEnergy(maxIterations=1000)

# --- NVT equilibration 100 ps ---
print("NVT equilibration (100 ps)...")
sim.context.setVelocitiesToTemperature(300*kelvin)
sim.reporters.append(StateDataReporter('nvt.log', 500,
    step=True, potentialEnergy=True, temperature=True))
sim.step(25000)

# --- NPT equilibration 100 ps ---
print("NPT equilibration (100 ps)...")
system.addForce(MonteCarloBarostat(1*atmosphere, 300*kelvin))
sim.context.reinitialize(preserveState=True)
sim.reporters.append(StateDataReporter('npt.log', 500,
    step=True, potentialEnergy=True, volume=True, temperature=True))
sim.step(25000)

# --- Production 10 ns ---
print("Production MD (10 ns)...")
sim.reporters.append(DCDReporter('traj.dcd', 5000))
sim.reporters.append(StateDataReporter('md.log', 5000,
    step=True, time=True, potentialEnergy=True,
    temperature=True, progress=True, remainingTime=True,
    speed=True, totalSteps=2500000))
sim.step(2500000)
print("Done.")
"""


def _collect_ff_files(ff_choices: dict) -> tuple[list, list, list, str]:
    """Return (ff_files, extra_xml_files, water_files, ecosystem).

    extra_xml_files are absolute-path XML files (e.g. custom 8OG parameters)
    that must be passed to ForceField() in addition to the standard ff_files.
    """
    ff_files, extra_xml, eco = [], [], None
    if ff_choices.get("protein") in FF_PROTEIN:
        entry = FF_PROTEIN[ff_choices["protein"]]
        ff_files += entry["files"]; eco = entry["eco"]
    if ff_choices.get("dna") in FF_DNA:
        ff_files += FF_DNA[ff_choices["dna"]]["files"]
    if ff_choices.get("rna") in FF_RNA:
        entry = FF_RNA[ff_choices["rna"]]
        ff_files += entry["files"]
        extra_xml += entry.get("extra_xml", [])
    seen = set()
    ff_files = [f for f in ff_files if not (f in seen or seen.add(f))]
    water_files = FF_WATER.get(eco or "charmm", FF_WATER["charmm"])["files"]
    return ff_files, extra_xml, water_files, eco or "charmm"


def generate_files_job(prepared_path: Path, ff_choices: dict, box: dict,
                        out_dir: Path, engine: str):
    global job_running
    job_running = True
    try:
        ff_files, extra_xml, water_files, eco = _collect_ff_files(ff_choices)
        out_dir.mkdir(parents=True, exist_ok=True)
        import shutil as _sh

        if engine == "openmm":
            # Fast path: just write the script, no system build needed here.
            # Solvation happens when the script is RUN on the server.
            emit("Writing OpenMM simulation script ...")
            script = make_openmm_script(
                pdb_filename=prepared_path.name,
                ff_files=ff_files, water_files=water_files, box=box,
                extra_xml=extra_xml,
            )
            (out_dir / "run_openmm.py").write_text(script, encoding="utf-8")
            _sh.copy(prepared_path, out_dir / prepared_path.name)
            for xf in extra_xml:
                _sh.copy(xf, out_dir / Path(xf).name)
            emit("  run_openmm.py — OK", "ok")

            all_ff_labels = ff_files + water_files + [Path(x).name for x in extra_xml]
            readme = (
                f"all_atom MD Setup — {datetime.now().strftime('%Y-%m-%d %H:%M')}\n"
                f"Engine    : OpenMM\n"
                f"Structure : {prepared_path.name}\n"
                f"FF        : {', '.join(all_ff_labels)}\n\n"
                f"Run on server:\n"
                f"  conda activate allatom\n"
                f"  python run_openmm.py\n\n"
                f"NOTE: Small molecule ligands are NOT supported by this script.\n"
                f"      They require separate parameterization (GAFF2 or CGenFF).\n"
            )
            (out_dir / "README.txt").write_text(readme, encoding="utf-8")

        elif engine == "gromacs":
            import openmm.app as mm
            from openmm import unit

            emit("Building solvated system for GROMACS export ...")
            emit("  (30-120 s depending on structure size)", "warn")
            pdb     = mm.PDBFile(str(prepared_path))
            ff      = mm.ForceField(*(ff_files + water_files + extra_xml))
            padding = float(box.get("padding", 1.0))
            ionic   = float(box.get("ionic_strength", 0.15))
            pos_ion = box.get("pos_ion", "Na+")
            neg_ion = box.get("neg_ion", "Cl-")

            modeller = mm.Modeller(pdb.topology, pdb.positions)
            emit(f"  Adding solvent (pad={padding} nm, ionic={ionic} M) ...")
            modeller.addSolvent(ff,
                padding=padding * unit.nanometer,
                ionicStrength=ionic * unit.molar,
                positiveIon=pos_ion, negativeIon=neg_ion)
            emit(f"  Solvated: {modeller.topology.getNumAtoms():,} atoms", "ok")

            emit("Creating OpenMM system object ...")
            system = ff.createSystem(
                modeller.topology,
                nonbondedMethod=mm.PME,
                nonbondedCutoff=1.2 * unit.nanometer,
                constraints=mm.HBonds,
                hydrogenMass=1.5 * unit.amu,
            )
            emit("  System created", "ok")

            emit("Writing GROMACS MDP files ...")
            for name, content in [("em.mdp", MDP_EM), ("nvt.mdp", MDP_NVT),
                                   ("npt.mdp", MDP_NPT), ("md.mdp", MDP_MD)]:
                (out_dir / name).write_text(content, encoding="utf-8")
            emit("  em/nvt/npt/md.mdp — OK", "ok")

            emit("Exporting topology via ParmEd ...")
            emit("  Verify topol.top before production runs", "warn")
            try:
                import parmed as pmd
                structure = pmd.openmm.load_topology(
                    modeller.topology, system, xyz=modeller.positions)
                structure.save(str(out_dir / "topol.top"), overwrite=True)
                structure.save(str(out_dir / "system.gro"), overwrite=True)
                emit("  topol.top  system.gro — OK", "ok")
            except Exception:
                emit("  ParmEd export failed:", "error")
                emit(traceback.format_exc(), "error")

            readme = (
                f"all_atom MD Setup — {datetime.now().strftime('%Y-%m-%d %H:%M')}\n"
                f"Engine    : GROMACS\n"
                f"Structure : {prepared_path.name}\n"
                f"FF        : {', '.join(ff_files + water_files + [Path(x).name for x in extra_xml])}\n\n"
                f"Run on server:\n"
                f"  gmx grompp -f em.mdp  -c system.gro -p topol.top -o em.tpr\n"
                f"  gmx mdrun  -v -deffnm em\n"
                f"  gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr -r em.gro\n"
                f"  gmx mdrun  -v -deffnm nvt\n"
                f"  gmx grompp -f npt.mdp -c nvt.gro -p topol.top -o npt.tpr -r nvt.gro\n"
                f"  gmx mdrun  -v -deffnm npt\n"
                f"  gmx grompp -f md.mdp  -c npt.gro -p topol.top -o md.tpr\n"
                f"  gmx mdrun  -v -deffnm md\n"
            )
            (out_dir / "README.txt").write_text(readme, encoding="utf-8")

        # --- ZIP ---
        emit("Packaging ZIP ...")
        zip_path = out_dir / "md_files.zip"
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
            for f in out_dir.iterdir():
                if f.suffix != ".zip":
                    zf.write(f, f.name)
        emit(f"  md_files.zip ({zip_path.stat().st_size // 1024} KB)", "ok")
        emit("DONE", "done")
        SESSION["zip_path"] = str(zip_path.relative_to(PROJECT_ROOT))

    except Exception:
        emit(f"ERROR:\n{traceback.format_exc()}", "error")
        emit("DONE", "done")
    finally:
        job_running = False


def run_pdbfixer_job(pdb_path: Path, chains: list, options: dict,
                     mutations: list, out_path: Path):
    global job_running
    job_running = True
    try:
        from pdbfixer import PDBFixer
        import openmm.app as mm

        emit(f"Loading {pdb_path.name} ...")
        fixer = PDBFixer(filename=str(pdb_path))

        # Chain selection
        all_ids   = [c.id for c in fixer.topology.chains()]
        to_remove = [c for c in all_ids if c not in chains]
        if to_remove:
            emit(f"Removing chains: {to_remove}")
            fixer.removeChains(chainIds=to_remove)
        emit(f"Keeping chains: {chains}", "ok")

        # Mutations
        if mutations:
            emit(f"Applying {len(mutations)} mutation(s) ...")
            by_chain: dict[str, list] = {}
            for m in mutations:
                by_chain.setdefault(m["chain"], []).append(
                    f"{m['from_res']}-{m['resnum']}-{m['to_res']}"
                )
            for chain_id, mut_list in by_chain.items():
                fixer.applyMutations(mut_list, chain_id)
            emit(f"  Done", "ok")

        # Missing residues
        if options.get("add_missing_residues"):
            emit("Finding missing residues ...")
            fixer.findMissingResidues()
            n = sum(len(v) for v in fixer.missingResidues.values())
            emit(f"  → {n} missing residue(s)", "warn" if n else "ok")
            if n:
                fixer.addMissingAtoms()

        # Heterogens
        if options.get("remove_heterogens"):
            keep = options.get("keep_water", False)
            emit("Removing heterogens" + (" (keeping water)" if keep else "") + " ...")
            fixer.removeHeterogens(keepWater=keep)
            emit("  Done", "ok")

        # Missing heavy atoms
        if options.get("add_missing_atoms"):
            emit("Adding missing heavy atoms ...")
            fixer.findMissingAtoms()
            n_atoms = sum(len(v) for v in fixer.missingAtoms.values())
            n_term  = len(fixer.missingTerminals)
            emit(f"  → {n_atoms} atom(s), {n_term} terminal(s)",
                 "warn" if n_atoms else "ok")
            fixer.addMissingAtoms()
            emit("  Done", "ok")

        # Hydrogens
        if options.get("add_hydrogens"):
            ph = float(options.get("ph", 7.0))
            emit(f"Adding hydrogens (pH {ph}) ...")
            fixer.addMissingHydrogens(pH=ph)
            emit("  Done", "ok")

        # Write
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w") as f:
            mm.PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
        n_out = fixer.topology.getNumAtoms()
        emit(f"Saved → {out_path.relative_to(PROJECT_ROOT)}  ({n_out:,} atoms)", "ok")
        SESSION["prepared_path"] = str(out_path.relative_to(PROJECT_ROOT))
        emit("DONE", "done")

    except Exception:
        emit(f"ERROR:\n{traceback.format_exc()}", "error")
        emit("DONE", "done")
    finally:
        job_running = False


# ── Routes ─────────────────────────────────────────────────────────────────

@app.route("/")
def index():
    return render_template("index.html")


@app.route("/api/pdbs")
def api_pdbs():
    return jsonify(find_pdbs())


@app.route("/api/upload", methods=["POST"])
def api_upload():
    f = request.files.get("file")
    if not f:
        return jsonify({"error": "No file received"}), 400
    ext = Path(f.filename).suffix.lower()
    if ext not in (".pdb", ".cif", ".mmcif"):
        return jsonify({"error": "Only .pdb or .cif files accepted"}), 400

    if ext == ".pdb":
        dest = UPLOAD_DIR / f.filename
        f.save(dest)
        return jsonify({"path": str(dest.relative_to(PROJECT_ROOT))})

    # CIF → PDB conversion via PDBFixer
    import tempfile, os as _os
    tmp = tempfile.NamedTemporaryFile(suffix=ext, delete=False)
    tmp_path = Path(tmp.name)
    try:
        f.save(tmp_path)
        from pdbfixer import PDBFixer
        import openmm.app as mm
        fixer = PDBFixer(filename=str(tmp_path))
        pdb_name = Path(f.filename).stem + "_from_cif.pdb"
        dest = UPLOAD_DIR / pdb_name
        with open(dest, "w") as out:
            mm.PDBFile.writeFile(fixer.topology, fixer.positions, out, keepIds=True)
        return jsonify({
            "path": str(dest.relative_to(PROJECT_ROOT)),
            "converted_from": f.filename,
        })
    except Exception:
        return jsonify({"error": traceback.format_exc()}), 500
    finally:
        try: tmp_path.unlink()
        except: pass


@app.route("/api/analyze", methods=["POST"])
def api_analyze():
    data = request.json
    pdb_path = PROJECT_ROOT / data["pdb"]
    if not pdb_path.exists():
        return jsonify({"error": "File not found"}), 404
    try:
        return jsonify(analyze_pdb(pdb_path))
    except Exception:
        return jsonify({"error": traceback.format_exc()}), 500


@app.route("/api/run_pdbfixer", methods=["POST"])
def api_run_pdbfixer():
    global job_running
    if job_running:
        return jsonify({"error": "A job is already running"}), 409
    data      = request.json
    pdb_path  = PROJECT_ROOT / data["pdb"]
    chains    = data["chains"]
    options   = data["options"]
    mutations = data.get("mutations", [])
    suffix    = "_".join(chains).lower()
    out_path  = PREPARED_DIR / f"{pdb_path.stem}_{suffix}_fixed.pdb"
    clear_queue()
    threading.Thread(
        target=run_pdbfixer_job,
        args=(pdb_path, chains, options, mutations, out_path),
        daemon=True,
    ).start()
    return jsonify({"status": "started", "output": str(out_path.relative_to(PROJECT_ROOT))})


def _run_mutations_job(pdb_path: Path, mutations: list, out_path: Path):
    global job_running
    job_running = True
    try:
        from pdbfixer import PDBFixer
        import openmm.app as mm
        emit(f"Loading {pdb_path.name} ...")
        fixer = PDBFixer(filename=str(pdb_path))
        emit(f"Applying {len(mutations)} mutation(s) ...")
        by_chain: dict = {}
        for m in mutations:
            by_chain.setdefault(m["chain"], []).append(
                f"{m['from_res']}-{m['resnum']}-{m['to_res']}"
            )
        for chain_id, mut_list in by_chain.items():
            fixer.applyMutations(mut_list, chain_id)
            emit(f"  Chain {chain_id}: {', '.join(mut_list)}", "ok")
        # findMissingResidues must be called before findMissingAtoms
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w") as f:
            mm.PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
        emit(f"Saved → {out_path.relative_to(PROJECT_ROOT)}  ({fixer.topology.getNumAtoms():,} atoms)", "ok")
        SESSION["prepared_path"] = str(out_path.relative_to(PROJECT_ROOT))
        emit("DONE", "done")
    except Exception:
        emit(traceback.format_exc(), "error")
        emit("DONE", "done")
    finally:
        job_running = False


@app.route("/api/run_mutations", methods=["POST"])
def api_run_mutations():
    global job_running
    if job_running:
        return jsonify({"error": "A job is already running"}), 409
    data      = request.json
    pdb_path  = PROJECT_ROOT / data["pdb"]
    mutations = data["mutations"]
    out_path  = PREPARED_DIR / f"{pdb_path.stem}_mut.pdb"
    clear_queue()
    threading.Thread(
        target=_run_mutations_job,
        args=(pdb_path, mutations, out_path),
        daemon=True,
    ).start()
    return jsonify({"status": "started", "output": str(out_path.relative_to(PROJECT_ROOT))})


@app.route("/api/classify", methods=["POST"])
def api_classify():
    data = request.json
    pdb_path = PROJECT_ROOT / data["pdb"]
    if not pdb_path.exists():
        return jsonify({"error": "File not found"}), 404
    try:
        comp = classify_components(pdb_path)
        # Return FF options alongside
        return jsonify({
            "components": comp,
            "ff_options": {
                "protein": {k: v["label"] for k, v in FF_PROTEIN.items()},
                "dna":     {k: v["label"] for k, v in FF_DNA.items()},
                "rna":     {k: v["label"] for k, v in FF_RNA.items()},
            },
        })
    except Exception:
        return jsonify({"error": traceback.format_exc()}), 500


@app.route("/api/check_ff", methods=["POST"])
def api_check_ff():
    return jsonify(check_ff_compatibility(request.json))


@app.route("/api/generate", methods=["POST"])
def api_generate():
    global job_running
    if job_running:
        return jsonify({"error": "A job is already running"}), 409
    data     = request.json
    pdb_path = PROJECT_ROOT / data["pdb"]
    ff_ch    = data["ff_choices"]
    box      = data["box"]
    engine   = data.get("engine", "openmm")   # "openmm" | "gromacs"
    ts       = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir  = OUTPUT_DIR / f"{engine}_{pdb_path.stem}_{ts}"
    clear_queue()
    threading.Thread(
        target=generate_files_job,
        args=(pdb_path, ff_ch, box, out_dir, engine),
        daemon=True,
    ).start()
    return jsonify({"status": "started", "engine": engine})


@app.route("/api/download")
def api_download():
    zip_rel = SESSION.get("zip_path")
    if not zip_rel:
        return jsonify({"error": "No file ready"}), 404
    zip_path = PROJECT_ROOT / zip_rel
    if not zip_path.exists():
        return jsonify({"error": "File not found"}), 404
    return send_file(zip_path, as_attachment=True)


@app.route("/api/log_stream")
def api_log_stream():
    def generate():
        while True:
            try:
                msg  = log_queue.get(timeout=30)
                yield f"data: {msg}\n\n"
                data = json.loads(msg)
                if data.get("level") == "done":
                    break
            except queue.Empty:
                yield 'data: {"msg":"...","level":"ping"}\n\n'
    return Response(generate(), mimetype="text/event-stream",
                    headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"})


if __name__ == "__main__":
    print(f"Project root : {PROJECT_ROOT}")
    print("Open browser : http://localhost:5000")
    app.run(debug=False, port=5000, threaded=True)
