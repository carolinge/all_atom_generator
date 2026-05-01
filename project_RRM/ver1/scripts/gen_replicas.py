# -*- coding: utf-8 -*-
"""Generate 3 replicas each for WT and 8OG (G3) production runs.

Each replica:
  - lives in its own directory under output/ (for staging/inspection)
  - uses a unique random seed (different velocities + ion placement)
  - runs 100 ns of production
  - bundles its own sub.sh ready for sbatch

Deployment to the cluster
-------------------------
After running this script, upload directories to:

    bio:/data/biophys/carolinge/clawork/37_OXR/replicas/

(the same path is also reachable as ~/work/37_OXR/replicas/ via symlink).
The MD output (~3 GB DCD per replica) lives on /data, NOT /home — the
home filesystem has an 80 GB quota and was filled by the first replica
batch in April 2026.
"""
from pathlib import Path
from datetime import datetime
import shutil

PROJECT = Path(__file__).resolve().parent.parent
# Local output for staging (small, easy to inspect locally).
# Heavy MD output (trajectories) goes to /data on the cluster — see sub.sh.
OUT     = PROJECT / "output"
# Cluster destination for replica directories (mounted via symlink ~/work/37_OXR
# on the bio host; absolute path is /data/biophys/carolinge/clawork/37_OXR).
CLUSTER_DEST = "/data/biophys/carolinge/clawork/37_OXR/replicas"
XML_8OG = PROJECT / "app" / "ff_custom" / "8OG_RNA_amber.xml"

SYSTEMS = {
    "WT":  {
        "src_pdb": OUT / "4BS2_WT_MD" / "4BS2_prepared.pdb",
        "extra_xml": [],
        "label": "4BS2 WT",
    },
    "8OG": {
        "src_pdb": OUT / "4BS2_8OG_G3_MD" / "4BS2_8OG_G3.pdb",
        "extra_xml": [str(XML_8OG)],
        "label": "4BS2 8OG (G3)",
    },
}

PROD_NS    = 100                     # production length (ns)
DT_FS      = 4                       # 4 fs HMR
PROD_STEPS = int(PROD_NS * 1000 / (DT_FS / 1000.0))   # = 25_000_000 for 100 ns @ 4 fs
TRAJ_EVERY = 5000                    # frames every 20 ps
LOG_EVERY  = 5000

FF_PROTEIN = ["amber14/protein.ff14SB.xml"]
FF_RNA     = ["amber14/RNA.OL3.xml"]
FF_WATER   = ["amber14/tip3p.xml"]


def make_run_script(pdb_name: str, ff_files: list, water_files: list,
                     extra_xml_names: list, seed: int, prod_steps: int) -> str:
    ff_list = ", ".join(f"'{f}'" for f in ff_files + water_files + extra_xml_names)
    ts = datetime.now().strftime("%Y-%m-%d %H:%M")
    return f"""\
# -*- coding: utf-8 -*-
# Generated {ts}  --  replica seed={seed},  production={prod_steps} steps ({PROD_NS} ns)
import random
random.seed({seed})

from pathlib import Path
from openmm import *
from openmm.app import *
from openmm.unit import *

HERE = Path(__file__).parent

pdb = PDBFile(str(HERE / '{pdb_name}'))
ff  = ForceField({ff_list})

modeller = Modeller(pdb.topology, pdb.positions)
modeller.addSolvent(ff,
    padding=1.0*nanometer,
    ionicStrength=0.15*molar,
    positiveIon='Na+', negativeIon='Cl-')
print(f"System: {{modeller.topology.getNumAtoms()}} atoms")

# Save the solvated topology before simulation so analysis can use it
with open(str(HERE / 'solvated.pdb'), 'w') as f:
    PDBFile.writeFile(modeller.topology, modeller.positions, f)

system = ff.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1.2*nanometer,
    constraints=HBonds,
    hydrogenMass=1.5*amu,
)

integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, {DT_FS}*femtoseconds)
integrator.setRandomNumberSeed({seed})

# Platform with fallback
sim = None
for pname in ('CUDA', 'OpenCL', 'CPU'):
    try:
        platform = Platform.getPlatformByName(pname)
        props = {{'CudaPrecision': 'mixed'}} if pname == 'CUDA' else {{}}
        sim = Simulation(modeller.topology, system, integrator, platform, props)
        print(f"Using platform: {{pname}}")
        break
    except Exception as e:
        print(f"  {{pname}} failed: {{e}}")
        if pname != 'CPU':
            integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, {DT_FS}*femtoseconds)
            integrator.setRandomNumberSeed({seed})
        continue
if sim is None:
    raise RuntimeError("No usable OpenMM platform found")
sim.context.setPositions(modeller.positions)

print("Minimizing...")
sim.minimizeEnergy(maxIterations=1000)

print("NVT equilibration (100 ps)...")
sim.context.setVelocitiesToTemperature(300*kelvin, {seed})
sim.reporters.append(StateDataReporter('nvt.log', 500,
    step=True, potentialEnergy=True, temperature=True))
sim.step(25000)

print("NPT equilibration (100 ps)...")
barostat = MonteCarloBarostat(1*atmosphere, 300*kelvin)
barostat.setRandomNumberSeed({seed})
system.addForce(barostat)
sim.context.reinitialize(preserveState=True)
sim.reporters.append(StateDataReporter('npt.log', 500,
    step=True, potentialEnergy=True, volume=True, temperature=True))
sim.step(25000)

print("Production MD ({PROD_NS} ns)...")
sim.reporters.append(DCDReporter('traj.dcd', {TRAJ_EVERY}))
sim.reporters.append(StateDataReporter('md.log', {LOG_EVERY},
    step=True, time=True, potentialEnergy=True,
    temperature=True, progress=True, remainingTime=True,
    speed=True, totalSteps={prod_steps}))
sim.step({prod_steps})

# Final state checkpoint
sim.saveCheckpoint(str(HERE / 'final.chk'))
print("Done.")
"""


SUB_SH_TEMPLATE = """\
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --partition=graphic
#SBATCH --constraint="gpu"
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=10
#SBATCH --mem=225000
#SBATCH --mail-type=END
#SBATCH --job-name="{jobname}"
#SBATCH --output=./%j.out
#SBATCH --error=./%j.err

set -e
module purge
module load cuda/12.5
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate allatom

python ./run_openmm.py
"""


def main():
    summary = []
    seed = 0
    for sys_name, info in SYSTEMS.items():
        for rep in (1, 2, 3):
            seed += 1   # 1..6 unique seeds across all replicas
            tag = f"4BS2_{sys_name}_r{rep}"
            d   = OUT / tag
            d.mkdir(parents=True, exist_ok=True)

            # 1) Copy starting structure
            pdb_dst = d / info["src_pdb"].name
            shutil.copy(info["src_pdb"], pdb_dst)

            # 2) Copy any extra XML
            extra_names = []
            for xf in info["extra_xml"]:
                xn = Path(xf).name
                shutil.copy(xf, d / xn)
                extra_names.append(xn)

            # 3) Write run_openmm.py
            ff_files = FF_PROTEIN + FF_RNA
            run_py = make_run_script(
                pdb_name=info["src_pdb"].name,
                ff_files=ff_files,
                water_files=FF_WATER,
                extra_xml_names=extra_names,
                seed=seed,
                prod_steps=PROD_STEPS,
            )
            (d / "run_openmm.py").write_text(run_py, encoding="utf-8")

            # 4) Write sub.sh
            sub_sh = SUB_SH_TEMPLATE.format(jobname=tag)
            (d / "sub.sh").write_text(sub_sh, encoding="utf-8")

            # 5) README
            readme = (
                f"{info['label']} — replica {rep}\n"
                f"{'=' * 50}\n"
                f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n"
                f"Seed: {seed}\n"
                f"Production length: {PROD_NS} ns ({PROD_STEPS} steps @ {DT_FS} fs)\n"
                f"Trajectory output: every {TRAJ_EVERY} steps "
                f"({TRAJ_EVERY*DT_FS/1000:.1f} ps)\n\n"
                f"To run on the cluster:\n"
                f"  cd {tag}\n"
                f"  sbatch sub.sh\n"
            )
            (d / "README.txt").write_text(readme, encoding="utf-8")

            summary.append((tag, seed, str(d)))

    print("Generated:")
    for tag, s, p in summary:
        print(f"  {tag}  seed={s}  -> {p}")
    print(f"\nTotal: {len(summary)} replicas, {PROD_NS} ns each "
          f"= {PROD_NS * len(summary)} ns")


if __name__ == "__main__":
    main()
