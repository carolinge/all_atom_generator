# -*- coding: utf-8 -*-
"""
run_md.py  —  ver_2 production MD driver (prmtop route)
=========================================================
Replica-aware OpenMM script. Reads the assembled prmtop + inpcrd from
the parent ``system/`` directory (built by build_amber_system.sh) and
runs the standard 4BS2-class protocol:

    minimize 1000 -> NVT 100 ps -> NPT 100 ps -> Production 100 ns

The replica seed is read from the file name ``replica.txt`` in the CWD
(written by spawn_replicas.sh). Velocities, integrator RNG, and barostat
RNG are all seeded with this value.

Checkpoint every 500_000 steps (~2 ns at 4 fs); on rerun, resumes from
the latest checkpoint without redoing equilibration.

Cluster invocation: see sub.sh in this same directory.
"""

import sys
from pathlib import Path

from openmm import (LangevinMiddleIntegrator, MonteCarloBarostat, Platform,
                    unit)
from openmm.app import (AmberInpcrdFile, AmberPrmtopFile, CheckpointReporter,
                        DCDReporter, HBonds, PME, Simulation,
                        StateDataReporter)

# ── Paths ────────────────────────────────────────────────────────────────
HERE = Path(__file__).parent
SYSTEM_DIR = HERE.parent / "system"   # replicas_v2/<name>/system/

# Locate the prmtop/inpcrd (single solvated pair per system).
# Skip *.dry.prmtop (unsolvated diagnostic from build_free_rna.sh).
prmtop_files = sorted(p for p in SYSTEM_DIR.glob("*.prmtop")
                       if not p.name.endswith(".dry.prmtop"))
inpcrd_files = sorted(p for p in SYSTEM_DIR.glob("*.inpcrd")
                       if not p.name.endswith(".dry.inpcrd"))
if not prmtop_files or not inpcrd_files:
    sys.exit(f"ERROR: no prmtop/inpcrd in {SYSTEM_DIR}. "
             f"Run build_amber_system.sh first.")
PRMTOP = prmtop_files[0]
INPCRD = inpcrd_files[0]
print(f"[setup] prmtop = {PRMTOP.name}")
print(f"[setup] inpcrd = {INPCRD.name}")

# ── Replica seed ─────────────────────────────────────────────────────────
seed_file = HERE / "replica.txt"
if seed_file.exists():
    SEED = int(seed_file.read_text().strip())
else:
    # If running standalone, default to seed 1 — log loudly
    SEED = 1
    print("[setup] WARNING: replica.txt not found — defaulting to seed=1")
print(f"[setup] replica seed = {SEED}")

# ── Run protocol parameters ──────────────────────────────────────────────
# Override production length via env var MD_NS (default 100 ns).
import os as _os
TEMP            = 300 * unit.kelvin
PRESSURE        = 1 * unit.atmosphere
DT              = 4 * unit.femtoseconds
NVT_STEPS       = 25_000        # 100 ps
NPT_STEPS       = 25_000        # 100 ps
_MD_NS          = float(_os.environ.get("MD_NS", "100"))
PROD_STEPS      = int(_MD_NS * 250_000)   # ns -> steps @ 4 fs
REPORT_INTERVAL = 5_000         # 20 ps
CKPT_INTERVAL   = 500_000       # 2 ns
print(f"[setup] production length: {_MD_NS:g} ns ({PROD_STEPS:,} steps)")

# ── Output paths ─────────────────────────────────────────────────────────
TRAJ_DCD   = HERE / "traj.dcd"
NVT_LOG    = HERE / "nvt.log"
NPT_LOG    = HERE / "npt.log"
MD_LOG     = HERE / "md.log"
CHECKPOINT = HERE / "checkpoint.chk"
FINAL_CHK  = HERE / "final.chk"

# ── Build system ─────────────────────────────────────────────────────────
print("[1/5] Loading prmtop/inpcrd ...")
prmtop = AmberPrmtopFile(str(PRMTOP))
inpcrd = AmberInpcrdFile(str(INPCRD))
print(f"      System: {prmtop.topology.getNumAtoms():,} atoms")

# Note: HMR + 4 fs requires hydrogenMass=1.5*amu. tleap's prmtop already
# carries the mass; we override it here so OpenMM redistributes correctly.
print("[2/5] Creating system ...")
system = prmtop.createSystem(
    nonbondedMethod=PME,
    nonbondedCutoff=1.2 * unit.nanometer,
    constraints=HBonds,
    hydrogenMass=1.5 * unit.amu,
)

integrator = LangevinMiddleIntegrator(TEMP, 1 / unit.picosecond, DT)
integrator.setRandomNumberSeed(SEED)

# ── Platform: CUDA > OpenCL > CPU ────────────────────────────────────────
print("[3/5] Platform selection ...")
sim = None
for pname in ("CUDA", "OpenCL", "CPU"):
    try:
        platform = Platform.getPlatformByName(pname)
        props = {"CudaPrecision": "mixed"} if pname == "CUDA" else {}
        sim = Simulation(prmtop.topology, system, integrator, platform, props)
        print(f"      Using platform: {pname}")
        break
    except Exception as exc:
        print(f"      {pname} failed: {exc}")
        # rebuild integrator (it gets consumed by failed Simulation init)
        if pname != "CPU":
            integrator = LangevinMiddleIntegrator(TEMP, 1 / unit.picosecond, DT)
            integrator.setRandomNumberSeed(SEED)
if sim is None:
    sys.exit("ERROR: no usable OpenMM platform.")

# ── Restart logic ────────────────────────────────────────────────────────
restart = CHECKPOINT.exists()
print(f"[4/5] {'RESTART' if restart else 'FRESH'} run")

if restart:
    sim.loadCheckpoint(str(CHECKPOINT))
    # Skip min/equil; jump straight to production reporters
else:
    sim.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        sim.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    print("      Minimizing ...")
    sim.minimizeEnergy(maxIterations=1000)

    print("      NVT equilibration (100 ps) ...")
    sim.context.setVelocitiesToTemperature(TEMP, SEED)
    sim.reporters.append(StateDataReporter(
        str(NVT_LOG), REPORT_INTERVAL,
        step=True, potentialEnergy=True, temperature=True))
    sim.step(NVT_STEPS)
    sim.reporters.clear()

    print("      NPT equilibration (100 ps) ...")
    barostat = MonteCarloBarostat(PRESSURE, TEMP)
    barostat.setRandomNumberSeed(SEED)
    system.addForce(barostat)
    sim.context.reinitialize(preserveState=True)
    sim.reporters.append(StateDataReporter(
        str(NPT_LOG), REPORT_INTERVAL,
        step=True, potentialEnergy=True, volume=True, temperature=True))
    sim.step(NPT_STEPS)
    sim.reporters.clear()

# ── Production ───────────────────────────────────────────────────────────
print(f"[5/5] Production ({PROD_STEPS / 250_000:.0f} ns @ {DT.value_in_unit(unit.femtoseconds):.0f} fs) ...")
sim.reporters.append(DCDReporter(str(TRAJ_DCD), REPORT_INTERVAL, append=restart))
sim.reporters.append(StateDataReporter(
    str(MD_LOG), REPORT_INTERVAL,
    step=True, time=True, potentialEnergy=True, temperature=True,
    progress=True, remainingTime=True, speed=True,
    totalSteps=PROD_STEPS, append=restart))
sim.reporters.append(CheckpointReporter(str(CHECKPOINT), CKPT_INTERVAL))

sim.step(PROD_STEPS)

sim.saveCheckpoint(str(FINAL_CHK))
print("Done.")
