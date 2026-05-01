# -*- coding: utf-8 -*-
"""
OpenMM MD Pipeline — 2RRM
结构    : fold_2026_03_23_21_01_model_4_from_cif.pdb
力场    : CHARMM36 + TIP3P water
流程    : 能量最小化 → NVT 100 ps → NPT 100 ps → Production 10 ns
用法    : conda activate allatom && python run_openmm.py
断点续跑: 若 output/ 中存在 checkpoint.xml，自动从断点恢复 Production
"""

import time
from pathlib import Path
from datetime import timedelta

from openmm import *
from openmm.app import *
from openmm.unit import *

# ── 路径 ──────────────────────────────────────────────────────────────────────
HERE   = Path(__file__).parent
PDB    = HERE / "fold_2026_03_23_21_01_model_4_from_cif.pdb"
OUTDIR = HERE / "output"
OUTDIR.mkdir(exist_ok=True)

CHECKPOINT    = OUTDIR / "checkpoint.xml"
TRAJ_DCD      = OUTDIR / "traj.dcd"
SOLVATED_PDB  = OUTDIR / "solvated.pdb"   # 溶剂化结构，供后处理用
LOG_NVT       = OUTDIR / "nvt.log"
LOG_NPT       = OUTDIR / "npt.log"
LOG_PROD      = OUTDIR / "production.log"
FINAL_PDB     = OUTDIR / "final_structure.pdb"

# ── 参数 ──────────────────────────────────────────────────────────────────────
TEMP            = 300 * kelvin
PRESSURE        = 1 * atmosphere
NVT_STEPS       = 25_000        # 100 ps  (4 fs step)
NPT_STEPS       = 25_000        # 100 ps
PROD_STEPS      = 2_500_000     # 10 ns
REPORT_INTERVAL = 5_000         # 每 20 ps 记录一次
CKPT_INTERVAL   = 50_000        # 每 200 ps 保存断点

# ── 后处理：残基编号修正 + wrap + 可视化轨迹 ──────────────────────────────────
def read_original_resids(pdb_path):
    """从原始 PDB 读取残基编号列表（按出现顺序，去重）"""
    resids = []
    seen = set()
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            chain = line[21]
            resSeq = int(line[22:26].strip())
            resName = line[17:20].strip()
            key = (chain, resSeq, resName)
            if key not in seen:
                seen.add(key)
                resids.append(resSeq)
    return resids


def postprocess():
    import numpy as np
    import mdtraj as md
    import mdtraj.formats as fmt

    print("\n[Post] Post-processing trajectory...")

    # 1. 用 OpenMM 先写一份蛋白质-only 的 PDB 作为拓扑
    #    （从 solvated.pdb 中提取，或用 final_structure.pdb 兜底）
    tmp_pdb = OUTDIR / "_protein_topo.pdb"
    if SOLVATED_PDB.exists():
        import MDAnalysis as mda
        u = mda.Universe(str(SOLVATED_PDB))
        u.select_atoms("protein").write(str(tmp_pdb))
    else:
        # 直接用 OpenMM 写全原子蛋白质
        from openmm.app import PDBFile as PDBFileApp, ForceField as FFApp, Modeller as ModApp
        from openmm.unit import nanometer as nm_unit
        _pdb = PDBFileApp(str(PDB))
        _ff  = FFApp("charmm36.xml", "charmm36/water.xml")
        _mod = ModApp(_pdb.topology, _pdb.positions)
        _mod.addHydrogens(_ff, pH=7.0)
        with open(tmp_pdb, "w") as f:
            PDBFileApp.writeFile(_mod.topology, _mod.positions, f)

    ref = md.load(str(tmp_pdb))
    n_protein = ref.n_atoms
    print(f"[Post] Protein atoms: {n_protein}")

    # 2. 恢复原始残基编号
    if PDB.exists():
        original_resids = read_original_resids(PDB)
        residues = list(ref.topology.residues)
        if len(residues) == len(original_resids):
            for res, orig_id in zip(residues, original_resids):
                res.resSeq = orig_id
            print(f"[Post] Residue IDs restored: {original_resids[0]}..{original_resids[-1]}")
        else:
            print(f"[Post] Residue count mismatch, skipping ID restoration")

    # 3. 从 DCD 读取蛋白质坐标
    with fmt.DCDTrajectoryFile(str(TRAJ_DCD)) as f:
        all_xyz, cell_lengths, cell_angles = f.read()

    n_frames = all_xyz.shape[0]
    prot_xyz = (all_xyz[:, :n_protein, :] / 10.0).astype(np.float32)
    box_len  = cell_lengths / 10.0
    box_ang  = cell_angles

    traj = md.Trajectory(prot_xyz, ref.topology,
                         unitcell_lengths=box_len, unitcell_angles=box_ang)

    # 4. 保存最终帧
    traj[-1].save_pdb(str(FINAL_PDB))
    print(f"[Post] final_structure.pdb saved ({n_frames} frames processed)")

    # 5. PBC wrap：蛋白质质心居中
    box_half = traj.unitcell_lengths / 2
    com      = traj.xyz.mean(axis=1, keepdims=True)
    traj.xyz += (box_half[:, None, :] - com)

    # 6. 写出可视化文件
    viz = OUTDIR / "viz"
    viz.mkdir(exist_ok=True)
    ref_pdb  = viz / "reference.pdb"
    traj_out = viz / "traj_protein_wrapped.dcd"

    traj[0].save_pdb(str(ref_pdb))
    traj.save_dcd(str(traj_out))

    # 清理临时文件
    tmp_pdb.unlink(missing_ok=True)

    print(f"\n  Visualization files -> {viz}/")
    print(f"    reference.pdb            <- topology (load first)")
    print(f"    traj_protein_wrapped.dcd <- wrapped trajectory")
    print(f"  VMD:     vmd {ref_pdb.name} {traj_out.name}")
    print(f"  ChimeraX: open {ref_pdb.name} ; open {traj_out.name}")

# ── 主程序 ────────────────────────────────────────────────────────────────────
def main():
    t0 = time.time()
    restart = CHECKPOINT.exists()

    print("=" * 60)
    print(f"  OpenMM Pipeline — 2RRM")
    print(f"  Output : {OUTDIR}")
    print(f"  Mode   : {'RESTART from checkpoint' if restart else 'FRESH start'}")
    print("=" * 60)

    # 1. 读取结构 & 力场
    print("\n[1/5] Loading structure...")
    pdb = PDBFile(str(PDB))
    ff  = ForceField("charmm36.xml", "charmm36/water.xml")

    # 2. 溶剂化
    print("[2/5] Solvating...")
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(ff, pH=7.0)
    modeller.addSolvent(ff,
        padding=1 * nanometer,
        ionicStrength=0.15 * molar,
        positiveIon="Na+", negativeIon="Cl-")
    print(f"      System size: {modeller.topology.getNumAtoms()} atoms")

    # 溶剂化结构存盘（后处理用）
    if not SOLVATED_PDB.exists():
        with open(SOLVATED_PDB, "w") as f:
            PDBFile.writeFile(modeller.topology, modeller.positions, f)

    # 3. 建立系统（HMR → 4 fs timestep）
    print("[3/5] Creating system...")
    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.2 * nanometer,
        constraints=HBonds,
        hydrogenMass=1.5 * amu,
    )

    for pname in ("CUDA", "OpenCL", "CPU"):
        try:
            platform = Platform.getPlatformByName(pname)
            props = {"CudaPrecision": "mixed"} if pname == "CUDA" else {}
            integrator = LangevinMiddleIntegrator(TEMP, 1 / picosecond, 4 * femtoseconds)
            sim = Simulation(modeller.topology, system, integrator, platform, props)
            print(f"[Platform] {pname}")
            break
        except Exception as e:
            print(f"[Platform] {pname} failed: {e}")
            continue
    else:
        raise RuntimeError("No usable platform found")

    if restart:
        print(f"[4/5] Loading checkpoint: {CHECKPOINT.name}")
        sim.loadCheckpoint(str(CHECKPOINT))
    else:
        sim.context.setPositions(modeller.positions)

        print("[4/5] Minimizing energy...")
        sim.minimizeEnergy(maxIterations=2000)

        print("      NVT equilibration (100 ps)...")
        sim.context.setVelocitiesToTemperature(TEMP)
        sim.reporters.append(StateDataReporter(str(LOG_NVT), REPORT_INTERVAL,
            step=True, potentialEnergy=True, temperature=True))
        sim.step(NVT_STEPS)
        sim.reporters.clear()

        print("      NPT equilibration (100 ps)...")
        system.addForce(MonteCarloBarostat(PRESSURE, TEMP))
        sim.context.reinitialize(preserveState=True)
        sim.reporters.append(StateDataReporter(str(LOG_NPT), REPORT_INTERVAL,
            step=True, potentialEnergy=True, volume=True, temperature=True))
        sim.step(NPT_STEPS)
        sim.reporters.clear()

    # 4. Production MD
    print("[5/5] Production MD (10 ns)...")
    sim.reporters.append(DCDReporter(str(TRAJ_DCD), REPORT_INTERVAL,
        append=restart))
    sim.reporters.append(StateDataReporter(str(LOG_PROD), REPORT_INTERVAL,
        step=True, time=True, potentialEnergy=True, temperature=True,
        progress=True, remainingTime=True, speed=True,
        totalSteps=PROD_STEPS, append=restart))
    sim.reporters.append(CheckpointReporter(str(CHECKPOINT), CKPT_INTERVAL))

    sim.step(PROD_STEPS)

    # 5. 后处理：wrap + 可视化轨迹 + 正确残基编号
    postprocess()

    elapsed = timedelta(seconds=int(time.time() - t0))
    print(f"\n[Done] Total wall time: {elapsed}")
    print(f"       Raw trajectory  : {TRAJ_DCD}")
    print(f"       Final structure : {FINAL_PDB}")

if __name__ == "__main__":
    main()
