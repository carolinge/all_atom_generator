# -*- coding: utf-8 -*-
"""
轨迹后处理脚本 — 独立运行，不需要重跑模拟
功能：
  1. 从 DCD 直接提取蛋白质原子坐标（无需完整溶剂化拓扑）
  2. 恢复原始 PDB 残基编号（从输入结构读取映射）
  3. PBC wrap（蛋白质居中）
  4. 输出 VMD/Chimera 可读的 reference.pdb + traj_protein_wrapped.dcd
  5. 输出保留原始残基编号的 final_structure.pdb（最终帧）

用法：
  conda activate allatom
  python postprocess.py
  python postprocess.py --traj output/traj.dcd
"""

import argparse
from pathlib import Path

HERE         = Path(__file__).parent
INPUT_PDB    = HERE / "fold_2026_03_23_21_01_model_4_from_cif.pdb"
OUTDIR       = HERE / "output"
TRAJ_DCD     = OUTDIR / "traj.dcd"
PROTEIN_PDB  = OUTDIR / "final_structure.pdb"
FINAL_PDB    = OUTDIR / "final_structure.pdb"
VIZ_DIR      = OUTDIR / "viz"


def read_original_resids(pdb_path: Path):
    """从原始 PDB 读取残基编号列表（按残基出现顺序，去重）"""
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


def restore_resids(topology, original_resids):
    """将 MDTraj topology 中蛋白质残基编号修正为原始值"""
    residues = list(topology.residues)
    if len(residues) != len(original_resids):
        print(f"  [警告] 残基数量不匹配: topology {len(residues)} vs 原始 {len(original_resids)}")
        print(f"         跳过残基编号修正")
        return
    for res, orig_id in zip(residues, original_resids):
        res.resSeq = orig_id
    print(f"  残基编号已修正: {original_resids[0]}..{original_resids[-1]}")


def run(traj_path: Path, protein_pdb: Path, input_pdb: Path):
    import numpy as np
    import mdtraj as md
    import mdtraj.formats as fmt

    if not traj_path.exists():
        raise FileNotFoundError(f"轨迹文件不存在: {traj_path}")
    if not protein_pdb.exists():
        raise FileNotFoundError(f"蛋白质 PDB 不存在: {protein_pdb}")

    # 1. 加载蛋白质拓扑
    print(f"Protein topology : {protein_pdb}")
    ref = md.load(str(protein_pdb))
    n_protein = ref.n_atoms
    print(f"Protein atoms    : {n_protein}")

    # 2. 恢复原始残基编号
    if input_pdb.exists():
        print(f"Original PDB     : {input_pdb}")
        original_resids = read_original_resids(input_pdb)
        restore_resids(ref.topology, original_resids)
    else:
        print(f"[注意] 原始 PDB 不存在 ({input_pdb})，残基编号保持不变")

    # 3. 从 DCD 读取所有帧（仅取前 n_protein 列 = 蛋白质）
    #    OpenMM 建系统时蛋白质原子永远排在最前面
    print(f"\nReading DCD: {traj_path}")
    with fmt.DCDTrajectoryFile(str(traj_path)) as f:
        all_xyz, cell_lengths, cell_angles = f.read()   # Angstrom

    n_frames = all_xyz.shape[0]
    n_total  = all_xyz.shape[1]
    print(f"Frames: {n_frames}  |  Total atoms in DCD: {n_total}")

    # 蛋白质坐标，Angstrom -> nm
    prot_xyz = (all_xyz[:, :n_protein, :] / 10.0).astype(np.float32)
    box_len  = cell_lengths / 10.0
    box_ang  = cell_angles

    # 4. 构建 MDTraj Trajectory（使用修正后的拓扑）
    traj = md.Trajectory(
        prot_xyz,
        ref.topology,
        unitcell_lengths=box_len,
        unitcell_angles=box_ang,
    )

    # 5. 保存最终帧
    print("\n[1/3] Saving final_structure.pdb (last frame)...")
    traj[-1].save_pdb(str(FINAL_PDB))
    print(f"      -> {FINAL_PDB}")

    # 6. PBC wrap：蛋白质质心平移到盒子中心
    print("[2/3] Centering protein in box...")
    box_half = traj.unitcell_lengths / 2
    com      = traj.xyz.mean(axis=1, keepdims=True)
    traj.xyz += (box_half[:, None, :] - com)

    # 7. 写出可视化文件
    VIZ_DIR.mkdir(parents=True, exist_ok=True)
    ref_pdb  = VIZ_DIR / "reference.pdb"
    traj_out = VIZ_DIR / "traj_protein_wrapped.dcd"

    print(f"[3/3] Writing visualization files ({n_frames} frames)...")
    traj[0].save_pdb(str(ref_pdb))
    traj.save_dcd(str(traj_out))

    print(f"\nDone! Output files:")
    print(f"  {FINAL_PDB}")
    print(f"  {ref_pdb}")
    print(f"  {traj_out}")
    print(f"\nVMD:     vmd \"{ref_pdb}\" \"{traj_out}\"")
    print(f"ChimeraX: open \"{ref_pdb}\" ; open \"{traj_out}\"")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MD trajectory post-processing")
    parser.add_argument("--traj",        type=Path, default=TRAJ_DCD)
    parser.add_argument("--protein-pdb", type=Path, default=PROTEIN_PDB,
                        dest="protein_pdb")
    parser.add_argument("--input-pdb",   type=Path, default=INPUT_PDB,
                        dest="input_pdb",
                        help="Original input PDB (for restoring residue IDs)")
    args = parser.parse_args()
    run(args.traj, args.protein_pdb, args.input_pdb)
