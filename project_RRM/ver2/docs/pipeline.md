# Pipeline：CIF → 修饰 RNA → 4BS2 复合物 MD（ver_2）

> **架构关键决策**（2026-05-01）：放弃 OpenMM XML 转换路线，改为 **prmtop 路线**。
> 原因：modXNA 的 base mol2 charges 是 fragment-level RESP，不能直接拼到 OL3 sugar+backbone（电荷不自洽）。modxna.sh + tleap 组装才得到正确电荷。OpenMM 的 `AmberPrmtopFile` 直接读 prmtop，根本不需要 OpenMM XML。

## 总流程图

```
4BS2.cif
   │
   ├─▶ PDBFixer ───────▶ structures/prepared/4BS2_clean.pdb
   │
   ├─▶ NeRF 几何拼接 ──▶ structures/modified/4BS2_8OG_at_G3.pdb
   │   (port from ptm_builder.py)
   │
   │     ┌──── 服务器（Linux + AmberTools）────┐
   │     ▼                                       │
   ├─▶ build_modxna_residues.sh                  │
   │   └─▶ 8OGI.lib, PUUI.lib, M1AI.lib  ◀─── 一次性，跨实验复用
   │                                             │
   ├─▶ build_amber_system.sh <pdb> <prefix>      │
   │   └─▶ tleap → prmtop + inpcrd  ◀─── 每个体系一次
   │     │
   │     ▼ rsync 回 Windows
   │
   └─▶ OpenMM 用 AmberPrmtopFile 直接读取 ──▶ run_md.py → MD 轨迹
                                                          │
                                                          ▼
                                                    分析（MDAnalysis）
```

## 详细步骤

### 1. 力场准备（一次性）

```bash
cd project_RRM/ver2/force_fields/
git clone https://github.com/modxna/modxna
# 检查 modxna/dat/lib_base/ 下是否有 8OG.mol2, PSU.mol2, M1A.mol2, M6A.mol2
# 检查 modxna/dat/frcmod.modxna
```

下载论文 SI 到 `papers_SI/`：

- Sarzyńska/Lahiri 2022 SI（ψ 升级版）
- Bussi `ff-m6a-fit5_AC.rtp`（m6A 升级版）

### 2. 在服务器上组装 modXNA 残基（一次性）

```bash
# 服务器登录后
ssh bio.example
cd ~/work/all_atom    # 假设仓库已经 rsync 过去
bash project_RRM/ver2/pipeline/setup_modxna_server.sh    # 安装 ambertools
bash project_RRM/ver2/pipeline/build_modxna_residues.sh  # 组装 .lib 文件
```

输出 `.lib` 文件到 `force_fields/openmm_xml/lib_amber/`：
- `8OGI.lib`、`8OG3.lib`、`8OG5.lib`（internal / 3'-cap / 5'-cap）
- `PUUI.lib`、`PUU3.lib`、`PUU5.lib`
- `M1AI.lib`、`M1A3.lib`、`M1A5.lib`

把这些 .lib 文件 rsync 回 Windows 即可一直用。

### 3. 结构准备

```python
# pipeline/prepare_system.py
from pdbfixer import PDBFixer
fixer = PDBFixer(filename="structures/original/4BS2.cif")
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.0)
# 保存到 structures/prepared/
```

### 4. 加修饰

复用 [`/app/ptm_builder.py`](../../../app/ptm_builder.py) 的 NeRF 逻辑（**不修改原文件**，复制并扩展）：

```python
# pipeline/modify_rna.py
from ptm_builder_v2 import apply_ptm   # ver_2 自己的 fork
new_pdb = apply_ptm("structures/prepared/4BS2_clean.pdb",
                    chain="B", resnum=3, ptm_code="8OG")
```

ver_2 需要新增的 patches：

- `PSU` (U → ψ)：N1-C1' 改 C5-C1' 键
- `M1A` (A → m1A)：N1 加甲基（CH3 + 删 N1-H）
- `M6A` (A → m6A)：N6 单 H 改成甲基

### 5. 服务器：tleap 装配 prmtop（每个体系一次）

```bash
bash project_RRM/ver2/pipeline/build_amber_system.sh \
    project_RRM/ver2/structures/modified/4BS2_8OG_G3.pdb \
    project_RRM/ver2/runs/4BS2_8OG_G3
```

输出：
- `runs/4BS2_8OG_G3.prmtop`
- `runs/4BS2_8OG_G3.inpcrd`
- `runs/4BS2_8OG_G3.solvated.pdb`

### 6. OpenMM：直接读 prmtop（Windows 或服务器都行）

```python
# pipeline/run_md.py
from openmm.app import AmberPrmtopFile, AmberInpcrdFile, Simulation, ...
prm = AmberPrmtopFile('runs/4BS2_8OG_G3.prmtop')
inp = AmberInpcrdFile('runs/4BS2_8OG_G3.inpcrd')
system = prm.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.2*nm, ...)
# minimize, equilibrate, produce
```

**不需要 OpenMM ForceField XML**。prmtop 已包含全部参数。

### 6. 验证标准（必跑）

任何新修饰参数集进入生产前：

1. **总电荷整数检查**（每残基）
2. **createSystem 不抛异常**
3. **free 单核苷 100 ns**：检查 χ 分布、O3'-C3' 距离、糖 pucker

只有全部通过才进入 4BS2 复合物。

## 集群部署

集群目录依然是 `bio:/data/biophys/carolinge/clawork/37_OXR/`。
ver_2 用 `replicas_v2/` 子目录，**不要覆盖** `replicas/`（ver_1 的产物，作为对照）。
