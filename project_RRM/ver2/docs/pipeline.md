# Pipeline：CIF → 修饰 RNA → 4BS2 复合物 MD

## 总流程图

```
4BS2.cif
   │
   ├─▶ PDBFixer ───────▶ structures/prepared/4BS2_clean.pdb
   │
   ├─▶ ptm_builder ────▶ structures/modified/4BS2_8OG_at_G3.pdb
   │   (NeRF 加 O8/H7)            structures/modified/4BS2_PSU_at_U2.pdb
   │
   └─▶ tleap (载入 modxna) ──▶ AMBER prmtop/inpcrd
       或
       parmed (mol2/frcmod → OpenMM XML) ──▶ run_openmm.py
                  │
                  ▼
            生产 MD（OpenMM 8.x，HMR 4 fs，OPC 水）
                  │
                  ▼
            轨迹分析（MDAnalysis/MDTraj）
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

### 2. mol2 + frcmod → OpenMM XML（待写）

`pipeline/convert_modxna_to_openmm.py`：

```python
import parmed as pmd
# 加载 modxna 的 mol2 + frcmod
mol = pmd.load_file('force_fields/modxna/dat/lib_base/8OG.mol2')
mol.load_parameters('force_fields/modxna/dat/frcmod.modxna')
# 写成 OpenMM-loadable XML
ff_xml = pmd.openmm.OpenMMParameterSet.from_structure(mol)
ff_xml.write('force_fields/openmm_xml/8OG_RNA.xml')
```

（精确 API 待验证；备选用 `openmmforcefields.SystemGenerator` 的转换路径。）

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

### 5. 系统装配 + 运行

```python
# pipeline/run_md.py（OpenMM 路线）
from openmm.app import *
ff = ForceField(
    'amber14/protein.ff14SB.xml',
    'amber14/RNA.OL3.xml',
    'amber14/tip3p.xml',                       # 或 OPC
    'force_fields/openmm_xml/8OG_RNA.xml',
)
# Modeller, addSolvent, createSystem, ...
```

### 6. 验证标准（必跑）

任何新修饰参数集进入生产前：

1. **总电荷整数检查**（每残基）
2. **createSystem 不抛异常**
3. **free 单核苷 100 ns**：检查 χ 分布、O3'-C3' 距离、糖 pucker

只有全部通过才进入 4BS2 复合物。

## 集群部署

集群目录依然是 `bio:/data/biophys/carolinge/clawork/37_OXR/`。
ver_2 用 `replicas_v2/` 子目录，**不要覆盖** `replicas/`（ver_1 的产物，作为对照）。
