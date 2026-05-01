# ver_1 — 归档（2026-03 至 2026-04）

## 当时做了什么

1. 用 [`scripts/gen_4bs2_md.py`](scripts/gen_4bs2_md.py) 从 [`4BS2.cif`](4BS2.cif) 准备了 WT 和 G3→8OG 两套结构
2. 自建力场：[`/app/ff_custom/8OG_RNA_amber.xml`](../../app/ff_custom/8OG_RNA_amber.xml)（**仍在 app/ 不要动**）
   - 电荷：照 G 类比迁移，O8/H7 估算后手工 redistribute
   - bonded：复用 G 的所有 bond/angle/dihedral
3. 用 [`scripts/gen_replicas.py`](scripts/gen_replicas.py) 生成 6 replicas（WT × 3 + 8OG × 3，独立 seed），每个 100 ns 生产
4. 集群部署：`bio:/data/biophys/carolinge/clawork/37_OXR/replicas/`
5. 分析栈：[`scripts/analyze_binding.py`](scripts/analyze_binding.py)、[`scripts/analyze_replicas.py`](scripts/analyze_replicas.py)、[`scripts/analyze_rmsf.py`](scripts/analyze_rmsf.py)、[`scripts/compare_wt_8og.py`](scripts/compare_wt_8og.py) 等

## 失败与教训

详见 [`/project_RRM/ver2/docs/lessons_from_ver1.md`](../ver2/docs/lessons_from_ver1.md)。简而言之：**自建 XML 是错误路线**。结构（NeRF 几何拼接）本身没问题，问题在力场——

1. C8=O8 / N7-H7 引入的新 bond/angle/improper 没补全
2. 重新分配电荷后总和大概率非整数（手工估算精度 ~0.01-0.1 e）
3. 8-oxoG 的 6,8-双酮互变异构体相邻原子电荷应该一起重 RESP，不是只改 4 个原子
4. χ_OL3 二面角对 N9 类型敏感（这一项保留 `RNA-N*` 实际是对的，问题不在这里）

## 不要在 ver_1 上继续工作

所有新工作在 [`/project_RRM/ver2/`](../ver2/) 下进行。这里保留是为了：
- 保留 6 replicas 的 100 ns 轨迹（[`output/4BS2_*_r*/`](output/)）作为对照（后续看 ver_2 是否能复现/不复现某些 contact pattern）
- 保留分析脚本的实现，ver_2 的分析重写可参考

## 文件清单

```
ver1/
├── 4BS2.cif                  原始晶体结构（CIF 格式）
├── 4juy.pdb                  另一个早期测试结构
├── output/                   所有 MD 运行产物（gitignored，本地 staging）
│   ├── 4BS2_WT_MD/           初版 WT 单 replica
│   ├── 4BS2_8OG_G3_MD/       初版 8OG 单 replica
│   ├── 4BS2_WT_r1..r3/       WT 三 replicas（gitignored）
│   ├── 4BS2_8OG_r1..r3/      8OG 三 replicas（gitignored）
│   └── analysis/             contact map / RMSF 输出
├── prepared/                 PDBFixer 输出（gitignored）
├── uploaded/                 上传的 PDB（gitignored）
├── scripts/                  生成 + 分析脚本
└── example_2RRM/             更早的 2RRM 单蛋白 MD 测试
```
