# 新成员上手指南

## 项目背景一句话

我们在用全原子 MD 研究 RNA 修饰（8oxoG / ψ / m1A / m6A）如何改变 TDP-43 RRM1+2（PDB 4BS2）对 RNA 的结合。详见 [`../README.md`](../README.md)。

## 第 0 步：明白哪些不是这个项目

- 根目录 [`/app/`](../../../app/)、[`/launch_ui.bat`](../../../launch_ui.bat) 是另一个独立工具（PDB 准备 Web UI），与本科学项目无关。
- 根目录 [`/docs/workflow_notes.md`](../../../docs/workflow_notes.md) 是更早的引擎选型讨论，可读但已过时。
- [`/project_RRM/ver1/`](../../ver1/) 是失败尝试的归档。**不要在 ver_1 上做新工作**。

## 第 1 步：环境

```bash
# Windows / Linux 通用
conda env create -f /project_RRM/environment_v2.yml
conda activate allatom_v2
```

服务器（A100 + CUDA 12.5）的环境构建脚本待补到 [`../setup_server_v2.sh`](../setup_server_v2.sh)。

## 第 2 步：拿到力场

```bash
cd /project_RRM/ver2/force_fields/
git clone https://github.com/modxna/modxna
```

记录 commit hash 到 `VERSIONS.md`。

论文 SI 文件需要手动下载：

- [Sarzyńska/Lahiri 2022 (DOI 10.1007/s10822-022-00447-4)](https://link.springer.com/article/10.1007/s10822-022-00447-4) → `papers_SI/sarzynska_2022_psu/`
- [Bussi m6A repo](https://github.com/bussilab/m6a-charge-fitting) → `papers_SI/bussi_m6a/`

## 第 3 步：跑第一个 replica

待 [`pipeline/`](../pipeline/) 脚本就绪后补本节。

预期一次完整运行：
1. `prepare_system.py 4BS2.cif` → 干净 PDB
2. `modify_rna.py --ptm 8OG --chain B --resnum 3` → 修饰 PDB
3. `run_md.py --replica 1 --time 100ns` → 生产
4. `analyze.py runs/4BS2_8OG_r1` → 报告

## 第 4 步：不能犯的错误

参考 [`lessons_from_ver1.md`](lessons_from_ver1.md)。一句话：**不要自己写力场 XML**，用已发表参数。

## 联系

`carolinge73@gmail.com`
