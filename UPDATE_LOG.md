# Project Update Log
# 蛋白构象采样与药物设计平台

> 格式说明：每次更新在顶部追加一个新块。
> Git commit message 直接截取最新块的 `## [日期]` 至下一个 `---` 之间的内容。

---

## [2026-03-19] 环境搭建 — conda env allatom

### 完成内容
- 创建 conda 环境 `allatom`（Python 3.11），基于 `environment.yml`
- 安装核心包：OpenMM 8.4.0、MDAnalysis 2.10.0、MDTraj 1.11.1、PDBFixer、ParmEd、Paramiko、Fabric
- 验证脚本：`test_env.py`，全部通过

### 技术决策
- 主力 MD 引擎确定为 **OpenMM**（全原子阶段）
- Windows 本地：OpenCL GPU 平台可用（GTX 1650）；CUDA 平台保留给 Linux 服务器
- `environment.yml` 可由合作者直接 `conda env create` 复现，后期持续维护

### 下一步
- [ ] 测试 paramiko 连接到服务器
- [ ] 搭建第一个 Ala dipeptide 测试 pipeline（验证本地→服务器全流程）

---

## [2026-03-19] 项目初始化

### 科学目标
- 利用粗粒化（Coarse-Grained）和全原子（All-Atom）模拟采样蛋白构象空间
- 将代表性构象输入对接软件或蛋白设计软件，用于药物发现

### 本次完成
- 建立项目目录结构 `d:/all_atom`
- 配置 VSCode 紫色主题（`.vscode/settings.json`）
- 确立更新日志规范（本文件）
- 讨论 Windows→Linux 工作流及 OpenMM vs GROMACS 选型（见 `docs/workflow_notes.md`）

### 技术决策记录
- 主力 MD 引擎：待定（GROMACS vs OpenMM，见讨论）
- 粗粒化力场候选：Martini 3 / SIRAH / OpenCG
- 对接软件候选：AutoDock Vina / Glide / DiffDock

### 下一步
- [ ] 确定 Windows→Linux 自动化传输方案
- [ ] 确定 MD 引擎选型
- [ ] 搭建第一个测试系统

---

<!-- 新条目加在上方 --- 之前 -->
