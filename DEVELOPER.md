# all_atom — Developer Guide

> 这是一个本地 MD 模拟建模平台，运行在 Windows 上，生成可直接在 Linux 服务器提交的模拟文件。

---

## 快速启动

```bash
# 1. 创建环境（首次）
conda env create -f environment.yml
conda activate allatom

# 2. 启动 UI
python app/server.py
# 浏览器打开 http://localhost:5000

# 或双击
launch_ui.bat
```

---

## 项目结构

```
d:/all_atom/
├── app/
│   ├── server.py          # Flask 后端，全部 API
│   ├── templates/
│   │   └── index.html     # 前端：4步向导，全部 JS/CSS 内联
│   └── static/            # 预留静态资源目录（暂空）
│
├── prepared/              # PDBFixer 输出（自动生成）
├── output/                # 生成的模拟文件包（自动生成）
├── uploaded/              # 用户上传的 PDB（自动生成）
│
├── environment.yml        # conda 环境定义（维护此文件）
├── launch_ui.bat          # Windows 一键启动脚本
├── DEVELOPER.md           # 本文件
└── UPDATE_LOG.md          # 更新记录
```

---

## 技术栈

| 层 | 工具 | 说明 |
|---|-----|------|
| 模拟引擎 | OpenMM 8.x | Python-native，CUDA 加速 |
| 结构处理 | PDBFixer | 缺失残基/原子/氢原子 |
| 格式转换 | ParmEd | OpenMM → GROMACS 拓扑 |
| 力场 | CHARMM36m / AMBER ff14SB | OpenMM 内置 XML |
| 后端 | Flask | 本地单用户，无需认证 |
| 前端 | Vanilla JS | 无框架，无构建工具 |

---

## API 端点

| 方法 | 路由 | 功能 |
|-----|------|------|
| GET  | `/api/pdbs` | 列出项目根目录下的 PDB 文件（排除 prepared/output/uploaded） |
| POST | `/api/upload` | 上传新 PDB 文件到 `uploaded/` |
| POST | `/api/analyze` | 用 PDBFixer 分析 PDB，返回链信息、缺失原子数 |
| POST | `/api/run_pdbfixer` | 后台运行 PDBFixer（线程），结果写入 `prepared/` |
| POST | `/api/classify` | 识别结构中的分子类型（蛋白/DNA/RNA/PTM/离子/未知） |
| POST | `/api/check_ff` | 检查所选力场组合是否兼容 |
| POST | `/api/generate` | 后台生成全套模拟文件，打包 ZIP |
| GET  | `/api/download` | 下载最近生成的 md_files.zip |
| GET  | `/api/log_stream` | SSE 流，实时推送后台任务日志 |

---

## 4步向导逻辑

```
Step 1  选择/上传 PDB → /api/analyze → 显示链和缺失信息
Step 2  选链 + 突变 + PDBFixer 选项 → /api/run_pdbfixer → 输出 prepared/*.pdb
Step 3  /api/classify 识别组件 → 用户选力场 → /api/check_ff 验证
Step 4  /api/generate → 生成 run_openmm.py + GROMACS 全套 + md_files.zip
```

步骤状态：`locked` → `active` → `done`，由前端 JS 管理，后端无状态。
后台任务日志通过 `SESSION` dict（模块级）传递 `zip_path`，其余步骤数据全在前端 `S` 对象中。

---

## 力场体系

### 当前支持（内置于 OpenMM）

| 类型 | CHARMM 路线 | AMBER 路线 |
|-----|------------|-----------|
| 蛋白 | `charmm36.xml` | `amber14/protein.ff14SB.xml` |
| DNA  | `charmm36/na.xml` | `amber14/DNA.OL15.xml` |
| RNA  | `charmm36/na.xml` | `amber14/RNA.OL3.xml` |
| 水   | `charmm36/water.xml` | `amber14/tip3p.xml` |

**规则：同一模拟必须使用同一生态（CHARMM 或 AMBER），不可混用。**

### 添加新力场

1. 把 OpenMM XML 格式的力场文件放到 `app/ff_custom/` 目录
2. 在 `server.py` 的 `FF_PROTEIN` / `FF_DNA` / `FF_RNA` 字典里加入新条目：
   ```python
   FF_RNA["my_modified_na"] = {
       "label": "My modified RNA FF",
       "eco": "charmm",
       "files": ["charmm36/na.xml", "ff_custom/my_mod.xml"]
   }
   ```
3. 重启服务器，Step 3 的下拉菜单自动出现新选项

### 添加新修饰残基识别

在 `server.py` 顶部的 `KNOWN_PTM` 字典添加：
```python
KNOWN_PTM["8OG"] = "8-oxoguanosine (RNA)"
```

---

## 分子分类逻辑

`classify_components()` 按残基名称查表：

- `STANDARD_AA`：20种标准氨基酸 + ACE/NME 封端
- `DNA_RES` / `RNA_RES`：标准核苷酸缩写
- `KNOWN_PTM`：已知修饰残基，归入蛋白链但单独列出
- `IONS` / `WATER`：特殊处理
- 其余：归入 `unknown`（未来接配体参数化）

---

## GROMACS 输出说明

通过 ParmEd 从 OpenMM System 对象转换生成：
```python
structure = pmd.openmm.load_topology(topology, system, xyz=positions)
structure.save("topol.top")
structure.save("system.gro")
```

**已知限制：**
- ParmEd 转换对 CHARMM36m CMAP 修正项支持不完整
- 建议在 GROMACS 运行前用 `gmx grompp` 验证 topol.top
- OpenMM 脚本（`run_openmm.py`）生成的结果更可靠

---

## 常见问题

**`addMissingHydrogens` 不存在**
PDBFixer 版本 ≥1.10，方法名是 `addMissingHydrogens`，不是 `addHydrogens`。

**按钮无法点击（禁止符号）**
- Step N 的按钮在前一步完成前保持 `disabled`
- 如果 PDB 分析报错，检查浏览器控制台（F12）和右侧日志面板

**服务器找不到窗口**
用 `launch_ui.bat` 启动，它会开一个独立终端窗口，关闭窗口即停止。

**ParmEd GROMACS 导出失败**
OpenMM 模拟脚本仍然有效，GROMACS 文件暂时跳过，等待手动处理。

---

## 未来扩展路线

| 优先级 | 功能 | 技术路线 |
|--------|------|---------|
| 高 | 小分子配体参数化 | CGenFF API 或 AmberTools GAFF2 |
| 高 | 修饰核苷（8-oxoG RNA）| 自建 XML 参数文件 + ParmEd |
| 中 | 糖基化 | CHARMM Glycan 参数集 |
| 中 | 粗粒度（Martini） | 单独 GROMACS 工作流 |
| 低 | 增强采样 | PLUMED 接口 |
| 低 | 远程服务器提交 | Paramiko + Fabric（已在环境中） |

---

## 环境维护

```bash
# 添加新包
conda install -n allatom -c conda-forge <package>
# 然后更新 environment.yml 中的 dependencies 列表

# 团队成员复现环境
conda env create -f environment.yml
conda activate allatom
```
