# Windows → Linux 工作流 & 引擎选型讨论

## 问题：Windows本地 → Linux服务器 自动化提交

### 核心痛点
- GROMACS 在 Windows 上没有原生支持，必须在 Linux 端运行
- 打包上传整个模拟环境（top/gro/mdp/forcefield）繁琐且容易出错
- 服务器权限受限，无 sudo

### 推荐方案：Python 脚本生成 + Paramiko 自动传输

```
本地 Windows
│
├── Python 脚本（生成所有输入文件 + bash提交脚本）
│       ↓  paramiko / fabric
└── 自动SCP传输到服务器 → 自动 sbatch / qsub 提交
```

**优势**：
- 所有参数在本地 Python 里管理，版本控制友好
- 一行命令完成"生成→传输→提交"全流程
- 不需要服务器上有任何特殊权限

**核心库**：
- `paramiko`：Python SSH/SCP 库，纯 Python，pip 即可装
- `fabric`：对 paramiko 的高层封装，适合批量任务
- （可选）`rich`：本地终端进度显示

### 典型流程示意

```python
# local_submit.py（本地运行）
from pathlib import Path
import paramiko

def generate_mdp(output_dir, nsteps=50000, temp=300):
    """在本地生成 GROMACS mdp 文件"""
    content = f"""
integrator  = md
nsteps      = {nsteps}
dt          = 0.002
ref_t       = {temp}
...
"""
    (output_dir / "md.mdp").write_text(content)

def upload_and_submit(local_dir, remote_dir, server, user, key_path):
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(server, username=user, key_filename=key_path)

    sftp = ssh.open_sftp()
    # 上传整个目录
    for f in local_dir.iterdir():
        sftp.put(str(f), f"{remote_dir}/{f.name}")

    # 提交作业
    _, stdout, _ = ssh.exec_command(f"cd {remote_dir} && sbatch submit.sh")
    print(stdout.read().decode())
    ssh.close()
```

---

## OpenMM 能完全替代 GROMACS 吗？

### 结论（先说）
**不能完全替代，但对我们的项目有显著优势，值得作为主力或并行引擎。**

### 详细对比

| 维度 | GROMACS | OpenMM |
|------|---------|--------|
| 性能（纯速度） | ✅ 极优，高度优化的C++ | 🟡 稍慢，但GPU上差距很小 |
| Windows友好度 | ❌ 不支持，必须Linux | ✅ pip安装，Windows/Linux通用 |
| 脚本自动化 | 🟡 CLI工具，需shell脚本 | ✅ 纯Python API，极易参数化 |
| 力场支持 | ✅ AMBER/CHARMM/OPLS/Martini | ✅ 通过OpenFF/ParmEd可读GROMACS力场 |
| 粗粒化（Martini） | ✅ 原生支持 | 🟡 有插件(OpenMartini)，不够成熟 |
| 社区生态 | ✅ 庞大，文献多 | ✅ 增长快，Python生态友好 |
| 集群提交 | 🟡 需要手动配置脚本 | ✅ Python脚本直接生成提交 |
| 可复现性 | 🟡 版本敏感 | ✅ pip锁定版本，conda env |
| 自由能计算 | ✅ 成熟（FEP、AWH等） | 🟡 有但生态不如GROMACS |

### 对我们项目的建议

1. **全原子采样（conformational sampling）**
   → OpenMM + OpenFF/CHARMM36 完全够用，Python化的优势压倒一切

2. **粗粒化（Martini 3）**
   → 目前仍推荐 GROMACS，Martini生态在GROMACS上最成熟
   → 但脚本生成+paramiko传输解决工作流问题

3. **混合策略（推荐）**
   - 粗粒化：GROMACS（脚本自动生成+传输）
   - 全原子精化：OpenMM（本地Python脚本，直接运行或传输）
   - 统一接口：MDAnalysis / MDTraj 做轨迹分析，两者输出都支持

### OpenMM 的关键优势场景

```python
# 这在GROMACS里需要写shell脚本+多个命令
# OpenMM里是纯Python，可以for循环跑多个系统

for protein in protein_list:
    system = create_system(protein, forcefield="charmm36")
    simulation = run_nvt(system, temp=300, steps=500000)
    coords = extract_representative_frames(simulation, n=10)
    save_for_docking(coords, output_dir / protein.name)
```

---

## 下一步行动项
- [ ] 测试 paramiko 连接到服务器
- [ ] 确认服务器上 GROMACS/OpenMM 版本
- [ ] 搭建第一个 Ala dipeptide 测试系统（作为pipeline验证）
