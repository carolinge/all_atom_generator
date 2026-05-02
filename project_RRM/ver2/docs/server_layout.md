# 服务器目录布局 + 工作流

## 集群信息

- 主机：`newton.pks.mpg.de`（Max Planck Institute for the Physics of Complex Systems, Dresden）
- 调度：Slurm
- 用户：`carolinge`
- 登录：`ssh bio`（`~/.ssh/config` 已配 alias）
- 项目根：`/data/biophys/carolinge/clawork/37_OXR/`（也可走 `~/work/37_OXR` 软链）
- /home 仅 80 GB —— 任何产物**不能**放 /home，必须在 /data

## 目录布局

```
/data/biophys/carolinge/clawork/37_OXR/
├── repo/                       ← rsync 自 d:/all_atom（sync_to_server.sh）
│   ├── app/                       (不动)
│   └── project_RRM/
│       ├── ver1/                  (归档)
│       └── ver2/
│           ├── force_fields/modxna/   ← 不动；modxna.sh 来源
│           ├── pipeline/              ← 所有 server-side 脚本
│           └── ...
│
├── replicas/                   ← ver_1 历史 replicas（保留）
│   └── 4BS2_{WT,8OG}_r{1..3}/
│
├── replicas_v2/                ← ver_2 运行目录（新）
│   └── <NAME>/                    e.g. 4BS2_8OG_G3
│       ├── system/
│       │   ├── <NAME>.prmtop
│       │   ├── <NAME>.inpcrd
│       │   └── <NAME>.solvated.pdb
│       ├── r1/
│       │   ├── run_md.py            （从 pipeline/ 复制）
│       │   ├── sub.sh               （从模板生成，job-name=<NAME>_r1）
│       │   ├── replica.txt          （seed 编号）
│       │   ├── traj.dcd             （产物）
│       │   ├── md.log / nvt.log / npt.log
│       │   └── checkpoint.chk / final.chk
│       ├── r2/, r3/                 （同上，不同 seed）
│       └── system/                  （prmtop+inpcrd 共享）
│
├── force_fields/               ← ver_2 共用力场产物
│   └── lib_amber/                 build_modxna_residues.sh 输出
│       ├── 8OGI.lib  (internal)
│       ├── 8OG3.lib  (3'-cap)
│       ├── 8OG5.lib  (5'-cap)
│       ├── PUUI.lib  PUU3.lib  PUU5.lib
│       └── M1AI.lib  M1A3.lib  M1A5.lib
│
├── structures/                 ← ver_2 共用准备好的 PDB（暂未填充）
│   ├── original/4BS2.cif
│   ├── prepared/                  PDBFixer 输出
│   └── modified/                  加好修饰的 PDB
│
└── pipeline_logs/              ← 长跑脚本的 log
```

## Conda 环境架构（双 env）

| Env | 来源 | 内容 | 用途 |
|------|------|------|------|
| `allatom_v2` | clone of ver_1 `allatom` | OpenMM 8.2 + ParmEd + PDBFixer + MDAnalysis + ... | OpenMM MD 跑 + 分析 |
| `amber24` | micromamba 新建 | ambertools 24.8 + python 3.11 | tleap / cpptraj / sander |

**为什么分开**：把 ambertools 装进 allatom_v2 会触发 conda 2022.10 classic solver 的 phantom `__glibc` conflict（30 分钟 solve 失败）。micromamba (libsolv) 在新 env 上 1-3 分钟完成，干净。`server_setup.sh` 自动两个都建。

## 一次性设置（首次或有更新）

```bash
# 1) 在本机 d:/all_atom 上：
bash project_RRM/ver2/pipeline/sync_to_server.sh

# 2) 登录服务器：
ssh bio
cd /data/biophys/carolinge/clawork/37_OXR/repo

# 3) 一次性建两个 env + 装 micromamba（约 5 分钟）：
bash project_RRM/ver2/pipeline/server_setup.sh

# 4) 组装 modXNA 修饰核苷的 .lib 文件（一次性，跨 replica 复用）：
bash project_RRM/ver2/pipeline/build_modxna_residues.sh
# 脚本自动 prepend amber24 PATH，不需要 activate
# 应输出 9 个 .lib: 8OG{I,3,5}.lib  PUU{I,3,5}.lib  M1A{I,3,5}.lib
```

## 每个新体系（如 4BS2 + 8OG@G3）

```bash
# 在本机准备好已经加好 8OG 几何的 PDB（待 ptm_builder_v2 完成）：
#   project_RRM/ver2/structures/modified/4BS2_8OG_G3.pdb

# 同步到服务器：
bash project_RRM/ver2/pipeline/sync_to_server.sh

# 登录服务器：
ssh bio
cd /data/biophys/carolinge/clawork/37_OXR/repo
conda activate allatom_v2

# 建 prmtop（每个体系一次）：
bash project_RRM/ver2/pipeline/build_amber_system.sh \
    4BS2_8OG_G3 \
    project_RRM/ver2/structures/modified/4BS2_8OG_G3.pdb

# 衍生 r1/r2/r3 + 提交：
bash project_RRM/ver2/pipeline/spawn_replicas.sh 4BS2_8OG_G3
for r in /data/biophys/carolinge/clawork/37_OXR/replicas_v2/4BS2_8OG_G3/r*; do
    (cd "$r" && sbatch sub.sh)
done
squeue -u carolinge
```

## 监控

```bash
squeue -u carolinge                                            # 队列
sacct -u carolinge --starttime=2026-05-01                      # 历史
tail -f replicas_v2/4BS2_8OG_G3/r1/md.log                      # 进度
ls -la replicas_v2/4BS2_8OG_G3/r1/traj.dcd                     # 轨迹大小
```

## 已知坑

1. **/home 80 GB 配额**——所有产物必须在 `/data/biophys/carolinge/clawork/37_OXR/` 下。脚本默认就是这个。
2. **CRLF**——从 Windows 上传的 *.sh 必须 LF。`sync_to_server.sh` + `spawn_replicas.sh` 都自动 `sed -i 's/\r$//'` 处理。
3. **CONECT records**——非标准残基的 PDB 要带。我们的 prmtop 路线是从 tleap 走 .lib 库的，不依赖 CONECT；但 build_amber_system.sh 仍读 PDB，所以保留 CONECT 是好习惯。
4. **CUDA 版本失配**——run_md.py 已经有 CUDA→OpenCL→CPU fallback。
5. **mid-run 崩溃**——ver_2 的 run_md.py **支持** restart：检测到 `checkpoint.chk` 就跳过 min/equil 直接续 production。重跑只需再 `sbatch sub.sh`。
