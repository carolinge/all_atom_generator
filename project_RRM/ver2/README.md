# ver_2 — TDP-43 RRM + 修饰 RNA MD（主线）

## 与 ver_1 的关键区别

| 维度 | ver_1 | ver_2 |
|------|-------|-------|
| 力场来源 | 自建 XML（电荷类比，无 RESP） | 已发表参数（modXNA 主，Sarzyńska/Bussi 增强） |
| 8oxoG | 自写 OpenMM XML | modXNA `8OG.mol2` + `frcmod.modxna` |
| ψ | — | Sarzyńska/Lahiri 2022（DOI 10.1007/s10822-022-00447-4）SI |
| m1A / m6A | — | Bussi `ff-m6a-fit5_AC.rtp` + Aduri 2007 (m1A)；**只做 free-RNA 验证** |
| 几何修饰 | NeRF（[`/app/ptm_builder.py`](../../app/ptm_builder.py)） | 同款 NeRF（已验证）+ 扩展 PSU/M1A/M6A patches |
| 文档 | 散落在 scripts/ 注释 | 集中在 [`docs/`](docs/) |

## 目录

```
ver2/
├── README.md                  本文件
├── docs/
│   ├── force_field_sources.md   每个修饰的来源、DOI、下载位置、版本固定
│   ├── lessons_from_ver1.md     ver_1 失败的 4 类原因
│   ├── pipeline.md              CIF → 准备 → 加修饰 → 力场 → MD → 分析
│   └── onboarding.md            新成员从 0 到第一个 replica 的步骤
├── force_fields/
│   ├── modxna/                  git submodule / clone（Bergonzo 2024）
│   └── openmm_xml/              modxna mol2/frcmod → OpenMM XML 转换结果
├── structures/
│   ├── original/                4BS2.cif（从 ver1 复制）
│   ├── prepared/                PDBFixer 输出
│   └── modified/                加了 8OG / ψ 的 PDB
├── pipeline/                    生成 + 运行 + 分析脚本（待写）
├── runs/                        本地 MD staging
└── results/                     分析输出
```

## 当前进度

### 完成

- [x] git clone modXNA → [`force_fields/modxna/`](force_fields/modxna/)（commit `595cff1b`，GPL-3.0）
- [x] 检查 modXNA 实际覆盖 — 见 [`results/inspect_modxna.md`](results/inspect_modxna.md)
- [x] 架构决策：**prmtop 路线**（见 [`docs/pipeline.md`](docs/pipeline.md)）；OpenMM 通过 `AmberPrmtopFile` 直接读 tleap 输出
- [x] 服务器端 pipeline 完整骨架：
  - [`pipeline/sync_to_server.sh`](pipeline/sync_to_server.sh) — Windows → bio rsync（CRLF→LF 自动处理）
  - [`pipeline/server_setup.sh`](pipeline/server_setup.sh) — 一次性建 conda env `allatom_v2` + 装 AmberTools
  - [`pipeline/build_modxna_residues.sh`](pipeline/build_modxna_residues.sh) — modxna.sh 组装 9 个 .lib（3 mods × 3 variants）
  - [`pipeline/build_amber_system.sh`](pipeline/build_amber_system.sh) — tleap → prmtop+inpcrd
  - [`pipeline/spawn_replicas.sh`](pipeline/spawn_replicas.sh) — 衍生 r1/r2/r3 + 渲染 sub.sh
  - [`pipeline/run_md.py`](pipeline/run_md.py) — production driver（含 checkpoint restart）
  - [`pipeline/sub.sh.template`](pipeline/sub.sh.template) — Slurm 模板（A100, 24h, allatom_v2）
- [x] [`docs/server_layout.md`](docs/server_layout.md) — 服务器目录 + 完整工作流文档

### 下一步

- [ ] **服务器验证 pass 1**：你来跑 sync → server_setup → build_modxna_residues，把 9 个 .lib 文件的 ls 输出给我
- [ ] 根据 .lib 实际命名（特别是 8OG 是否还有 O6/O8 重名），写 `pipeline/ptm_builder_v2.py`：
  - 移植 [`/app/ptm_builder.py`](../../app/ptm_builder.py) 的 8OG NeRF 逻辑（需要看 lib 后调整 atom 名）
  - 加 PUU patch（U → ψ：拆 N1-C1' 键、装 C5-C1' 键、加 N1-H、删 C5-H5）
  - 加 M1A patch（A → m1A：N1 上加甲基 C11/H1×3）
- [ ] 写 `pipeline/prepare_system.py`（PDBFixer + ptm_builder_v2 一键流程）
- [ ] 在 4BS2 跑第一个 8OG@G3 replica 验证整条 pipeline
- [ ] m1A / m6A free-RNA 验证（无蛋白）—— m6A 这一步要单做（不在 4BS2 主线）
