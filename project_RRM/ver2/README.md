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

- [x] git clone modXNA → [`force_fields/modxna/`](force_fields/modxna/)（commit `595cff1b`，GPL-3.0）
- [x] 验证 LICENSE：GPL-3.0 ✅
- [x] 检查 modXNA 实际覆盖 — 详见 [`results/inspect_modxna.md`](results/inspect_modxna.md)：
  - **8OG** ✅ 标准
  - **PUU** ✅ 标准 ψ（**不是 PSU** —— PSU 是设计衍生物）
  - **M1A** ✅ 标准
  - **m6A** ❌ modXNA 不含；外部使用 Bussi 2022
- [x] **架构修正**：放弃 OpenMM XML 路线，改用 prmtop 路线（详见 [`docs/pipeline.md`](docs/pipeline.md)）
- [x] 服务器端组装脚本骨架完成：
  - [`pipeline/setup_modxna_server.sh`](pipeline/setup_modxna_server.sh)
  - [`pipeline/build_modxna_residues.sh`](pipeline/build_modxna_residues.sh)
  - [`pipeline/build_amber_system.sh`](pipeline/build_amber_system.sh)
- [ ] **下一步**：在服务器上跑通 setup → build_modxna → 得到 8OGI.lib / PUUI.lib / M1AI.lib
- [ ] 移植 NeRF 几何拼接到 `pipeline/`（注意 8OG.mol2 的 atom 6 `O6` 应重命名为 `O8`）
- [ ] 写 `pipeline/run_md.py`（OpenMM 用 AmberPrmtopFile 读 prmtop）
- [ ] 在 4BS2 上跑第一个 8OG replica 验证 pipeline
- [ ] m1A / m6A free-RNA 验证（无蛋白）
