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

- [ ] git clone modXNA → [`force_fields/modxna/`](force_fields/modxna/)
- [ ] 验证 modXNA LICENSE
- [ ] 检查 modXNA 实际覆盖（8OG / PSU / M1A / M6A 是否齐全）
- [ ] 写 modxna mol2/frcmod → OpenMM XML 的转换脚本（用 ParmEd）
- [ ] 移植 [`ptm_builder.py`](../../app/ptm_builder.py) 的 NeRF 逻辑到 [`pipeline/`](pipeline/)，扩展 PSU/M1A/M6A patches
- [ ] 在 4BS2 上跑第一个 8OG replica 验证 modXNA pipeline
- [ ] m1A / m6A free-RNA 验证（无蛋白）
