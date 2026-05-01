# Force Field Source Versions

每次更新某个力场源时，追加新行到对应表。**永远不删除旧行**——保留版本可追溯性。

## modXNA

| 日期 | Commit hash | 来源 URL | 备注 |
|------|------------|----------|------|
| 2026-05-01 | `595cff1b4cf7c6260d9f482c951a79b66b43c550` | github.com/modxna/modxna | 当前主线；`Merge PR #8 from cb-fix-lib_backbone`（2026-04-14） |

**License**: GPL-3.0（拷贝在 `modxna/LICENSE`）。
**含义**：我们使用、修改 OK；如果发布修改版需要保持 GPL-3.0。学术内部使用无限制。

**镜像 URL（备用）**: github.com/drroe/modXNA（modXNA 官网首页指向，[modxna.chpc.utah.edu](https://modxna.chpc.utah.edu/) 显示这是 Daniel R. Roe 的镜像）

**目录关键路径**：
- `dat/lib_base/{8OG,PSU,M1A,M6A}.mol2` — 我们 4 个目标修饰的 base 片段
- `dat/lib_backbone/` — backbone 片段（连接 base 与磷酸）
- `dat/lib_sugar/` — sugar 片段（ribose 标准 + 各种修饰糖）
- `dat/frcmod.modxna` — 共享的 bond/angle/torsion 参数（覆盖 modrna08 + 扩展）
- `modxna.sh` — bash 主脚本，把 base+sugar+backbone 模块组装成完整核苷
- `Catalog – ModXNA.2026_03_05.pdf` — 完整目录（2 MB，待人工查阅）

**已知问题（待验证）**：
- `8OG.mol2` 中存在两个原子名为 `O6` —— atom 6（C8 旁，应该是 8-位 O8）和 atom 10（C6 旁，标准 G O6）。命名歧义可能影响 tleap 加载或 OpenMM 模板匹配。在 [`pipeline/convert_modxna_to_openmm.py`](../pipeline/convert_modxna_to_openmm.py)（待写）里需要重命名 atom 6 → `O8` 或验证 modXNA 自己的内部约定。

---

## Sarzyńska/Lahiri 2022 — pseudouridine 升级版

| 日期 | 来源 | 文件 | 备注 |
|------|------|------|------|
| 待下载 | DOI 10.1007/s10822-022-00447-4 SI | — | 仅在 modXNA Ψ 在 4BS2 上验证不通过时启用 |

---

## Bussi 2022 (Piomponi) — m6A 升级版

| 日期 | 来源 | 文件 | 备注 |
|------|------|------|------|
| 待下载 | github.com/bussilab/m6a-charge-fitting | `ff-m6a-fit5_AC.rtp` | GROMACS RTP 格式，需手工转 AMBER；仅在 modXNA m6A 验证不通过时启用 |
