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
- `dat/lib_base/{8OG,PUU,M1A}.mol2` — 我们 3 个目标修饰的 base 片段（修正后）
- `dat/lib_backbone/RPO.mol2` — RNA 磷酸 backbone
- `dat/lib_sugar/RC3.mol2` — RNA ribose sugar
- `dat/frcmod.modxna` — 共享 bond/angle/torsion 参数
- `modxna.sh` — bash 主脚本（需要 cpptraj+tleap+sander，Linux only）
- `Catalog – ModXNA.2026_03_05.pdf` — 完整目录（2 MB）

**已发现的命名陷阱（吃过亏的，写在前面）**：

| modXNA 代码 | 实际是 | 我们项目里**该用**的 | 验证方法 |
|------|------|-----|-----|
| `8OG` | 8-oxoguanosine ✅ | `8OG` | 16 个 base 原子，O8（被命名为 `O6`）+ N7-H7 |
| `PSU` | 2-thio-5-isobutyl-pseudouridine（设计衍生物） | `PUU`（HEAD01=C5） | PSU.mol2 含 S8 + C16/C10 烷基 → 不是天然 ψ |
| `PUU` | **天然 pseudouridine** ✅ | `PUU` | HEAD01=C5（C-糖苷是 ψ 的定义特征） |
| `M1A` | N1-methyladenosine ✅ | `M1A` | 18 个 base 原子，N1-C11 甲基 |
| `M6A` | N6,N6-dimethyladenosine | **不要用** | M6A.mol2 含 C11+C12 双甲基挂 N6，无 N6-H |
| `DMA` | 2,8-dimethyladenosine | **不要用** | C2-Me（C15）+ C8-Me（C6） |

**8OG 内部命名陷阱**：`8OG.mol2` 中两个原子都叫 `O6` —— atom 6（C8 旁，实际是 8-位 O8）和 atom 10（C6 旁，标准 G 的 O6）。modxna.sh 组装时会按内部约定处理，但任何**不经过 modxna.sh** 直接读 mol2 的脚本都需要先把 atom 6 重命名为 `O8`。

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
