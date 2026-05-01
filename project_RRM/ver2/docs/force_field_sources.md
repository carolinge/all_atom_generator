# 力场来源记录

> 每次更新某个修饰的参数时，更新对应表格行 + commit hash。

## 总览

**好消息**：modXNA 一站式覆盖全部 4 个修饰（已 clone + 验证 `dat/lib_base/{8OG,PSU,M1A,M6A}.mol2` 全部存在）。
后两栏的"升级路线"仅在 modXNA 在 4BS2 上验证不通过时才启用。

| 修饰 | 主用 | 升级路线（仅当 modXNA 不收敛时） | 状态 |
|------|------|-----|------|
| 8oxoG | modXNA `8OG.mol2` + `frcmod.modxna` | — | clone OK，待 OpenMM XML 转换 |
| Ψ (PSU) | modXNA `PSU.mol2` + `frcmod.modxna` | Sarzyńska/Lahiri 2022 (DOI 10.1007/s10822-022-00447-4) SI | clone OK，待 OpenMM XML 转换 |
| m1A | modXNA `M1A.mol2` + `frcmod.modxna` | Xu/MacKerell 2016 (CHARMM 路线，DOI 10.1002/jcc.24307) | clone OK，待 OpenMM XML 转换 |
| m6A | modXNA `M6A.mol2` + `frcmod.modxna` | Bussi `ff-m6a-fit5_AC.rtp` (github.com/bussilab/m6a-charge-fitting) | clone OK，待 OpenMM XML 转换 |

modXNA 的 **commit hash 和 license** 见 [`../force_fields/VERSIONS.md`](../force_fields/VERSIONS.md)。

## 历史基线

**Aduri et al. 2007** (J. Chem. Theory Comput. 3, 1464; DOI 10.1021/ct600329w)
覆盖 MODOMICS 上 107 个 RNA 修饰，包括本项目全部 4 个。
RESP/HF/6-31G(d) 电荷合法，但 bonded terms 主要靠类比 GAFF —— 这是后人重做的原因。
SI 文件历史发布途径：ACS 论文 SI；现亦打包在 modXNA 中。

## CHARMM 备用

**Xu/MacKerell 2016** (J. Comput. Chem. 37, 896; DOI 10.1002/jcc.24307)
112 个修饰核苷的 CHARMM36 参数集，对 m1A 是金标准（AMBER 端没有同水平的）。
若 AMBER 路线在 m1A 上不收敛，考虑切到 CHARMM。

## 版本固定原则

每次 git clone 力场仓库后，**记录 commit hash 到 [`force_fields/VERSIONS.md`](../force_fields/VERSIONS.md)**（待建）。
论文 SI 文件下载后存到 [`force_fields/papers_SI/`](../force_fields/papers_SI/)，并附上下载日期 + 论文 DOI。
