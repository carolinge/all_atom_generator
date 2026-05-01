# 力场来源记录

> 每次更新某个修饰的参数时，更新对应表格行 + commit hash。

## 总览

**修正后状态**（详见 [`../results/inspect_modxna.md`](../results/inspect_modxna.md)）：

| 修饰 | 主用 | 备注 |
|------|------|------|
| 8oxoG | **modXNA `8OG.mol2`** + `frcmod.modxna` | ✅ 标准结构（16 个 base 原子，含 O8 + N7-H7） |
| Ψ (pseudouridine) | **modXNA `PUU.mol2`** ⚠️ 不是 PSU | ✅ 标准结构（HEAD01=C5，C-糖苷连接） |
| m1A | **modXNA `M1A.mol2`** | ✅ 标准结构（N1 上接甲基 C11） |
| m6A | **Bussi 2022** `ff-m6a-fit5_AC.rtp` | ❌ modXNA 不含标准 m6A，必须外部 |

⚠️ **关键陷阱**：modXNA 是为药物化学（ASO/siRNA）设计的，命名约定与天然 RNA 修饰不一致：
- modXNA `PSU` = 2-thio-5-isobutyl 的设计 ψ 衍生物（**不是天然 ψ**）
- modXNA `M6A` = N6,N6-dimethyladenosine（**不是天然 m6A**）
- modXNA `DMA` = 2,8-dimethyladenosine
- 天然 ψ 在 modXNA 里叫 **PUU**

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
