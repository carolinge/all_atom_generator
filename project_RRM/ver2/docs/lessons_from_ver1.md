# ver_1 失败原因分析

> 来源：研究 agent 文献综述 + 直接审视 [`/app/ff_custom/8OG_RNA_amber.xml`](../../../app/ff_custom/8OG_RNA_amber.xml)
> 日期：2026-05-01

ver_1 的核心错误：**自己手写 OpenMM XML 来定义 8oxoG**，而不是使用已发表的 RESP-derived 参数集。

下面按可能性排序。

## 1. 新键对应的 bonded terms 大概率没补全（最可能）

H8 → O8 后引入了：

- 新键 **C8=O8**（sp2 羰基，~1.21 Å, k≈570 kcal/mol/Å²）
- 新键角 **N7-C8-O8** 和 **N9-C8-O8**
- 新 **improper torsion**（C8 sp2 平面性约束）
- 新二面角 **N9-C8=O8 / N7-C8=O8** 周围的 dihedral terms

XML 里你写了 `<Bond>` 但 OpenMM 的 angle / torsion / improper 是从底层 amber14 库的 `RNA-C / RNA-O / RNA-NA` **类型组合** 里抽。这些类型组合在标准 G 里不会同时出现（因为 G 没有 C8=O8），OL3 库里多半根本没有这条路径，结果是要么报错，要么静默用错值。

**症状**：minimization 阶段直接爆炸，或前 100 ps 内 NaN。

## 2. 总电荷大概率非整数

AMBER 残基的总电荷必须严格 −1.0（含磷酸）或 0.0（cap-neutral）。RESP 拟合时强制约束总电荷整数；**手工估算 + redistribute 几乎不可能精确**。

快速核算 ver_1 的 8OG 残基（电荷之和）：原 G 的 H8+C8+N7 ≈ 0.16 + 0.07 + (−0.61) = −0.38；ver_1 改成 O8+C8+N7+H7 = (−0.56) + 0.40 + (−0.46) + 0.35 = −0.27。**已多出 ~0.11 e**。

OpenMM 不会强制整数化，所以**不会立刻崩溃**，但会有长时间尺度下的离子云漂移、和 PME 的隐性误差。

**症状**：跑得动但 ns 之后行为奇怪。

## 3. 8-oxoG 的 6,8-双酮互变异构体处理不完整

真实 8oxoG 是 6,8-双酮形式（N1-H **和** N7-H 都存在）。ver_1 的 XML 在结构上确实加了 N7-H 和保留了 N1-H —— 这一步是对的。

**但电荷上**：相邻 N1/N3/C2/N2 的电荷应该跟着 RESP 一起重新拟合，因为 8-oxo 修饰改变了整个环的电子分布。ver_1 直接复用了 G 的旧值，会有微偏差（典型 0.05-0.15 e per atom）。

**症状**：base-pair 几何缓慢偏离实验结构；syn/anti 平衡比例错。

## 4. χ_OL3 二面角依赖 N9 原子类型

OL3 修正的 χ 二面角参数是按 N9 类型（`N*`）匹配的。ver_1 保留了 `RNA-N*`，所以 χ_OL3 还能 match 上 —— **这一项实际上没问题**。

记录在这里只是为了下次写新修饰时不要踩坑：如果给 N9 换了新类型，χ_OL3 会被静默丢弃，必须手动补 χ 项。

## 解决方案

**不要继续改自建 XML**。用 modXNA 提供的 RESP 参数 + tleap 组装路线（不是 OpenMM XML 路线）。架构和命名陷阱见 [`pipeline.md`](pipeline.md) 和 [`../force_fields/VERSIONS.md`](../force_fields/VERSIONS.md)。

## 第二轮发现的陷阱（2026-05-01 探索期间）

ver_2 启动后第一周发现的**新一类失败**——和 ver_1 失败模式无关，但同等致命：

### 5. 误把 modXNA 的 PSU 当成天然 pseudouridine

modXNA 是为药物化学（antisense oligonucleotides）设计的，命名不遵循天然修饰的 PDB 约定：
- `PSU.mol2` = **2-thio-5-isobutyl-pseudouridine**（设计衍生物，含 S8 + 烷基侧链），不是天然 ψ
- 天然 ψ 在 modXNA 里叫 **`PUU.mol2`**（HEAD01=C5，C-糖苷是 ψ 的定义特征）
- `M6A.mol2` = N6,N6-dimethyl，**不是 m6A**
- `DMA.mol2` = 2,8-dimethyl，跟 N6,N6-dimethyl-A（俗称 DMA）不一样

**症状**（如果没发现就跑了）：模拟会运行起来，但生物学结论完全错误——你以为在研究天然修饰，实际在跑设计药物分子。

**避免方法**：永远先看 base mol2 的 ATOM 段（特别是 HEAD01 和"额外原子"），而不是只信文件名。

### 6. modXNA base mol2 的电荷不能直接拼到 OL3 sugar/backbone

base mol2 里的 USER_CHARGES 是 fragment-level RESP——按 modXNA 的 base+sugar+backbone 模块化协议设计的。把这些电荷直接和 amber14/RNA.OL3.xml 的标准 sugar+phosphate 拼在一起，**总电荷不会是整数**（重复了 ver_1 第 2 项失败）。

**正确做法**：必须运行 `modxna.sh` 在 Linux + AmberTools 环境组装出 `.lib`，让 modXNA 的内部协议把连接处的电荷重分配到位。组装好的 .lib 总电荷 = -1.0（含磷酸）。

详见 [`../results/inspect_modxna.md`](../results/inspect_modxna.md) 的"After stripping"列——4 个修饰里没有一个直接拼出整数总电荷。

## 验证标准

任何新参数集投入生产前，必须通过：

1. **总电荷检查**：每个残基（含磷酸 −1.0；末端 0.0）
2. **bonded terms 完整性**：用 [`/app/ff_custom/8OG_RNA_amber.xml`](../../../app/ff_custom/8OG_RNA_amber.xml) 边的 [`/project_RRM/ver1/scripts/validate_8OG_ff.py`](../../ver1/scripts/) 风格的脚本（能 createSystem 不抛异常）
3. **free-base 100 ns**：单核苷在水中跑 100 ns，χ 二面角分布与 NMR/QM 参考一致（8oxoG 应该有 ~80% anti, ~20% syn）

ver_2 不重复 ver_1 的"先跑了再说"模式 —— **每个修饰先做 free-RNA 验证，再上 4BS2 复合物**。
