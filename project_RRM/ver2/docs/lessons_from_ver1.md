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

**不要继续改自建 XML**。直接用 modXNA：

```bash
git clone https://github.com/modxna/modxna project_RRM/ver2/force_fields/modxna
```

modXNA 提供 `8OG.mol2` + `frcmod.modxna`，这是 RESP/HF/6-31G\* 拟合 + 完整 bonded terms 的 OL3-RNA 兼容参数集。

转换为 OpenMM XML 的步骤见 [`pipeline.md`](pipeline.md)。

## 验证标准

任何新参数集投入生产前，必须通过：

1. **总电荷检查**：每个残基（含磷酸 −1.0；末端 0.0）
2. **bonded terms 完整性**：用 [`/app/ff_custom/8OG_RNA_amber.xml`](../../../app/ff_custom/8OG_RNA_amber.xml) 边的 [`/project_RRM/ver1/scripts/validate_8OG_ff.py`](../../ver1/scripts/) 风格的脚本（能 createSystem 不抛异常）
3. **free-base 100 ns**：单核苷在水中跑 100 ns，χ 二面角分布与 NMR/QM 参考一致（8oxoG 应该有 ~80% anti, ~20% syn）

ver_2 不重复 ver_1 的"先跑了再说"模式 —— **每个修饰先做 free-RNA 验证，再上 4BS2 复合物**。
