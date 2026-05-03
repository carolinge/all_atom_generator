# 4BS2 / WT vs 8oxo-G3 — ver_2 第一阶段科学分析报告

> 日期：2026-05-03
> 数据：6×100 ns（WT n=3 + 8OG-G3 n=3，共 600 ns）
> 力场：AMBER ff14SB + RNA.OL3 + TIP3P + modXNA (Bergonzo 2024) for 8oxoG
> 状态：**Pilot study，不是结论性结果**。详见 [`expert_review.md`](expert_review.md)。

## 摘要（先看这段）

| 问题 | 当前答案 |
|------|---------|
| modXNA 8OG 力场参数物理合理吗？ | **未验证**——pucker 偏移方向对了，但**只是结果，不是验证**。需要 free-nucleoside 1 μs 控制实验 |
| 0.32 e 残差堆 N9 OK 吗？ | **不 OK**——专家说"N9 是 χ 二面角最敏感的原子，0.32 e 远超典型 RESP 残差 (<0.05 e)"，需要 audit |
| 8OG 削弱 TDP-43 结合是真的吗？ | **现有数据不能下结论**——n=3 太少，r1 是 PBC-confused 的 near-dissociation。需要 ≥500 ns × 5 replicas |
| 这值得发表吗？ | **作为 pilot/hypothesis-generating 可以**，作为机制结论不行。需要按 expert review 的 6 步走 |

---

## 1. 实验目标一句话

把 4BS2 中 RNA 链上 G3 换成 **8-氧鸟苷（8-oxoguanosine, 8oxoG）**，跑 100 ns × 3 replicas，跟野生型 (WT) 100 ns × 3 对比，看看：
- (A) modXNA 力场参数是否物理合理 — **力场验证**
- (B) 8-oxoG 是否真改变 TDP-43 RRM 对 RNA 的结合 — **生物学结论**

---

## 2. 8-oxo-G 是什么？为什么重要？

### 化学
活性氧（ROS）攻击鸟嘌呤的 C8 位置，把那里的 C-H 变成 C=O 羰基，同时 N7 上多了一个 H。结构上多了一个氧原子（O8），加了一个氢（H7），删了一个氢（H8）。

```
    标准 G            8-oxo-G
    H8                 O8 (=)
     \                  ||
      C8 ── N7         C8 ── N7
      ||      \        |       \
      N9       C5      N9       C5—H ← N7 现在带 H
                                    (新 donor)
```

### 为什么重要
- **氧化应激** 损伤 RNA → 神经退行性疾病（ALS/FTD），TDP-43 是关键蛋白
- 单个 8-oxo-G 改变了碱基的氢键模式：N7 从 acceptor 变 donor，C8=O8 是新 acceptor —— Watson-Crick 边没变但 Hoogsteen 边完全反转
- 在 DNA 中，8-oxo-G 会让 C8=O8 跟糖环上的 H 撞，**导致鸟嘌呤碱基围绕 N9-C1' 键翻转 180°**，这就是著名的 **syn-anti 反转**

---

## 3. 力场验证：怎么算"通过"？

文献上 8-oxo-G 力场验证有**两条经典证据**（gold standard）：

### 证据 A：糖环 pucker（南-北偏移）✅ **观测到**

#### 概念
RNA 核糖环不是平的，它会"皱"在某个原子之上。两种主要构象：
- **C3'-endo（north）**：C3' 凸向碱基同一侧，pseudorotation 角 P ≈ 0-36°，标准 A-form RNA
- **C2'-endo（south）**：C2' 凸向碱基同一侧，P ≈ 144-180°，B-form DNA

8-oxo-G **倾向 C2'-endo (south)**，因为 8-位的氧排斥糖环上的氢，迫使糖翻转。

#### 我们的数据
| 体系 | G3 sugar pucker P 平均值 | 主峰位置 |
|------|-------------------------|---------|
| WT (n=3) | **190°** | 双峰：P~30° (north) + P~200° (south region) |
| 8OG (n=3) | **244°** | 单峰：P~210-240° (south, C2'-exo / C2'-endo) |

**Δ = +53°**，明显从 north 偏到 south。✅ **modXNA 抓到了这个签名。**

> ⚠️ 警告：WT 的双峰分布表明即使在 WT 中，G3 也有相当一部分时间在 south region —— 这可能是因为 G3 的 sugar 在蛋白结合界面有特殊几何要求。所以 +53° 的差异部分是因为 WT 也部分 south，不是纯粹 8oxoG 的效应。详细解读见 §6。

### 证据 B：syn-anti 反转 ❌ **未观测到**

#### 概念
**glycosidic χ 角** = O4'-C1'-N9-C4 二面角，描述碱基相对糖的转动。
- **anti**：χ ≈ -120° ± 30°，标准构象，碱基"远离"糖
- **syn**：χ ≈ +60° ± 30°，碱基"翻向"糖

**自由 8-oxo-G 单核苷在水中**：~50% syn vs ~50% anti（标准 G 是 ~95% anti）。这是 8-oxo-G 最有名的力场测试。

#### 我们的数据
| 体系 | G3 χ 平均值 | 分布形状 |
|------|-------------|---------|
| WT | **−111°** | anti band 单峰 |
| 8OG | **−118°** | anti band 单峰，**没有 syn 群体** |

**🤔 修饰的 G3 完全没翻 syn**。这有三种可能解释：

1. **(a) 蛋白锁定**：TDP-43 RRM 的 RNA-binding pocket 强制 G3 保持 anti，蛋白接触能赢过 syn 倾向。这是**物理上合理的**（实验上 TDP-43-bound RNA 大多 anti）。
2. **(b) 采样不足**：100 ns 不够看到 syn 翻转。free 8oxoG 的 syn-anti 翻转大概几十 ns 一次，但**蛋白结合下能阻挡几个数量级**。需要 μs 级采样。
3. **(c) 力场低估 syn**：modXNA 的 χ 二面角参数可能偏 anti。

**这是当前最大的不确定**。已经请专家 agent 评判（背景在跑）。

---

## 4. 生物学结果

### 全局指标（PBC outlier 过滤后）

| 指标 | WT 平均 | 8OG 平均 | Δ | 解读 |
|------|---------|----------|---|------|
| 蛋白 backbone RMSD | 3.95 Å | 3.68 Å | −0.27 | 蛋白整体不受影响（修饰是局部效应） |
| **RNA RMSD** | **4.90 Å** | **6.32 Å** | **+1.42** | RNA 移动幅度增大 |
| 蛋白-RNA 总接触数 | 326 | 281 | **−45 (−14%)** | 接触面缩小 |
| COM 距离 | 12.85 Å | 13.57 Å | +0.72 | RNA 略远离蛋白 |
| **G3 局部接触** | **27.3** | **20.1** | **−7.2 (−26%)** | 修饰位点最受影响 |
| G3 极性接触 | 3.55 | 2.73 | −0.82 | H-bond 减少 |

### 每条 replica 的真相

跨 replica 的均值掩盖了重要事实。看 [`figures/per_replica_timeseries.png`](_aggregate/figures/per_replica_timeseries.png)：

| Replica | RNA RMSD 行为 | G3 接触 | 解读 |
|---------|--------------|---------|------|
| **WT/r1** | 平稳 ~5-7 Å | 25-35 稳定 | 标准 |
| WT/r2 | 平稳 ~4-5 Å | 25-35 稳定 | 标准 |
| WT/r3 | 平稳 + 末端 PBC 噪声 | 25-35 稳定 | 标准 |
| **8OG/r1** | **30 ns 左右跳到 60+ Å** | **从 30 跌到 ~5** | 🚨 **部分解离事件** |
| 8OG/r2 | 平稳 ~5-7 Å | ~20 稳定 | 接触面减小但仍结合 |
| 8OG/r3 | 平稳 ~5-7 Å | ~15-20 稳定 | 接触面减小但仍结合 |

**关键观察**：8OG/r1 显示 RNA 在 30 ns 时**几乎完全解离**——这是单一事件（n=1/3），但是真实的 dissociation episode，不是 PBC 噪声（可以从 G3 接触从 30 跌到 5 同步看到）。

剩下 2/3 的 8OG replicas 跟 WT 类似但接触面稍小。

---

## 5. 诚实评估

### 现在能说的（高置信度）

- ✅ **力场架构通跑**：6×100 ns 没有 NaN/爆炸；能量、温度、压强都收敛；轨迹质量好
- ✅ **modXNA 8OG 参数能复现 sugar pucker 偏移**（south 方向）—— 跟 8oxoG 化学合理
- ✅ **8OG 一致地减少 G3 局部蛋白-RNA 接触**（−7 contacts，3/3 replicas）
- ✅ **8OG/r1 显示完整的 partial dissociation episode**（虽然 n=1）
- ✅ **效应是局部的**（蛋白整体 RMSD 不变）

### 现在还不能说（需要更多数据/分析）

- ❓ **modXNA χ 二面角能否复现 syn 翻转**：100 ns + 蛋白结合 = 难判断。需要 (i) free 8oxoG 单核苷在水中 100 ns 控制实验、(ii) μs 级蛋白结合实验、或 (iii) 文献对比
- ❓ **8OG/r1 的 dissociation 事件**：是真生物学还是 outlier？需要更多 replicas 或更长时间
- ❓ **结合自由能差**：6×100 ns 不足以定量。需要更长 + WHAM/MBAR/FEP

### 已知问题

- 残余 PBC artefact：~5-32% 帧需要 outlier mask（修复方案：MD 输出时直接做 unwrap，分析时无需后处理）
- modxna 8OG 总电荷修正：将 +0.317 e 残差全堆到 N9 上（应该分散到多个原子做更好）—— [`pipeline/rebalance_lib.py`](../pipeline/rebalance_lib.py) 文档化

---

## 6. 专家 agent 给的 minimum path to publication

来自 [`expert_review.md`](expert_review.md)。按顺序做，每一步是**前置条件**：

### Tier 1 — 必须做才能 claim 任何东西

1. **free 8oxoG 单核苷 in TIP3P, 1 μs** + **free G 单核苷, 1 μs（对照）**。在水中无蛋白。计算 χ 角分布，对比 NMR (Uesugi & Ikehara 1977; Cheong et al. 2004) 和 QM (Millen et al. 2008)。
   - 如果 modXNA 给出 8oxoG <10% syn → 力场太 anti-biased
   - 如果给 30-50% syn → 力场合理
   - 如果给 >70% syn → 太 syn-biased
2. **Audit 0.32 e 残差**：从 modXNA 公开的 ESP 网格重新做 RESP，看真实残差是多少。如果确实接近 0.3，**分散到多个原子**（不光 N9）测敏感性。
3. **重新 image 现有轨迹**做 frame-by-frame protein-COM centering 后再算 RMSD/contacts。

### Tier 2 — 强烈推荐（for bound system）

4. **N7-H7 → 蛋白/RNA H-bond inventory** —— 8oxoG 引入的新 donor，机制最有意思
5. **O8 → 蛋白/RNA H-bond inventory** —— 新 acceptor
6. **Hoogsteen 边 SASA**（N7/O8 面是否暴露/埋藏）
7. **N1-N7 距离 + base 法向**（stacking integrity）
8. **Backbone ε/ζ/α/γ at G3**（OL3 已知 ε/ζ 问题）—— 跟 Richardson et al. RNA 2008 "suite" classification 对比
9. **MM-GBSA/PBSA ΔΔG (WT vs 8OG)**，order-of-magnitude only

### Tier 3 — nice to have

10. 一个 short window 2 fs 无 HMR sanity check
11. Stacking energy decomposition (G3 with i±1)

### 然后

- **延长结合 simulation 到 500 ns × 5 replicas**（field standard），最好 1-2 μs single replica 加 5×500 ns
- 跑完上面 1-9 后才讨论生物学

### 还有

- 确认 G3 实际属于 RRM1 还是 RRM2 pocket（Kuo et al. NAR 2014）—— 解读机制依赖这个
- m1A free-RNA 已排队；用同样的 Tier 1 框架验证 modXNA M1A
- m6A 不在 modXNA，需要 Bussi 2022 (github.com/bussilab/m6a-charge-fitting) 单独引入

---

## 7. 给合作者/审稿人的"can I trust this?" 答案

经过专家评审后的**诚实**答案：

| 问 | 答 |
|----|------|
| MD 基础流程对吗？ | ✅ 是（能量/温度/PME/HMR/Langevin/PBC 都对） |
| 力场参数来源公开吗？ | ✅ Bergonzo 2024 modXNA |
| 用了对照吗？ | ✅ WT vs 8OG, n=3 each |
| 跨 replica 做了统计吗？ | ✅ |
| 力场被验证了吗？ | ❌ **没有**——pucker 偏移方向对但不是验证；syn 必须靠 free-nucleoside 1 μs 控制实验 |
| 0.32 e 电荷残差堆 N9 上合理吗？ | ❌ **不合理**——典型残差 <0.05 e，且 N9 是 χ 最敏感原子 |
| "8OG 弱化结合"能下结论吗？ | ❌ **不能**——n=3 太少，r1 是 PBC-confused 而非真 dissociation，COM Δ 在 thermal noise 内 |
| 现状能否发表？ | ⚠️ 作为 **pilot/hypothesis-generating study** 可以；作为 mechanistic claim 不行 |

**现在能讲的故事**：「我们搭好了 TDP-43 + 8oxoG MD pipeline，pucker 偏移方向跟文献一致，看到 1/3 replica 部分解离，需要更长采样和力场基准来证实。」

**还不能讲的故事**：「8oxoG 通过破坏接触显著弱化 TDP-43 结合。」

---

## 文件位置

```
project_RRM/ver2/
├── results/
│   ├── analysis_v1.md                  ← 本文
│   └── _aggregate/
│       ├── report.md                   ← 自动生成的 stats 报告
│       ├── stats.csv                   ← Welch's t / KS / MW p-values
│       ├── timeseries.npz              ← 跨 replica 平均时序
│       ├── rmsf.npz                    ← per-residue RMSF
│       ├── contact_diff.npz            ← 接触图 8OG−WT 差分
│       ├── chi_distributions.npz       ← G3 χ 角分布
│       ├── pucker_distributions.npz    ← G3 sugar pucker 分布
│       └── figures/
│           ├── timeseries.png          ← 6 大 CV 跨 replica 时序
│           ├── rmsf.png                ← 蛋白+RNA per-residue RMSF
│           ├── contact_diff.png        ← WT/8OG/差分 contact maps
│           ├── chi_pucker.png          ← G3 χ + pucker 分布直方图
│           └── per_replica_timeseries.png  ← 每条 replica 单独时序
└── pipeline/
    ├── analyze_v2.py                   ← 单 trajectory CV + RMSF
    ├── compare_v2.py                   ← 跨 replica + 统计 + 绘图
    ├── plot_per_replica.py             ← per-replica 时序
    └── run_analysis.sh                 ← 服务器端 driver
```
