# project_RRM — TDP-43 RRM + 修饰 RNA 全原子 MD

科学项目：研究 RNA 修饰对 TDP-43 RRM1+2（PDB **4BS2**）-RNA 结合的影响。

## 目标修饰

| 修饰 | 缩写 | 起始残基 | 状态 |
|------|------|---------|------|
| 8-oxoguanosine | 8OG / 8oxoG | G | 4BS2 主线（G→8OG） |
| Pseudouridine | PSU / Ψ | U | 4BS2 主线（U→Ψ） |
| N1-methyladenosine | M1A / m1A | A | **free-RNA 验证**（4BS2 无 A） |
| N6-methyladenosine | M6A / m6A | A | **free-RNA 验证**（4BS2 无 A） |

## 目录

- [`ver1/`](ver1/) — 2026-03 至 2026-04 的初次尝试（自建 8OG XML，6×100 ns replicas）。**已归档**。
- [`ver2/`](ver2/) — 当前主线：基于已发表力场（modXNA + Sarzyńska + Bussi）重做。
- [`environment_v2.yml`](environment_v2.yml) — ver_2 的 conda 环境定义（独立于根目录 `environment.yml`）。

## 与项目其他部分的关系

根目录的 [`app/`](../app/) 是独立的 PDB 准备 Web UI 工具，**不属于本科学项目**。两者不应互相依赖。
