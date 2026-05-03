# Charge audit — modXNA fragment stitching residuals

Goal: explain the +0.32 e residual on internal 8OG that rebalance_lib.py was dumping onto N9. Per expert review:

> Typical RESP integer-rounding residuals are <0.05 e; 0.32 e suggests something was off in the modXNA fragment-stitching or capping scheme

## Backbone (RPO) and sugar (RC3) baseline

### RPO  (`P` head, `O5'` tail)
  raw_q = -0.999000 e   kept_q = -0.884700 e   (4 kept / 9 stripped)
  HEAD01STRIP atoms: ["O3'", 'C3', 'H31', 'H32', 'H33']
  TAIL01STRIP atoms: ["C5'", "H5'1", "H5'2", "H5'3"]

### RC3  (`C5'` head, `O3'` tail)
  raw_q = +0.109450 e   kept_q = -0.004700 e   (15 kept / 7 stripped)
  HEAD01STRIP: ['O1', 'H6']
  HEAD03STRIP: ['C6', 'H10', 'H11', 'H9']
  TAIL01STRIP: ['H7']

## Base fragments + assembly arithmetic

If the modXNA model is `internal residue charge = kept_backbone + kept_sugar + kept_base`, then the expected internal residue charge depends only on these three numbers. Target for ANY internal nucleotide = **-1.000000 e**.

| code | base raw_q | base kept_q | BB.kept + SU.kept + base.kept | residual vs −1.0 |
|------|------:|------:|------:|------:|
| `8OG` | -0.4875 | -0.7495 | -1.6389 | -0.6389 |
| `PUU` | -0.3598 | -0.4947 | -1.3841 | -0.3841 |
| `M1A` | -0.1708 | -0.3453 | -1.2347 | -0.2347 |

## 8OG detail

- Raw 8OG.mol2: 24 atoms, sum charge = -0.487472 e
- HEAD01STRIP removes 8 atoms, removed charge +0.261986 e
- Kept (base of 16 atoms): -0.749458 e

### Per-atom (kept) charges in 8OG base fragment:

| name | type | charge |
|------|------|------:|
| `N9` | `N*` | -0.311958 |
| `C8` | `CK` | +0.514195 |
| `O6` | `O` | -0.635575 |
| `N7` | `NA` | -0.537134 |
| `C5` | `CB` | -0.199516 |
| `C6` | `C` | +0.616353 |
| `O6` | `O` | -0.614804 |
| `N1` | `NA` | -0.639957 |
| `H1` | `H` | +0.432650 |
| `C2` | `CA` | +0.772412 |
| `N2` | `N2` | -0.986096 |
| `H21` | `H` | +0.420959 |
| `H22` | `H` | +0.420959 |
| `N3` | `NC` | -0.784574 |
| `C4` | `CB` | +0.368687 |
| `H7` | `H` | +0.413941 |

## Interpretation

If `BB.kept + SU.kept + base.kept` ≈ −1.0 e for 8OG, then modxna.sh's stripped charge corrections (`charge -0.8832`, `charge -0.01191`, `charge -0.10489`) would NOT be needed — the kept-charge sum is already close to integer. Our patch removed those corrections, which is correct IF the corrections were stale code.

If the residual is ~0.3 e for 8OG but ~0 for standard G (modrna08), then the 8OG base mol2's raw charges have been adjusted by modxna's authors to compensate for the removed `charge VAL` lines — but the compensation is on the WRONG atom, OR the mol2 has not been updated to match the broken modxna.sh.

Check this hypothesis: compare 8OG residual vs PUU and M1A residuals above. If all 3 are ~0.3 e, it's a consistent modxna issue. If only 8OG is off, it's specific to that fragment.
