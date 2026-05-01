# modXNA Base-Fragment Inspection Report

Source: `project_RRM\ver2\force_fields\modxna` (commit recorded in [`VERSIONS.md`](../force_fields/VERSIONS.md))

Targets: 8OG, PUU, M1A

## Verification questions

1. Do the 4 base mol2 files exist and parse cleanly?
2. What's the **net charge contribution** of the base (after STRIP atoms removed)? In a properly assembled OL3-RNA nucleotide the base should contribute **~+0.00 to +0.40 e** (the rest of the −1.0 sits on the phosphate + sugar).
3. Are there any **duplicate atom names** that would break tleap or OpenMM template matching?
4. Are there any **atom types outside standard OL3** that would need new definitions in the OpenMM XML?

## 8OG  —  `8OG.mol2`

- Header: `modXNA-fragment:8OG:HEAD01:@N9:HEAD01STRIP:@C1',H1',O4',HHO4,C2',H2'1,H2'3,H2'2`
- Total atoms: **24**
- Total charge (all atoms): **-0.487472 e**
- HEAD01 atom (connects to sugar): `N9`
- HEAD01STRIP atoms (dropped on assembly): `C1', H1', O4', HHO4, C2', H2'1, H2'3, H2'2`
- After stripping (8 atoms removed): **16 atoms** remain, sum charge = **-0.749458 e**
- ⚠️ **Duplicate atom names**: {'O6': 2}
    - `O6`: idx=6 type=O q=-0.6356; idx=10 type=O q=-0.6148
- Atom types used: `C, CA, CB, CK, CT, H, H2, HC, HO, N*, N2, NA, NC, O, OH`
- All types are in standard OL3 set ✅

| idx | name | type | charge | strip? |
|----:|------|------|-------:|--------|
|   1 | `O4'` | `OH` | -0.755712 | STRIP |
|   2 | `C1'` | `CT` | +0.503913 | STRIP |
|   3 | `H1'` | `H2` | +0.043368 | STRIP |
|   4 | `N9` | `N*` | -0.311958 |  |
|   5 | `C8` | `CK` | +0.514195 |  |
|   6 | `O6` | `O` | -0.635575 |  |
|   7 | `N7` | `NA` | -0.537134 |  |
|   8 | `C5` | `CB` | -0.199516 |  |
|   9 | `C6` | `C` | +0.616353 |  |
|  10 | `O6` | `O` | -0.614804 |  |
|  11 | `N1` | `NA` | -0.639957 |  |
|  12 | `H1` | `H` | +0.432650 |  |
|  13 | `C2` | `CA` | +0.772412 |  |
|  14 | `N2` | `N2` | -0.986096 |  |
|  15 | `H21` | `H` | +0.420959 |  |
|  16 | `H22` | `H` | +0.420959 |  |
|  17 | `N3` | `NC` | -0.784574 |  |
|  18 | `C4` | `CB` | +0.368687 |  |
|  19 | `C2'` | `CT` | -0.392425 | STRIP |
|  20 | `H2'1` | `HC` | +0.131147 | STRIP |
|  21 | `HHO4` | `HO` | +0.469401 | STRIP |
|  22 | `H2'2` | `HC` | +0.131147 | STRIP |
|  23 | `H2'3` | `HC` | +0.131147 | STRIP |
|  24 | `H7` | `H` | +0.413941 |  |

## PUU  —  `PUU.mol2`

- Header: `modXNA-fragment:PUU:HEAD01:@C5:HEAD01STRIP:@C1',H1',O4',HHO4,C2',H2'1,H2'3,H2'2`
- Total atoms: **19**
- Total charge (all atoms): **-0.359769 e**
- HEAD01 atom (connects to sugar): `C5`
- HEAD01STRIP atoms (dropped on assembly): `C1', H1', O4', HHO4, C2', H2'1, H2'3, H2'2`
- After stripping (8 atoms removed): **11 atoms** remain, sum charge = **-0.494733 e**
- Atom name uniqueness: ✅ all 19 unique
- Atom types used: `C, CM, CT, H, H1, H4, HC, HO, NA, O, OH`
- All types are in standard OL3 set ✅

| idx | name | type | charge | strip? |
|----:|------|------|-------:|--------|
|   1 | `O4'` | `OH` | -0.695231 | STRIP |
|   2 | `C1'` | `CT` | +0.392990 | STRIP |
|   3 | `H1'` | `H1` | -0.003632 | STRIP |
|   4 | `C5` | `CM` | -0.330057 |  |
|   5 | `C6` | `CM` | +0.094238 |  |
|   6 | `H6` | `H4` | +0.172312 |  |
|   7 | `N1` | `NA` | -0.547392 |  |
|   8 | `H1` | `H` | +0.364333 |  |
|   9 | `C2` | `C` | +0.628142 |  |
|  10 | `O2` | `O` | -0.622147 |  |
|  11 | `N3` | `NA` | -0.622056 |  |
|  12 | `H3` | `H` | +0.361837 |  |
|  13 | `C4` | `C` | +0.607616 |  |
|  14 | `O4` | `O` | -0.601559 |  |
|  15 | `C2'` | `CT` | -0.253660 | STRIP |
|  16 | `H2'1` | `HC` | +0.096528 | STRIP |
|  17 | `HHO4` | `HO` | +0.404913 | STRIP |
|  18 | `H2'2` | `HC` | +0.096528 | STRIP |
|  19 | `H2'3` | `HC` | +0.096528 | STRIP |

## M1A  —  `M1A.mol2`

- Header: `modXNA-fragment:M1A:HEAD01:@N9:HEAD01STRIP:@C1',H1',O4',HHO4,C2',H2'1,H2'3,H2'2`
- Total atoms: **26**
- Total charge (all atoms): **-0.170762 e**
- HEAD01 atom (connects to sugar): `N9`
- HEAD01STRIP atoms (dropped on assembly): `C1', H1', O4', HHO4, C2', H2'1, H2'3, H2'2`
- After stripping (8 atoms removed): **18 atoms** remain, sum charge = **-0.345302 e**
- Atom name uniqueness: ✅ all 26 unique
- Atom types used: `CA, CB, CK, CM, CT, H, H1, H2, H5, HC, HO, N*, N2, NB, NC, OH`
- All types are in standard OL3 set ✅

| idx | name | type | charge | strip? |
|----:|------|------|-------:|--------|
|   1 | `O4'` | `OH` | -0.694110 | STRIP |
|   2 | `C1'` | `CT` | +0.500029 | STRIP |
|   3 | `H1'` | `H2` | -0.007277 | STRIP |
|   4 | `N9` | `N*` | -0.058508 |  |
|   5 | `C8` | `CK` | +0.033658 |  |
|   6 | `H8` | `H5` | +0.161421 |  |
|   7 | `N7` | `NB` | -0.568384 |  |
|   8 | `C5` | `CB` | +0.113204 |  |
|   9 | `C6` | `CM` | +0.286708 |  |
|  10 | `N6` | `N2` | -1.019819 |  |
|  11 | `H61` | `H` | +0.441095 |  |
|  12 | `H62` | `H` | +0.441095 |  |
|  13 | `N1` | `N*` | -0.050133 |  |
|  14 | `C2` | `CA` | +0.013696 |  |
|  15 | `H2` | `H5` | +0.123001 |  |
|  16 | `N3` | `NC` | -0.651883 |  |
|  17 | `C4` | `CB` | +0.187216 |  |
|  18 | `C2'` | `CT` | -0.486099 | STRIP |
|  19 | `H2'1` | `HC` | +0.153265 | STRIP |
|  20 | `HHO4` | `HO` | +0.402202 | STRIP |
|  21 | `H2'2` | `HC` | +0.153265 | STRIP |
|  22 | `H2'3` | `HC` | +0.153265 | STRIP |
|  23 | `C11` | `CT` | -0.107662 |  |
|  24 | `1H11` | `H1` | +0.103331 |  |
|  25 | `2H11` | `H1` | +0.103331 |  |
|  26 | `3H11` | `H1` | +0.103331 |  |

---

# Reference: misleadingly-named modXNA codes

These are **NOT** what their name suggests in modXNA's library. Inspecting them here so the gotchas are documented in one place.

## PSU  —  `PSU.mol2`

- Header: `modXNA-fragment:PSU:HEAD01:@N1:HEAD01STRIP:@C1',H1',C2',H2'1,H2'2,H2'3,O4',HHO4`
- Total atoms: **25**
- Total charge (all atoms): **-0.004443 e**
- HEAD01 atom (connects to sugar): `N1`
- HEAD01STRIP atoms (dropped on assembly): `C1', H1', C2', H2'1, H2'2, H2'3, O4', HHO4`
- After stripping (8 atoms removed): **17 atoms** remain, sum charge = **-0.104937 e**
- Atom name uniqueness: ✅ all 25 unique
- Atom types used: `C, CA, CM, CT, H, H2, HA, HC, HO, N*, NA, O, OH, S`
- All types are in standard OL3 set ✅

| idx | name | type | charge | strip? |
|----:|------|------|-------:|--------|
|   1 | `C1'` | `CT` | +0.308415 | STRIP |
|   2 | `C2'` | `CT` | -0.183129 | STRIP |
|   3 | `O4'` | `OH` | -0.674134 | STRIP |
|   4 | `O4` | `O` | -0.579212 |  |
|   5 | `C4` | `C` | +0.494158 |  |
|   6 | `N3` | `NA` | -0.105296 |  |
|   7 | `C2` | `CA` | -0.153171 |  |
|   8 | `S8` | `S` | -0.378105 |  |
|   9 | `N1` | `N*` | +0.375959 |  |
|  10 | `C6` | `CM` | -0.086027 |  |
|  11 | `C5` | `CM` | -0.302524 |  |
|  12 | `H1'` | `H2` | +0.008190 | STRIP |
|  13 | `H2'1` | `HC` | +0.073239 | STRIP |
|  14 | `H2'2` | `HC` | +0.073239 | STRIP |
|  15 | `H3` | `H` | +0.246952 |  |
|  16 | `C16` | `CT` | -0.039382 |  |
|  17 | `H5` | `HA` | +0.208170 |  |
|  18 | `H2'3` | `HC` | +0.073239 | STRIP |
|  19 | `HHO4` | `HO` | +0.421435 | STRIP |
|  20 | `H161` | `HC` | +0.083272 |  |
|  21 | `H162` | `HC` | +0.083272 |  |
|  22 | `C10` | `CT` | -0.024265 |  |
|  23 | `H101` | `HC` | +0.023754 |  |
|  24 | `H102` | `HC` | +0.023754 |  |
|  25 | `H103` | `HC` | +0.023754 |  |

## M6A  —  `M6A.mol2`

- Header: `modXNA-fragment:M6A:HEAD01:@N9:HEAD01STRIP:@C1',H1',O4',HHO4,C2',H2'1,H2'3,H2'2`
- Total atoms: **28**
- Total charge (all atoms): **-0.370738 e**
- HEAD01 atom (connects to sugar): `N9`
- HEAD01STRIP atoms (dropped on assembly): `C1', H1', O4', HHO4, C2', H2'1, H2'3, H2'2`
- After stripping (8 atoms removed): **20 atoms** remain, sum charge = **-0.517566 e**
- Atom name uniqueness: ✅ all 28 unique
- Atom types used: `CA, CB, CK, CQ, CT, H1, H2, H5, HC, HO, N*, N2, NB, NC, OH`
- All types are in standard OL3 set ✅

| idx | name | type | charge | strip? |
|----:|------|------|-------:|--------|
|   1 | `O4'` | `OH` | -0.733550 | STRIP |
|   2 | `C1'` | `CT` | +0.312926 | STRIP |
|   3 | `H1'` | `H2` | +0.056955 | STRIP |
|   4 | `N9` | `N*` | -0.096204 |  |
|   5 | `C8` | `CK` | +0.127564 |  |
|   6 | `H8` | `H5` | +0.172977 |  |
|   7 | `N7` | `NB` | -0.528158 |  |
|   8 | `C5` | `CB` | -0.079141 |  |
|   9 | `C6` | `CA` | +0.489502 |  |
|  10 | `N6` | `N2` | -0.193297 |  |
|  11 | `C11` | `CT` | -0.106112 |  |
|  12 | `C12` | `CT` | -0.106112 |  |
|  13 | `N1` | `NC` | -0.796893 |  |
|  14 | `C2` | `CQ` | +0.500405 |  |
|  15 | `H2` | `H5` | +0.036578 |  |
|  16 | `N3` | `NC` | -0.803603 |  |
|  17 | `C4` | `CB` | +0.380194 |  |
|  18 | `C2'` | `CT` | -0.345216 | STRIP |
|  19 | `H2'1` | `HC` | +0.127451 | STRIP |
|  20 | `HHO4` | `HO` | +0.473360 | STRIP |
|  21 | `H2'2` | `HC` | +0.127451 | STRIP |
|  22 | `H2'3` | `HC` | +0.127451 | STRIP |
|  23 | `1H11` | `H1` | +0.080789 |  |
|  24 | `2H11` | `H1` | +0.080789 |  |
|  25 | `3H11` | `H1` | +0.080789 |  |
|  26 | `1H12` | `H1` | +0.080789 |  |
|  27 | `2H12` | `H1` | +0.080789 |  |
|  28 | `3H12` | `H1` | +0.080789 |  |

## DMA  —  `DMA.mol2`

- Header: `modXNA-fragment:DMA:HEAD01:@N9:HEAD01STRIP:@C1',H1',O4',HHO4,C2',H2'1,H2'3,H2'2`
- Total atoms: **28**
- Total charge (all atoms): **-0.401368 e**
- HEAD01 atom (connects to sugar): `N9`
- HEAD01STRIP atoms (dropped on assembly): `C1', H1', O4', HHO4, C2', H2'1, H2'3, H2'2`
- After stripping (8 atoms removed): **20 atoms** remain, sum charge = **-0.559160 e**
- ⚠️ **Duplicate atom names**: {'C6': 2, 'H61': 2, 'H62': 2}
    - `C6`: idx=6 type=CT q=-0.3035; idx=9 type=CA q=+0.5850
    - `H61`: idx=11 type=H q=+0.4186; idx=23 type=HC q=+0.1269
    - `H62`: idx=12 type=H q=+0.4186; idx=24 type=HC q=+0.1269
- Atom types used: `CA, CB, CQ, CT, H, H2, HC, HO, N*, N2, NB, NC, OH`
- All types are in standard OL3 set ✅

| idx | name | type | charge | strip? |
|----:|------|------|-------:|--------|
|   1 | `O4'` | `OH` | -0.742185 | STRIP |
|   2 | `C1'` | `CT` | +0.393143 | STRIP |
|   3 | `H1'` | `H2` | +0.040067 | STRIP |
|   4 | `N9` | `N*` | -0.180926 |  |
|   5 | `C8` | `CA` | +0.420766 |  |
|   6 | `C6` | `CT` | -0.303492 |  |
|   7 | `N7` | `NB` | -0.635399 |  |
|   8 | `C5` | `CB` | -0.051779 |  |
|   9 | `C6` | `CA` | +0.584952 |  |
|  10 | `N6` | `N2` | -0.948703 |  |
|  11 | `H61` | `H` | +0.418569 |  |
|  12 | `H62` | `H` | +0.418569 |  |
|  13 | `N1` | `NC` | -0.841661 |  |
|  14 | `C2` | `CQ` | +0.668380 |  |
|  15 | `C15` | `CT` | -0.338232 |  |
|  16 | `N3` | `NC` | -0.769874 |  |
|  17 | `C4` | `CB` | +0.265903 |  |
|  18 | `C2'` | `CT` | -0.368648 | STRIP |
|  19 | `H2'1` | `HC` | +0.124958 | STRIP |
|  20 | `HHO4` | `HO` | +0.460541 | STRIP |
|  21 | `H2'2` | `HC` | +0.124958 | STRIP |
|  22 | `H2'3` | `HC` | +0.124958 | STRIP |
|  23 | `H61` | `HC` | +0.126888 |  |
|  24 | `H62` | `HC` | +0.126888 |  |
|  25 | `H63` | `HC` | +0.126888 |  |
|  26 | `1H15` | `HC` | +0.117701 |  |
|  27 | `2H15` | `HC` | +0.117701 |  |
|  28 | `3H15` | `HC` | +0.117701 |  |

