# External reference data

> **Provenance for the Cunning 2024 CBASS ED50 data** · Updated 2026-09-03 · Index: [`README.md`](../../README.md) · used by `code/sensitivity/26_thermal_context.R`.

## `cunning2024_apulchra_ed50.csv`

Per-genet acute thermal-tolerance thresholds (Fv/Fm **ED50**, °C) for *Acropora
pulchra* from **Mahana, Mo'orea**, measured by CBASS rapid acute heat-stress
assay (~18 h: ramp → 3 h hold → ramp down). This file is external context only:
do **not** join it to LTH source patches A/C/D until SNP matching confirms which
Cunning genet, if any, each LTH source patch matches.

**Source:** Cunning R, Matsuda SB, Bartels E, D'Alessandro M, Detmer AR,
Harnay P, Levy J, Lirman D, Moeller HV, Muller EM, Nedimyer K, Pfab F,
Putnam HM (2024). *On the use of rapid acute heat tolerance assays to resolve
ecologically relevant differences among corals.* **Coral Reefs**.
doi:10.1007/s00338-024-02577-7. Data: github.com/jrcunning/CBASS_methods
(`data/classic_cbass/processed/ed50.csv`; collection site "mahana" per
`collection_metadata.csv`, collected 2022-12-03, CBASS 2022-12-04).

**Columns copied:** `geno` genet ID (1–20); `ed50`/`std.error` ED50 (°C) and SE from
the primary dose-response fit; `ed50.f`/`std.error.f` ED50 (°C) and SE from the
alternate refined fit.

**Why it's here (LTH project #17):** Same species and island as the LTH
heat × wound experiment. Used by `code/sensitivity/26_thermal_context.R` to (a)
show how far the LTH chronic treatments (28 °C, 31 °C) sit below acute Fv/Fm
ED50 values and (b) compare the amount of among-genotype thermal-tolerance
variation detected by acute (CBASS) and chronic (LTH) methods.

**Caveats.** ED50 is an **acute** (18 h) threshold, while the LTH experiment is
**chronic** (weeks at +3 °C). Use ED50 only as a reference for how far 28 °C and
31 °C are below acute limits; it is not the temperature scale of the chronic
experiment. The LTH source-patch labels (A, C, D) are **field labels** and are not
genotype-matched to Cunning's numbered genets. Correlating individual genets
requires matching LTH RNA-seq SNPs to Cunning's `genet_map` first.
