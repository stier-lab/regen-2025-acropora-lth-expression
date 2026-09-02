# LTH RNA-seq — analysis notes

**Project:** LTH #17 — heat × wound × *Acropora pulchra*, Mahana/Tiahura, Mo'orea 2025 (*Expression by Temperature*).
This is the consolidated brief; the repo root `README.md` is the index for running things.

---

## 1. Purpose & status

The RNA-seq analysis has not been run yet. The **phenotype analysis** (Methods + Results — physiology, morphology, growth, genet variation, thermal context) is finished and reproduces from one command; the narrative (Introduction, Discussion, Abstract) is out of scope here. Treat this file as background for the analyst, not a prescribed pipeline.

- **RNA-seq status:** libraries await sequencing — this is the open genomics piece.
- **Phenotype pipeline:** `Rscript code/_run_all.R` (~4 min) rebuilds every figure and table.
- **Single source of truth for every statistic:** `output/tables/20_master_results.csv` (each row: effect size + test stat + df + p + CI; `_paper_ready.csv` is the formatted version). Cite the table; don't hand-copy.
- **Audit:** `code/30_manuscript_audit.R` recomputes the phenotype numbers and flags drift (advisory, phenotype only; 15/15 checks pass). It does not check the Intro/Discussion/Abstract or transcriptomics.
- Phenotype Results are already written in `manuscript/Manuscript_LTH.docx`.

### Where things live

| Need | Path |
|---|---|
| Orientation / how to run | `README.md` (repo root) |
| Full phenotype results narrative (incl. §10 limitations) | `RESULTS.docx` |
| Every statistic (source of truth) | `output/tables/20_master_results.csv` |
| RNA-seq plate layout / sequencing plan | `notes/archive/sequencing-plan-keck-LTH.md`, `notes/archive/Experimental_Plan_Gene_Expression.md`, `code/plate_fig.R` |
| Plate maps / IDs | `notes/LTH_PlateLayout_with_IDs`, `data/raw/plate_layout/` |
| Per-library phenotype covariate table | `output/tables/31_rnaseq_phenotype_covariates.csv` (built by `code/31_rnaseq_covariate_table.R`) |
| Raw un-recoded library lookup | `output/tables/31_rnaseq_library_lookup_raw.csv` |
| Tentative expression-phenotype integration plan | `docs/rnaseq/expression_integration_analysis_plan.md` |
| Coral expression literature synthesis | `docs/rnaseq/coral_expression_literature_synthesis.md` |
| Candidate gene evidence table | `docs/rnaseq/candidate_genes_reference.csv` |
| Literature search log | `docs/rnaseq/literature_search_log.md` |
| Per-genet resilience scores | `output/tables/19_genet_resilience_summary.csv` |
| Cached Cunning 2024 ED50 | `data/external/cunning2024_apulchra_ed50.csv` |
| Thicket GPS / metadata | `data/raw/metadata/metadata.csv` |
| Verified references + cite→PDF index | `manuscript/references.bib`, `literature/LITERATURE.docx` |
| Literature library (101 PDFs, mapped) | `literature/LITERATURE.docx` |

---

## 2. Sequencing / library design (as built)

**144 libraries = 4 tanks × 3 genets × 2 wound states × 3 days × 2 temperatures.**

- **Temperature:** 28 °C (ambient) vs 31 °C (chronic heat, weeks). Tanks — 28 °C: 3, 6, 9, 12; 31 °C: 4, 5, 10, 11 (4 tanks/temp, 8 total).
- **Wound:** Wounded (about 1 cm clipped from the growing tip) vs no-wound control. Verify the exact no-wound handling before calling this a sham control in the manuscript.
- **Thicket/genet label:** 3 field-collected parent thickets A, C, D (field labels used as genet proxies; not yet matched to Cunning's numbered genets).
- **Tissue:** wound-margin (M) biopsy only, vs matched unwounded margin.
- **Timepoints:** Days 1 (24 h post-wound), 3, 10. Day 15 samples exist but are **not** in the 144-library analysis design.
- **Counts:** per temp × day × wound = 4 tanks × 3 genets = 12; per temp × day (W+U) = 24; per day (both temps) = 48; × 3 days = **144 libraries**.

**Plating / batch design.** Two 96-well plates (72 primary coral libraries + 8 fixed controls + up to 8 optional anchors each). Each Temp × Day × Wound set of 12 splits 6/6 across the plates, so temperature, wound, day, tank, and genotype are all orthogonal to plate/lane; genotype→plate assignment was rotated/flipped across days. Both plates pooled on **one NovaSeq run**. Fixed controls per plate: 2 extraction NTCs, 2 library blanks, 2 ERCC/spike-in standards, 2 cross-plate technical-replicate anchors.

**Wet-lab.** EZNA Total RNA kit (from DNA/RNA-shield-preserved samples) → NEBNext Ultra Directional RNA Library prep (¼ reactions), on a single **NovaSeq 25B lane** (3.1 billion read pairs; expected read depth depends on the final number of libraries and controls). Provisional QC targets to confirm with the Bay lab: RNA input 50–100 ng, RIN ≥ 7 or DV200 > 70%, ≥ 20 M paired-end reads/library, ≥ 50% of reads mapping to host; record the host:symbiont ratio. State in the final metadata whether symbiont reads were removed, modeled separately, or retained for a parallel symbiont analysis.

**Sampling timeline (2025):** collect/glue May 18–19; acclimate 28 °C May 19–26; ramp to 32 °C May 27–Jun 2; pre-wound PAM Jun 3; Day 0 wounding + buoyant weight Jun 4; Day 1 Jun 5; Day 3 Jun 7; Day 10 Jun 14; Day 15 Jun 19.

---

## 3. Analysis proposal & genet-matching plan

The analyst chooses the expression tool and normalization, but the minimum design constraints are fixed: include temperature, wound state, day, thicket/genet label, tank/blocking, plate, extraction batch, and library-prep batch in the design or in documented QC checks. The phenotype results raise the questions below; the expression data can **test, extend, or revise** each. We name predicted *processes*, not gene symbols. Start genome-wide, then use the gene sets for interpretation.

**DE questions (phenotype-anchored):**
1. **Does heat suppress the regeneration program while sparing healing?** Wounds close equally at both temperatures, but new-corallite formation is delayed or suppressed under heat (interval-censored time ratio = 1.32, p = 1.4e-7; first-observed Cox HR = 0.22). Test whether early-healing processes (re-epithelialization, immune, ECM remodeling, proliferation) turn on similarly at both temperatures while skeletal/biomineralization + corallite-patterning are lower under heat — or whether heat suppresses healing too, which would revise the phase-decoupling story. Main comparison: the temperature-by-day contrast in wounded margins, with Day 10 prioritized because new-corallite formation diverges most strongly then.
2. **What distinguishes resilient genet C from sensitive A, D?** C holds onto photochemistry, pigmentation, and symbionts far better (heat effect 2.7–3.5× weaker than A/D). Test for a smaller heat-induced shift and/or higher expression of stress-response genes before or early in the heat response. Because the current 144-library design is D1/D3/D10, do not call this a D0 baseline test unless D0 libraries are added. Main comparison: genet C vs genets A/D by temperature, prioritizing unwounded margins where the genet spread is sharpest. A WGCNA module tracking `19_genet_resilience_summary.csv` is one route.
3. **Sustained response, not an acute heat-shock spike?** 31 °C sits ~4.4 °C below the acute Fv/Fm ED50 (35.4 °C; Cunning et al. 2024), and wounding came after 7 days at temperature, so the short-lived HSP burst may have faded. Expect expression that reflects maintained stress or acclimation after chronic exposure, not only a quick HSP70 burst.
4. **Does wounding change genet differences?** The genet heat-sensitivity spread remains strong in both scopes: unwounded margins (A = 0.99, D = 0.90, C = 0.43) and wounded margins/morphology (A = 0.56, D = 0.44, C = -0.05). Test whether the genet effect and genet-by-temperature effect differ by wound state. Symbiont (*Symbiodiniaceae*) reads, if retained after host/symbiont mapping, could corroborate the symbiont-density loss and C's retention.

**Genet-matching (A/C/D ↔ Cunning genets) — high-value if the SNP match works.** Cunning et al. 2024 (*Coral Reefs*, doi:10.1007/s00338-024-02577-7) measured **acute CBASS Fv/Fm ED50 for 20 genotyped *A. pulchra* genets from Mahana** (range 34.4–36.6 °C; ED50 predicts bleaching, R = 0.74; collected Dec 2022, "mahana"). We measured **chronic** resilience for 3 thickets from the same site. The test: call genotype-distinguishing SNPs from the host RNA-seq reads, match A/C/D to Cunning's reference, then ask whether **acute CBASS ED50 predicts our chronic wound-context ranking (C > D > A)** — a cross-method validation and a link to the genet-by-temperature expression difference.

- **Supporting GPS** (`data/raw/metadata/metadata.csv`, `coord_lat`/`coord_long`): A = 17.49735 °S, 149.91557 °W (72 frags); C = 17.49808, 149.91595 (72); D = 17.49726, 149.91581 (64). They sit ~40–90 m apart, all in the Mahana/Tiahura stand Cunning sampled. Proximity is **suggestive, not conclusive** (*A. pulchra* forms clonal thickets).
- **External ask we can chase:** Cunning's per-genet host SNP genotypes (not just ED50 + genet number) — from the CBASS_methods repo (`github.com/jrcunning/CBASS_methods`, `data/reproducibility/genet_map.xlsx`) or by asking Cunning/Putnam directly (co-authors Detmer & Moeller are in the UCSB/Mo'orea network). *How* (reference genome, variant caller, identity metric) is yours.
- **Fallback:** if genotypes can't be matched, the population-level statement still holds (already in the manuscript): both acute (CBASS) and chronic (LTH) methods independently detect substantial heritable thermal-tolerance variation in Mahana *A. pulchra*.

---

## 4. Phenotype ↔ expression integration map

Every link runs both ways (confirm / extend / revise). Timecourse alignment:

| Day | Phenotype state | Expression snapshot likely captures |
|---|---|---|
| **1** | wound fresh; healing starting | early wound response (ECM remodeling, immune, proliferation) |
| **3** | coenosarc coverage advancing (≥90% by D5–7) | healing program at peak; regeneration not yet engaged |
| **10** | new-corallite formation diverges by temp | healing→regeneration transition — where heat bites |

| # | Phenotype anchor | Suggested contrast |
|---|---|---|
| 1 | Healing proceeds, regeneration stalls under heat (time ratio 1.32, 95% CI 1.19–1.47; Cox HR 0.22, 95% CI 0.07–0.69; 67% closed-but-never-regenerated vs 0%); diverges by Day 10 | **Temp × timepoint, wounded margins** — genes/modules up at Day 10 in 28 °C not 31 °C. Refuted if healing (D1–3) is also suppressed |
| 2 | Genet C resilience (multivariate displacement 1.02 vs 3.71 A / 3.38 D); spread largest unwounded | **Genet (C vs A, D) × temperature, unwounded margins**; WGCNA vs resilience score. Barshis 2013 frontloading |
| 3 | Wounding homogenizes genet response (spread compresses when wounded) | **Genet × wound** — genet + genet×temp effect smaller in wounded margins |
| 4 | Wound response is local (whole-colony physiology tracks temp, barely wound) | **Wound main effect at margin, per day** — positive control that margin sampling worked |
| 5 | Chronic/constitutive, not acute HSP burst (31 °C is ~4.4 °C below ED50; wound applied after 7 d at temp) | interpretation cue: expect sustained/frontloaded signatures over transient HSP70 |
| 6 | A/C/D not yet matched to Cunning CBASS genets (ED50 34.3–36.6 °C, range 2.2 °C) | SNP-call → match → test acute ED50 vs chronic ranking C > D > A |

**Suggested entry points:** sample QC and PCA/MDS batch checks → wound main effect as a positive control (Anchor 4) → temperature-by-day in wounded margins (Anchor 1) → genet C vs A/D by temperature in unwounded margins (Anchor 2) → WGCNA modules vs phenotype axes → SNP calling + Cunning matching (Anchor 6, in parallel). Per-coral covariates are in `output/tables/31_rnaseq_phenotype_covariates.csv` (joined by `Fragment_ID`).

---

## 5. Candidate genes and literature synthesis

Use the literature as prior context, not as a target panel. Do not filter the
DE analysis to these genes.

Current evidence package:

- Concise synthesis: `docs/rnaseq/coral_expression_literature_synthesis.md`
- Machine-readable evidence table: `docs/rnaseq/candidate_genes_reference.csv`
- Search provenance and follow-up gaps: `docs/rnaseq/literature_search_log.md`

Highest-value candidate processes for LTH:

| LTH question | Candidate processes | Most useful contrast |
|---|---|---|
| Did the sampled margin capture a wound signal? | c-Fos/AP-1, JNK/MAPK, TLR/NOD/NF-kB, ECM remodeling, oxidative response | wounded vs unwounded at Day 1 and Day 3 |
| Does heat spare closure but block regeneration? | Wnt, FGF, calicoblast fate, galaxin/SAARP/mucin, carbonic anhydrase, Ca2+ handling, skeletal matrix proteins | 31 C vs 28 C in wounded Day 10 margins |
| Is the heat response chronic rather than an acute spike? | HSPs, ubiquitin-proteasome, redox, mitochondria, metabolism, cell adhesion, immune/apoptosis | 31 C vs 28 C in unwounded margins across Day 1, 3, 10 |
| Does source C differ from A/D? | frontloaded or dampened stress modules, redox/proteostasis, symbiosis, cell adhesion, ion transport | temperature x source-thicket label, with SNP clusters added when available |
| Do organismal traits map to expression? | WGCNA modules for stress, symbiosis, pigment, calcification, growth, regeneration | module eigengenes vs PAM, color, symbionts, growth, and morphology milestones |

Evidence tiers in the CSV:

| Tier | Meaning |
|---|---|
| A | Direct coral expression evidence from Acropora or an LTH-like heat/injury design |
| B | Direct coral expression evidence from another genus, life stage, or related stress design |
| C | Review, meta-analysis, genomic association, or phenotype-only support |
| D | Lab-only, contradicted, or not yet source-verified |

Key guardrail: HSP70/HSP90 are useful but not sufficient. In chronic,
sublethal heat, the signal may appear as dampening, altered recovery, redox or
proteostasis load, metabolism, symbiosis, or cell-adhesion changes rather than a
simple acute heat-shock spike.

---

## 6. Optional Intro/Discussion framing angles

Phenotype-first raw material, written before your expression results exist — **not the spine** (the recommended spine leads with transcriptomics). Take one, blend them, or discard them.

- **A. Energetic-triage / phase-decoupling (phenotype-first).** Heat spares tissue-healing (coenosarc closure) but stalls regeneration (skeletal regrowth). Hook: survival is a poor proxy for structural recovery — a restoration angle (picking genets on bleaching survival may propagate corals that persist without regenerating). *Tension:* this is an allocation inference the phenotype alone can't prove.
- **B. Molecular basis of the healing→regeneration transition (transcriptomics-led).** The phenotype locates *where* recovery fails; the expression data ask *why*. Leads with mechanism, not allocation — likely the stronger spine.
- **C. Genotype-dependent thermal tolerance of regeneration (cross-method).** Regenerative capacity is heritable (C vs A/D); it pairs with the independent acute CBASS assay (Cunning et al. 2024) and the candidate thermal-tolerance genes. Can sit under A or B, or stand alone as a Discussion section.

**Phenotype numbers these rest on** (→ `output/tables/20_master_results.csv`):
- Coenosarc closure indistinguishable 28/31 °C; new-corallite regeneration delayed under heat (time ratio 1.32, 95% CI 1.19–1.47, p = 1.4e-7; first-observed Cox HR 0.22, 95% CI 0.07–0.69).
- Closed-but-never-regenerated: 67% at 31 °C vs 0% at 28 °C (median healing→regeneration lag 10 vs 8 d).
- Whole-colony physiology responds to temperature, barely to wounding; skeletal growth (% mass change) **34% lower** at 31 °C.
- Genet C multivariate displacement 1.02 vs 3.71 (A) / 3.38 (D); most likely to regenerate at 31 °C.
- 31 °C sits ~4.4 °C below the population acute ED50 (35.4 °C; Cunning et al. 2024) — chronic/sublethal.

**Fragility to carry forward:** the regeneration result is strongest for interval-censored new-corallite onset (p = 1.4e-7); tip-exist is also delayed, and tip-extension points the same way but is n.s. The first-observed Cox model sits near the PH diagnostic boundary, so the interval model + censored fraction are the cleaner anchors. Three genets resolve variation but not architecture. Apical-tip excision ≠ surface wound bed (Munk 2024). See `RESULTS.docx` §10. Chlorophyll-a was **not** run (the metadata slot is kept for provenance; analysis uses PAM, color-card scores, and symbiont counts).

**Verified citation bank** (DOIs in `manuscript/references.bib`; index `literature/LITERATURE.docx`):
- *Biphasic healing↔regeneration; growth & injury ecology of *Acropora*:* Henry & Hart 2005 (closure precedes regeneration); Yap & Gomez 1984 (*A. pulchra* extension 13–16 cm yr⁻¹); Highsmith 1982; Madin et al. 2014.
- *Thermal sensitivity of *Acropora* / *A. pulchra*:* Hoegh-Guldberg 1999; Hughes et al. 2017; Berg et al. 2020 (persistent photosystem damage under sustained heat).
- *Heat drains the energy budget:* Warner et al. 1999 (PSII damage); Hoegh-Guldberg 1999; Jokiel & Coles 1990; Jokiel & Coles 1977; Comeau et al. 2014 (calcification suppression).
- *Prior heat × injury (recovery as a single rate — the gap):* Meesters & Bak 1993; Bonesso et al. 2017; Traylor-Knowles et al. 2016.
- *Heritable genotype-level tolerance & cross-method link:* Dixon et al. 2015; Shaw et al. 2016; Cunning et al. 2024 (CBASS Fv/Fm ED50, the acute assay to match).
- *Wound-geometry caveat:* Munk 2024 (surface wound-bed polyp reappearance ≠ apical-tip excision).
